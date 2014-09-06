#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# package import
import permutation2020.python.permutation as pm
import permutation2020.python.utils as utils
from permutation2020.python.gene_sequence import GeneSequence
from permutation2020.python.sequence_context import SequenceContext

# external imports
import pandas as pd
import pysam
from multiprocessing import Pool
import argparse
import datetime
import logging

logger = logging.getLogger(__name__)  # module logger

def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def multiprocess_permutation(bed_dict, mut_df, opts):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys())
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    num_permutations = opts['num_permutations']
    result_list = [[0, 0] for k in range(num_permutations)]
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for chrom_result in process_results:
                    for j in range(num_permutations):
                        result_list[j][0] += chrom_result[j][0]
                        result_list[j][1] += chrom_result[j][1]
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            info = (bed_dict[chroms[i]], mut_df, opts)
            chrom_result = singleprocess_permutation(info)
            for j in range(num_permutations):
                result_list[j][0] += chrom_result[j][0]
                result_list[j][1] += chrom_result[j][1]

    return result_list


@utils.log_error_decorator
def singleprocess_permutation(info):
    bed_list, mut_df, opts = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    num_permutations = opts['num_permutations']
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    result = [[0, 0] for k in range(num_permutations)]
    for bed in bed_list:
        # prepare info for running permutation test
        gene_mut = mut_df[mut_df['Gene']==bed.gene_name]
        cols = ['Chromosome', 'Start_Position', 'Reference_Allele',
                'Tumor_Allele', 'Variant_Classification', 'Protein_Change']
        mut_info = gene_mut[cols]
        gs.set_gene(bed)
        sc = SequenceContext(gs)

        # count total mutations in gene
        total_mut = len(mut_info)

        # fix nucleotide letter if gene is on - strand
        if bed.strand == '-':
            mut_info['Tumor_Allele'].map(lambda x: utils.rev_comp(x))

        # get coding positions, mutations unmapped to the reference tx will have
        # NA for a coding position
        pos_list = []
        for ix, row in mut_info.iterrows():
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], row['Start_Position'])
            pos_list.append(coding_pos)
        mut_info['Coding Position'] = pos_list

        # recover mutations that could not be mapped to the reference transcript
        # for a gene before being dropped (next step)
        unmapped_mut_info = utils.recover_unmapped_mut_info(mut_info, bed, sc, opts)

        # drop mutations wich do not map to reference tx
        mut_info = mut_info.dropna(subset=['Coding Position'])  # mutations need to map to tx
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        unmapped_muts = total_mut - len(mut_info)

        if len(mut_info) > 0:
            mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
            mut_info['Context'] = mut_info['Coding Position'].apply(lambda x: sc.pos2context[x])

            # group mutations by context
            cols = ['Context', 'Tumor_Allele']
            unmapped_mut_df = pd.DataFrame(unmapped_mut_info)
            tmp_df = pd.concat([mut_info[cols], unmapped_mut_df[cols]])
            context_cts = tmp_df['Context'].value_counts()
            context_to_mutations = dict((name, group['Tumor_Allele'])
                                        for name, group in tmp_df.groupby('Context'))

            # calculate non silent count
            tmp_result = pm.non_silent_ratio_permutation(context_cts,
                                                         context_to_mutations,
                                                         sc,  # sequence context obj
                                                         gs,  # gene sequence obj
                                                         num_permutations)
        else:
            tmp_result = [[0, 0] for k in range(num_permutations)]

        # increment the non-silent/silent counts for each permutation
        for j in range(num_permutations):
            result[j][0] += tmp_result[j][0]
            result[j][1] += tmp_result[j][1]

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def parse_arguments():
    # make a parser
    info = 'Simulates the non-silent mutation ratio by randomly permuting mutations'
    parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts "stdout")')

    # program arguments
    help_str = 'gene FASTA file from extract_gene_seq.py script'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'DNA mutations file'
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help=help_str)
    help_str = 'BED file annotation of genes'
    parser.add_argument('-b', '--bed',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Number of processes to use. 0 indicates using a single '
                'process without using a multiprocessing pool '
                '(more means Faster, default: 0).')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Number of permutations for null model. p-value precision '
                'increases with more permutations (Default: 10000).')
    parser.add_argument('-n', '--num-permutations',
                        type=int, default=10000,
                        help=help_str)
    help_str = ('Number of DNA bases to use as context. 0 indicates no context. '
                '1 indicates only use the mutated base.  1.5 indicates using '
                'the base context used in CHASM '
                '(http://wiki.chasmsoftware.org/index.php/CHASM_Overview). '
                '2 indicates using the mutated base and the upstream base. '
                '3 indicates using the mutated base and both the upstream '
                'and downstream bases. (Default: 1.5)')
    parser.add_argument('-c', '--context',
                        type=float, default=1.5,
                        help=help_str)
    help_str = ('Use mutations that are not mapped to the the single reference '
                'transcript for a gene specified in the bed file indicated by '
                'the -b option.')
    parser.add_argument('-u', '--use-unmapped',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Path to the genome fasta file. Required if --use-unmapped flag '
                'is used. (Default: None)')
    parser.add_argument('-g', '--genome',
                        type=str, default='',
                        help=help_str)
    help_str = 'Output text file of results'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    log_level = args.log_level
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging

    opts = vars(args)
    if opts['use_unmapped'] and not opts['genome']:
        print('You must specify a genome fasta with -g if you set the '
              '--use-unmapped flag to true.')
        sys.exit(1)

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts):
    # hack to index the FASTA file
    gene_fa = pysam.Fastafile(opts['input'])
    gene_fa.close()

    # Get Mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df)

    # perform permutation test
    bed_dict = utils.read_bed(opts['bed'], [])
    permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)

    # convert to dataframe to save to file
    non_silent_ratio_df = pd.DataFrame(permutation_result,
                                       columns=['non-silent count',
                                                'silent count'])
    # save output
    non_silent_ratio_df.to_csv(opts['output'], sep='\t', index=False)

    return non_silent_ratio_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
