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
import permutation2020.cython.cutils as cutils
import permutation2020.python.mutation_context as mc
import permutation2020.python.indel as indel

# external imports
import numpy as np
import pandas as pd
import pysam
import csv
from multiprocessing import Pool
import argparse
import datetime
import logging
import copy
import itertools as it
import IPython

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


def multiprocess_permutation(bed_dict, mut_df, opts, indel_df=None):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys(), key=lambda x: len(bed_dict[x]), reverse=True)
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    file_handle = open(opts['output'], 'w')
    mywriter = csv.writer(file_handle, delimiter='\t', lineterminator='\n')
    if opts['maf']:
        header = ['gene', 'strand', 'Chromosome', 'Start_Position',
                  'End_Position', 'Reference_Allele', 'Tumor_Allele',
                  'Context', 'DNA_Change', 'Protein_Change', 'Variant_Classification']
    else:
        header = ['gene', 'ID', 'non-silent snv', 'silent snv', 'nonsense', 'lost stop',
                  'splice site', 'lost start', 'missense', 'recurrent missense',
                  'normalized missense position entropy', 'frameshift indel',
                  'inframe indel']
    mywriter.writerow(header)
    num_permutations = opts['num_permutations']

    # simulate indel counts
    if not opts['maf']:
        fs_cts, inframe_cts, gene_names = indel.simulate_indel_counts(indel_df,
                                                                      bed_dict,
                                                                      num_permutations)
        name2ix = {gene_names[z]: z for z in range(len(gene_names))}

    # simulate snvs
    obs_result = []
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                # iterate through each chromosome result
                for chrom_result in process_results:
                    # add columns for indels
                    if not opts['maf']:
                        tmp_chrom_result = []
                        for gname, grp in it.groupby(chrom_result, lambda x: x[0]):
                            for l, row in enumerate(grp):
                                gene_ix = name2ix[gname]
                                fs_count = fs_cts[l, gene_ix]
                                inframe_count = inframe_cts[l, gene_ix]
                                tmp_chrom_result.append(row+[fs_count, inframe_count])
                        chrom_result = tmp_chrom_result

                    # write output to file
                    mywriter.writerows(chrom_result)
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            # perform simulation
            info = (bed_dict[chroms[i]], mut_df, opts)
            chrom_results = singleprocess_permutation(info)

            # add indel columns
            if not opts['maf']:
                tmp_chrom_result = []
                for gname, grp in it.groupby(chrom_results, lambda x: x[0]):
                    for l, row in enumerate(grp):
                        gene_ix = name2ix[gname]
                        fs_count = fs_cts[l, gene_ix]
                        inframe_count = inframe_cts[l, gene_ix]
                        tmp_chrom_result.append(row+[fs_count, inframe_count])
                chrom_results = tmp_chrom_result

            # write to file
            mywriter.writerows(chrom_results)
    file_handle.close()


@utils.log_error_decorator
def singleprocess_permutation(info):
    bed_list, mut_df, opts = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    num_permutations = opts['num_permutations']
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    # go through each gene to permform simulation
    result = []
    for bed in bed_list:
        # compute context counts and somatic bases for each context
        gene_tuple = mc.compute_mutation_context(bed, gs, mut_df, opts)
        context_cts, context_to_mutations, mutations_df, gs, sc = gene_tuple

        if context_to_mutations:
            ## get information about observed non-silent counts
            # get info about mutations
            tmp_mut_info = mc.get_aa_mut_info(mutations_df['Coding Position'],
                                              mutations_df['Tumor_Allele'].tolist(),
                                              gs)
            # calc mutation info summarizing observed mutations
            tmp_non_silent = cutils.calc_summary_info(tmp_mut_info['Reference AA'],
                                                      tmp_mut_info['Somatic AA'],
                                                      tmp_mut_info['Codon Pos'])
            ## Do permutations
            if opts['maf']:
                # if user specified MAF format then output all mutations in
                # MAF format
                tmp_result = pm.maf_permutation(context_cts,
                                                context_to_mutations,
                                                sc,
                                                gs,
                                                num_permutations)

            else:
                # Summarized results for feature for each simulation for each
                # gene
                tmp_result = pm.summary_permutation(context_cts,
                                                    context_to_mutations,
                                                    sc,  # sequence context obj
                                                    gs,  # gene sequence obj
                                                    num_permutations)
            result += tmp_result

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def parse_arguments():
    # make a parser
    info = 'Simulates by randomly permuting mutation positions. Saves results to file'
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
                'increases with more permutations (Default: 1).')
    parser.add_argument('-n', '--num-permutations',
                        type=int, default=1,
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
    parser_grouper = parser.add_mutually_exclusive_group()
    parser_grouper.add_argument('--summary',
                                action='store_true',
                                default=True,
                                help='Flag for saving results as summarized '
                                'features (Default: True).')
    parser_grouper.add_argument('--maf',
                                action='store_true',
                                default=False,
                                help='Flag for saving results in MAF format '
                                '(Default: False).')
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
    help_str = ('Specify the seed for the pseudo random number generator. '
                'By default, the seed is randomly chosen based. The seed will '
                'be used for the monte carlo simulations.')
    parser.add_argument('-seed', '--seed',
                        type=int, default=None,
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
    indel_df = indel.keep_indels(mut_df)  # return indels only
    indel_df['Start_Position'] = indel_df['Start_Position'] - 1  # convert to 0-based
    indel_df['indel len'] += 1
    logger.info('There were {0} indels identified.'.format(len(indel_df)))
    mut_df = mut_df.dropna(subset=['Tumor_Allele', 'Start_Position', 'Chromosome'])
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df)

    # read in bed info
    bed_dict = utils.read_bed(opts['bed'], [])

    # perform permutation
    multiprocess_permutation(bed_dict, mut_df, opts, indel_df)

    # save indels
    if opts['maf']:
        with open(opts['output'], 'a') as handle:
            mywriter = csv.writer(handle, delimiter='\t')
            for maf_lines in indel.simulate_indel_maf(indel_df, bed_dict,
                                                      opts['num_permutations']):
                mywriter.writerows(maf_lines)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
