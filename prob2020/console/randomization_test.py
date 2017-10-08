#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

# package imports
import prob2020.python.utils as utils
from prob2020.python.gene_sequence import GeneSequence
from prob2020.python.sequence_context import SequenceContext
import prob2020.python.mutation_context as mc
import prob2020.python.count_frameshifts as cf
import prob2020.python.process_result as pr
import prob2020.python.p_value as mypval

# external imports
import argparse
import pysam
import pandas as pd
import numpy as np
from multiprocessing import Pool
import logging

logger = logging.getLogger(__name__)  # module logger


@utils.log_error_decorator
def singleprocess_permutation(info):
    # initialize input
    bed_list, mut_df, opts, fs_cts_df, p_inactivating = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    # list of columns that are needed
    cols = ['Chromosome', 'Start_Position', 'Reference_Allele',
            'Tumor_Allele', 'Variant_Classification',]
    # conditionally add protein_change column if exists
    if 'Protein_Change' in mut_df.columns:
        cols += ['Protein_Change']

    # figure out which genes actually have a mutation
    genes_with_mut = set(mut_df['Gene'].unique())

    # iterate through each gene
    result = []
    for bed in bed_list:
        if bed.gene_name not in genes_with_mut:
            # skip genes with no mutations
            continue

        # prepare info for running permutation test
        mut_info = mut_df.loc[mut_df['Gene']==bed.gene_name, cols]
        gs.set_gene(bed)
        sc = SequenceContext(gs, seed=opts['seed'])

        # count total mutations in gene
        total_mut = len(mut_info)

        # fix nucleotide letter if gene is on - strand
        if bed.strand == '-':
            rc = mut_info['Tumor_Allele'].map(lambda x: utils.rev_comp(x))
            mut_info.loc[:, 'Tumor_Allele'] = rc

        # get coding positions, mutations unmapped to the reference tx will have
        # NA for a coding position
        pos_list = []
        for ix, row in mut_info.iterrows():
            coding_pos = bed.query_position(bed.strand, row['Chromosome'], row['Start_Position'])
            pos_list.append(coding_pos)
        mut_info.loc[:, 'Coding Position'] = pos_list

        # recover mutations that could not be mapped to the reference transcript
        # for a gene before being dropped (next step)
        unmapped_mut_info = mc.recover_unmapped_mut_info(mut_info, bed, sc, opts)

        # drop mutations wich do not map to reference tx
        mut_info = mut_info.dropna(subset=['Coding Position'])  # mutations need to map to tx
        mut_info['Coding Position'] = mut_info['Coding Position'].astype(int)
        num_mapped_muts = len(mut_info)
        unmapped_muts = total_mut - num_mapped_muts

        # construct sequence context
        #gs.add_germline_variants(mut_info['Reference_Allele'].tolist(),
        #                         mut_info['Coding Position'].tolist())

        # calculate results of permutation test
        if opts['kind'] == 'oncogene':
            # calculate position based permutation results
            tmp_result = mypval.calc_position_p_value(mut_info, unmapped_mut_info, sc,
                                                      gs, bed, opts['score_dir'],
                                                      opts['num_iterations'],
                                                      opts['stop_criteria'],
                                                      0,  # no recurrent mutation pseudo count
                                                      opts['recurrent'],
                                                      opts['fraction'])
            result.append(tmp_result + [total_mut, unmapped_muts])
        elif opts['kind'] == 'tsg':
            # calculate results for deleterious mutation permutation test
            #fs_ct = fs_cts_df['total'][bed.gene_name]
            #fs_unmapped = fs_cts_df['unmapped'][bed.gene_name]
            # replaced fs_ct with zero to stop using the frameshifts in
            # simulation
            tmp_result = mypval.calc_deleterious_p_value(mut_info, unmapped_mut_info,
                                                         sc, gs, bed,
                                                         opts['num_iterations'],
                                                         opts['stop_criteria'],
                                                         opts['deleterious'],
                                                         0,  # no deleterious mutation pseudo count
                                                         opts['seed'])
            result.append(tmp_result + [num_mapped_muts, unmapped_muts])
                                        #fs_ct, fs_unmapped])
        elif opts['kind'] == 'hotmaps1d':
            # save null distribution if user option specified
            if opts['null_distr_dir']:
                if not os.path.exists(opts['null_distr_dir']): os.mkdir(opts['null_distr_dir'])
                save_path = os.path.join(opts['null_distr_dir'], bed.gene_name + '.{0}.txt')
            else:
                save_path = None
            # calculate position based permutation results
            mywindow = list(map(int, opts['window'].split(',')))
            tmp_result = mypval.calc_hotmaps_p_value(mut_info, unmapped_mut_info, sc,
                                                     gs, bed,
                                                     mywindow,
                                                     opts['num_iterations'],
                                                     opts['stop_criteria'],
                                                     opts['report_index'],
                                                     null_save_path=save_path)
            result.extend(tmp_result)
        elif opts['kind'] == 'protein':
            tmp_result = mypval.calc_protein_p_value(mut_info, unmapped_mut_info,
                                                     sc, gs, bed,
                                                     opts['neighbor_graph_dir'],
                                                     opts['num_iterations'],
                                                     opts['stop_criteria'],
                                                     opts['recurrent'],
                                                     opts['fraction'])
            result.append(tmp_result + [total_mut, unmapped_muts])
        else:
            # calc results for entropy-on-effect permutation test
            tmp_result = mypval.calc_effect_p_value(mut_info, unmapped_mut_info,
                                                    sc, gs, bed,
                                                    opts['num_iterations'],
                                                    0, #  no recurrent mutation pseudo count
                                                    opts['recurrent'],
                                                    opts['fraction'])
            result.append(tmp_result + [total_mut, unmapped_muts])

    gene_fa.close()
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result


def multiprocess_permutation(bed_dict, mut_df, opts,
                             fs_cts_df=None, p_inactivating=None):
    """Handles parallelization of permutations by splitting work
    by chromosome.
    """
    chroms = sorted(bed_dict.keys(), key=lambda x: len(bed_dict[x]), reverse=True)
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    result_list = []
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts, fs_cts_df, p_inactivating)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for chrom_result in process_results:
                    result_list += chrom_result
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            info = (bed_dict[chroms[i]], mut_df, opts, fs_cts_df, p_inactivating)
            result_list += singleprocess_permutation(info)

    return result_list


def parse_arguments():
    # make a parser
    info = 'Performs a randomization-based test on the oncogene and TSG score'
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
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

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
    help_str = 'Directory containing score information in pickle files (Default: None).'
    parser.add_argument('-s', '--score-dir',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Directory containing neighbor graph information in pickle files (Default: None).'
    parser.add_argument('-ng', '--neighbor-graph-dir',
                        type=str, default=None,
                        help=help_str)
    help_str = ('Number of processes to use. 0 indicates using a single '
                'process without using a multiprocessing pool '
                '(more means Faster, default: 0).')
    parser.add_argument('-p', '--processes',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Number of iterations for null model. p-value precision '
                'increases with more iterations, however this will also '
                'increase the run time (Default: 10000).')
    parser.add_argument('-n', '--num-iterations',
                        type=int, default=10000,
                        help=help_str)
    help_str = ('Number of iterations more significant then the observed statistic '
                'to stop further computations. This decreases compute time spent in resolving '
                'p-values for non-significant genes. (Default: 1000).')
    parser.add_argument('-sc', '--stop-criteria',
                        type=int, default=1000,
                        help=help_str)
    help_str = ('Kind of permutation test to perform ("oncogene" or "tsg"). "position-based" permutation '
                'test is intended to find oncogenes using position based statistics. '
                'The "deleterious" permutation test is intended to find tumor '
                'suppressor genes. (Default: oncogene)')
    parser.add_argument('-k', '--kind',
                        type=str, default='oncogene',
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
    help_str = ('Only keep unique mutations for each tumor sample.'
                'Mutations reproted from heterogeneous sources may contain'
                ' duplicates, e.g. a tumor sample was sequenced twice.')
    parser.add_argument('--unique',
                        action='store_true',
                        default=False,
                        help=help_str)
    help_str = ('Minimum number of mutations at a position for it to be '
                'considered a recurrently mutated position (Default: 3).')
    parser.add_argument('-r', '--recurrent',
                        type=int, default=3,
                        help=help_str)
    help_str = ('Fraction of total mutations in a gene. This define the '
                'minimumm number of mutations for a position to be defined '
                'as recurrently mutated (Defaul: .02).')
    parser.add_argument('-f', '--fraction',
                        type=float, default=.02,
                        help=help_str)
    help_str = ('Perform tsg permutation test if gene has '
                'at least a user specified number of deleterious mutations (default: 1)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=1,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. Values greater than one indicate all '
                'genes will be tested (Default: 1.01).')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=1.01,
                        help=help_str)
    help_str = ('Deleterious mutation pseudo-count for null distribution '
                'statistics. (Default: 0)')
    parser.add_argument('-dp', '--deleterious-pseudo-count',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Recurrent missense mutation pseudo-count for null distribution '
                'statistics. (Default: 0)')
    parser.add_argument('-rp', '--recurrent-pseudo-count',
                        type=int, default=0,
                        help=help_str)
    help_str = ('Specify the seed for the pseudo random number generator. '
                'By default, the seed is randomly chosen based. The seed will '
                'be used for the permutation test monte carlo simulations.')
    parser.add_argument('-seed', '--seed',
                        type=int, default=None,
                        help=help_str)
    help_str = 'Output of probabilistic 20/20 results'
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
    utils.start_logging(log_file=log_file,
                        log_level=log_level,
                        verbose=args.verbose)  # start logging

    opts = vars(args)
    if opts['use_unmapped'] and not opts['genome']:
        print('You must specify a genome fasta with -g if you set the '
              '--use-unmapped flag to true.')
        sys.exit(1)

    # log user entered command
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts, mut_df=None, frameshift_df=None):
    # hack to index the FASTA file
    gene_fa = pysam.Fastafile(opts['input'])
    gene_fa.close()

    # Get Mutations
    if mut_df is None:
        mut_df = pd.read_csv(opts['mutations'], sep='\t')
    orig_num_mut = len(mut_df)

    # rename columns to fit my internal column names
    rename_dict = {
        'Hugo_Symbol': 'Gene',
        'Tumor_Sample_Barcode': 'Tumor_Sample',
        'Tumor_Seq_Allele2' : 'Tumor_Allele'
    }
    mut_df.rename(columns=rename_dict, inplace=True)

    # drop rows with missing info
    na_cols = ['Gene', 'Tumor_Allele', 'Start_Position', 'Chromosome']
    mut_df = mut_df.dropna(subset=na_cols)
    logger.info('Kept {0} mutations after droping mutations with missing '
                'information (Droped: {1})'.format(len(mut_df), orig_num_mut - len(mut_df)))

    # count frameshifts
    if opts['kind'] == 'tsg':
        if frameshift_df is None:
            # read in mutations
            if mut_df is None:
                mut_df = pd.read_csv(opts['mutations'], sep='\t')

            # count number of frameshifts
            frameshift_df = cf.count_frameshift_total(mut_df, opts['bed'],
                                                      opts['use_unmapped'])

        # calculate the proportion of inactivating
        #num_inact = len(mut_df[mut_df['Variant_Classification'].isin(utils.variant_inactivating)])
        #num_non_inact = len(mut_df[mut_df['Variant_Classification'].isin(utils.variant_non_inactivating)])
        num_fs = len(mut_df[mut_df['Variant_Classification'].isin(utils.variant_frameshift)])
        num_all = len(mut_df[mut_df['Variant_Classification'].isin(utils.all_variants)])
        #p_inactivating = float(num_inact) / (num_inact + num_non_inact)
        p_inactivating = float(num_fs) / num_all

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df, opts['unique'])

    # log random number seed choice if provided
    if opts['seed'] is not None:
        logger.info('Pseudo Random Number Generator Seed: {0}'.format(opts['seed']))

    # read BED file
    bed_dict = utils.read_bed(opts['bed'])

    # Perform BH p-value adjustment and tidy up data for output
    if opts['kind'] == 'oncogene':
        permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)
        permutation_df = pr.handle_oncogene_results(permutation_result,
                                                    opts['num_iterations'])
    elif opts['kind'] == 'tsg':
        permutation_result = multiprocess_permutation(bed_dict, mut_df, opts,
                                                      frameshift_df, p_inactivating)
        permutation_df = pr.handle_tsg_results(permutation_result)
    elif opts['kind'] == 'hotmaps1d':
        permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)
                                                      #frameshift_df, p_inactivating)
        permutation_df = pr.handle_hotmaps_results(permutation_result)
    elif opts['kind'] == 'protein':
        permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)
        permutation_df = pr.handle_protein_results(permutation_result)
    elif opts['kind'] == 'effect':
        permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)
        permutation_df = pr.handle_effect_results(permutation_result)

    # save output
    if opts['output']:
        permutation_df.to_csv(opts['output'], sep='\t', index=False)

    return permutation_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
