#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

# package import
import prob2020.python.permutation as pm
import prob2020.python.utils as utils
from prob2020.python.gene_sequence import GeneSequence
import prob2020.cython.cutils as cutils
import prob2020.python.mutation_context as mc

# external imports
import numpy as np
import pandas as pd
import pysam
from multiprocessing import Pool
import argparse
import logging
import copy

logger = logging.getLogger(__name__)  # module logger

# column header for output
cols = ['non-silent count', 'silent count', 'nonsense count',
        'lost stop count', 'splice site count', 'lost start count',
        'missense count']

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
    if not opts['by_sample']:
        obs_result = []
    else:
        uniq_samp = mut_df['Tumor_Sample'].unique()
        obs_result = pd.DataFrame(np.zeros((len(uniq_samp), len(cols))),
                                  index=uniq_samp, columns=cols)

    # initialize list containing output
    if not opts['score_dir']:
        result_list = [[0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]
    else:
        result_list = [[0, 0, 0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]

    # iterate over each chromosome
    for i in range(0, len(chroms), num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            tmp_num_proc = len(chroms) - i if i + num_processes > len(chroms) else num_processes
            info_repeat = ((bed_dict[chroms[tmp_ix]], mut_df, opts)
                            for tmp_ix in range(i, i+tmp_num_proc))
            process_results = pool.imap(singleprocess_permutation, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for chrom_result, obs_mutations in process_results:
                    for j in range(num_permutations):
                        result_list[j][0] += chrom_result[j][0]
                        result_list[j][1] += chrom_result[j][1]
                        result_list[j][2] += chrom_result[j][2]
                        result_list[j][3] += chrom_result[j][3]
                        result_list[j][4] += chrom_result[j][4]
                        result_list[j][5] += chrom_result[j][5]
                        result_list[j][6] += chrom_result[j][6]
                        if opts['score_dir']:
                            result_list[j][7] += chrom_result[j][7]
                            result_list[j][8] += chrom_result[j][8]

                    if not opts['by_sample']:
                        obs_result.append(obs_mutations)
                    else:
                        obs_result = obs_result + obs_mutations
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            info = (bed_dict[chroms[i]], mut_df, opts)
            chrom_result, obs_mutations = singleprocess_permutation(info)
            for j in range(num_permutations):
                result_list[j][0] += chrom_result[j][0]
                result_list[j][1] += chrom_result[j][1]
                result_list[j][2] += chrom_result[j][2]
                result_list[j][3] += chrom_result[j][3]
                result_list[j][4] += chrom_result[j][4]
                result_list[j][5] += chrom_result[j][5]
                result_list[j][6] += chrom_result[j][6]
                if opts['score_dir']:
                    result_list[j][7] += chrom_result[j][7]
                    result_list[j][8] += chrom_result[j][8]
            if not opts['by_sample']:
                obs_result.append(obs_mutations)
            else:
                obs_result = obs_result + obs_mutations

    return result_list, obs_result


@utils.log_error_decorator
def singleprocess_permutation(info):
    bed_list, mut_df, opts = info
    current_chrom = bed_list[0].chrom
    logger.info('Working on chromosome: {0} . . .'.format(current_chrom))
    num_permutations = opts['num_permutations']
    gene_fa = pysam.Fastafile(opts['input'])
    gs = GeneSequence(gene_fa, nuc_context=opts['context'])

    # variables for recording the actual observed number of non-silent
    # vs. silent mutations
    if not opts['by_sample']:
        obs_silent = 0
        obs_non_silent = 0
        obs_nonsense = 0
        obs_loststop = 0
        obs_splice_site = 0
        obs_loststart = 0
        obs_missense = 0
        obs_vest = 0
        obs_mga_entropy = 0
    else:
        uniq_samp = mut_df['Tumor_Sample'].unique()
        obs_df = pd.DataFrame(np.zeros((len(uniq_samp), len(cols))),
                              index=uniq_samp, columns=cols)

    # go through each gene to permform simulation
    if opts['score_dir']:
        result = [[0, 0, 0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]
    else:
        result = [[0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]
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
            # update the observed count
            if not opts['by_sample']:
                # calc deleterious mutation info
                #tmp_non_silent = cutils.calc_non_silent_info(tmp_mut_info['Reference AA'],
                                                             #tmp_mut_info['Somatic AA'],
                                                             #tmp_mut_info['Codon Pos'])
                # calc mutation info summarizing observed mutations
                tmp_result = cutils.calc_summary_info(tmp_mut_info['Reference AA'],
                                                      tmp_mut_info['Somatic AA'],
                                                      tmp_mut_info['Codon Pos'],
                                                      bed.gene_name,
                                                      opts['score_dir'],
                                                      #min_frac=opts['fraction'],
                                                      min_frac=0.0,
                                                      #min_recur=opts['recurrent']
                                                      min_recur=3
                                                      )
                obs_non_silent += tmp_result[0]
                obs_silent += tmp_result[1]
                obs_nonsense += tmp_result[2]
                obs_loststop += tmp_result[3]
                obs_splice_site += tmp_result[4]
                obs_loststart += tmp_result[5]
                obs_missense += tmp_result[6]
                if opts['score_dir']:
                    obs_vest += tmp_result[-2]
                    obs_mga_entropy += tmp_result[-3]
            else:
                for tsamp in mutations_df['Tumor_Sample'].unique():
                    ixs = np.where(mutations_df['Tumor_Sample']==tsamp)[0]
                    ref_aa = [r for i, r in enumerate(tmp_mut_info['Reference AA']) if i in ixs]
                    somatic_aa = [s for i, s in enumerate(tmp_mut_info['Somatic AA']) if i in ixs]
                    codon_pos = [c for i, c in enumerate(tmp_mut_info['Codon Pos']) if i in ixs]
                    #tmp_non_silent = cutils.calc_non_silent_info(ref_aa,
                                                                 #somatic_aa,
                                                                 #codon_pos)
                    # get summary info
                    tmp_result = cutils.calc_summary_info(ref_aa,
                                                          somatic_aa,
                                                          codon_pos,
                                                          bed.gene_name,
                                                          opts['score_dir'],
                                                          min_frac=0.0,
                                                          min_recur=3)
                    if opts['score_dir']:
                        tmp_result.pop(-4)
                        tmp_result.pop(-4)
                        tmp_result.pop(-1)
                    # update df
                    #obs_df.loc[tsamp,:] = obs_df.loc[tsamp,:] + np.array(tmp_non_silent)
                    obs_df.loc[tsamp,:] = obs_df.loc[tsamp,:] + np.array(tmp_result)

            ## Do permutations
            # calculate non silent count
            #tmp_result = pm.non_silent_ratio_permutation(context_cts,
                                                         #context_to_mutations,
                                                         #sc,  # sequence context obj
                                                         #gs,  # gene sequence obj
                                                         #num_permutations)
            tmp_result = pm.summary_permutation(context_cts,
                                                context_to_mutations,
                                                sc,  # sequence context obj
                                                gs,  # gene sequence obj
                                                opts['score_dir'],
                                                num_permutations)
        else:
            if opts['score_dir']:
                tmp_result = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]
            else:
                tmp_result = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0] for k in range(num_permutations)]

        # increment the non-silent/silent counts for each permutation
        offset = 3
        for j in range(num_permutations):
            result[j][0] += tmp_result[j][0+offset]
            result[j][1] += tmp_result[j][1+offset]
            result[j][2] += tmp_result[j][2+offset]
            result[j][3] += tmp_result[j][3+offset]
            result[j][4] += tmp_result[j][4+offset]
            result[j][5] += tmp_result[j][5+offset]
            result[j][6] += tmp_result[j][6+offset]
            if opts['score_dir']:
                result[j][7] += tmp_result[j][9+offset]
                result[j][8] += tmp_result[j][10+offset]

    gene_fa.close()
    if not opts['by_sample']:
        obs_result = [obs_non_silent, obs_silent, obs_nonsense,
                      obs_loststop, obs_splice_site, obs_loststart, obs_missense]
        if opts['score_dir']:
            obs_result.extend([obs_mga_entropy, obs_vest])
    else:
        obs_result = obs_df
    logger.info('Finished working on chromosome: {0}.'.format(current_chrom))
    return result, obs_result


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
    help_str = 'Directory containing score information in pickle files (Default: None).'
    parser.add_argument('-s', '--score-dir',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Report counts for observed mutations stratified by the tumor sample'
    parser.add_argument('-bs', '--by-sample',
                        action='store_true',
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
    help_str = 'Output text file of observed results (optional).'
    parser.add_argument('-oo', '--observed-output',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Output text file of simulation results'
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
    global cols
    if opts['score_dir']:
        cols.extend(['Total MGAEntropy', 'Total Missense VEST'])

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

    # read in bed info
    bed_dict = utils.read_bed(opts['bed'])

    # perform permutation test
    #permutation_result = multiprocess_permutation(bed_dict, mut_df, opts)
    sim_result, obs_result = multiprocess_permutation(bed_dict, mut_df, opts)

    # report number of observed non-silent and silent mutations
    #obs_result = [x[1] for x in permutation_result]  # actually observed num mutations
    #obs_result = permutation_result[1]  # actually observed num mutations
    if not opts['by_sample']:
        total_non_silent = sum(o[0] for o in obs_result)
        total_silent = sum(o[1] for o in obs_result)
        total_nonsense = sum(o[2] for o in obs_result)
        total_loststop = sum(o[3] for o in obs_result)
        total_splice_site = sum(o[4] for o in obs_result)
        total_loststart = sum(o[5] for o in obs_result)
        total_missense = sum(o[6] for o in obs_result)
        if opts['score_dir']:
            total_mgaentropy = sum(o[7] for o in obs_result)
            total_vest = sum(o[8] for o in obs_result)
        logger.info('There were {0} non-silent SNVs and {1} silent SNVs actually '
                    'observed from the provided mutations.'.format(total_non_silent,
                                                                    total_silent))
        logger.info('There were {0} missense SNVs, {1} nonsense SNVs, {2} lost stop SNVs, '
                    ', {3} lost start, and {4} splice site SNVs'.format(total_missense,
                                                                        total_nonsense,
                                                                        total_loststop,
                                                                        total_loststart,
                                                                        total_splice_site))
    else:
        obs_non_silent_df = obs_result

    #sim_result = [s[0] for s in permutation_result]  # results with permutation
    #sim_result = permutation_result[0]

    # convert to dataframe to save to file
    non_silent_ratio_df = pd.DataFrame(sim_result,
                                       columns=cols)
    # save simulation output
    non_silent_ratio_df.to_csv(opts['output'], sep='\t', index=False)

    # save observed values if file provided
    if opts['observed_output']:
        if not opts['by_sample']:
            obs_result = [total_non_silent, total_silent, total_nonsense,
                        total_loststop, total_splice_site, total_loststart,
                        total_missense]
            if opts['score_dir']:
                obs_result.extend([total_mgaentropy, total_vest])
            obs_non_silent_df = pd.DataFrame([obs_result], columns=cols)
            obs_non_silent_df.to_csv(opts['observed_output'], sep='\t', index=False)
        else:
            obs_non_silent_df.to_csv(opts['observed_output'], sep='\t')

    return non_silent_ratio_df


def cli_main():
    opts = parse_arguments()
    main(opts)


if __name__ == "__main__":
    cli_main()
