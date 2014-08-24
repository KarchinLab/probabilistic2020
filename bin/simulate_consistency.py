#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# package imports
import permutation2020.python.utils as utils
import permutation2020.python.simulation_plots as plot_data
from permutation2020.python.random_split import RandomSplit
import permutation2020.python.simulation as sim

import permutation_test as pt
import pandas as pd
from scipy import stats
import pysam
import datetime
import argparse
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


def consistency_comparison(df1, df2, bed_dict, non_tested_genes, opts):
    """Function called by multiprocessing to run predictions.

    """
    permutation_result1 = pt.multiprocess_permutation(bed_dict, df1, opts)
    permutation_result2 = pt.multiprocess_permutation(bed_dict, df2, opts)
    if opts['kind'] == 'oncogene':
        permutation_df1 = pt.handle_oncogene_results(permutation_result1, non_tested_genes)
        permutation_df2 = pt.handle_oncogene_results(permutation_result2, non_tested_genes)

        # calculate jaccard similarity
        recurrent_jaccard = sim.jaccard_index(permutation_df1['recurrent BH q-value'],
                                              permutation_df2['recurrent BH q-value'])
        entropy_jaccard = sim.jaccard_index(permutation_df1['entropy BH q-value'],
                                            permutation_df2['entropy BH q-value'])
        delta_entropy_jaccard = sim.jaccard_index(permutation_df1['delta entropy BH q-value'],
                                                  permutation_df2['delta entropy BH q-value'])

        # rank genes
        recurrent_rank1, recurrent_rank2 = sim.rank_genes(permutation_df1['recurrent p-value'].copy(),
                                                          permutation_df2['recurrent p-value'].copy(),
                                                          permutation_df1['recurrent BH q-value'].copy(),
                                                          permutation_df2['recurrent BH q-value'].copy())
        entropy_rank1, entropy_rank2 = sim.rank_genes(permutation_df1['entropy p-value'].copy(),
                                                      permutation_df2['entropy p-value'].copy(),
                                                      permutation_df1['entropy BH q-value'].copy(),
                                                      permutation_df2['entropy BH q-value'].copy())
        delta_entropy_rank1, delta_entropy_rank2 = sim.rank_genes(permutation_df1['delta entropy p-value'].copy(),
                                                                  permutation_df2['delta entropy p-value'].copy(),
                                                                  permutation_df1['delta entropy BH q-value'].copy(),
                                                                  permutation_df2['delta entropy BH q-value'].copy())

        # calc spearman rank correlation
        sp_rho_recurrent, sp_pval_recurrent = stats.pearsonr(recurrent_rank1,
                                                             recurrent_rank2)
        sp_rho_entropy, sp_pval_entropy = stats.pearsonr(entropy_rank1, entropy_rank2)
        sp_rho_delta_entropy, sp_pval_delta_entropy = stats.pearsonr(delta_entropy_rank1,
                                                                     delta_entropy_rank2)

        # calc kendall tau correlation
        if recurrent_rank1.shape[0]:
            kt_rho_recurrent, kt_pval_recurrent = stats.kendalltau(recurrent_rank1,
                                                                   recurrent_rank2)
        else:
            kt_rho_recurrent, kt_pval_recurrent = 0.0, 1.0
        if entropy_rank1.shape[0]:
            kt_rho_entropy, kt_pval_entropy = stats.kendalltau(entropy_rank1,
                                                               entropy_rank2)
        else:
            kt_rho_entropy, kt_pval_entropy = 0.0, 1.0
        if delta_entropy_rank1.shape[0]:
            kt_rho_delta_entropy, kt_pval_delta_entropy = stats.kendalltau(delta_entropy_rank1,
                                                                           delta_entropy_rank2)
        else:
            kt_rho_delta_entropy, kt_pval_delta_entropy = 0.0, 1.0


        results = pd.DataFrame({'jaccard index': [recurrent_jaccard,
                                                  entropy_jaccard,
                                                  delta_entropy_jaccard],
                                'spearman correlation': [sp_rho_recurrent,
                                                         sp_rho_entropy,
                                                         sp_rho_delta_entropy],
                                'kendall tau correlation': [kt_rho_recurrent,
                                                            kt_rho_entropy,
                                                            kt_rho_delta_entropy],
                                'spearman p-value': [sp_pval_recurrent,
                                                     sp_pval_entropy,
                                                     sp_pval_delta_entropy],
                                'kendall tau p-value': [kt_pval_recurrent,
                                                        kt_pval_entropy,
                                                        kt_pval_delta_entropy]},
                                index=['{0} recurrent'.format(opts['kind']),
                                       '{0} entropy'.format(opts['kind']),
                                       '{0} delta entropy'.format(opts['kind'])])
    else:
        permutation_df1 = pt.handle_tsg_results(permutation_result1)
        permutation_df2 = pt.handle_tsg_results(permutation_result2)

        # calculate jaccard similarity
        deleterious_jaccard = sim.jaccard_index(permutation_df1['deleterious BH q-value'],
                                                permutation_df2['deleterious BH q-value'])

        # rank genes
        deleterious_rank1, deleterious_rank2 = sim.rank_genes(permutation_df1['deleterious p-value'].copy(),
                                                              permutation_df2['deleterious p-value'].copy(),
                                                              permutation_df1['deleterious BH q-value'].copy(),
                                                              permutation_df2['deleterious BH q-value'].copy())

        # calc spearman rank correlation
        sp_rho_deleterious, sp_pval_deleterious = stats.pearsonr(deleterious_rank1,
                                                                 deleterious_rank2)

        # calc kendall tau correlation
        if deleterious_rank1.shape[0]:
            kt_rho_deleterious, kt_pval_deleterious = stats.kendalltau(deleterious_rank1,
                                                                       deleterious_rank2)
        else:
            kt_rho_deleterious, kt_pval_deleterious = 0.0, 1.0

        results = pd.DataFrame({'jaccard index': [deleterious_jaccard],
                                'spearman correlation': [sp_rho_deleterious],
                                'kendall tau correlation': [kt_rho_deleterious],
                                'spearman p-value': [sp_pval_deleterious],
                                'kendall tau p-value': [kt_pval_deleterious]},
                                index=['{0} deleterious'.format(opts['kind'])])

    return results


@utils.log_error_decorator
def singleprocess_simulate(info):
    dfg, bed_dict, non_tested_genes, opts = info  # unpack tuple
    df_left, df_right = next(dfg.dataframe_generator())
    single_result = consistency_comparison(df_left, df_right, bed_dict, non_tested_genes, opts)
    return single_result


def parse_arguments():
    # make a parser
    info = 'Performs a permutation test on the oncogene and TSG score'
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
    parser.add_argument('-iter', '--iterations',
                        type=int,
                        action='store',
                        default=10,
                        help='Number of iterations for each sample rate in the simulation')
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
    help_str = ('Perform recurrent mutation permutation test if gene has '
                'atleast a user specified number of recurrent mutations (default: 2)')
    parser.add_argument('-r', '--recurrent',
                        type=int, default=2,
                        help=help_str)
    help_str = ('Perform tsg permutation test if gene has '
                'at least a user specified number of deleterious mutations (default: 5)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=5,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. (Default: .10)')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=.10,
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

    # specify genes to skip
    if opts['kind'] == 'oncogene':
        # find genes with tsg score above threshold to filter out for oncogene
        # permutation test
        non_tested_genes = utils._get_high_tsg_score(mut_df, opts['tsg_score'])
    else:
        # don't filter out genes for tsg permutation test
        non_tested_genes = []

    # select valid single nucleotide variants only
    mut_df = utils._fix_mutation_df(mut_df)

    # read in bed file
    bed_dict = utils.read_bed(opts['bed'], non_tested_genes)

    # object that generates features from randomly choosen sample names while
    # still respecting the stratification of tumor types
    sample_rate = .5  # always do half splits
    dfg = RandomSplit(mut_df.copy(),
                      sub_sample=sample_rate,
                      num_iter=opts['iterations'])

    # perform simulation
    multiproces_output = sim.multiprocess_simulate(dfg, bed_dict, non_tested_genes,
                                                   opts.copy(), singleprocess_simulate)
    sim_results = {i: mydf for i, mydf in enumerate(multiproces_output)}

    # record result for a specific sample rate
    tmp_results = sim.calculate_stats(sim_results)
    result = {}
    result[sample_rate] = tmp_results
    result[sample_rate].to_csv(opts['output'], sep='\t')

    # make pandas panel objects out of summarized
    # results from simulations
    panel_result = pd.Panel(result)

    # sim.save_simulation_result(panel_result, opts['output'])

    # aggregate results for plotting
    plot_results = {'Permutation Test': panel_result}

    return panel_result


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
