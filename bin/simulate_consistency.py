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
import probabilistic2020 as prob
import count_frameshifts as cf
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


def consistency_comparison(df1, df2, opts):
    """Function called by multiprocessing to run predictions.

    """
    # count frameshifts
    if opts['kind'] != 'oncogene':
        fs1 = cf.count_frameshifts(df1, opts['bed'], opts['bins'],
                                   opts['sample_number'], opts['use_unmapped'])
        fs2 = cf.count_frameshifts(df2, opts['bed'], opts['bins'],
                                   opts['sample_number'], opts['use_unmapped'])
    else:
        fs1 = None
        fs2 = None

    # get prediction results
    permutation_df1 = prob.main(opts, mutation_df=df1, frameshift_df=fs1)
    permutation_df2 = prob.main(opts, mutation_df=df2, frameshift_df=fs2)
    if opts['kind'] == 'oncogene':
        ## calculate jaccard similarity
        # calculate recurrent info
        sort_cols = ['recurrent p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        recurrent_jaccard = sim.overlap(permutation_df1['recurrent BH q-value'],
                                              permutation_df2['recurrent BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['recurrent BH q-value'],
                                                   permutation_df2['recurrent BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        recurrent_ji, recurrent_mean_ji, recurrent_weighted_mean_ji = ji_mean_tuple
        # calculate entropy consistency
        sort_cols = ['entropy p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        entropy_jaccard = sim.overlap(permutation_df1['entropy BH q-value'],
                                            permutation_df2['entropy BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['entropy BH q-value'],
                                                   permutation_df2['entropy BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        entropy_ji, entropy_mean_ji, entropy_weighted_mean_ji = ji_mean_tuple
        # calculate delta entropy consistency
        sort_cols = ['delta entropy p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        delta_entropy_jaccard = sim.overlap(permutation_df1['delta entropy BH q-value'],
                                                  permutation_df2['delta entropy BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['delta entropy BH q-value'],
                                                   permutation_df2['delta entropy BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        delta_entropy_ji, delta_entropy_mean_ji, delta_entropy_weighted_mean_ji = ji_mean_tuple

        results = pd.DataFrame({'jaccard index': [recurrent_jaccard,
                                                  entropy_jaccard,
                                                  delta_entropy_jaccard],
                                'mean jaccard index': [recurrent_mean_ji,
                                                       entropy_mean_ji,
                                                       delta_entropy_mean_ji],
                                'weighted mean jaccard index': [recurrent_weighted_mean_ji,
                                                                entropy_weighted_mean_ji,
                                                                delta_entropy_weighted_mean_ji],

                                },
                                index=['{0} recurrent'.format(opts['kind']),
                                       '{0} entropy'.format(opts['kind']),
                                       '{0} delta entropy'.format(opts['kind'])])

        # record jaccard index at intervals
        ji_intervals = range(opts['step_size'], opts['depth']+1, opts['step_size'])
        ji_df = pd.DataFrame({'recurrent ji': recurrent_ji,
                              'entropy ji': entropy_ji,
                              'delta entropy ji': delta_entropy_ji},
                             index=ji_intervals)
    elif opts['kind']=='tsg':
        # calculate jaccard similarity
        sort_cols = ['deleterious p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        deleterious_jaccard = sim.overlap(permutation_df1['deleterious BH q-value'],
                                                permutation_df2['deleterious BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['deleterious BH q-value'],
                                                   permutation_df2['deleterious BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        deleterious_ji, deleterious_mean_ji, deleterious_weighted_mean_ji = ji_mean_tuple

        # frameshift p-values
        sort_cols = ['frameshift p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        fs_jaccard = sim.overlap(permutation_df1['frameshift BH q-value'],
                                                permutation_df2['frameshift BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['frameshift BH q-value'],
                                                   permutation_df2['frameshift BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        fs_ji, fs_mean_ji, fs_weighted_mean_ji = ji_mean_tuple

        # combined p-value
        sort_cols = ['combined p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        comb_jaccard = sim.overlap(permutation_df1['combined BH q-value'],
                                                permutation_df2['combined BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['combined BH q-value'],
                                                   permutation_df2['combined BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        comb_ji, comb_mean_ji, comb_weighted_mean_ji = ji_mean_tuple

        results = pd.DataFrame({'jaccard index': [deleterious_jaccard,
                                                  fs_jaccard,
                                                  comb_jaccard],
                                'mean jaccard index': [deleterious_mean_ji,
                                                       fs_mean_ji,
                                                       comb_mean_ji],
                                'weighted mean jaccard index': [deleterious_weighted_mean_ji,
                                                                fs_weighted_mean_ji,
                                                                comb_weighted_mean_ji],
                                },
                                index=['{0} deleterious'.format(opts['kind']),
                                       'frameshift', 'combined'])

        # record jaccard index at intervals
        ji_intervals = range(opts['step_size'], opts['depth']+1, opts['step_size'])
        ji_df = pd.DataFrame({'deleterious ji': deleterious_ji,
                              'frameshift ji': fs_ji,
                              'combined ji': comb_ji},
                             index=ji_intervals)
    else:
        # calculate jaccard similarity
        sort_cols = ['entropy-on-effect p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        on_effect_jaccard = sim.overlap(permutation_df1['entropy-on-effect BH q-value'],
                                              permutation_df2['entropy-on-effect BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['entropy-on-effect BH q-value'],
                                                   permutation_df2['entropy-on-effect BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        on_effect_ji, on_effect_mean_ji, on_effect_weighted_mean_ji = ji_mean_tuple

        # calc frameshift
        sort_cols = ['frameshift p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        fs_jaccard = sim.overlap(permutation_df1['frameshift BH q-value'],
                                              permutation_df2['frameshift BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['frameshift BH q-value'],
                                                   permutation_df2['frameshift BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        fs_ji, fs_mean_ji, fs_weighted_mean_ji = ji_mean_tuple

        # calc combined
        sort_cols = ['combined p-value', 'gene']
        permutation_df1 = permutation_df1.sort(columns=sort_cols)
        permutation_df2 = permutation_df2.sort(columns=sort_cols)
        comb_jaccard = sim.overlap(permutation_df1['combined BH q-value'],
                                              permutation_df2['combined BH q-value'])
        ji_mean_tuple = sim.weighted_overlap(permutation_df1['combined BH q-value'],
                                                   permutation_df2['combined BH q-value'],
                                                   max_depth=opts['depth'],
                                                   step_size=opts['step_size'],
                                                   weight_factor=opts['weight'])
        comb_ji, comb_mean_ji, comb_weighted_mean_ji = ji_mean_tuple

        results = pd.DataFrame({'jaccard index': [on_effect_jaccard,
                                                  fs_jaccard,
                                                  comb_jaccard],
                                'mean jaccard index': [on_effect_mean_ji,
                                                       fs_mean_ji,
                                                       comb_mean_ji],
                                'weighted mean jaccard index': [on_effect_weighted_mean_ji,
                                                                fs_weighted_mean_ji,
                                                                comb_weighted_mean_ji],
                                },
                                index=['entropy-on-effect', 'frameshift', 'combined'])

        # record jaccard index at intervals
        ji_intervals = range(opts['step_size'], opts['depth']+1, opts['step_size'])
        ji_df = pd.DataFrame({'entropy-on-effect ji': on_effect_ji,
                              'frameshift ji': fs_ji,
                              'combined ji': comb_ji},
                             index=ji_intervals)
    return results, ji_df


@utils.log_error_decorator
def singleprocess_simulate(info):
    dfg, opts = info  # unpack tuple
    df_left, df_right = next(dfg.dataframe_generator())
    single_result = consistency_comparison(df_left, df_right, opts)
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
    help_str = 'Step size for progessing further down list of top genes'
    parser.add_argument('-step', '--step-size',
                        type=int,
                        default=40,
                        help=help_str)
    help_str = 'Weight factor for ranked biased consistency measure'
    parser.add_argument('-weight', '--weight',
                        type=float,
                        default=.1,
                        help=help_str)
    help_str = 'Maximum depth of genes from top of list to consider for consistency'
    parser.add_argument('-depth', '--depth',
                        type=int,
                        default=200,
                        help=help_str)
    help_str = ('Output the jaccard index at intervals specified by the "-depth" '
                'and the "-step" options. Specify the file path to save output.')
    parser.add_argument('-jc', '--jaccard-curve',
                        type=str, required=True,
                        help=help_str)
    parser.add_argument('-iter', '--iterations',
                        type=int,
                        action='store',
                        default=10,
                        help='Number of iterations for each sample rate in the simulation')
    help_str = ('Sample rate (i.e. fraction of total samples) for sampling with '
                'replacement. If 0, then sampling without replacement is performed '
                'at a sample rate of .5 (i.e. split data in half). (Default: 0)')
    parser.add_argument('-w', '--with-replacement',
                        type=float, default=0.0,
                        help=help_str)
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
    help_str = ('Background non-coding rate of INDELs with lengths matching '
                'frameshifts in --frameshift-counts option. Enter path to file '
                'generated by calc_non_coding_frameshift_rate.py.')
    parser.add_argument('-non-coding', '--non-coding-background',
                        type=str,
                        help=help_str)
    help_str = 'Number of bins to categorize framshift lengths by'
    parser.add_argument('-bins', '--bins',
                        type=int, default=5,
                        help=help_str)
    help_str = 'Number of sequenced samples'
    parser.add_argument('-sn', '--sample-number',
                        type=int,
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
    help_str = ('Perform tsg permutation test if gene has '
                'at least a user specified number of deleterious mutations (default: 2)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=2,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. (Default: .10)')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=.10,
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
    # read in mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')

    # pre-filter out high TSG score genes if using oncogene method
    if opts['kind'] == 'oncogene':
        non_tested_genes = utils._get_high_tsg_score(mut_df, opts['tsg_score'])
        mut_df = mut_df[~mut_df['Gene'].isin(non_tested_genes)]
        opts['tsg_score'] = 1.0

    # object that generates features from randomly choosen sample names while
    # still respecting the stratification of tumor types
    if not opts['with_replacement']:
        sample_rate = .5  # always do half splits
        dfg = RandomSplit(mut_df.copy(),
                          sub_sample=sample_rate,
                          num_iter=opts['iterations'])
    else:
        sample_rate = opts['with_replacement']
        dfg = RandomSplit(mut_df.copy(),
                          sub_sample=sample_rate,
                          num_iter=opts['iterations'],
                          with_replacement=True)

    # correct sample numbers for down sampling
    if opts['sample_number']:
        opts['sample_number'] = int(sample_rate * opts['sample_number'])

    # perform simulation
    multiproces_output = sim.multiprocess_simulate(dfg, opts.copy(),
                                                   singleprocess_simulate)
    sim_results = {i: output[0] for i, output in enumerate(multiproces_output)}
    sim_ji_results = {i: output[1] for i, output in enumerate(multiproces_output)}

    # record result for a specific sample rate
    tmp_results = sim.calculate_stats(sim_results)
    result = {}
    result[sample_rate] = tmp_results
    result[sample_rate].to_csv(opts['output'], sep='\t')

    # ouptput jaccard index curve information
    tmp_results = sim.calculate_stats(sim_ji_results)
    tmp_results.to_csv(opts['jaccard_curve'], sep='\t')

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
