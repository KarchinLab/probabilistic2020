#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

# package imports
import permutation2020.python.utils as utils
import permutation_test as pt
import frameshift_binomial_test as fs
import count_frameshifts as cf

import argparse
import pandas as pd
import numpy as np
import logging
import datetime

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


def parse_arguments():
    # make a parser
    info = 'Performs a statistical test for oncogene, TSG, or driver gene'
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
                'increases with more permutations, however this will also '
                'increase the run time (Default: 10000).')
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
    help_str = 'Frameshift counts from count_frameshifts.py'
    parser.add_argument('-fc', '--frameshift-counts',
                        type=str,
                        default=None,
                        help=help_str)
    help_str = ('Background non-coding rate of INDELs with lengths matching '
                'frameshifts in --frameshift-counts option. Enter path to file '
                'generated by calc_non_coding_frameshift_rate.py.')
    parser.add_argument('-non-coding', '--non-coding-background',
                        type=str,
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
                'at least a user specified number of deleterious mutations (default: 2)')
    parser.add_argument('-d', '--deleterious',
                        type=int, default=2,
                        help=help_str)
    help_str = ('Maximum TSG score to allow gene to be tested for oncogene '
                'permutation test. (Default: .1)')
    parser.add_argument('-t', '--tsg-score',
                        type=float, default=.1,
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
    help_str = ('Number of bins to categorize framshift lengths by. Only needed'
                ' if path to frameshift counts (--frameshift-counts) is not '
                'specified for predicting tumor suppressor genes. (Default: 3)')
    parser.add_argument('-bins', '--bins',
                        type=int, default=3,
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


def main(opts,
         mutation_df=None,
         frameshift_df=None):
    # get output file
    myoutput_path = opts['output']
    opts['output'] = ''

    # perform permutation test
    result_df = pt.main(opts, mutation_df)

    # clean up p-values for combined p-value calculation
    if opts['kind'] == 'tsg':
        p_val_col = 'deleterious p-value'
        q_val_col = 'deleterious BH q-value'
    elif opts['kind'] == 'effect':
        p_val_col = 'entropy-on-effect p-value'
        q_val_col = 'entropy-on-effect BH q-value'
    elif opts['kind'] == 'oncogene':
        p_val_col = 'entropy p-value'
        q_val_col = 'entropy BH q-value'
    #smallest_p_val = 1. / (opts['num_permutations'] * 10.)
    smallest_p_val = 1. / (opts['num_permutations']*1.001)
    result_df[p_val_col] = result_df[p_val_col].where(result_df[p_val_col]!=0,
                                                      smallest_p_val)
    result_df[p_val_col] = result_df[p_val_col].fillna(1)
    result_df[q_val_col] = result_df[q_val_col].fillna(1)

    if opts['kind'] != 'oncogene':
        # count frameshifts for user if they did not specify a file
        # that was output from count_frameshifts.py script.
        if frameshift_df is None and opts['frameshift_counts'] is None:
            if mutation_df is None:
                mutation_df = pd.read_csv(opts['mutations'], sep='\t')
            samp_num = len(mutation_df['Tumor_Sample'].unique())
            frameshift_df = cf.count_frameshifts(mutation_df, opts['bed'], opts['bins'],
                                                 samp_num, opts['use_unmapped'])

        # perform frameshift test if not trying to identify oncogenes
        frameshift_result = fs.main(opts, fs_cts=frameshift_df)
        result_df = pd.merge(result_df, frameshift_result, how='outer',
                             left_on='gene', right_on='gene')

        # drop genes that never occur
        if opts['kind'] == 'tsg' or opts['kind'] == 'effect':
            no_ssvs = (result_df['Total Mutations']==0) & (result_df['total frameshifts']==0)
            result_df = result_df[~no_ssvs]

        # calculate combined results
        result_df['combined p-value'] = result_df[[p_val_col, 'frameshift p-value']].apply(utils.fishers_method, axis=1)
        result_df['combined BH q-value'] = utils.bh_fdr(result_df['combined p-value'])
        result_df = result_df.sort(columns='combined p-value')
    elif opts['kind'] == 'oncogene':
        result_df = result_df[result_df['Total Mutations']>0]
        result_df['entropy BH q-value'] = utils.bh_fdr(result_df['entropy p-value'])
        result_df['delta entropy BH q-value'] = utils.bh_fdr(result_df['delta entropy p-value'])
        result_df['recurrent BH q-value'] = utils.bh_fdr(result_df['recurrent p-value'])

    if myoutput_path:
        # write output if specified
        result_df.to_csv(myoutput_path, sep='\t', index=False)

    result_df = result_df.set_index('gene', drop=False)

    return result_df


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
