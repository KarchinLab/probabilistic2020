#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))
sys.path.append(os.path.join(file_dir, '../../'))

# package imports
import prob2020
import prob2020.python.utils as utils
import prob2020.python.p_value as mypval
import prob2020.python.indel as indel
import prob2020.console.randomization_test as rt

import argparse
import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)  # module logger

def parse_arguments():
    # make a parser
    info = 'Performs a statistical test for oncogene, TSG, or driver gene'
    parent_parser = argparse.ArgumentParser(description=info)

    # logging arguments
    parent_parser.add_argument('-ll', '--log-level',
                               type=str,
                               action='store',
                               default='',
                               help='Write a log file (--log-level=DEBUG for debug mode, '
                               '--log-level=INFO for info mode)')
    parent_parser.add_argument('-l', '--log',
                               type=str,
                               action='store',
                               default='stdout',
                               help='Path to log file. (accepts "stdout")')
    parent_parser.add_argument('-v', '--verbose',
                               action='store_true',
                               default=False,
                               help='Flag for more verbose log output')

    # add subparsers
    subparsers = parent_parser.add_subparsers(title='Driver Gene Type', dest='kind')
    parser_og = subparsers.add_parser('oncogene',
                                      help='Find statistically significant oncogene-like genes.',
                                      description='Find statsitically significant oncogene-like genes. '
                                      'Evaluates clustering of missense mutations and high in '
                                      'silico pathogenicity scores for missense mutations.')
    help_info = 'Find statistically significant Tumor Suppressor-like genes.'
    parser_tsg = subparsers.add_parser('tsg',
                                       help=help_info,
                                       description=help_info + ' Evaluates for a higher proportion '
                                       'of inactivating mutations than expected.')
    help_info = 'Find codons with significant clustering of missense mutations in sequence.'
    parser_hotmaps = subparsers.add_parser('hotmaps1d',
                                           help=help_info,
                                           description=help_info + ' Evaluates for a higher ammount of '
                                           'clustering of missense mutations.')
    #parser_protein = subparsers.add_parser('protein', help='Find statistically significant '
                                           #'3D clustering in genes based on protein structure.')

    # program arguments
    for i, parser in enumerate([parser_og, parser_tsg, parser_hotmaps]):
        # group of parameters
        major_parser = parser.add_argument_group(title='Major options')
        advance_parser = parser.add_argument_group(title='Advanced options')

        # set the CLI params
        help_str = 'gene FASTA file from extract_gene_seq.py script'
        major_parser.add_argument('-i', '--input',
                                  type=str, required=True,
                                  help=help_str)
        help_str = ('DNA mutations file (MAF file). Columns can be in any order, '
                    'but should contain the correct column header names.')
        major_parser.add_argument('-m', '--mutations',
                                  type=str, required=True,
                                  help=help_str)
        help_str = 'BED file annotation of genes'
        major_parser.add_argument('-b', '--bed',
                                  type=str, required=True,
                                  help=help_str)
        help_str = ('Number of processes to use for parallelization. 0 indicates using a single '
                    'process without using a multiprocessing pool '
                    '(more means Faster, default: 0).')
        major_parser.add_argument('-p', '--processes',
                                  type=int, default=0,
                                  help=help_str)
        help_str = ('Number of iterations for null model. p-value precision '
                    'increases with more iterations, however this will also '
                    'increase the run time (Default: 100,000).')
        major_parser.add_argument('-n', '--num-iterations',
                                  type=int, default=100000,
                                  help=help_str)
        help_str = ('Number of iterations more significant then the observed statistic '
                    'to stop further computations. This decreases compute time spent in resolving '
                    'p-values for non-significant genes. (Default: 1000).')
        advance_parser.add_argument('-sc', '--stop-criteria',
                                    type=int, default=1000,
                                    help=help_str)
        help_str = ('Number of DNA bases to use as context. 0 indicates no context. '
                    '1 indicates only use the mutated base.  1.5 indicates using '
                    'the base context used in CHASM '
                    '(http://wiki.chasmsoftware.org/index.php/CHASM_Overview). '
                    '2 indicates using the mutated base and the upstream base. '
                    '3 indicates using the mutated base and both the upstream '
                    'and downstream bases. (Default: 1.5)')
        major_parser.add_argument('-c', '--context',
                                  type=float, default=1.5,
                                  help=help_str)
        if i == 0:
            help_str = 'Directory containing VEST score information in pickle files (Default: None).'
            major_parser.add_argument('-s', '--score-dir',
                                      type=str, default=None,
                                      help=help_str)
            help_str = ('Minimum number of mutations at a position for it to be '
                        'considered a recurrently mutated position (Default: 3).')
            advance_parser.add_argument('-r', '--recurrent',
                                        type=int, default=3,
                                        help=help_str)
            help_str = ('Fraction of total mutations in a gene. This define the '
                        'minimumm number of mutations for a position to be defined '
                        'as recurrently mutated (Defaul: .02).')
            advance_parser.add_argument('-f', '--fraction',
                                        type=float, default=.02,
                                        help=help_str)
        elif i == 1:
            help_str = ('Perform tsg randomization-based test if gene has '
                        'at least a user specified number of deleterious mutations (default: 1)')
            advance_parser.add_argument('-d', '--deleterious',
                                        type=int, default=1,
                                        help=help_str)
        elif i == 2:
            help_str = ('Sequence window size for HotMAPS 1D algorithm '
                        'by number of codons (Default: 3)')
            advance_parser.add_argument('-w', '--window',
                                        type=str, default='3',
                                        help=help_str)
            help_str = ('Flag for reporting index (row number, starts at zero) in associated mutation file')
            advance_parser.add_argument('-r', '--report-index',
                                        action='store_true', default=False,
                                        help=help_str)
            help_str = ('Path to directory to save empirical null distribution')
            advance_parser.add_argument('-nd', '--null-distr-dir',
                                        type=str,
                                        help=help_str)
        elif i == 3:
            help_str = 'Directory containing codon neighbor graph information in pickle files (Default: None).'
            major_parser.add_argument('-ng', '--neighbor-graph-dir',
                                      type=str, required=True,
                                      help=help_str)
            help_str = ('Minimum number of mutations at a position for it to be '
                        'considered a recurrently mutated position (Default: 3).')
            advance_parser.add_argument('-r', '--recurrent',
                                type=int, default=3,
                                help=help_str)
            help_str = ('Fraction of total mutations in a gene. This define the '
                        'minimumm number of mutations for a position to be defined '
                        'as recurrently mutated (Default: .02).')
            advance_parser.add_argument('-f', '--fraction',
                                        type=float, default=.02,
                                        help=help_str)
        help_str = ('Only keep unique mutations for each tumor sample. '
                    'Mutations reported from heterogeneous sources may contain'
                    ' duplicates, e.g. a tumor sample was sequenced twice.')
        advance_parser.add_argument('--unique',
                                    action='store_true',
                                    default=False,
                                    help=help_str)
        help_str = ('Use mutations that are not mapped to the the single reference '
                    'transcript for a gene specified in the bed file indicated by '
                    'the -b option.')
        advance_parser.add_argument('-u', '--use-unmapped',
                                    action='store_true',
                                    default=False,
                                    help=help_str)
        help_str = ('Path to the genome fasta file. Required if --use-unmapped flag '
                    'is used. (Default: None)')
        advance_parser.add_argument('-g', '--genome',
                                    type=str, default='',
                                    help=help_str)
        help_str = ('Specify the seed for the pseudo random number generator. '
                    'By default, the seed is randomly chosen. The seed will '
                    'be used for the monte carlo simulations (Default: 101).')
        advance_parser.add_argument('-seed', '--seed',
                                    type=int, default=101,
                                    help=help_str)
        help_str = 'Output text file of probabilistic 20/20 results'
        major_parser.add_argument('-o', '--output',
                                  type=str, required=True,
                                  help=help_str)
    args = parent_parser.parse_args()

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
    logger.info('Version: {0}'.format(prob2020.__version__))
    logger.info('Command: {0}'.format(' '.join(sys.argv)))
    return opts


def main(opts,
         mutation_df=None,
         frameshift_df=None):
    # get output file
    myoutput_path = opts['output']
    opts['output'] = ''

    # perform randomization-based test
    result_df = rt.main(opts, mutation_df)

    # clean up p-values for combined p-value calculation
    if opts['kind'] == 'tsg':
        p_val_col = 'inactivating p-value'
        q_val_col = 'inactivating BH q-value'
    elif opts['kind'] == 'effect':
        p_val_col = 'entropy-on-effect p-value'
        q_val_col = 'entropy-on-effect BH q-value'
    elif opts['kind'] == 'oncogene':
        p_val_col = 'entropy p-value'
        q_val_col = 'entropy BH q-value'
    elif opts['kind'] == 'protein':
        p_val_col = 'normalized graph-smoothed position entropy p-value'
        q_val_col = 'normalized graph-smoothed position entropy BH q-value'
    elif opts['kind'] == 'hotmaps1d':
        p_val_col = 'p-value'
        q_val_col = 'q-value'
    result_df[p_val_col] = result_df[p_val_col].fillna(1)
    result_df[q_val_col] = result_df[q_val_col].fillna(1)

    if opts['kind'] == 'tsg':
        # drop genes that never occur
        if opts['kind'] == 'tsg' or opts['kind'] == 'effect':
            no_ssvs = (result_df['Total SNV Mutations']==0)
            result_df = result_df[~no_ssvs]

        result_df = result_df.sort_values(by=p_val_col)
    elif opts['kind'] == 'oncogene':
        # get FDR
        result_df = result_df[result_df['Total Mutations']>0]
        result_df['entropy BH q-value'] = mypval.bh_fdr(result_df['entropy p-value'])

        # combine p-values
        result_df['tmp entropy p-value'] = result_df['entropy p-value']
        result_df['tmp vest p-value'] = result_df['vest p-value']
        result_df.loc[result_df['entropy p-value']==0, 'tmp entropy p-value'] = 1. / opts['num_iterations']
        result_df.loc[result_df['vest p-value']==0, 'tmp vest p-value'] = 1. / opts['num_iterations']
        result_df['combined p-value'] = result_df[['tmp entropy p-value', 'tmp vest p-value']].apply(mypval.fishers_method, axis=1)
        result_df['combined BH q-value'] = mypval.bh_fdr(result_df['combined p-value'])
        del result_df['tmp vest p-value']
        del result_df['tmp entropy p-value']

    if myoutput_path:
        # write output if specified
        result_df.to_csv(myoutput_path, sep='\t', index=False)

    result_df = result_df.set_index('gene', drop=False)

    return result_df


def cli_main():
    # run main with CLI options
    opts = parse_arguments()
    main(opts)


if __name__ == "__main__":
    cli_main()
