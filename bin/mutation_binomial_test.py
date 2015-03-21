#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import pandas as pd
import numpy as np
from scipy.stats import binom
import prob2020.python.utils as utils
import argparse

# logging
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


def frequency_test(mut, bases_at_risk):
    """Perform a binomial test on the frequency of frameshifts while accounting
    for differing rates depending on frameshift length.

    Parameters
    ----------
    fs : pd.DataFrame
        dataframe with genes as index and columns as frameshift lengths.
    bases_at_risk : pd.Series/np.array
        contains the number of bases at risk for a frameshift. It equals the
        number of samples times the gene length.
    noncoding_bg : pd.DataFrame
        If using non-coding background for frameshifts provide a table from
        the read_noncoding_background_rate function. If using a gene based
        background, then pass in the python object None.

    Returns
    -------
    fs_result : pd.DataFrame
        p-value/q-value for each gene for binomial frameshift test
    """
    # initialize p-values
    p_values = pd.Series(np.zeros(len(mut)),
                         index=mut.index)
    non_silent_ct = pd.Series(np.zeros(len(mut)),
                              index=mut.index)

    # get rid of the ' bases at risk' name in the headers
    bases_at_risk = bases_at_risk.rename(columns=lambda x: x.replace(' bases at risk', ''))

    coding_fs_cts = mut.sum()
    bg = coding_fs_cts.astype(float) / bases_at_risk.sum()

    # column names for SNV contexts at risk
    snv_cols = [col for col in bases_at_risk.columns
                if 'indel' not in col]

    # iterate through each gene to calculate p-value
    logger.info('Calculating binomial test p-values . . .')
    for k in range(len(mut)):
        g_obs = mut.iloc[k,:].sum()
        w = mut.iloc[k,:].astype(float) / g_obs

        # get weighting factor
        Pg = 0
        for i in range(len(mut.columns)):
            Pg += w[i] * bg[i]

        num_bases_at_risk = bases_at_risk.ix[k, snv_cols].sum()
        p_val = binomial_test(g_obs, num_bases_at_risk, Pg)
        p_values[k] = p_val
        non_silent_ct[k] = g_obs

    logger.info('Finished calculating binomial test p-values.')

    # format results
    qval = utils.bh_fdr(p_values)
    cols = ['gene', 'non-silent count', 'mutation frequency p-value', 'mutation frequency BH q-value']
    mut_result = pd.DataFrame({cols[0]: mut.index,
                               cols[1]: non_silent_ct,
                               cols[2]: p_values,
                               cols[3]: qval},
                              index=mut.index)[cols]  # make sure cols in right order
    mut_result = mut_result.sort(columns='mutation frequency p-value', ascending=True)
    return mut_result


def binomial_test(n, N, P):
    """Perform binomial test on the observed n being higher than expected.

    Specifically, N bases are at risk and of those there are n that a frameshift
    occurred at. Given the background probability of a frameshift at a specific
    base, the p-value is calculated as the probability of observing n or greater
    frameshifts. Since N is large and n is small, it is computationally more
    efficient to take 1 - Pr(i<=n).

    Parameters
    ----------
    n : int
        number of frameshifts observed (i.e. num bases where frameshifts occurred)
    N : int
        number of bases at risk (num samples * gene length)
    P : float
        background probability that a frameshift would occur at a single base

    Returns
    -------
    pval : np.array
        p-value for binomial test
    """
    if n <= 0:
        return 1.0
    pval = binom.sf(n, N, P)
    return pval


def parse_arguments():
    info = ('Peform a Binomial test on the frequency of non-silent mutations '
            'compared to the background.')
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

    # script arguments
    help_str = 'Mutation counts from count_mutations.py'
    parser.add_argument('-mc', '--mutation-counts',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Output file'
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
    logger.info('Command: {0}'.format(' '.join(sys.argv)))

    return vars(args)


def main(opts):
    # read in gene frameshift counts
    mut_cts = pd.read_csv(opts['mutation_counts'], sep='\t', index_col=0)

    # drop extra info
    # mut_cts = mut_cts.drop(['indel unmapped', 'gene length'], axis=1)

    # separate out bases at risk columns
    bases_at_risk_cols = [col for col in mut_cts.columns if 'bases at risk' in col]
    bases_at_risk = mut_cts[bases_at_risk_cols]
    mut_cts = mut_cts.drop(bases_at_risk_cols, axis=1)

    # perform binomial test
    result = frequency_test(mut_cts, bases_at_risk)

    # save results
    if opts['output']:
        result.to_csv(opts['output'], sep='\t', index=False)
    logger.info('Finished!')

    return result


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
