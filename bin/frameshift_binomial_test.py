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
import IPython

# logging
import logging
logger = logging.getLogger(__name__)  # module logger


def frameshift_test(fs, bases_at_risk, noncoding_bg):
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
    # use specific background type
    if noncoding_bg:
        logger.info('Reading non-coding frameshift background rate . . .')
        bg = noncoding_bg.ix['non-coding frameshift', :]
    else:
        coding_fs_cts = fs.sum()
        bg = coding_fs_cts.astype(float) / bases_at_risk.sum()

    # iterate through each gene to calculate p-value
    # initialize p-values
    p_values = pd.Series(np.zeros(len(fs)),
                         index=fs.index)
    logger.info('Calculating binomial test p-values . . .')
    for k in range(len(fs)):
        g_obs = fs.iloc[k,:].sum()
        if g_obs == 0:
            p_values[k] = 1.0
            continue

        # get weighting factor
        w = fs.iloc[k,:].astype(float) / g_obs
        Pg = 0
        for i in range(len(fs.columns)):
            Pg += w[i] * bg[i]

        p_val = binomial_test(g_obs, bases_at_risk[k], Pg)
        p_values[k] = p_val
    logger.info('Finished calculating binomial test p-values.')

    # format results
    qval = utils.bh_fdr(p_values)

    cols = ['gene', 'frameshift p-value', 'frameshift BH q-value']
    fs_result = pd.DataFrame({cols[0]: fs.index,
                              cols[1]: p_values,
                              cols[2]: qval},
                              index=fs.index)[cols]  # make sure cols in right order

    fs_result = fs_result.sort(columns='frameshift p-value', ascending=True)
    return fs_result


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


def read_noncoding_background_rate(path):
    """Reads in the non-coding background rate for frameshift text file.

    Needed for binomial test.
    """
    background_df = pd.read_csv(path, sep='\t', index_col=0)
    drop_list = ['Genome Length', 'Black List Length',
                 'Non-coding Length', 'Number of Samples']
    background_df = background_df.drop(drop_list, axis=1)
    bases_at_risk = background_df.pop('Bases at Risk').iloc[0]
    background_df.ix['non-coding frameshift'] = background_df.ix['non-coding frameshift'].astype(float) / bases_at_risk
    return background_df


def parse_arguments():
    info = ('Peform a Binomial test on the frequency of frameshift mutations '
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
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')

    # script arguments
    help_str = 'Frameshift counts from count_frameshifts.py'
    parser.add_argument('-fc', '--frameshift-counts',
                        type=str, required=True,
                        help=help_str)
    help_str = ('Background non-coding rate of INDELs with lengths matching '
                'frameshifts in --frameshift-counts option. Enter path to file '
                'generated by calc_non_coding_frameshift_rate.py.')
    parser.add_argument('-non-coding', '--non-coding-background',
                        type=str,
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
    utils.start_logging(log_file=log_file,
                        log_level=log_level,
                        verbose=args.verbose)  # start logging
    logger.info('Command: {0}'.format(' '.join(sys.argv)))

    return vars(args)


def main(opts,
         fs_cts=None):
    # read in gene frameshift counts
    if fs_cts is None:
        fs_cts = pd.read_csv(opts['frameshift_counts'], sep='\t', index_col=0)
    total_fs = fs_cts.pop('total')
    fs_cts = fs_cts.drop(['unmapped', 'gene length'], axis=1)
    gene_bases_at_risk = fs_cts.pop('bases at risk')

    # read non-coding background rate file
    if opts['non_coding_background']:
        f = read_noncoding_background_rate(opts['non_coding_background'])
    else:
        f = None

    # perform binomial test
    result = frameshift_test(fs_cts, gene_bases_at_risk, f)
    result['total frameshifts'] = total_fs

    # save results
    if opts['output']:
        result.to_csv(opts['output'], sep='\t', index=False)
    logger.info('Finished!')

    return result


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
