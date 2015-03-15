#!/usr/bin/env python
# fix problems with pythons terrible import system
import sys
import os
file_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(file_dir, '../'))

import pandas as pd
import numpy as np
from scipy.stats import binom
import permutation2020.python.utils as utils
import argparse
import IPython

# logging
import logging
logger = logging.getLogger(__name__)  # module logger


def frameshift_test(fs, bases_at_risk, noncoding_bg,
                    overdisperion=False):
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
    overdisperion : bool
        flag indicating if should use overdispersion model

    Returns
    -------
    fs_result : pd.DataFrame
        p-value/q-value for each gene for binomial frameshift test
    """
    # use specific background type
    if noncoding_bg:
        logger.info('Reading non-coding frameshift background rate . . .')
        bg = noncoding_bg.ix['non-coding frameshift', :]
        expected_rate = bg.sum()
        #cv = np.mean(noncoding_bg.ix['cv']) # np.mean(noncoding_bg.ix['bootstrap std']/noncoding_bg.ix['non-coding frameshift'])
        cv = noncoding_bg.ix['cv']
        logger.info('The coefficient of variation is {0}'.format(cv.tolist()))
    else:
        coding_fs_cts = fs.sum()
        bg = coding_fs_cts.astype(float) / bases_at_risk.sum()

    if overdisperion:
        import rpy2.robjects as ro
        import rpy2.robjects.numpy2ri
        # Code for R's VGAM package using rpy2
        ro.r("suppressPackageStartupMessages(library(VGAM))")  # load randomForest library
        ro.r('''betabinom.test <- function(n, N, alpha, beta){
            return(1 - pbetabinom.ab(n, N, alpha, beta))
        }''')
        ro.r('''betabinom.lrt <- function(n, N, alpha.null, beta.null, alpha.alt, beta.alt){
            N.vec <- rep(N, by=length(n))
            loglik.alt <- sum(dbetabinom.ab(n, N.vec, alpha.alt, beta.alt, log=T))
            loglik.null <- sum(dbetabinom.ab(n, N.vec, alpha.null, beta.null, log=T))
            lrtStatistic <- -2 * (loglik.alt - loglik.null)
            result <- c(1 - pchisq(lrtStatistic, 1), lrtStatistic)
            return(result)
        }''')
        betabinom_test = ro.r['betabinom.test']
        betabinom_lrt = ro.r['betabinom.lrt']

    # iterate through each gene to calculate p-value
    # initialize p-values
    p_values = pd.Series(np.zeros(len(fs)),
                         index=fs.index)
    lrt_stats = pd.Series(np.zeros(len(fs)),
                          index=fs.index)
    effect_size = pd.Series(np.zeros(len(fs)),
                            index=fs.index)
    logger.info('Calculating binomial test p-values . . .')
    for k in range(len(fs)):
        g_obs = fs.iloc[k,:].sum()
        if g_obs == 0:
            p_values[k] = 1.0
            continue

        # get the relative ratio of observed frameshift rate to expected
        obs_rate = g_obs.astype(float) / bases_at_risk[k]
        relRate = obs_rate / expected_rate

        # get weighting factor
        w = fs.iloc[k,:].astype(float) / g_obs
        Pg = 0
        for i in range(len(fs.columns)):
            Pg += w[i] * bg[i]


        if overdisperion:
            # just redine the background rates as Pg
            Pg = bg

            # calc shape params
            ab_null = Pg * (1-Pg) / (cv*Pg)**2 - 1
            alpha_null = (Pg * ab_null)
            beta_null = ((1-Pg) * ab_null)
            Pg = relRate*Pg
            ab_alt = Pg * (1-Pg) / (cv*Pg)**2 - 1
            alpha_alt = (Pg * ab_alt)
            beta_alt = ((1-Pg) * ab_alt)
            #p_val = betabinom_test(g_obs, bases_at_risk[k], alpha, beta)  # returns a float vector
            obs_fs = ro.IntVector(fs.iloc[k,:])
            alpha_null_vec = ro.FloatVector(alpha_null)
            beta_null_vec = ro.FloatVector(beta_null)
            alpha_alt_vec = ro.FloatVector(alpha_alt)
            beta_alt_vec = ro.FloatVector(beta_alt)
            p_val = betabinom_lrt(obs_fs, bases_at_risk[k],
                                  alpha_null_vec, beta_null_vec,
                                  alpha_alt_vec, beta_alt_vec)  # returns a float vector
            p_values[k] = p_val[0]
            lrt_stats[k] = p_val[1]
            effect_size[k] = relRate
        else:
            p_val = binomial_test(g_obs, bases_at_risk[k], Pg)
            p_values[k] = p_val
    logger.info('Finished calculating binomial test p-values.')

    # format results
    qval = utils.bh_fdr(p_values)

    if not overdisperion:
        cols = ['gene', 'frameshift p-value', 'frameshift BH q-value']
        fs_result = pd.DataFrame({cols[0]: fs.index,
                                  cols[1]: p_values,
                                  cols[2]: qval},
                                  index=fs.index)[cols]  # make sure cols in right order
    else:
        cols = ['gene', 'LRT statistic', 'Effect Size', 'frameshift p-value', 'frameshift BH q-value']
        fs_result = pd.DataFrame({cols[0]: fs.index,
                                  cols[1]: lrt_stats,
                                  cols[2]: effect_size,
                                  cols[3]: p_values,
                                  cols[4]: qval},
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
    help_str = 'Flag indicating using overdispersed model'
    parser.add_argument('--overdispersion',
                        action='store_true',
                        default=False,
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
    result = frameshift_test(fs_cts, gene_bases_at_risk, f, overdisperion=opts['overdispersion'])
    result['total frameshifts'] = total_fs

    # save results
    if opts['output']:
        result.to_csv(opts['output'], sep='\t', index=False)
    logger.info('Finished!')

    return result


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
