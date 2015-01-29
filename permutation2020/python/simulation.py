import numpy as np
import pandas as pd
from scipy import stats
import time
import sys
import utils
from multiprocessing import Pool
import itertools as it
import logging

logger = logging.getLogger(__name__)  # module logger

def multiprocess_simulate(dfg, opts,
                          singleprocess_func):
    # handle the number of processes to use, should it even use the
    # multiprocessing module?
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    opts['processes'] = 0  # do not use multi-processing within permutation test

    # handle multiprocessing of simulation if necessary
    process_results = None
    result_list = []
    for i in range(0, dfg.num_iter, num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            del process_results  # possibly help free up more memory
            time.sleep(5)  # wait 5 seconds, might help make sure memory is free
            tmp_num_pred = dfg.num_iter - i if  i + num_processes > dfg.num_iter else num_processes
            # df_generator = dfg.dataframe_generator()
            info_repeat = it.repeat((dfg, opts), tmp_num_pred)
            #pool = Pool(processes=tmp_num_pred)
            process_results = pool.imap(singleprocess_func, info_repeat)
            process_results.next = utils.keyboard_exit_wrapper(process_results.next)
            try:
                for tmp_result in process_results:
                    result_list.append(tmp_result)
            except KeyboardInterrupt:
                pool.close()
                pool.join()
                logger.info('Exited by user. ctrl-c')
                sys.exit(0)
            pool.close()
            pool.join()
        else:
            info = (dfg, opts)
            tmp_result = singleprocess_func(info)
            result_list.append(tmp_result)
    return result_list


def calculate_sem(wp):
    """Calculates the standard error of the mean for a pd.Panel object.

    **Note:** The pd.Panel.apply method seems to have a bug preventing
    me from using it. So instead I am using the numpy apply function
    for calculating sem.

    Parameters
    ----------
    wp : pd.Panel
        panel that stratifies samples

    Returns
    -------
    tmp_sem : pd.DataFrame
        standard error of the mean calculated along the sample axis
    """
    tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)  # hack because pandas apply method has a bug
    tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           columns=wp.minor_axis,
                           index=wp.major_axis)
    return tmp_sem


def calculate_stats(result_dict,
                    metrics=['precision', 'recall', 'ROC AUC',
                             'PR AUC', 'count', 'average num drivers']):
    """Computes mean and sem of classification performance metrics.

    Parameters
    ----------
    result_dict : dict
        dictionary with the i'th sample as the key and data frames
        with "oncogene"/"tsg" (row) classification performance metrics
        (columns) as values

    Returns
    -------
    result_df : pd.DataFrame
        Data frame with mean and sem of classification performance
        metrics. (rows: "oncogene"/"tsg", columns: summarized metrics)
    """
    wp = pd.Panel(result_dict)
    tmp_means = wp.mean(axis=0)
    tmp_sem = calculate_sem(wp)
    result_df = pd.merge(tmp_means, tmp_sem,
                         left_index=True, right_index=True,
                         suffixes=(' mean', ' sem'))
    return result_df


def save_simulation_result(mypanel, mypath):
    # make 'oncogenes'/'tsg' as the 'items' axis
    mypan = mypanel.swapaxes('items', 'major')

    # collapse pd.panel into a data frame
    myitems = mypan.items
    num_items = len(myitems)
    mydf = mypan[myitems[0]]
    rename_dict = {'count mean': 'count mean {0}'.format(myitems[0]),
                   'count sem': 'count sem {0}'.format(myitems[0])}
    mydf = mydf.rename(columns=rename_dict)
    if num_items > 1:
        for i in range(1, num_items):
            tmp_df = mypan[myitems[i]]
            rename_dict = {'count mean': 'count mean {0}'.format(myitems[i]),
                           'count sem': 'count sem {0}'.format(myitems[i])}
            mydf = pd.merge(mydf, tmp_df.rename(columns=rename_dict),
                            left_index=True, right_index=True)

    # save data frame to specified paths
    mydf.to_csv(mypath, sep='\t')


def rank_genes(p1, p2,
               fdr1, fdr2,
               thresh=.1,
               na_fill=1.):
    """Rank the p-values of statistically significant genes in either list.

    The thresh variable is the FDR threshold for statistically significant in either
    list. The union of statistically significant genes in both lists is taken.
    This union of genes is then ranked in both lists by p-value.

    Parameters
    ----------
    p1 : pd.Series
        P-values of first test
    p2 : pd.Series
        P-values of second test
    fdr1 : pd.Series
        q-values of first test
    fdr2 : pd.Series
        q-values of second test
    thresh : float
        FDR threshold for statistical significance
    na_fill : float
        value to fill missing p-values (NA's)

    Returns
    -------
    top_rank1 : pd.Series
        gene ranks for first test
    top_rank2 : pd.Series
        gene ranks for second test
    """
    all_ixs = list(set(fdr1[fdr1<thresh].index) | set(fdr2[fdr2<thresh].index))
    top_s1 = p1[all_ixs]
    top_s2 = p2[all_ixs]
    top_rank1 = top_s1.fillna(na_fill).rank()
    top_rank2 = top_s2.fillna(na_fill).rank()
    top_rank2 = top_rank2[top_rank1.index]  # match ordering of index
    return top_rank1, top_rank2


def overlap(fdr1, fdr2,
            thresh=.1,
            depth=None,
            keep_same=False):
    """Calculates the fraction overlap for statistically significant genes.

    The thresh parameter determines statistical significance for each list.
    It represents the FDR threshold. If the depth parameter is provided, then
    it overides the thresh parameter. Simply, the depth parameter specifies
    how many of the top genes should be compared.

    Parameters
    ----------
    fdr1 : pd.Series
        q-values for first test, index should be gene names
    fdr2 : pd.Series
        q-values for second test, index should be gene names
    thresh : float
        FDR threshold for statistical signficance
    depth : int
        Number of top genes to evaluate the jaccard index for
    keep_same : bools
        Whether an FDR thrshold or depth limit should NOT be
        applied. If true, this ignores the user input for
        parameters thresh and depth.

    Returns
    -------
    overlap_sim : float
        overlap measuring simularity of statistically significant
        genes in both tests
    """
    if not keep_same:
        if not depth:
            s1_genes = set(fdr1[fdr1<thresh].index)
            s2_genes = set(fdr2[fdr2<thresh].index)
        else:
            s1_genes = set(fdr1[:depth].index)
            s2_genes = set(fdr2[:depth].index)
    else:
        s1_genes = set(fdr1.index)
        s2_genes = set(fdr2.index)

    num_intersect = len(s1_genes & s2_genes)
    num_total = len(s1_genes)
    if num_total:
        # there are significant genes
        overlap_sim = num_intersect / float(num_total)
    else:
        # no significant gene case
        overlap_sim = 0
    return overlap_sim


def weighted_overlap(fdr1, fdr2,
                     max_depth,
                     step_size,
                     weight_factor):
    # calculate jaccard index at specified intervals
    num_depths = (max_depth) // step_size
    num_depths_total = len(fdr2) // step_size
    ov = np.zeros(num_depths)
    ov_all = np.zeros(num_depths_total)
    for i, depth in enumerate(range(step_size, num_depths_total+1, step_size)):
        if depth <= max_depth:
            ov_tmp = overlap(fdr1.iloc[:depth].copy(), fdr2.iloc[:depth+100].copy(),
                             keep_same=True, depth=max_depth)
            ov[i] = ov_tmp
            ov_all[i] = ov_tmp
        else:
            # ov_all[i] = overlap(fdr1.iloc[:max_depth].copy(), fdr2, depth=depth)
            ov_all[i] = overlap(fdr1.iloc[:max_depth].copy(), fdr2.iloc[:depth+100].copy())

    # calculate the weighting for jaccard index
    p = weight_factor ** (1./(num_depths-1))
    w = p*np.ones(num_depths_total)
    w[0] = 1
    w = np.cumprod(w)  # calculate geometric weights
    w = w / w.sum()  # normalize weights to 1

    weighted_mean_ov = np.dot(w, ov_all)
    mean_ov = np.mean(ov)

    return ov, mean_ov, weighted_mean_ov


def jaccard_index(fdr1, fdr2,
                  thresh=.1,
                  depth=None):
    """Calculates the Jaccard Index for statistically significant genes.

    The thresh parameter determines statistical significance for each list.
    It represents the FDR threshold. If the depth parameter is provided, then
    it overides the thresh parameter. Simply, the depth parameter specifies
    how many of the top genes should be compared.

    Parameters
    ----------
    fdr1 : pd.Series
        q-values for first test, index should be gene names
    fdr2 : pd.Series
        q-values for second test, index should be gene names
    thresh : float
        FDR threshold for statistical signficance
    depth : int
        Number of top genes to evaluate the jaccard index for

    Returns
    -------
    jaccard_sim : float
        Jaccard index measuring simularity of statistically significant
        genes in both tests
    """
    if not depth:
        s1_genes = set(fdr1[fdr1<thresh].index)
        s2_genes = set(fdr2[fdr2<thresh].index)
    else:
        s1_genes = set(fdr1[:depth].index)
        s2_genes = set(fdr2[:depth].index)

    num_intersect = len(s1_genes & s2_genes)
    num_union = len(s1_genes | s2_genes)
    if num_union:
        # there are significant genes
        jaccard_sim = num_intersect / float(num_union)
    else:
        # no significant gene case
        jaccard_sim = 0
    return jaccard_sim


def weighted_jaccard_index(fdr1, fdr2,
                           max_depth,
                           step_size,
                           weight_factor):
    # calculate jaccard index at specified intervals
    num_depths = (max_depth) // step_size
    num_depths_total = (len(fdr2)) // step_size
    ji = np.zeros(num_depths)
    ji_all = np.zeros(num_depths_total)
    for i, depth in enumerate(range(step_size, num_depths_total+1, step_size)):
        if depth <= max_depth:
            ji_tmp = jaccard_index(fdr1.iloc[:depth].copy(), fdr2, depth=max_depth)
            ji[i] = ji_tmp
            ji_all[i] = ji_tmp
        else:
            ji_all[i] = jaccard_index(fdr1.iloc[:max_depth].copy(), fdr2, depth=depth)

    # calculate the weighting for jaccard index
    p = weight_factor ** (1./(num_depths-1))
    w = p*np.ones(num_depths_total)
    w[0] = 1
    w = np.cumprod(w)  # calculate geometric weights
    w = w / w.sum()  # normalize weights to 1

    weighted_mean_ji = np.dot(w, ji_all)
    mean_ji = np.mean(ji)

    return ji, mean_ji, weighted_mean_ji
