cimport cython
import numpy as np
cimport numpy as np
from libc.math cimport sqrt, log, exp
from libcpp.vector cimport vector

import mymath

# define data types for kde function
DTYPE_INT = np.int
DTYPE_UINT = np.uint
DTYPE_FLOAT = np.double
DTYPE_LONG_FLOAT = np.longdouble
# define compile time data types
ctypedef np.int_t DTYPE_INT_t
ctypedef np.uint_t DTYPE_UINT_t
ctypedef np.double_t DTYPE_FLOAT_t
ctypedef np.longdouble_t DTYPE_LONG_FLOAT_t

# define constants
cdef:
    DTYPE_FLOAT_t HALF = .5  # statically type this constant
    DTYPE_FLOAT_t PI = np.pi
    DTYPE_FLOAT_t ROOT_2PI = sqrt(2*PI)
    DTYPE_FLOAT_t NEGINF = -np.inf


def uniform_kde_entropy(DTYPE_INT_t[::1] pos, bw):
    cdef:
        int used_bw, n=pos.shape[0], p, prev_pos = -1, num_unique_pos = 0, i
        #np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] pos_array
        np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] kernel_cts
        np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] unique_kernel_cts
        np.ndarray[DTYPE_FLOAT_t, ndim=1, mode='c'] norm_kernel_cts
        vector[int] unique_pos_cts

    # handle low count cases
    if n < 1:
        # return immediately if no mutations
        return 1.0, 3

    # define bandwidth
    if bw is None:
        used_bw = _best_uniform_kernel_bandwidth(pos)
    else:
        used_bw = bw

    # compute the counts for each position
    kernel_cts = np.ascontiguousarray(discrete_uniform_kernel(pos, used_bw),
                                      dtype=DTYPE_INT)


    # remove duplicate positions
    for i in range(n):
        if prev_pos != <int> pos[i]:
            unique_pos_cts.push_back(kernel_cts[i])
            prev_pos = pos[i]
            num_unique_pos += 1
    unique_kernel_cts = np.empty(num_unique_pos, dtype=DTYPE_INT)
    for i in range(num_unique_pos):
        unique_kernel_cts[i] = unique_pos_cts[i]

    #pos_array = np.ascontiguousarray(pos, dtype=DTYPE_INT)
    #unique_pos = np.unique(pos)
    #unique_cts = []
    #for p in unique_pos:
        #unique_cts.append(kernel_cts[pos_array==p].max())
    #unique_kernel_cts = np.ascontiguousarray(unique_cts, dtype=DTYPE_INT)

    # calculate entropy
    norm_kernel_cts = unique_kernel_cts / float(unique_kernel_cts.sum())
    ent = mymath.shannon_entropy(norm_kernel_cts) / np.log2(n)
    return ent, used_bw


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline DTYPE_INT_t[::1] discrete_uniform_kernel(DTYPE_INT_t[::1] pos,
                                                     DTYPE_INT_t bandwidth):
    cdef:
        int n = pos.shape[0], first_ix=0, second_ix=0, current_pos, tmp_ct = 0, i
        DTYPE_INT_t[::1] kernel_cts = np.empty(n, dtype=DTYPE_INT)

    for i in range(n):
        current_pos = pos[i]
        while second_ix < n and (pos[second_ix] <= current_pos + bandwidth):
            tmp_ct += 1
            second_ix += 1
        while pos[first_ix] < current_pos - bandwidth:
            tmp_ct -= 1
            first_ix += 1
        kernel_cts[i] = tmp_ct

    return kernel_cts


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline DTYPE_INT_t[::1] masked_discrete_uniform_kernel(DTYPE_INT_t[::1] pos,
                                                            DTYPE_INT_t[::1] mask,
                                                            DTYPE_INT_t bandwidth):
    cdef:
        int n = pos.shape[0], first_ix=0, second_ix=0, current_pos, tmp_ct = 0, i
        DTYPE_INT_t[::1] kernel_cts = np.empty(n, dtype=DTYPE_INT)

    for i in range(n):
        current_pos = pos[i]
        while second_ix < n and (pos[second_ix] <= (current_pos + bandwidth)):
            tmp_ct += 1 * (1 - mask[second_ix])
            second_ix += 1
        while pos[first_ix] < (current_pos - bandwidth):
            tmp_ct -= 1 * (1 - mask[first_ix])
            first_ix += 1
        kernel_cts[i] = tmp_ct

    return kernel_cts


@cython.boundscheck(False)
@cython.wraparound(False)
cdef DTYPE_INT_t _best_uniform_kernel_bandwidth(DTYPE_INT_t[::1] mypos,
                                                int kfolds=5):
    cdef:
        int n = mypos.shape[0]
        np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] rand_order
        np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] mask

        int test_start=0, test_end=0
        DTYPE_FLOAT_t fold_inc, fold_pos

        # default KDE bandwiths to try with cross-validation
        np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] DEFAULT_BANDWIDTHS =  np.array([3, 5, 10, 15,
                                                                                  25, 40, 50, 75],
                                                                                 dtype=DTYPE_INT)
        # likelihood matrices
        np.ndarray[DTYPE_FLOAT_t, ndim=2, mode='c'] prob_matrix = np.empty((n, DEFAULT_BANDWIDTHS.shape[0]),
                                                                           dtype=DTYPE_FLOAT)
        np.ndarray[DTYPE_FLOAT_t, ndim=2, mode='c'] log_prob_matrix

    # handle low count cases
    if n <= 1:
        return DEFAULT_BANDWIDTHS[0]
    elif n < 5:
        kfolds = n

    rand_order = np.random.choice(n, n, replace=False)
    fold_inc = 1./kfolds
    mask = np.zeros(n, dtype=DTYPE_INT)
    for fold_pos in np.linspace(fold_inc, 1, kfolds):
        # mask the test set
        mask[rand_order[test_start:test_end]] = 0
        test_start = <int>(n*(fold_pos-fold_inc))
        test_end = <int>(n*fold_pos)
        mask[rand_order[test_start:test_end]] = 1

        prob_matrix = calc_uniform_kde_likelihood(prob_matrix,
                                                  mypos,
                                                  mask,
                                                  DEFAULT_BANDWIDTHS)

    # find best bandwidth after cross-validation. The "best" bandwidth is
    # defined as having the maximum log likelihood
    prob_matrix += .01 / ((2*DEFAULT_BANDWIDTHS+1) * <DTYPE_FLOAT_t> n)
    log_prob_matrix = np.log(prob_matrix)
    bws_log_prob = log_prob_matrix.sum(axis=0)
    best_bw_ix = bws_log_prob.argmax()
    return DEFAULT_BANDWIDTHS[best_bw_ix]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[DTYPE_FLOAT_t, ndim=2, mode='c'] calc_uniform_kde_likelihood(np.ndarray[DTYPE_FLOAT_t, ndim=2, mode='c'] lh,
                                                                             DTYPE_INT_t[::1] pos,
                                                                             DTYPE_INT_t[::1] mask,
                                                                             DTYPE_INT_t[::1] bws):
    cdef:
        int num_bws = bws.shape[0]
        int num_masked = 0
        int j, i
        int n = pos.shape[0]
        DTYPE_INT_t bw
        DTYPE_INT_t[::1] bw_cts
    for i in range(num_bws):
        bw = bws[i]
        bw_cts = masked_discrete_uniform_kernel(pos, mask, bw)
        for j in range(n):
            if mask[j] == 1:
                lh[j, i] = (<DTYPE_FLOAT_t> bw_cts[j]) / <DTYPE_FLOAT_t> (n * (2*bw+1))
    return lh
