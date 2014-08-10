cimport cython
from libc.math cimport sqrt, log, exp
import numpy as np
cimport numpy as np
import scipy.misc as misc
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
ctypedef fused one_or_two_dim_array:
    np.ndarray[DTYPE_FLOAT_t, ndim=1]
    np.ndarray[DTYPE_FLOAT_t, ndim=2]

# define constants
cdef:
    DTYPE_FLOAT_t HALF = .5  # statically type this constant
    DTYPE_FLOAT_t PI = np.pi
    DTYPE_FLOAT_t ROOT_2PI = sqrt(2*PI)
    DTYPE_FLOAT_t NEGINF = -np.inf


def kde_entropy(np.ndarray[DTYPE_INT_t, ndim=1, mode='c'] x,
                bandwidth_param=None,
                int kfolds=5):
    cdef:
        # np.ndarray[DTYPE_FLOAT_t, ndim=2] log_prob
        # np.ndarray[DTYPE_FLOAT_t, ndim=2] norm_log_prob
        np.ndarray[DTYPE_FLOAT_t, ndim=1] log_prob
        DTYPE_FLOAT_t[::1] log_prob_mv
        np.ndarray[DTYPE_FLOAT_t, ndim=1] norm_log_prob
        DTYPE_INT_t[::1] x_mv = x
        # np.ndarray[DTYPE_FLOAT_t, ndim=1] single_bandwidth_array
        # np.ndarray[DTYPE_FLOAT_t, ndim=1] cv_bws
        DTYPE_FLOAT_t sum_log_prob
        DTYPE_FLOAT_t bandwidth
        DTYPE_INT_t n
    # only make sense to perform this if there is actually more than
    # one mutation
    n = x.shape[0]
    if n <= 1:
        return 1.0, None

    # check if bandwidth selection is necessary
    # if isinstance(bandwidth, numbers.Number):
    if bandwidth_param:
        bandwidth = bandwidth_param
    else:
        # lower fold number if few mutations
        if len(x) < kfolds:
            kfolds = len(x)
        # x.sort()
        # if user provided a list of bandwidths for cross-validation
        # then perform cross-validation with those instead of the default.
        # Otherwise, fall back to DEFAULT_BANDWIDTHS.
        # if isinstance(bandwidth, list):
        #    cv_bws = np.array(list, dtype=DTYPE_FLOAT)
        #    bandwidth = _best_kde_bandwidth(x, bws=cv_bws, kfolds=kfolds)
        # else:
        bandwidth = _best_kde_bandwidth(x, kfolds=kfolds)

    # put the bandwidth in a format accepted by the update_kde_log_prob_matrix
    # function. the bandwidth array only has a single element which is either
    # specified by the user or found from cross-validation
    # single_bandwidth_array = np.array([bandwidth], dtype=DTYPE_FLOAT)

    # populate log probabilities
    # log_prob = np.zeros((n, 1), dtype=DTYPE_FLOAT)
    # log_prob = np.zeros(n, dtype=DTYPE_FLOAT)
    # x = np.ma.masked_array(x)  #
    #update_kde_log_prob_matrix(x, log_prob, bws=bandwidth)
    # log_prob = update_kde_log_prob_matrix(x, log_prob, bws=single_bandwidth_array)
    log_prob_mv = update_kde_log_prob_matrix(x_mv,
                                             #log_prob,
                                             bw=bandwidth)

    # sum_log_prob = masked_logsumexp(log_prob)
    sum_log_prob = logsumexp(log_prob_mv, n)
    log_prob = np.asarray(log_prob_mv)
    norm_log_prob = log_prob - sum_log_prob
    # print norm_log_prob

    # perform calculations
    #kde_skl = KernelDensity(bandwidth=bandwidth)
    #kde_skl.fit(x[:, np.newaxis])
    #log_pdf = kde_skl.score_samples(x[:, np.newaxis])  # return log-likelihood
    #log_sum = logsumexp(log_pdf)
    #norm_log_pdf = log_pdf - log_sum

    ent = mymath.shannon_entropy(np.exp(norm_log_prob)) / np.log2(n)
    # ent = np.exp(mymath.log_shannon_entropy(norm_log_prob)) / float(np.log(n))
    return ent, bandwidth


@cython.boundscheck(False)
@cython.wraparound(False)
cdef DTYPE_FLOAT_t _best_kde_bandwidth(np.ndarray[DTYPE_INT_t, ndim=1] y,
                                       int kfolds=5):
    cdef:
        int n, num_bws, best_bw_ix
        np.ndarray[DTYPE_INT_t, ndim=1] rand_order  # random index order for cv
        float fold_inc  # fraction of total data that is in the cv test set
        float fold_pos  # end position of cv test set range (float: .1, .2, ...)
        np.ndarray[DTYPE_INT_t, ndim=1] mask
        np.ndarray[DTYPE_INT_t, ndim=1] masked_y
        np.ndarray[DTYPE_FLOAT_t, ndim=2] log_prob_matrix
        np.ndarray[DTYPE_FLOAT_t, ndim=1] bws_log_prob
        np.ndarray[DTYPE_FLOAT_t, ndim=1] bws=np.array([], dtype=DTYPE_FLOAT)

        # default KDE bandwiths to try with cross-validation
        np.ndarray[DTYPE_FLOAT_t, ndim=1] DEFAULT_BANDWIDTHS =  np.array([3.0, 5.0, 10.0, 15.0,
                                                                          20.0, 25.0, 40.0, 50.0],
                                                                          dtype=DTYPE_FLOAT)

    # setup dimension info
    n = y.shape[0]  # length of array
    num_bws = bws.shape[0]
    if num_bws == 0:
        bws = DEFAULT_BANDWIDTHS
        num_bws = bws.shape[0]
    log_prob_matrix = np.zeros((n, num_bws), dtype=DTYPE_FLOAT)

    # put indices in random order to make CV more robust
    rand_order = np.random.choice(n, n, replace=False)

    # perform k-fold cross-validation
    fold_inc = 1./kfolds
    mask = np.zeros(n, dtype=DTYPE_INT)
    masked_y = np.ma.array(y, mask=mask)
    for fold_pos in np.linspace(fold_inc, 1, kfolds):
        # mask the test set
        test_start = <int>(n*(fold_pos-fold_inc))
        test_end = <int>(n*fold_pos)
        masked_y.mask = np.ma.nomask
        masked_y[rand_order[test_start:test_end]].mask = np.ma.masked
        log_prob_matrix = masked_update_kde_log_prob_matrix(y, log_prob_matrix, bws)

    # find best bandwidth after cross-validation. The "best" bandwidth is
    # defined as having the maximum log likelihood
    bws_log_prob = log_prob_matrix.sum(axis=0)
    best_bw_ix = bws_log_prob.argmax()
    return bws[best_bw_ix]


@cython.boundscheck(False)
@cython.wraparound(False)
cdef np.ndarray[DTYPE_FLOAT_t, ndim=2] masked_update_kde_log_prob_matrix(np.ndarray[DTYPE_INT_t, ndim=1] x,
                                                                         np.ndarray[DTYPE_FLOAT_t, ndim=2] log_prob,
                                                                         np.ndarray[DTYPE_FLOAT_t, ndim=1] bws):
    cdef:
        np.ndarray[DTYPE_UINT_t, ndim=1] masked_ixs
        # DTYPE_INT_t n
        DTYPE_INT_t masked_ix
        int i, j, num_bws=bws.shape[0], n, num_masked_ixs
        # DTYPE_FLOAT_t tmp_val

    if not isinstance(x, np.ma.MaskedArray):
        x = np.ma.masked_array(x)

    masked_ixs = np.where(x.mask==1)[0].view(DTYPE_UINT)
    num_masked_ixs = masked_ixs.shape[0]
    n = x.shape[0] - num_masked_ixs  # number of not masked elements
    for i in range(num_bws):
        bw = bws[i]  # get current bandwidth value
        for j in range(num_masked_ixs):
            masked_ix = masked_ixs[j]  # get current masked index
            log_prob[masked_ix, i] = masked_kde_log_prob(x.data[masked_ix],
                                                         x, n, bw)
            # log_prob[masked_ix, i] = tmp_val
    return log_prob


@cython.boundscheck(False)
@cython.wraparound(False)
cdef DTYPE_FLOAT_t[::1] update_kde_log_prob_matrix(DTYPE_INT_t[::1] x,
                                                                # np.ndarray[DTYPE_FLOAT_t, ndim=1] log_prob,
                                                                DTYPE_FLOAT_t bw):
                                                                # np.ndarray[DTYPE_FLOAT_t, ndim=1] bws):
    cdef:
        # DTYPE_INT_t n, i, j, , num_bws = bws.shape[0]
        # int i, num_bws=bws.shape[0]
        int j, n, i
        # DTYPE_FLOAT_t[::1] log_prob_mv = log_prob
        # DTYPE_FLOAT_t tmp_val
        DTYPE_FLOAT_t u
        # log sum of gaussian kernels (with out constant factors)
        DTYPE_FLOAT_t log_kern_sum, lp

    n = x.shape[0]  # number of elements
    # cdef DTYPE_INT_t[::1] x_mv = x  # cython memory view for x
    cdef DTYPE_FLOAT_t[::1] log_prob_mv = np.empty(n),
    cdef DTYPE_FLOAT_t[:, ::1] z = np.empty((n, n))
    # cdef np.ndarray[DTYPE_FLOAT_t, ndim=2, mode='c'] z_array
    # for i in range(num_bws):
        # bw = bws[i]  # get current bandwidth value
    u = lp = log_kern_sum = 0.0
    for j in range(n):
        for i in range(j, n):
            u = (x[j] - x[i]) / bw
            z[i, j] = z[i, j] = -HALF * u * u
        # log_kern_sum = logsumexp(z[j], n)
        # lp = -log(ROOT_2PI) - log(n*bw) + log_kern_sum
        # log_prob_mv[j] = lp
        pass
    # z_array = np.asarray(z)
    log_kern_sum, lp = 0.0, 0.0
    log_prob_mv = matrix_logsumexp(z, n, n)
    for j in range(n):
        log_prob_mv[j] -= log(ROOT_2PI) + log(n*bw)


            #log_prob_mv[j] = kde_log_prob(x_mv[j],
                                      #x_mv, n, bw)
            # log_prob[j, i] = tmp_val
    return log_prob_mv


cdef DTYPE_FLOAT_t masked_kde_log_prob(DTYPE_INT_t val_of_interest,
                                       np.ndarray[DTYPE_INT_t, ndim=1] x,
                                       int num_not_masked,
                                       DTYPE_FLOAT_t bw):
    cdef:
        # temporary arrays
        np.ndarray[DTYPE_FLOAT_t, ndim=1] u
        np.ndarray[DTYPE_FLOAT_t, ndim=1] z
        # log sum of gaussian kernels (with out constant factors)
        DTYPE_FLOAT_t log_kern_sum
        # log probability of value of interest
        DTYPE_FLOAT_t log_prob

    u = (val_of_interest - x) / bw
    # z = -HALF*np.ma.dot(u, u)
    z = -HALF * np.ma.power(u, 2)
    log_kern_sum = masked_logsumexp(z)
    log_prob = -log(ROOT_2PI) - log(num_not_masked*bw) + log_kern_sum
    return log_prob


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline DTYPE_FLOAT_t kde_log_prob(DTYPE_INT_t val_of_interest,
                                       DTYPE_INT_t[::1] x,
                                       DTYPE_INT_t n,
                                       DTYPE_FLOAT_t bw):
    cdef:
        # temporary arrays
        DTYPE_FLOAT_t u
        DTYPE_FLOAT_t[::1] z = np.empty(n)
        # log sum of gaussian kernels (with out constant factors)
        DTYPE_FLOAT_t log_kern_sum
        # log probability of value of interest
        DTYPE_FLOAT_t log_prob

    for i in range(n):
        u = (val_of_interest - x[i]) / bw
        z[i] = -HALF * u * u
    log_kern_sum = logsumexp(z, n)
    log_prob = -log(ROOT_2PI) - log(n*bw) + log_kern_sum
    return log_prob


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline DTYPE_FLOAT_t[::1] matrix_logsumexp(DTYPE_FLOAT_t[:, ::1] a,
                                                DTYPE_INT_t nrows,
                                                DTYPE_INT_t ncols):
    cdef:
        DTYPE_FLOAT_t shifted_log_sum, shifted_sum = 0.0
        DTYPE_INT_t i, j
        np.ndarray[DTYPE_FLOAT_t, ndim=1, mode='c'] tmp = np.empty(ncols)
        DTYPE_FLOAT_t[::1] a_max
        DTYPE_FLOAT_t[::1] out = np.empty(ncols)
    tmp.fill(NEGINF)
    a_max = tmp

    # find max
    for i in range(ncols):
        for j in range(nrows):
            if a[j, i] > a_max[i]:
                a_max[i] = a[j, i]

    # use logsumexp trick
    for i in range(ncols):
        shifted_sum = 0.0
        for j in range(nrows):
            shifted_sum += exp(a[j, i] - a_max[i])
        shifted_log_sum = log(shifted_sum)
        out[i] = a_max[i] + shifted_log_sum
    return out


cdef inline DTYPE_FLOAT_t logsumexp(DTYPE_FLOAT_t[::1] a,
                                    DTYPE_INT_t n):
    cdef:
        DTYPE_FLOAT_t a_max = NEGINF, shifted_log_sum, shifted_sum = 0.0, out
        DTYPE_INT_t i

    # find max
    for i in range(n):
        if a[i] > a_max:
            a_max = a[i]

    # use logsumexp trick
    for i in range(n):
        shifted_sum += exp(a[i] - a_max)
    shifted_log_sum = log(shifted_sum)
    out = a_max + shifted_log_sum
    return out


cdef DTYPE_FLOAT_t masked_logsumexp(one_or_two_dim_array a):
    cdef:
        DTYPE_FLOAT_t a_max
        DTYPE_FLOAT_t shifted_log_sum
        # np.ndarray[DTYPE_FLOAT_t, ndim=1] tmp_exp
        one_or_two_dim_array tmp_exp
        DTYPE_FLOAT_t out

    a_max = a.max()
    tmp_exp = np.ma.exp(a - a_max)
    shifted_log_sum = np.log(np.ma.sum(tmp_exp))
    out = a_max + shifted_log_sum
    return out
