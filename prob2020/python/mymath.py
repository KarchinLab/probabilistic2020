import numpy as np
from scipy.misc import logsumexp
# from sklearn.neighbors import KernelDensity
# from sklearn.grid_search import GridSearchCV


def shannon_entropy(p):
    """Calculates shannon entropy in bits.

    Parameters
    ----------
    p : np.array
        array of probabilities

    Returns
    -------
    shannon entropy in bits
    """
    return -np.sum(np.where(p!=0, p * np.log2(p), 0))


def log_shannon_entropy(log_p):
    """Calculates shannon entropy in bits.

    Parameters
    ----------
    p : np.array
        array of probabilities

    Returns
    -------
    shannon entropy in bits
    """
    out = -logsumexp(log_p + np.log(log_p))
    return out


def max_shannon_entropy(n):
    """Returns max possible entropy given "n" mutations.

    The maximum possible entropy is the entropy of the
    uniform distribution. The uniform distribution has
    entropy equal to log(n) (I will use base 2).

    Parameters
    ----------
    n : int
        total mutation counts

    Returns
    -------
    max possible shannon entropy in bits
    """
    if n <= 0:
        return 0.
    return float(np.log2(n))


def normalized_mutation_entropy(counts, total_cts=None):
    """Calculate the normalized mutation entropy based on a list/array
    of mutation counts.

    Note: Any grouping of mutation counts together should be done before hand

    Parameters
    ----------
    counts : np.array_like
        array/list of mutation counts

    Returns
    -------
    norm_ent : float
        normalized entropy of mutation count distribution.
    """
    cts = np.asarray(counts, dtype=float)
    if total_cts is None:
        total_cts = np.sum(cts)
    if total_cts > 1:
        p = cts / total_cts
        ent = shannon_entropy(p)
        max_ent = max_shannon_entropy(total_cts)
        norm_ent = ent / max_ent
    else:
        norm_ent = 1.0
    return norm_ent


def kl_divergence(p, q):
    """Compute the Kullback-Leibler (KL) divergence for discrete distributions.

    Parameters
    ----------
    p : np.array
        "Ideal"/"true" Probability distribution
    q : np.array
        Approximation of probability distribution p

    Returns
    -------
    kl : float
        KL divergence of approximating p with the distribution q
    """
    # make sure numpy arrays are floats
    p = p.astype(float)
    q = q.astype(float)

    # compute kl divergence
    kl = np.sum(np.where(p!=0, p*np.log2(p/q), 0))
    return kl


def js_divergence(p, q):
    """Compute the Jensen-Shannon Divergence between two discrete distributions.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_div : float
        js divergence between the two distrubtions
    """
    m = .5 * (p+q)
    js_div = .5*kl_divergence(p, m) + .5*kl_divergence(q, m)
    return js_div


def js_distance(p, q):
    """Compute the Jensen-Shannon distance between two discrete distributions.

    NOTE: JS divergence is not a metric but the sqrt of JS divergence is a
    metric and is called the JS distance.

    Parameters
    ----------
    p : np.array
        probability mass array (sums to 1)
    q : np.array
        probability mass array (sums to 1)

    Returns
    -------
    js_dist : float
        Jensen-Shannon distance between two discrete distributions
    """
    js_dist = np.sqrt(js_divergence(p, q))
    return js_dist


"""
def kde_entropy(x, bandwidth=None, folds=5):
    # lower fold number if few mutations
    n = len(x)
    if n < folds:
        folds = n

    # check if bandwidth selection is necessary
    if bandwidth is None:
        grid = GridSearchCV(KernelDensity(),
                            {'bandwidth': np.array([3.0, 5.0, 10.0, 15.0,
                                                    20.0, 25.0, 40.0, 50.0])},
                            cv=folds) # 5-fold cross-validation
        grid.fit(x[:, None])
        bandwidth = grid.best_params_['bandwidth']

    # perform calculations
    kde_skl = KernelDensity(bandwidth=bandwidth)
    kde_skl.fit(x[:, np.newaxis])
    log_pdf = kde_skl.score_samples(x[:, np.newaxis])  # return log-likelihood
    log_sum = logsumexp(log_pdf)
    norm_log_pdf = log_pdf - log_sum
    ent = shannon_entropy(np.exp(norm_log_pdf)) / float(np.log2(n))
    return ent, bandwidth
"""
