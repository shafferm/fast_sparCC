import numpy as np
import sparcc_functions
from scipy.spatial.distance import squareform
from functools import partial

__author__ = 'shafferm'


def permute_w_replacement(frame):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey and adapted***
    Permute the frame values across the given axis.
    Create simulated dataset were the counts of each component (column)
    in each sample (row), are randomly sampled from the all the
    counts of that component in all samples.

    Parameters
    ----------
    frame : DataFrame
        Frame to permute.

    Returns
    -------
    Permuted DataFrame (new instance).
    '''
    from numpy.random import randint
    s = frame.shape[0]
    fun = lambda x: x.values[randint(0,s,(1,s))][0]
    perm = frame.apply(fun, axis=0)
    return perm


def bootstrapped_correlation(bootstrap, cor):
    in_cor = sparcc_functions.basis_corr(bootstrap)[0]
    in_cor = squareform(in_cor, checks=False)
    return np.abs(in_cor) >= cor


def bootstrap_correlations_single(df, cor, bootstraps=100):
    """"""
    # calculate p-values
    abs_cor = np.abs(squareform(cor, checks=False))
    n_sig = np.zeros(abs_cor.shape)
    for i in xrange(bootstraps):
        bootstrap = permute_w_replacement(df)
        n_sig[bootstrapped_correlation(bootstrap, abs_cor)] += 1
    # expand back out to square
    p_val_square = squareform(1.*n_sig/bootstraps, checks=False)

    # get p_values out
    p_vals = []
    for i in xrange(p_val_square.shape[0]):
        for j in xrange(i+1, p_val_square.shape[0]):
            p_vals.append(p_val_square[i, j])
    return p_vals


def bootstrap_correlations(df, cor, bootstraps=100, procs=None):
    """"""
    # setup
    import multiprocessing

    if procs is None:
        if multiprocessing.cpu_count() == 1:
            procs=1
        else:
            procs = multiprocessing.cpu_count()-1

    if procs == 1:
        bootstrap_correlations_single(df, cor, bootstraps)

    pool = multiprocessing.Pool(procs)
    print "Number of processors used: " + str(procs)

    # take absolute value of all values in cor for calculating two-sided p-value
    abs_cor = np.abs(squareform(cor, checks=False))
    # create an empty array of significant value counts in same shape as abs_cor
    n_sig = np.zeros(abs_cor.shape)

    # make partial function for use in multiprocessing
    pfun = partial(bootstrapped_correlation, cor=abs_cor, df=df)
    # run multiprocessing
    multi_results = pool.map(pfun, range(bootstraps))
    pool.close()
    pool.join()

    # find number of significant results across all bootstraps
    for i in multi_results:
        n_sig[i] += 1
    p_val_square = squareform(1.*n_sig/bootstraps, checks=False)

    # get p_values out
    p_vals = []
    for i in xrange(p_val_square.shape[0]):
        for j in xrange(i + 1, p_val_square.shape[0]):
            p_vals.append(p_val_square[i, j])
    return p_vals
