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


def bootstrapped_correlation(bootstrap, df, cor):
    bootstrap = permute_w_replacement(df)
    in_cor = sparcc_functions.basis_corr(bootstrap)
    return np.abs(in_cor) >= cor


def bootstrap_correlations(df, cor, bootstraps=100, procs=1):
    """"""
    # take absolute value of all values in cor for calculating two-sided p-value
    abs_cor = np.abs(squareform(cor, checks=False))
    # create an empty array of significant value counts in same shape as abs_cor
    n_sig = np.zeros(abs_cor.shape)

    if procs == 1:
        for i in xrange(bootstraps):
            n_sig += bootstrapped_correlation(i, df, abs_cor)
    else:
        import multiprocessing
        pool = multiprocessing.Pool(procs)
        print "Number of processors used: " + str(procs)

        # make partial function for use in multiprocessing
        pfun = partial(bootstrapped_correlation, cor=abs_cor, df=df)
        # run multiprocessing
        multi_results = pool.map(pfun, xrange(bootstraps))
        pool.close()
        pool.join()

        # find number of significant results across all bootstraps
        n_sig = np.sum(multi_results, axis=0)

    # get p_values out
    p_val_square = squareform(1. * n_sig / bootstraps, checks=False)
    p_vals = []
    for i in xrange(p_val_square.shape[0]):
        for j in xrange(i + 1, p_val_square.shape[0]):
            p_vals.append(p_val_square[i, j])
    return p_vals
