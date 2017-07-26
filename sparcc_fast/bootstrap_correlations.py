import numpy as np
import sparcc_functions
from scipy.spatial.distance import squareform
from functools import partial
from sparcc_functions import permute_w_replacement

__author__ = 'shafferm'


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
