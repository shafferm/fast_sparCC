import numpy as np
import pandas as pd
import sparcc_functions
from scipy.spatial.distance import squareform
from functools import partial

__author__ = 'shafferm'


def sparcc_correlations(df):
    """"""
    # calculate correlations
    cor, cov = sparcc_functions.basis_corr(df)
    return cor


def boostrapped_correlation(bootstrap, cor):
    in_cor = sparcc_functions.basis_corr(bootstrap, oprint=False)[0]
    in_cor = squareform(in_cor, checks=False)
    return np.abs(in_cor) >= cor


def calculate_sparcc_pvalues(df, cor, bootstraps=100):
    # calculate p-values
    abs_cor = np.abs(squareform(cor, checks=False))
    n_sig = np.zeros(abs_cor.shape)
    bootstrap_frames = sparcc_functions.make_bootstraps(df, bootstraps)
    for i in bootstrap_frames:
        n_sig[np.abs(boostrapped_correlation(i, abs_cor)) >= abs_cor] += 1
    p_vals = squareform(1.*n_sig/bootstraps, checks=False)
    return pd.DataFrame(p_vals, columns=cor.index.values, index=cor.index.values)


def boostrapped_correlation_multi(bootstrap, cor):
    in_cor = sparcc_functions.basis_corr(bootstrap)[0]
    in_cor = squareform(in_cor, checks=False)
    return np.abs(in_cor) >= cor


def calculate_sparcc_pvalues_multi(df, cor, bootstraps=100, procs=None):
    print df.shape; print cor.shape; print bootstraps; print procs
    """"""
    # setup
    import multiprocessing
    if multiprocessing.cpu_count() == 1:
        print "only one processor, will run with one process"
        return calculate_sparcc_pvalues(df, cor, bootstraps)
    elif procs is None:
        pool = multiprocessing.Pool(multiprocessing.cpu_count()-1)
    elif multiprocessing.cpu_count()-1 < procs:
        print "more processors requested than avaliable, will run with all avaliable"
        pool = multiprocessing.Pool(multiprocessing.cpu_count())
    else:
        pool = multiprocessing.Pool(procs)

    # take absolute value of all values in cor for calculating two-sided p-value
    abs_cor = np.abs(squareform(cor, checks=False))
    # create an empty array of significant value counts in same shape as abs_cor
    n_sig = np.zeros(abs_cor.shape)
    print type(n_sig); print n_sig.shape
    print type(abs_cor); print abs_cor.shape
    # make bootstraps
    bootstrap_frames = sparcc_functions.make_bootstraps(df, bootstraps)
    print len(bootstrap_frames); print type(bootstrap_frames); print type(bootstrap_frames[0]); print bootstrap_frames[0].shape
    # make partial function for use in multiprocessing
    pfun = partial(boostrapped_correlation_multi, cor=abs_cor)
    print pfun
    # run multiprocessing
    multi_results = pool.map(pfun, bootstrap_frames)
    pool.close()
    pool.join()

    # find number of significant results across all bootstraps
    for i in multi_results:
        n_sig[i] += 1
    p_vals = squareform(1.*n_sig/bootstraps, checks=False)

    return pd.DataFrame(p_vals, columns=cor.index.values, index=cor.index.values)
