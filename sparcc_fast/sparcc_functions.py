import pandas as pd
import numpy as np
from numpy.random.mtrand import dirichlet
from functools import partial
from scipy.spatial.distance import squareform

__author__ = 'shafferm'


def variation_mat(frame):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j."""

    x = 1.*np.asarray(frame)
    n, m = x.shape
    if m > 1000:
        return variation_mat_slow(frame)
    else:
        xx = np.tile(x.reshape((n, m, 1)), (1, 1, m))
        xx_t = xx.transpose(0, 2, 1)
        try:
            l = np.log(1.*xx/xx_t)
            v_mat = l.var(axis=0, ddof=1)
            return v_mat
        except MemoryError:
            return variation_mat_slow(frame)


def variation_mat_slow(frame):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    Slower version to be used in case the fast version runs out of memory."""
    frame_a = 1.*np.asarray(frame)
    k = frame_a.shape[1]
    v_mat = np.zeros((k, k))
    for i in range(k-1):
        for j in range(i+1, k):
            y = np.array(np.log(frame_a[:, i]/frame_a[:, j]))
            # set ddof to divide by (n-1), rather than n, thus getting an unbiased estimator (rather than the ML one).
            v = np.var(y, ddof=1)
            v_mat[i, j] = v
            v_mat[j, i] = v
    return v_mat


def to_fractions(frame, p_counts=1, axis=0):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/ and adapted***
    Covert counts to fraction using given method.

    Parameters
    ----------
    frame : pd.DataFrame
        2D array of counts. Columns are components, rows are samples.
    p_counts : int/float (default 1)
        The value of the pseudo counts to add to all counts.
        Used only if method is dirichlet
    axis : {0 | 1}
        0 : normalize each row.
        1 : normalize each column.

    Returns
    -------
    fracs: frame/array
        Estimated component fractions.
        Returns new instance of same class as input frame.
    """
    def dir_fun(x):
        a = x+p_counts
        f = dirichlet(a)
        return f
    if isinstance(frame, pd.DataFrame):
        fracs = frame.apply(dir_fun, 1-axis)
    else:
        fracs = np.apply_along_axis(dir_fun, 1-axis, frame)
    return fracs


def basis_var(var_mat, m, v_min=1e-10):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Estimate the variances of the basis of the compositional data x.
    Assumes that the correlations are sparse (mean correlation is small).
    The element of var_mat are refered to as t_ij in the SparCC paper.
    """
    # compute basis variances
    try:
        m_inv = np.linalg.inv(m)
    except:
        m_inv = np.linalg.pinv(m)
    v_vec = var_mat.sum(axis=1)  # elements are t_i's of SparCC paper
    v_base = np.dot(m_inv, v_vec)  # basis variances.
    # if any variances are <0 set them to V_min
    v_base[v_base <= 0] = v_min
    return v_base


def c_from_v(var_mat, v_base):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Given the estimated basis variances and observed fractions variation matrix,
    compute the basis correlation & covaraince matrices.
    """
    v_i, v_j = np.meshgrid(v_base, v_base)
    cov_base = 0.5*(v_i + v_j - var_mat)
    c_base = cov_base/np.sqrt(v_i)/np.sqrt(v_j)
    return c_base, cov_base


def new_excluded_pair(c_mat, th=0.1, previously_excluded=None):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Find component pair with highest correlation among pairs that
    weren't previously excluded.
    Return the i,j of pair if it's correlaiton >= than th.
    Otherwise return None.
    """
    if previously_excluded is None:
        previously_excluded = []
    c_temp = np.triu(abs(c_mat), 1)  # work only on upper triangle, excluding diagonal
    c_temp[zip(*previously_excluded)] = 0
    i, j = np.unravel_index(np.argmax(c_temp), c_temp.shape)
    cmax = c_temp[i, j]
    if cmax > th:
        return i, j
    else:
        return None


def run_sparcc(f, th=.1, xiter=10):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Estimate the correlations of the basis of the compositional data f.
    Assumes that the correlations are sparse (mean correlation is small).
    """
    # observed log-ratio variances
    var_mat = variation_mat(f)
    var_mat_temp = var_mat.copy()
    # Make matrix from eqs. 13 of SparCC paper such that: t_i = M * Basis_Varainces
    d = var_mat.shape[0]  # number of components
    m = np.ones((d, d)) + np.diag([d-2]*d)
    # get approx. basis variances and from them basis covariances/correlations
    v_base = basis_var(var_mat_temp, m)
    c_base, cov_base = c_from_v(var_mat, v_base)
    # Refine by excluding strongly correlated pairs
    excluded_pairs = []
    excluded_comp = np.array([])
    for xi in range(xiter):
        # search for new pair to exclude
        to_exclude = new_excluded_pair(c_base, th, excluded_pairs)  # i,j pair, or None
        if to_exclude is None:  # terminate if no new pairs to exclude
            break
        # exclude pair
        excluded_pairs.append(to_exclude)
        i, j = to_exclude
        m[i, j] -= 1
        m[j, i] -= 1
        m[i, i] -= 1
        m[j, j] -= 1
        inds = zip(*excluded_pairs)
        var_mat_temp[inds] = 0
        var_mat_temp.T[inds] = 0
        # search for new components to exclude
        nexcluded = np.bincount(np.ravel(excluded_pairs))  # number of excluded pairs for each component
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp = np.where(nexcluded >= d-3)[0]
        excluded_comp_new = set(excluded_comp) - excluded_comp_prev
        if len(excluded_comp_new) > 0:
            print excluded_comp
            # check if enough components left
            if len(excluded_comp) > d-4:
                raise ValueError('Too many component excluded. SparCC will not complete.')
            for xcomp in excluded_comp_new:
                var_mat_temp[xcomp, :] = 0
                var_mat_temp[:, xcomp] = 0
                m[xcomp, :] = 0
                m[:, xcomp] = 0
                m[xcomp, xcomp] = 1
        # run another sparcc iteration
        v_base = basis_var(var_mat_temp, m)
        c_base, cov_base = c_from_v(var_mat, v_base)
        # set excluded components infered values to nans
        for xcomp in excluded_comp:
            v_base[xcomp] = np.nan
            c_base[xcomp, :] = np.nan
            c_base[:, xcomp] = np.nan
            cov_base[xcomp, :] = np.nan
            cov_base[:, xcomp] = np.nan
    return v_base, c_base, cov_base


def calc_sparcc(frame, th, xiter, tol):
    fracs = to_fractions(frame)
    k = fracs.shape[1]
    # compute basis variances & correlations
    if k < 4:
        raise ValueError('Can not detect correlations between compositions of <4 components (%d given)' % k)
    v_sparse, cor_sparse, cov_sparse = run_sparcc(fracs, th=th, xiter=xiter)
    if np.max(np.abs(cor_sparse)) > 1 + tol:
        raise ValueError('Sparsity assumption violated. SparCC will not run.')
    return cor_sparse, np.diag(cov_sparse)


def repeater(x, reps):
    for i in xrange(reps):
        yield x


def basis_corr(frame, iters=20, th=.1, xiter=10, tol=1e-3, procs=1):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/ and adapted***
    ##Merged basis_corr from analysis_methods.py and get_correlations.py and main from get_correlations.py
    Compute correlations between all columns of a counts frame.
    This is a wrapper around pysurvey.analysis.basis_correlations.main

    Parameters
    ----------
    frame : pd.Dataframe
        2D array of counts. Columns are components, rows are samples.
    iters : int
        default 20, number of estimation iteration to average over.
    th : float
        0<th<1, default 0.1, exclusion threshold for SparCC.
    xiter : int
        default 10, number of exclusion iterations for sparcc.
    tol: float
        default 1e-3, tolerance for correlation range
    procs: int
        number of processors to use

    Returns
    -------
    cor_med: frame
        Estimated correlation matrix.
        Labels are column labels of input frame.
    cov_med: frame/None
        If method in {SparCC, clr} : Estimated covariance matrix.
        Labels are column labels of input frame.
        Otherwise: None.
    """
    cor_list = []  # list of cor matrices from different random fractions
    var_list = []  # list of cov matrices from different random fractions
    if procs == 1:
        for i in xrange(iters):
            cor_sparse, cov_sparse = calc_sparcc(frame, th, xiter, tol)
            cor_list.append(cor_sparse)
            var_list.append(cov_sparse)
    else:
        import multiprocessing
        pool = multiprocessing.Pool(procs)
        pfun = partial(calc_sparcc, th=th, xiter=xiter, tol=tol)
        results = pool.map(pfun, repeater(frame, iters))
        cor_list = [i[0] for i in results]
        var_list = [i[1] for i in results]
    cor_array = np.array(cor_list)
    var_med = np.nanmedian(var_list, axis=0)  # median variances
    cor_med = np.nanmedian(cor_array, axis=0)  # median correlations
    x, y = np.meshgrid(var_med, var_med)
    cov_med = cor_med * x**0.5 * y**0.5

    return squareform(cov_med, checks=False)


def sparcc(frame, iters=20, th=.1, xiter=10, tol=1e-3, procs=1):
    """
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/ and adapted***
    ##Merged basis_corr from analysis_methods.py and get_correlations.py and main from get_correlations.py
    Compute correlations between all columns of a counts frame.
    This is a wrapper around pysurvey.analysis.basis_correlations.main

    Parameters
    ----------
    frame : pd.Dataframe
        2D array of counts. Columns are components, rows are samples.
    iters : int
        default 20, number of estimation iteration to average over.
    th : float
        0<th<1, default 0.1, exclusion threshold for SparCC.
    xiter : int
        default 10, number of exclusion iterations for sparcc.
    tol: float
        default 1e-3, tolerance for correlation range
    procs: int
        number of processors to use

    Returns
    -------
    cor_med: frame
        Estimated correlation matrix.
        Labels are column labels of input frame.
    cov_med: frame/None
        If method in {SparCC, clr} : Estimated covariance matrix.
        Labels are column labels of input frame.
        Otherwise: None.
    """
    comps = frame.columns
    cor_list = []  # list of cor matrices from different random fractions
    var_list = []  # list of cov matrices from different random fractions
    if procs == 1:
        for i in xrange(iters):
            cor_sparse, cov_sparse = calc_sparcc(frame, th, xiter, tol)
            cor_list.append(cor_sparse)
            var_list.append(cov_sparse)
    else:
        import multiprocessing
        pool = multiprocessing.Pool(procs)
        pfun = partial(calc_sparcc, th=th, xiter=xiter, tol=tol)
        results = pool.map(pfun, repeater(frame, iters))
        cor_list = [i[0] for i in results]
        var_list = [i[1] for i in results]
    cor_array = np.array(cor_list)
    var_med = np.nanmedian(var_list, axis=0)  # median variances
    cor_med = np.nanmedian(cor_array, axis=0)  # median correlations
    x, y = np.meshgrid(var_med, var_med)
    cov_med = cor_med * x**0.5 * y**0.5

    cor = pd.DataFrame(cor_med, index=comps, columns=comps)
    if cov_med is None:
        cov = None
    else:
        cov = pd.DataFrame(cov_med, index=comps, columns=comps)
    return cor, cov
