import pandas as pd
import numpy as np
from numpy.random.mtrand import dirichlet

__author__ = 'shafferm'


def variation_mat(frame):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    '''
    x = 1.*np.asarray(frame)
    n,m = x.shape
    if m > 1000:
        return variation_mat_slow(frame)
    else:
        xx = np.tile(x.reshape((n,m,1)) ,(1,1,m))
        xx_t = xx.transpose(0,2,1)
        try:
            l = np.log(1.*xx/xx_t)
            V = l.var(axis=0, ddof=1)
            return V
        except MemoryError:
            return variation_mat_slow(frame)


def variation_mat_slow(frame, shrink=False):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Return the variation matrix of frame.
    Element i,j is the variance of the log ratio of components i and j.
    Slower version to be used in case the fast version runs out of memeory.
    '''
    frame_a = 1.*np.asarray(frame)
    k    = frame_a.shape[1]
    V      = np.zeros((k,k))
    for i in range(k-1):
        for j in range(i+1,k):
            y     = np.array(np.log(frame_a[:,i]/frame_a[:,j]))
            v = np.var(y, ddof=1) # set ddof to divide by (n-1), rather than n, thus getting an unbiased estimator (rather than the ML one).
            V[i,j] = v
            V[j,i] = v
    return V


def to_fractions(frame, p_counts=1, axis=0):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/ and adapted***
    Covert counts to fraction using given method.

    Parameters
    ----------
    method : string {'dirichlet' (default) | 'normalize' | 'pseudo'}
        dirichlet - randomly draw from the corresponding posterior
                    Dirichlet distribution with a uniform prior.
                    That is, for a vector of counts C,
                    draw the fractions from Dirichlet(C+1).
        normalize - simply divide each row by its sum.
        pseudo    - add given pseudo count (defualt 1) to each count and
                    do simple normalization.
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
    '''
    def dir_fun(x):
        a = x+p_counts
        f = dirichlet(a)
        return f
    if isinstance(frame, pd.DataFrame):
        fracs = frame.apply(dir_fun, 1-axis)
    else:
        fracs = np.apply_along_axis(dir_fun, 1-axis, frame)
    return fracs


def basis_var(f, Var_mat, M, V_min = 1e-10):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Estimate the variances of the basis of the compositional data x.
    Assumes that the correlations are sparse (mean correlation is small).
    The element of V_mat are refered to as t_ij in the SparCC paper.
    '''
    # compute basis variances
    try:    M_inv = np.linalg.inv(M)
    except: M_inv = np.linalg.pinv(M)
    V_vec = Var_mat.sum(axis=1)  # elements are t_i's of SparCC paper
    V_base = np.dot(M_inv, V_vec)  # basis variances.
    # if any variances are <0 set them to V_min
    V_base[V_base <= 0] = V_min
    return V_base


def C_from_V(Var_mat, V_base):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Given the estimated basis variances and observed fractions variation matrix,
    compute the basis correlation & covaraince matrices.
    '''
    Vi, Vj = np.meshgrid(V_base, V_base)
    Cov_base = 0.5*(Vi + Vj - Var_mat)
    C_base = Cov_base/ np.sqrt(Vi) / np.sqrt(Vj)
    return C_base, Cov_base


def new_excluded_pair(C, th=0.1, previously_excluded=[]):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Find component pair with highest correlation among pairs that
    weren't previously excluded.
    Return the i,j of pair if it's correlaiton >= than th.
    Otherwise return None.
    '''
    C_temp = np.triu(abs(C),1) # work only on upper triangle, excluding diagonal
    C_temp[zip(*previously_excluded)] = 0
    i,j = np.unravel_index(np.argmax(C_temp), C_temp.shape)
    cmax = C_temp[i,j]
    if cmax > th:
        return i,j
    else:
        return None


def run_sparcc(f, th=.1, xiter=10):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/***
    Estimate the correlations of the basis of the compositional data f.
    Assumes that the correlations are sparse (mean correlation is small).
    '''
    ## observed log-ratio variances
    Var_mat = variation_mat(f)
    Var_mat_temp = Var_mat.copy()
    ## Make matrix from eqs. 13 of SparCC paper such that: t_i = M * Basis_Varainces
    D = Var_mat.shape[0] # number of components
    M = np.ones((D,D)) + np.diag([D-2]*D)
    ## get approx. basis variances and from them basis covariances/correlations
    V_base = basis_var(f, Var_mat_temp, M)
    C_base, Cov_base = C_from_V(Var_mat, V_base)
    ## Refine by excluding strongly correlated pairs
    excluded_pairs = []
    excluded_comp  = np.array([])
    for xi in range(xiter):
        # search for new pair to exclude
        to_exclude = new_excluded_pair(C_base, th, excluded_pairs) #i,j pair, or None
        if to_exclude is None: #terminate if no new pairs to exclude
            break
        # exclude pair
        excluded_pairs.append(to_exclude)
        i,j = to_exclude
        M[i,j] -= 1
        M[j,i] -= 1
        M[i,i] -= 1
        M[j,j] -= 1
        inds = zip(*excluded_pairs)
        Var_mat_temp[inds]   = 0
        Var_mat_temp.T[inds] = 0
        # search for new components to exclude
        nexcluded = np.bincount(np.ravel(excluded_pairs)) #number of excluded pairs for each component
        excluded_comp_prev = set(excluded_comp.copy())
        excluded_comp      = np.where(nexcluded>=D-3)[0]
        excluded_comp_new  = set(excluded_comp) - excluded_comp_prev
        if len(excluded_comp_new)>0:
            print excluded_comp
            # check if enough components left
            if len(excluded_comp) > D-4:
                raise ValueError('Too many component excluded. SparCC will not complete.')
            for xcomp in excluded_comp_new:
                Var_mat_temp[xcomp,:] = 0
                Var_mat_temp[:,xcomp] = 0
                M[xcomp,:] = 0
                M[:,xcomp] = 0
                M[xcomp,xcomp] = 1
        # run another sparcc iteration
        V_base = basis_var(f, Var_mat_temp, M)
        C_base, Cov_base = C_from_V(Var_mat, V_base)
        # set excluded components infered values to nans
        for xcomp in excluded_comp:
            V_base[xcomp] = np.nan
            C_base[xcomp,:] = np.nan
            C_base[:,xcomp] = np.nan
            Cov_base[xcomp,:] = np.nan
            Cov_base[:,xcomp] = np.nan
    return V_base, C_base, Cov_base


def basis_corr(frame, iter=20, th=.1, xiter=10):
    '''
    ***STOLEN FROM https://bitbucket.org/yonatanf/pysurvey/ and adapted***
    ##Merged basis_corr from analysis_methods.py and get_correlations.py and main from get_correlations.py
    Compute correlations between all columns of a counts frame.
    This is a wrapper around pysurvey.analysis.basis_correlations.main

    Parameters
    ----------
    counts : array_like
        2D array of counts. Columns are components, rows are samples.

    Returns
    -------
    cor_med: frame
        Estimated correlation matrix.
        Labels are column labels of input frame.
    cov_med: frame/None
        If method in {SparCC, clr} : Estimated covariance matrix.
        Labels are column labels of input frame.
        Otherwise: None.

    =======   ============ =======   ================================================
    kwarg     Accepts      Default   Desctiption
    =======   ============ =======   ================================================
    iter      int          20        number of estimation iteration to average over.
    th        0<th<1       0.1       exclusion threshold for SparCC.
    xiter     int          10        number of exclusion iterations for sparcc.
    =======   ============ ========= ================================================
    '''
    comps = frame.columns
    cor_list = []  # list of cor matrices from different random fractions
    var_list = []  # list of cov matrices from different random fractions
    for i in xrange(iter):
        fracs = to_fractions(frame)
        k = fracs.shape[1]
        # compute basis variances & correlations
        if k<4:
            raise ValueError, 'Can not detect correlations between compositions of <4 components (%d given)' %k
        v_sparse, cor_sparse, cov_sparse = run_sparcc(fracs, th=th, xiter=xiter)
        tol = 1e-3 # tolerance for correlation range
        if np.max(np.abs(cor_sparse)) > 1 + tol:
            raise ValueError('Sparsity assumption violated. SparCC will not run.')
        var_list.append(np.diag(cov_sparse))
        cor_list.append(cor_sparse)
    cor_array = np.array(cor_list)
    var_med = np.nanmedian(var_list,axis=0) #median variances
    cor_med = np.nanmedian(cor_array,axis=0) #median correlations
    x,y = np.meshgrid(var_med,var_med)
    cov_med = cor_med * x**0.5 * y**0.5

    cor = pd.DataFrame(cor_med, index=comps, columns=comps)
    if cov_med is None:
        cov = None
    else:
        cov = pd.DataFrame(cov_med, index=comps, columns=comps)
    return cor, cov
