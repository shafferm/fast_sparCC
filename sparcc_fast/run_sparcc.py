import utils
from sparcc_functions import basis_corr
from bootstrap_correlations import bootstrap_correlations


def sparcc_correlation(table):
    df = utils.biom_to_pandas(table)
    cor, cov = basis_corr(df)
    correls = utils.df_to_correls(cor)
    return correls


def sparcc_correlation_w_bootstraps(table, procs, bootstraps):
    df = utils.biom_to_pandas(table)
    cor, cov = basis_corr(df)
    correls = utils.df_to_correls(cor)
    correls['p_value'] = bootstrap_correlations(df, cor, bootstraps, procs)
    return correls
