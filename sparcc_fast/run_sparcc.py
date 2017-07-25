import utils
from sparcc_functions import sparcc
from bootstrap_correlations import bootstrap_correlations


def sparcc_correlation(table, procs):
    df = utils.biom_to_pandas(table)
    cor, _ = sparcc(df, procs=procs)
    correls = utils.df_to_correls(cor)
    return correls


def sparcc_correlation_w_bootstraps(table, procs, bootstraps):
    df = utils.biom_to_pandas(table)
    cor, _ = sparcc(df)
    correls = utils.df_to_correls(cor)
    correls['p_value'] = bootstrap_correlations(df, cor, bootstraps, procs)
    return correls
