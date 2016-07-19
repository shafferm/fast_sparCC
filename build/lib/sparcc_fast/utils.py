import numpy as np
from scipy.stats import rankdata
import pandas as pd
import networkx as nx


__author__ = 'shafferm'


def bh_adjust(p_vals):
    """benjamini-hochberg p-value adjustment"""
    p_vals = np.array(p_vals)
    return p_vals*len(p_vals)/rankdata(p_vals, method='dense')


def bonferroni_adjust(p_vals):
    """bonferroni p-value adjustment"""
    return [i*len(p_vals) for i in p_vals]


def print_delimited(out_fp, lines, header=None):
    """print a tab delimited file with optional header"""
    out = open(out_fp, 'w')
    if header is not None:
        out.write('\t'.join([str(i) for i in header])+'\n')
    for line in lines:
        out.write('\t'.join([str(i) for i in line])+'\n')
    out.close()


def biom_to_pandas(table):
    # convert to pandas dataframe
    return pd.DataFrame(np.transpose(table.matrix_data.todense()), index=table.ids(), columns=table.ids(axis="observation"))


def df_to_correls(cor, p_vals=None):
    # TODO: change correls to dataframe
    header = ['feature1', 'feature2', 'r']
    correls = list()
    for i in xrange(len(cor.index)):
        for j in xrange(i+1, len(cor.index)):
            correls.append([str(cor.index[i]), str(cor.index[j]), cor.iat[i, j]])

    if p_vals is not None:
        # add calculated p-value ot correls
        header.append('p-value')
        for i, correl in enumerate(correls):
            correls[i].append(p_vals.loc[correl[0],correl[1]])
    return correls, header


def adjust_pvals(p_adjust, correls, header):
    # adjust p-value if desired
    header.append('adjusted_p')
    p_adjusted = p_adjust([i[3] for i in correls])
    for i in xrange(len(correls)):
        correls[i].append(p_adjusted[i])
    return correls, header
