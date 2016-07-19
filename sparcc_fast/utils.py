import numpy as np
import pandas as pd


__author__ = 'shafferm'


def bh_adjust(pvalues):
    """benjamini-hochberg p-value adjustment stolen from
    http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
    """
    pvalues = np.array(pvalues)
    n = pvalues.shape[0]
    new_pvalues = np.empty(n)
    values = [(pvalue, i) for i, pvalue in enumerate(pvalues)]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in xrange(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues


def bonferroni_adjust(pvalues):
    pvalues = np.array(pvalues)
    n = float(pvalues.shape[0])
    new_pvalues = n * pvalues
    return new_pvalues


def biom_to_pandas(table):
    # convert to pandas dataframe
    return pd.DataFrame(np.transpose(table.matrix_data.todense()), index=table.ids(), columns=table.ids(axis="observation"))


def df_to_correls(cor):
    header = ['feature1', 'feature2', 'r']
    correls = list()
    for i in xrange(len(cor.index)):
        for j in xrange(i+1, len(cor.index)):
            correls.append([str(cor.index[i]), str(cor.index[j]), cor.iat[i, j]])
    return pd.DataFrame(correls, columns=header)
