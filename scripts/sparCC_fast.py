#!/usr/local/bin/python2

import argparse
from biom import load_table
from sparcc_fast.sparcc_functions import basis_corr
from sparcc_fast import utils
from sparcc_fast.bootstrap_correlations import bootstrap_correlations

__author__ = 'shafferm'


# TODO: implement sparse pandas dataframe


def sparcc_paper_filter(table):
    """if a observation averages more than 2 reads per sample then keep,
    if a sample has more than 500 reads then keep"""
    table = table[table.sum(axis=1) > 500]
    table = table.loc[:, table.mean(axis=0) > 2]
    return table


def min_sample_filter(table, min_samples):
    """remove observations not present in more than a minimum number of samples"""
    zeroes_per_column = (table > 0).sum(axis=0)
    return table.loc[:, zeroes_per_column > min_samples]


def main(args):
    print "reading in table"
    table = load_table(args.input)
    df = utils.biom_to_pandas(table)

    # filter
    if args.min_sample is not None:
        table = min_sample_filter(df, args.min_sample)
        print "Table filtered: " + str(table.shape[1]) + " observations"
        print ""
    elif args.sparcc_filter is True:
        table = sparcc_paper_filter(df)
        print "Table filtered: " + str(table.shape[1]) + " observations"
        print ""

    print "calculating correlations"
    cor, cov = basis_corr(df)
    correls = utils.df_to_correls(cor)

    if not args.corr_only:
        correls['p_value'] = bootstrap_correlations(df, cor, args.boots, args.procs)
        # adjust p-value if desired
        if args.padjust == "FDR":
            correls['p_adjusted'] = utils.bh_adjust(correls['p_value'])
        elif args.padjust == "bonferroni":
            correls, header = utils.bonferroni_adjust(correls['p_value'])
    correls.to_csv(args.output, sep='\t', index=False)


if __name__ == '__main__':
    """main, takes argparser"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of input biom file")
    parser.add_argument("-o", "--output", help="output file location", default="correls.txt")
    parser.add_argument("-p", "--procs", help="number of processors to use, only matters if calculating p-values",
                        type=int)
    parser.add_argument("-b", "--boots", help="number of bootstraps", type=int, default=100)
    parser.add_argument("--corr_only", help="only calculate correlations, don't calculate p-values",
                        action="store_true", default=False)
    parser.add_argument("--p_adjust", help="multiple testing corretion method: FDR, bonferroni or none", default="FDR")
    parser.add_argument("--sparcc_filter", help="filter input table according to parameters defined in Inferring"
                                                "Correlation Networks from Genomic Survey Data",
                        default=False,action="store_true")
    parser.add_argument("--min_samples", help="minimum number of samples a observation must be present in to be kept"
                                              "in anaylsis",
                        type=int)
    main(parser.parse_args())
