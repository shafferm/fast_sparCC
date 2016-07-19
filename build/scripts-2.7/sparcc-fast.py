import argparse
from biom import load_table
from sparcc_fast.bootstrap_correlations import sparcc_correlations, calculate_sparcc_pvalues, calculate_sparcc_pvalues_multi
from sparcc_fast import utils
from operator import itemgetter
import networkx as nx

__author__ = 'shafferm'

## TODO: Add in filter table by number of observations

def main(args):
    print "reading in table"
    table = load_table(args.input)
    df = utils.biom_to_pandas(table)
    print "calculating correlations"
    cor = sparcc_correlations(df)

    if args.corr_only:
        correls, header = utils.df_to_correls(cor)
        if args.network is not None:
            print "making network"
            if args.edge_cut is None:
                net = utils.correls_to_net_cor(correls, min_cor=.5, metadata=utils.get_metadata_from_table(table))
            else:
                net = utils.correls_to_net_cor(correls, min_cor=args.edge_cut,
                                               metadata=utils.get_metadata_from_table(table))
            nx.write_gml(net, args.network)
    else:
        # calculate p-values
        if args.procs == 1:
            print "calculating p-values"
            p_values = calculate_sparcc_pvalues(df, cor, bootstraps=args.boots)
        else:
            print "calculating p-values: multiprocessed"
            p_values = calculate_sparcc_pvalues_multi(df, cor, bootstraps=args.boots, procs=args.procs)
        # convert to correls list
        correls, header = utils.df_to_correls(cor, p_values)
        # adjust p-value if desired
        if args.padjust == "FDR":
            correls, header = utils.adjust_pvals(utils.bh_adjust, correls, header)
        elif args.padjust == "bonferroni":
            correls, header = utils.adjust_pvals(utils.bonferroni_adjust, correls, header)
        # make network if desired
        if args.network is not None:
            print "making network"
            if args.edge_cut is None:
                net = utils.correls_to_net(correls, min_p=.05, metadata=utils.get_metadata_from_table(table))
            else:
                net = utils.correls_to_net(correls, min_p=args.edge_cut, metadata=utils.get_metadata_from_table(table))
            nx.write_gml(net, args.network)
    correls.sort(key=itemgetter(-1))
    utils.print_delimited(args.output, correls, header)


if __name__ == '__main__':
    """main, takes argparser"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="location of input biom file")
    parser.add_argument("-o", "--output", help="output file location", default="correls.txt")
    parser.add_argument("--corr_only", help="only calculate correlations, don't calculate p-values",
                        action="store_true", default=False)
    parser.add_argument("-p", "--procs", help="number of processors to use", type=int)
    parser.add_argument("-b", "--boots", help="number of bootstraps", type=int, default=100)
    parser.add_argument("--padjust", help="multiple testing corretion method: FDR, bonferroni or none", default="FDR")
    parser.add_argument("-n", "--network", help="network output file name, network will only be made if this argument"
                                                "is provided")
    parser.add_argument("--edge_cut", help="edge cut off for network formation, if p-value adjustment is done that will"
                                           "be used, if p-value calculation with no adjustment raw p-value will be"
                                           "used, if no p-vaules are calcualted correlation coefficent will be used",
                        type=float)
    main(parser.parse_args())
