# fast_sparCC
A fast command line interface to find correlations in biom tables with sparCC.

A way to use sparCC on biom formatted tables such as those output by QIIME. Outputs a tab delimited file with pairs of features, correlation values and optional p-values.
Includes options for filtering input tables, multiprocessing for bootstraping correlations values to determine significance and p-value adjustment with FDR or Bonferroni correction.

Example usage:

Calculate correlations only:
```
fast_sparCC.py -i example.biom -o correls.txt --corr_only
```

Calculate correlations only on table filtered based on Friedman and Alm 2012:
```
fast_sparCC.py -i example.biom -o correls.txt --corr_only --sparcc_filter
```

Calculate correlations and p-values based on 1000 bootstraps:
```
fast_sparCC.py -i example.biom -o correls.txt -b 1000
```

Calculate correlations and p-values based on 1000 bootstraps and 10 processors:
```
fast_sparCC.py -i example.biom -o correls.txt -b 1000 --procs 10
```
