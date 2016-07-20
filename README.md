# fast_sparCC
A fast command line interface to find correlations in biom tables with sparCC.

A way to use sparCC on biom formatted tables such as those output by QIIME. Outputs a tab delimited file with pairs of features, correlation values and optional p-values.
Includes options for filtering input tables, multiprocessing for bootstraping correlations values to determine significance and p-value adjustment with FDR or Bonferroni correction.

##Installation Instructions
Install fast_sparCC by using navigating to the folder of your choice and using these commands:
```
git clone https://github.com/shafferm/fast_sparCC.git
cd fast_sparCC
python setup.py install
```

##Example usage:
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

##API for use in custom python scripts:
Additionally the functions can be imported for use in your own python scripts. This can be used with and without
bootstrapping to calculate p-values.

Example usage of sparcc correlations:
```
#!/usr/local/bin/python2

from biom import load_table
from sparcc_fast import sparcc_correlation
table = load_table("example.biom")
correls = sparcc_correlation(table)
correls.to_csv("correls.txt")
```

Example usage of sparcc correlations with bootstrapping:
```
#!/usr/local/bin/python2

from biom import load_table
from sparcc_fast import sparcc_correlation_w_bootstraps
table = load_table("example.biom")
correls = sparcc_correlation_with_bootstraps(table, procs=3, bootstraps=1000)
correls.to_csv("correls.txt")
```
