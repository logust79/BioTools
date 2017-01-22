# Useful Online Annotation Scripts for Variants, Genes and HPO
[![Build Status](https://travis-ci.org/logust79/logust79.svg?branch=master)](https://travis-ci.org/logust79/BioTools)
CommonFuncs.py has some common functions
Variants.py has some variant classes (automatically clean variant formats using CommonFuncs.py, help find online ExAC and kaviar [using http://exac.hms.harvard.edu/rest/] and sometimes cadd if given a cadd tsv file)
Genes.py has some gene classes (help find online ExAC pLI, mis-z, symbol, alias using mygene.org)
### CommonFuncs examples
`clean_variant()` is very useful in that it can trim redundant bases, and that it also can fill missing bases.
```
v='2-27665608-TC--'
print clean_variant(v)
```
### Variant classes examples
```
import sqlite3
import Variants
import json

db_conn=sqlite3.connect('irdc.db')
V = Variants.Variants(db_conn,['20-61523355-T-C','X-153694021-C-T'])
print json.dumps(V.exac, indent=4)
```
