[![Build Status](https://travis-ci.org/logust79/BioTools.svg?branch=master)](https://travis-ci.org/logust79/BioTools)
[![Coverage Status](https://coveralls.io/repos/github/logust79/BioTools/badge.svg?branch=master)](https://coveralls.io/github/logust79/BioTools?branch=master)
# Useful Online Annotation Scripts for Variants, Genes and HPO
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
```python
import sqlite3
import Variants
import json

db_conn=sqlite3.connect('irdc.db')
V = Variants.Variants(db_conn,['20-61523355-T-C','X-153694021-C-T'])
print json.dumps(V.exac, indent=4)
```
# Compare dataframes
In the field of genetic diagosis, end results are refreshed when new evidence is introduced. Manual inspection necessitates highlight of changes to avoid wasting time on unchanged data.
```python
import Compare
import pandas as pd

fields_to_check = {
    'field1_to_check':None,
    'field2_to_check':None,
}
Compare.compare_dfs('original_df', 'new_df', 'index', fields_to_check)
```
```
{
    '+':['index1','index3'],
    '-':['index2','index6'],
    '<>':{
        'index4':{
            'change':{
                'field1_to_check':{
                    'from':'foo',
                    'to':'bar'
                }
            }
        }
    }
}
```
You can also pass your customised compare functions as values of fields for alternative comparisons
```python
import Compare
import pandas as pd

def filter_cb(a,b):
    dc = {
        'PASS':3,
        'VQSR':2,
        'FAIL':1,
    }
    if dc[b] > dc[a]: return True
    else: return False

fields_to_check = {
    'FILTER':filter_cb,
    'field2_to_check':None,
}
Compare.compare_dfs('original_df', 'new_df', 'index', fields_to_check)
```
Or you can use closures to create some even more customised comparison methods (useful when you want to variate some cutoffs.

```python
import Compare
import pandas as pd

def exac_cb_factory(thrd):
    def exac_cb_inner(a,b):
        # convert np.nan to None
        # since np.nan != np.nan
        a = None if pd.isnull(a) else a
        b = None if pd.isnull(b) else b
        # equal?
        if a == b: return False
        # one > thrd, one < thrd?
        if sorted([a,b,thrd])[1] == thrd: return True
        # both > thrd?
        if min(a,b,thrd) == thrd: return False
        else: return True
    
    return exac_cb_inner

exac_cutoff = 0.01
exac_compare = exac_cb_factory(exac_cutoff)
fields_to_check = {
    'EXAC_AF':exac_compare,
    'field2_to_check':None,
}
Compare.compare_dfs('original_df', 'new_df', 'index', fields_to_check)
```
When `None` is given as field value, it uses good old `==` for comparison. It takes care of the case where `np.nan != np.nan`
