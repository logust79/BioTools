# Useful Online Annotation Scripts for Variants, Genes and HPO
CommonFuncs.py has some common functions
Variants.py has some variant classes (automatically clean variant formats using CommonFuncs.py, help find online ExAC and kaviar [using http://exac.hms.harvard.edu/rest/] and sometimes cadd if given a cadd tsv file)
Genes.py has some gene classes (help find online ExAC pLI, mis-z, symbol, alias using mygene.org)

### Variant classes examples
```
import sqlite3
import Variants
import json

db_conn=sqlite3.connect('irdc.db')
V = Variants.Variants(db_conn,['20-61523355-T-C','X-153694021-C-T'])
print json.dumps(V.exac, indent=4)
```
