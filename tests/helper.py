'''
helper functions for tests
'''
from __future__ import print_function, division
import pandas as pd
import numpy as np
from numpy import NaN
import sys
sys.path.append('..')
from Compare import *
import os
import errno    

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise
'''
parse exac info
'''
def parse_exac(exac, hom=False):
    if 'allele_freq' in exac['variant']:
        # is covered and exists on exac
        if not hom:
            return exac['variant']['allele_freq']
        else:
            return float(exac['variant']['hom_count']) / (sum(exac['variant']['pop_ans'].values()) * 2 or 1)
    if exac['base_coverage'] and exac['base_coverage'][0]['has_coverage']:
        # is covered but not exsits on exac
        return 0
    # not covered
    return None

'''
some customised method for comparing dfs
'''
def field4_cb(a,b):
    dc = {
        'PASS':3,
        'VQSR':2,
        'FAIL':1,
    }
    if dc[b] > dc[a]: return True
    else: return False

def field5_cb_factory(thrd):
    def field5_cb_inner(a,b):
        # convert np.nan to None
        # since np.nan != np.nan
        a = -1 if pd.isnull(a) else a
        b = -1 if pd.isnull(b) else b
        # equal?
        if a == b: return False
        # one > thrd, one < thrd?
        if sorted([a,b,thrd])[1] == thrd: return True
        # both > thrd?
        if min(a,b,thrd) == thrd: return False
        else: return True
    
    return field5_cb_inner

'''
read in excel file sheets,
parse with header
'''
import time
def compare_excel(file1, file2, sheet, key, fields):
    f1 = pd.ExcelFile(file1)
    f2 = pd.ExcelFile(file2)
    df1 = f1.parse(sheet)
    df2 = f2.parse(sheet)
    t1 = time.time()
    result = compare_dfs(df1,df2,key,fields)
    t2 = time.time()
    print('compare 1 took--')
    print(t2-t1)

    return result
