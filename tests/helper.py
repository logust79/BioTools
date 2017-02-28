'''
helper functions for tests
'''
from __future__ import print_function, division
import pandas as pd
import numpy as np
from numpy import NaN
import sys
sys.path.append('..')
from Compare import compare_dfs
import os

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
read in excel file sheets,
parse with header
'''
def compare_excel(file1, file2, sheet, key, fields):
    f1 = pd.ExcelFile(file1)
    f2 = pd.ExcelFile(file2)
    df1 = f1.parse(sheet)
    df2 = f2.parse(sheet)
    result = compare_dfs(df1,df2,key,fields)
    return result
