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
    # NaN != NaN, so replace NaN with None
    diff_cols = df1.irow(1)[fields].replace([NaN],None)!=df2.irow(1)[fields].replace([NaN],None)
    cols = [colname for colname, diff_cols in zip(fields, diff_cols.values) if diff_cols]
    print(cols)
    print(set(df1[key]))
    print(df1[df1[key].isin([1,2])])
