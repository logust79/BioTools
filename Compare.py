'''
When new evidence comes in, highlight changes only.
Highlight methods are given as a dictionary as follows
{
    '+':{'style':{'background-color': '#a8f215'}},
    '-':{'style':{'background-color': '#aeaeae'}},
    '<>':{
        'style':{'background-color': '#f7f02c'},
        'data':[3,4,5],
    }
}
'''
from __future__ import print_function, division
import copy
import pandas
from numpy import NaN

'''
comare between df1(original) and df2(new). 
df1 and df2 are pandas dataframes
use key to tell if new(+) or obsolete(-).
key should be a number (index)
each ary's element[key] has to be unique, i.e., 
len(set([i[key] for i in ary])) == len(ary)
fields is a list of indexes where changes are checked
result = {
    '+': list of keys,
    '-': list of keys,
    '<>': {
        key1:{
            change:[indices],
        }...
    }
}
'''
def compare_dfs(df1, df2, key, fields):
    result = {
            '+':None,
            '-':None,
            '<>':{},
        }
    old_keys = set(df1[key])
    new_keys = set(df2[key])
    # len not right?
    bad = []
    if len(old_keys) != len(df1.index):
        bad.append('first')
    if len(new_keys) != len(df2.index):
        bad.append('second')
    if bad:
        msg = "% array's key column's elements are not unique to each other" % ' and '.join(bad)
        raise ValueError(msg)
    # carry on
    # old and delete?
    die_keys = old_keys - new_keys
    result['-'] = die_keys
    # new?
    born_keys = new_keys - old_keys
    result['+'] = born_keys
    
    # get common keys
    common_keys = old_keys - die_keys
    # detect changes.
    for k in common_keys:
        # NaN != NaN, so replace NaN with ''.
        # after subsetting, replacing NaN with None stops working!
        diff_cols = df1[df1[key]==k].iloc[0][fields].replace([NaN],'')!=df2[df2[key]==k].iloc[0][fields].replace([NaN],'')
        changes = [colname for colname, diff_cols in zip(fields, diff_cols.values) if diff_cols]
        if changes:
            result['<>'][k] = {
                    'change':changes,
            }
    return result

'''
similar to compare_dfs, but accepts methods when checking fields
it should return `True` for no change, and `False` for change
fields = {
    field1: method1,
    field2: method2,
    ...
}
if method == None, use '==' as default.
'''
def compare_dfs_with_methods(df1, df2, key, fields):
    pass
