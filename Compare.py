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
fields accept methods when checking fields
it should return `False` for no change, and `True` for change
fields = {
    field1: method1,
    field2: method2,
    ...
}
if method == None, use '!=' as default.
result = {
    '+': list of keys,
    '-': list of keys,
    '<>': {
        key1:{
            changes:{
                column1:{
                    from:
                    to:
                },
            },
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
        msg = "%s array's key column's elements have duplicates" % ' and '.join(bad)
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
    # easy_fields just use plain == as method
    easy_fields = [i for i,v in fields.items() if not v]
    hard_fields = [i for i in fields if i not in easy_fields]
    for k in common_keys:
        # get row
        d1 = df1[df1[key]==k].iloc[0]
        d2 = df2[df2[key]==k].iloc[0]
        # NaN != NaN, so replace NaN with ''.
        # after subsetting, replacing NaN with None stops working!
        # easy:
        easy1 = d1[easy_fields].replace([NaN],'')
        easy2 = d2[easy_fields].replace([NaN],'')
        same_cols = (easy1==easy2)
        diff_cols = [colname for colname, col in zip(easy_fields, same_cols.values) if not col]
        # hard:
        hard1 = d1[hard_fields]
        hard2 = d2[hard_fields]
        diff_cols += [fd for fd in hard_fields if fields[fd](hard1[fd], hard2[fd])]
        if diff_cols:
            result['<>'][k] = {
                    'change':{
                        fd:{
                            'from':d1.get(fd,None),
                            'to':d2.get(fd,None),
                        } for fd in diff_cols
                    },
            }
    return result

