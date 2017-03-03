# some useful utilities to query/update sqlite db
import sqlite3

def _chop_array(arr, size=999):
    for i in range(0, len(arr), size):
        yield arr[i:i + size]

def dict_factory(cursor, row):
    # convert from sqlite tuple to dictionary
    # if row is None, return None
    if row == None: return None
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def update_db(db_conn,table,fields,value_dict):
    # update database
    # first get all columns except id
    db_c = db_conn.cursor()
    sql = 'PRAGMA table_info(%s)' % table
    db_c.execute(sql)
    columns = [i[1] for i in db_c.fetchall()]
    remain_fields = list(set(columns)-set(fields+['id']))
    # constructing sql
    #sql = 'INSERT OR REPLACE INTO %(table)s (id, %(field)s, ' % {'table':table,'field':field}
    sql = 'INSERT OR REPLACE INTO %s (id, ' % table
    sql += ', '.join(fields + remain_fields)
    sql += ') VALUES ( '
    sql += '?, '*len(fields)
    sql += '?'
    if remain_fields:
        sql += ',' + ','.join(['(SELECT %s FROM %s WHERE id = ?)' % (f, table) for f in remain_fields])
    sql += ')'
    # update using executemany.
    # might take a long time for large array, chop it into small chunks to commit
    values_to_insert = []
    for k,v in value_dict.items():
        values_to_insert.append((k,)+tuple(v)+tuple([k]*len(remain_fields)))
        #db_c.execute(sql, (k,)+tuple(v)+tuple([k]*len(remain_fields)))
    for a in _chop_array(values_to_insert,200):
        db_c.executemany(sql,a)
        db_conn.commit()

def batch_query(db_c,table,arr,pointer='id'):
    # batch query sqlite database
    # since there is a limit on number of variables at 999, chop them
    result = []
    for a in _chop_array(arr):
        sql = 'SELECT * FROM %s WHERE %s in (%s)' % (
                        table,
                        pointer,
                        ','.join(['?']*len(a)))
        temp = db_c.execute(sql,a)
        result.extend([i for i in temp])
    return result
