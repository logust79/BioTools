# some useful utilities to query/update sqlite db
import sqlite3

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
    sql += ','.join(['(SELECT %s FROM %s WHERE id = ?)' % (f, table) for f in remain_fields])
    sql += ')'
    # update
    print sql
    for k,v in value_dict.iteritems():
        db_c.execute(sql, (k,)+tuple(v)+tuple([k]*len(remain_fields)))
    db_conn.commit()

def batch_query(db_c,table,arr):
    # batch query sqlite database
    sql = 'SELECT * FROM %s WHERE id in (%s)' % (
                    table,
                    ','.join(['?']*len(arr)))
    result = db_c.execute(sql,arr)
    return result
