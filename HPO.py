'''
hpo objects, including doing analyses such as hpo similarity, and group of hpos similarities
hp.obo can be downloaded from: http://www.obofoundry.org/ontology/hp.html
ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt.gz can be downloaded from https://github.com/moonso/phizz/blob/master/phizz/resources/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt.gz
'''
import requests
import json
import os
import time
import Genes
from CommonFuncs import *
from collections import defaultdict
from sqlite_utils import *
import pandas # for reading csv and dumping to sqlite

def _initiate_db(db_conn):
    db_c = db_conn.cursor()
    db_c.execute('''CREATE TABLE IF NOT EXISTS hpo
        (id text PRIMARY KEY UNIQUE, name text, alt_id text, parents text, ancestors text, genes text)''')
    db_conn.commit()

def _flatten_array_of_arrays(arrs):
    # flatten array of arrays
    result = []
    [result.extend(i) for i in arrs]
    return result

def _fetch_one(self,field):
    db_c = self.db_conn.cursor()
    db_c.execute('SELECT * FROM hpo WHERE id=?',(self._id,))
    db_hpo = dict_factory(db_c,db_c.fetchone())
    if db_hpo == None or db_hpo[field] == None:
        # first check if hpo db has been constructed
        raise ValueError('no %s can be retrieved' % field)
    return db_hpo[field]

class Hpo:
    def __init__(self,db_conn,id=None):
        _initiate_db(db_conn)
        self.db_conn = db_conn
        self._check_db()
        self.id = id
        if id:
            # alt id?
            sql = 'SELECT id FROM hpo WHERE id=?'
            db_c = db_conn.cursor()
            db_c.execute(sql,(id,))
            data = db_c.fetchone()
            if data:
                self._id = self.id
            else:
                #look into alt_id
                sql = "SELECT id FROM hpo WHERE alt_id LIKE ?"
                db_c.execute(sql,('%'+id+'%',))
                data = db_c.fetchone()
                if not data:
                    msg = '%s cannot be recognised' % id
                    raise ValueError(msg)
                self._id = data[0]
    
    def _find_ancestors(self,id,anc,data):
        # construct_db helper function, to find all ancestors of a node
        is_a = data[id]['is_a']
        if not is_a:
            return anc
        else:
            anc.extend(is_a)
            for i in is_a:
                return self._find_ancestors(i,anc,data)

    def _check_db(self):
        # hpo has been constructed? if not, raise error. if yes, do nothing
        db_c = self.db_conn.cursor()
        db_c.execute('SELECT * FROM hpo WHERE id=?',('HP:0000001',))
        db_hpo = dict_factory(db_c,db_c.fetchone())
        if db_hpo == None:
            self.construct_db()
    def construct_db(self):
        # construct hpo database using the data/hpo.csv file
        csvfile = os.path.join(os.path.dirname(__file__),'data','hpo.csv')
        df = pandas.read_csv(csvfile)
        df.to_sql('hpo', self.db_conn, if_exists='replace', index=False)
    
    @property
    def name(self):
        if getattr(self,'_name',None) is None:
            self._name = _fetch_one(self,'name')
        return self._name

    @property
    def alt_ids(self):
        if getattr(self,'_alt_ids',None) is None:
            self._alt_ids = json.loads(_fetch_one(self,'alt_id'))
        return self._alt_ids

    @property
    def parents(self):
        if getattr(self,'_parents',None) is None:
            self._parents = json.loads(_fetch_one(self,'parents'))
        return self._parents

    @property
    def ancestors(self):
        if getattr(self,'_ancestors',None) is None:
            self._ancestors = json.loads(_fetch_one(self,'ancestors'))
        return self._ancestors

    @property
    def genes(self):
        if getattr(self,'_genes',None) is None:
            self._genes = json.loads(_fetch_one(self,'genes'))
        return self._genes

    def names_to_ids(self, names):
        # translate an array of names to a dictionary of name => keys
        db_c = self.db_conn.cursor()
        result = batch_query(db_c,'HPO',names,pointer='name')
        data = {}
        for i in result:
            temp = dict_factory(db_c, i)
            data[temp['name']] = temp['id']
        return data
