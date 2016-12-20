import logging
import sqlite3
from CommonFuncs import *
from sqlite_utils import *
import copy
import json
import re

def _initiate_db(db_conn):
    db_c = db_conn.cursor()
    db_c.execute('''CREATE TABLE IF NOT EXISTS genes
        (id text NOT NULL UNIQUE, gene_id text, pLI real, mis_z real, genomic_pos_hg19 text, genomic_pos text, symbol text, alias text, PRIMARY KEY (id, gene_id))''')
    db_conn.commit()

def _update_db(db_conn, mgs):
    # a wrapper of sqlite_utils.update_db for genes
    fields = ['gene_id','pLI','mis_z','genomic_pos_hg19','genomic_pos','symbol','alias']
    # transform mgs to a dict
    data = {}
    for i in mgs:
        genes = []
        if isinstance(i['ensembl'], list):
            # sometimes ensembl returns a list, each element corresponds to a gene_id
            genes = [j['gene'] for j in i['ensembl']]
        else:
            genes = [i['ensembl']['gene']]
        for g in genes:
            data[g] = [
                i['_id'],
                i['exac']['all']['p_li'] if 'exac' in i and 'all' in i['exac'] else -1, #pLI
                i['exac']['all']['mis_z'] if 'exac' in i and 'all' in i['exac'] else -1, #mis_z
                json.dumps(i['genomic_pos_hg19'],indent=4), #genomic_pos_hg19
                json.dumps(i['genomic_pos'],indent=4), #genomic_pos
                i['symbol'],
                json.dumps(i.get('alias',[]),indent=4), #alias
            ]
    # update
    update_db(
            db_conn,
            'genes',
            fields,
            data
            )

def _fetch_one(self,field):
    db_c = self.db_conn.cursor()
    db_c.execute('SELECT * FROM genes WHERE id=?',(self.id,))
    db_gene = dict_factory(db_c,db_c.fetchone())
    if db_gene == None or db_gene[field] == None:
        # query mygene
        print 'query mygene'
        mg = my_gene(self.id)
        # update db
        _update_db(self.db_conn, [mg])
        # refetch
        db_c.execute('SELECT * FROM genes WHERE id=?',(self.id,))
        db_gene = dict_factory(db_c,db_c.fetchone())
    return db_gene[field]

class Gene(object):
    def __init__(self, db_conn, id=None):
        # id ok?
        if id[:4] != 'ENSG':
            raise ValueError("can't recognise gene id. It has to be an ensembl id!")
        self.id = id
        _initiate_db(db_conn)
        self.db_conn = db_conn
    
    @property
    def gene_id(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gene_id', None) is None:
            self._gene_id = _fetch_one(self,'gene_id')
        return self._gene_id

    @property
    def pLI(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_pLI', None) is None:
            self._pLI = _fetch_one(self,'pLI')
        return self._pLI

    @property
    def mis_z(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_mis_z', None) is None:
            self._mis_z = _fetch_one(self,'mis_z')
        return self._mis_z

    @property
    def genomic_pos_hg19(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp19', None) is None:
            self._gp19 = json.loads(_fetch_one(self,'genomic_pos_hg19'))
        return self._gp19

    @property
    def genomic_pos(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp', None) is None:
            self._gp = json.loads(_fetch_one(self,'genomic_pos'))
        return self._gp
    
    @property
    def symbol(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_symbol', None) is None:
            self._symbol = _fetch_one(self,'symbol')
        return self._symbol

    @property
    def alias(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_alias', None) is None:
            self._alias = json.loads(_fetch_one(self,'alias'))
        return self._alias

def _fetch_many(self,field):
    db_c = self.db_conn.cursor()
    result = batch_query(db_c,'genes',self.ids)
    data = {}
    new_genes = []
    final = {}
    for i in result:
        temp = dict_factory(db_c, i)
        data[temp['id']] = temp[field]
    for g in self.ids:
        if g in data and data[g] != None:
            final[g] = data[g]
        else:
            new_genes.append(g)
    if new_genes:
        print new_genes
        print 'querying mygenes'
        new_result = my_genes(new_genes)
        # update database
        _update_db(self.db_conn,new_result)
        # query again
        new_result = batch_query(db_c,'genes',new_genes)
        for i in new_result:
            temp = dict_factory(db_c, i)
            final[temp['id']] = temp[field]
    return final

class Genes(object):
    def __init__(self, db_conn, ids=[]):
        # id ok?
        for id in ids:
            if id[:4] != 'ENSG': raise ValueError('id has to be an Ensembl id, such as ENSG00000050453. (%s)' % id)
        _initiate_db(db_conn)
        self.db_conn = db_conn
        self.ids = ids
    
    def gene_ids_to_ids(self,gene_ids=[]):
        # convert from gene_ids to ids
        db_c = self.db_conn.cursor()
        db_result = batch_query(db_c,'genes',gene_ids,'gene_id')
        new_genes = []
        data = {}
        final = {}
        for i in db_result:
            temp = dict_factory(db_c,i)
            data[temp['gene_id']] = data.get(temp['gene_id'],[])
            data[temp['gene_id']].append(temp['id'])
        for g in gene_ids:
            if g in data and data[g] != None:
                final[g] = data[g]
            else:
                new_genes.append(g)
        if new_genes:
            print new_genes
            print 'querying mygenes'
            new_result = my_genes(new_genes)
            # update database
            _update_db(self.db_conn,new_result)
            # query again
            new_result = batch_query(db_c,'genes',new_genes,'gene_id')
            for i in new_result:
                temp = dict_factory(db_c, i)
                final[temp['gene_id']] = final.get(temp['gene_id'],[])
                final[temp['gene_id']].append(temp['id'])
        return final

    
    @property
    def gene_id(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gene_id', None) is None:
            self._gene_id = _fetch_many(self,'gene_id')
        return self._gene_id
    
    @property
    def pLI(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_pLI', None) is None:
            self._pLI = _fetch_many(self,'pLI')
        return self._pLI

    @property
    def mis_z(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_mis_z', None) is None:
            self._mis_z = _fetch_many(self,'mis_z')
        return self._mis_z

    @property
    def genomic_pos_hg19(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp19', None) is None:
            self._gp19 = {k:json.loads(v) for k,v in _fetch_many(self,'genomic_pos_hg19').iteritems()}
        return self._gp19

    @property
    def genomic_pos(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp', None) is None:
            self._gp = {k:json.loads(v) for k,v in _fetch_many(self,'genomic_pos').iteritems()}
        return self._gp
    
    @property
    def symbol(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_symbol', None) is None:
            self._symbol = _fetch_many(self,'symbol')
        return self._symbol

    @property
    def alias(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_alias', None) is None:
            self._alias = {k:json.loads(v) for k,v in _fetch_many(self,'alias').iteritems()}
        return self._alias
