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
import gzip

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
        self._check_db()
        raise ValueError('no %s can be retrieved' % field)
    return db_hpo[field]

class Hpo:
    def __init__(self,id,db_conn):
        _initiate_db(db_conn)
        self.id = id
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
        self.db_conn = db_conn
    
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
            raise ValueError('No data seems to be in the hpo database. Did you forget to load it using HPO.construct_db("hp.obo")?')
    def construct_db(self,obofile,hpo_gene_file):
        # construct hpo database using the obo file
        data = obo_parser(obofile)
        # get ancestors
        for k,v in data.items():
            data[k]['ancestors'] = self._find_ancestors(k,[],data)
        
        # get hpo_gene
        hpo_gene = defaultdict(list)
        if hpo_gene_file.endswith('gz'):
            f = gzip.open(hpo_gene_file,'rb')
        else:
            f = open(hpo_gene_file,'r')
        for l in f:
            if l[0] == '#': continue
            l = l.rstrip().split('\t')
            hpo_gene[l[0]].append(l[2])
        #convert hpo_gene's gene_ids to ensembl ids
        G = Genes.Genes(self.db_conn)
        for k,v in hpo_gene.items():
            hpo_gene[k] = _flatten_array_of_arrays(G.entrezIds_to_ensemblIds(v).values())
        # convert to array of tuples
        values = []
        for k,v in data.items():
            values.append((
                k,
                v['name'],
                json.dumps(v['alt_id']),
                json.dumps(v['is_a']),
                json.dumps(v['ancestors']),
                json.dumps(hpo_gene[k])
            ))
        # write to database
        db_c = self.db_conn.cursor()
        sql = 'INSERT INTO hpo VALUES (?,?,?,?,?,?)'
        db_c.executemany(sql,values)
        self.db_conn.commit()
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
