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
        (id text NOT NULL UNIQUE, entrez_id text, pLI real, mis_z real, genomic_pos_hg19 text, genomic_pos text, symbol text, alias text, PRIMARY KEY (id, entrez_id))''')
    db_conn.commit()

def _update_db(self, mgs):
    # a wrapper of sqlite_utils.update_db for genes
    fields = ['entrez_id','pLI','mis_z','genomic_pos_hg19','genomic_pos','symbol','alias']
    # transform mgs to a dict
    good_result = []
    # remove not found
    bad_genes = []
    for i in mgs:
        if i.get('notfound',None):
            bad_genes.append(i['query'])
        else:
            good_result.append(i)
    if bad_genes:
        self._bad_genes.extend(bad_genes)
        print '-----some queries are not found----'
        print json.dumps(bad_genes)
    data = {}
    for i in good_result:
        gene = None
        genomic_pos = None
        # some genes miss ensembl ids, fill them manually for the time being
        if '_id' not in i: continue #not found
        if i['_id'] == '7012':
            i['ensembl'] = {'gene':'ENSG00000270141'}
        elif i['_id'] == '6315':
            i['ensembl'] = {'gene':'ENSG00000230223'}
        elif i['_id'] == '2657':
            i['ensembl'] = {'gene':'ENSG00000130283'}
        elif i['_id'] == '84876':
            i['ensembl'] = {'gene':'ENSG00000276045'}
        elif i['_id'] == '9103':
            i['ensembl'] = {'gene':'ENSG00000244682'}
        elif 'genomic_pos_hg19' not in i:
            i['genomic_pos_hg19'] = {}
        if 'ensembl' not in i:
            self._bad_genes.append(i['query'])
            logging.warning('Warning: %s is not registered in ensembl' % i['query'])
            continue
        '''
        debug
        '''
        if 'ensembl' not in i or 'genomic_pos' not in i:
            print json.dumps(i,indent=4)
        if isinstance(i['ensembl'], list):
            # sometimes ensembl returns a list, each element corresponds to an id
            # check which is the active ensembl id
            # genomic_pos has only one in valid chromosomes
            print 'use ensembl API to check ensemblid'
            gene = [j for j in i['ensembl'] if check_ensemblId(j['gene'])][0]['gene']
            for val in i['genomic_pos']:
                if val['chr'] in VALID_CHROMOSOMES:
                    genomic_pos = val
                    break
        else:
            gene = i['ensembl']['gene']
            genomic_pos = i['genomic_pos']

        data[gene] = [
            i['_id'],
            i['exac']['all']['p_li'] if 'exac' in i and 'all' in i['exac'] else -1, #pLI
            i['exac']['all']['mis_z'] if 'exac' in i and 'all' in i['exac'] else -1, #mis_z
            json.dumps(i['genomic_pos_hg19'],indent=4), #genomic_pos_hg19
            json.dumps(genomic_pos,indent=4), #genomic_pos
            i['symbol'],
            json.dumps(i.get('alias',[]),indent=4), #alias
        ]
    # update
    update_db(
            self.db_conn,
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
        _update_db(self, [mg])
        # refetch
        db_c.execute('SELECT * FROM genes WHERE id=?',(self.id,))
        db_gene = dict_factory(db_c,db_c.fetchone())
    return db_gene[field]

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
        elif g not in self._bad_genes:
            new_genes.append(g)
    if new_genes:
        print 'querying mygenes from fetch_many'
        new_result = my_genes(new_genes)
        # update database
        _update_db(self,new_result)
        # query again
        new_result = batch_query(db_c,'genes',new_genes)
        for i in new_result:
            temp = dict_factory(db_c, i)
            final[temp['id']] = temp[field]
    return final

class Gene(object):
    def __init__(self, db_conn, id=None):
        # id ok?
        if id[:4] != 'ENSG':
            raise ValueError("can't recognise gene id. It has to be an ensembl id!")
        self.id = id
        _initiate_db(db_conn)
        self.db_conn = db_conn
    
    @property
    def entrez_id(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_entrez_id', None) is None:
            self._entrez_id = _fetch_one(self,'entrez_id')
        return self._entrez_id

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

class Genes(object):
    def __init__(self, db_conn, ids=[]):
        # id ok?
        for id in ids:
            if id[:4] != 'ENSG': raise ValueError('id has to be an Ensembl id, such as ENSG00000050453. (%s)' % id)
        _initiate_db(db_conn)
        self.db_conn = db_conn
        self.ids = ids
        self._bad_genes = [] # this is for bad ids that have no ordinary entries in mygenes. Avoid repetitive queries.
    
    def entrezIds_to_ensemblIds(self,entrez_ids=[]):
        # convert from entrez ids to ensembl ids
        db_c = self.db_conn.cursor()
        entrez_ids = list(set(entrez_ids) - set(self._bad_genes))
        db_result = batch_query(db_c,'genes',entrez_ids,'entrez_id')
        new_genes = []
        data = {}
        final = {}
        for i in db_result:
            temp = dict_factory(db_c,i)
            data[temp['entrez_id']] = data.get(temp['entrez_id'],[])
            data[temp['entrez_id']].append(temp['id'])
        for g in entrez_ids:
            if g in data and data[g] != None:
                final[g] = data[g]
            else:
                new_genes.append(g)
        if new_genes:
            print 'querying mygenes from entrezIds_to_ensemblIds'
            new_result = my_genes(new_genes)
            # update database
            _update_db(self,new_result)
            # query again
            new_result = batch_query(db_c,'genes',new_genes,'entrez_id')
            for i in new_result:
                temp = dict_factory(db_c, i)
                final[temp['entrez_id']] = final.get(temp['entrez_id'],[])
                final[temp['entrez_id']].append(temp['id'])
        return final

    def symbols_to_ensemblIds(self,symbols=[]):
        # convert from symbols to ensembl ids
        db_c = self.db_conn.cursor()
        # remove bad symbols
        symbols = list(set(symbols) - set(self._bad_genes))
        # seek symbols
        db_result = batch_query(db_c,'genes',symbols,'symbol')
        new_genes = []
        data = {}
        final = {}
        for i in db_result:
            temp = dict_factory(db_c,i)
            data[temp['symbol']] = temp['id']
        for g in symbols:
            if g in data and data[g] != None:
                final[g] = data[g]
            else:
                new_genes.append(g)
        # seek aliases
        sql = '''SELECT * FROM genes WHERE alias like ? '''
        found = []
        for g in new_genes:
            temp = [j for j in db_c.execute(sql,('%"'+g+'"%',))]
            if temp:
                final[g] = temp[0][0]
                found.append(g)
        for g in found:
            new_genes.remove(g)
        if new_genes:
            print 'querying mygenes for symbols to ensemblIds'
            new_result = my_genes_by_symbol(new_genes,species='human')
            # update database
            _update_db(self,new_result)
            # query again
            new_result = batch_query(db_c,'genes',new_genes,'symbol')
            for i in new_result:
                temp = dict_factory(db_c, i)
                final[i] = temp['id']
        return final

    @property
    def entrez_id(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_entrez_id', None) is None:
            self._entrez_id = _fetch_many(self,'entrez_id')
        return self._entrez_id
    
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
