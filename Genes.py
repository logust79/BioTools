import logging
import sqlite3
from CommonFuncs import *
from sqlite_utils import *
import copy
import json

def _initiate_db(db_conn):
    db_c = db_conn.cursor()
    db_c.execute('''CREATE TABLE IF NOT EXISTS genes
        (id text PRIMARY KEY, pLI real, mis_z real, genomic_pos_hg19 text, genomic_pos text, symbol text, alias text)''')
    db_conn.commit()

def _update_db(db_conn, mgs):
    # a wrapper of sqlite_utils.update_db for genes
    fields = ['pLI','mis_z','genomic_pos_hg19','genomic_pos','symbol','alias']
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
                i['exac']['all']['p_li'] if 'exac' in i and 'all' in i['exac'] else -1, #pLI
                i['exac']['all']['mis_z'] if 'exac' in i and 'all' in i['exac'] else -1, #mis_z
                json.dumps(i['genomic_pos_hg19'],indent=4), #genomic_pos_hg19
                json.dumps(i['genomic_pos'],indent=4), #genomic_pos
                i['symbol'],
                json.dumps(i.get('alias',[]),indent=4), #alias
            ]
    # update
    print mgs
    update_db(
            db_conn,
            'genes',
            fields,
            data
            )

class Gene(object):
    def __init__(self, db_conn, gene_id=None):
        # gene_id ok?
        if gene_id[:4] != 'ENSG': raise ValueError('gene_id has to be an Ensembl id, such as ENSG00000050453')
        _initiate_db(db_conn)
        self.db_conn = db_conn
        self.gene_id = gene_id

    @property
    def pLI(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_pLI', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['pLI'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['exac']['all']['p_li']
            else:
                temp = db_gene['pLI']
            self._pLI = temp
        return self._pLI

    @property
    def mis_z(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_mis_z', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['mis_z'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['exac']['all']['mis_z']
            else:
                temp = db_gene['mis_z']
            self._mis_z = temp
        return self._mis_z

    @property
    def genomic_pos_hg19(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp19', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['genomic_pos_hg19'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['genomic_pos_hg19']
            else:
                temp = json.loads(db_gene['genomic_pos_hg19'])
            self._gp19 = temp
        return self._gp19

    @property
    def genomic_pos(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['genomic_pos'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['genomic_pos']
            else:
                temp = json.loads(db_gene['genomic_pos'])
            self._gp = temp
        return self._gp
    
    @property
    def symbol(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_symbol', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['symbol'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['symbol']
            else:
                temp = db_gene['symbol']
            self._symbol = temp
        return self._symbol

    @property
    def alias(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_alias', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM genes WHERE id=?',(self.gene_id,))
            db_gene = dict_factory(db_c,db_c.fetchone())
            if db_gene == None or db_gene['alias'] == None:
                # query mygene
                mg = my_gene(self.gene_id)
                # update db
                _update_db(self.db_conn, [mg])
                temp = mg['alias']
            else:
                temp = json.loads(db_gene['alias'])
            self._alias = temp
        return self._alias

class Genes(object):
    def __init__(self, db_conn, gene_ids=[]):
        # gene_id ok?
        for gene_id in gene_ids:
            if gene_id[:4] != 'ENSG': raise ValueError('gene_id has to be an Ensembl id, such as ENSG00000050453. (%s)' % gene_id)
        _initiate_db(db_conn)
        self.db_conn = db_conn
        self.gene_ids = gene_ids

    @property
    def pLI(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_pLI', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            pLI = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = temp['pLI']
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    pLI[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        pLI[g] = i['exac']['all']['p_li'] if 'exac' in i and 'all' in i['exac'] else -1
            self._pLI = pLI
        return self._pLI

    @property
    def mis_z(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_mis_z', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            mis_z = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = temp['mis_z']
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    pLI[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        mis_z[g] = i['exac']['all']['mis_z'] if 'exac' in i and 'all' in i['exac'] else -1
            self._mis_z = mis_z
        return self._mis_z

    @property
    def genomic_pos_hg19(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp19', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            gp19 = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = json.loads(temp['genomic_pos_hg19'])
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    gp19[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        gp19[g] = i['genomic_pos_hg19']
            self._gp19 = gp19
        return self._gp19

    @property
    def genomic_pos(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_gp', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            gp = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = json.loads(temp['genomic_pos'])
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    gp[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        gp[g] = i['genomic_pos']
            self._gp = gp
        return self._gp
    
    @property
    def symbol(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_symbol', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            symbol = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = temp['symbol']
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    symbol[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        symbol[g] = i['symbol']
            self._symbol = symbol
        return self._symbol

    @property
    def alias(self):
        # check local database first. if na, use CommonFuncs to annotate, then store in db
        if getattr(self, '_alias', None) is None:
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'genes',self.gene_ids)
            data = {}
            new_genes = []
            alias = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = json.loads(temp['alias'])
            for g in self.gene_ids:
                if g in data and data[g] != None:
                    alias[g] = data[g]
                else:
                    new_genes.append(g)
            if new_genes:
                print 'querying web'
                new_result = my_genes(new_genes)
                # update database
                _update_db(self.db_conn,new_result)
                for i in new_result:
                    genes = []
                    if isinstance(i['ensembl'], list):
                        genes = [j['gene'] for j in i['ensembl']]
                    else:
                        genes = [i['ensembl']['gene']]
                    for g in genes:
                        alias[g] = i.get('alias',[])
            self._alias = alias
        return self._alias
