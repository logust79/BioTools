import logging
import sqlite3
from CommonFuncs import *
from sqlite_utils import *
import gnomad_utils
import kaviar_utils
import copy

def _initiate_db(db_conn):
    db_c = db_conn.cursor()
    db_c.execute('''CREATE TABLE IF NOT EXISTS variants
        (id text PRIMARY KEY UNIQUE, exac text, gnomad text, kaviar_af real, cadd_phred real)''')
    db_conn.commit()

def _liftover(variant, frm, to):
    # a wrap over CommonFuncs liftover, add some warnings
    lo = liftover(variant,frm,to)
    if len(lo) != 1:
        # no result or more than 1 match. Very rare though
        logging.warning('Warning: %s liftover from %s to hg19 returned result length not equal to one' % (self.cleaned_variant_id, build))
    return lo[0]

def _read_cadd(infile):
    # read cadd_file and populate it as dict
    cadd_in = open(infile,'r')
    cadd_data = {}
    for line in cadd_in:
        if line[0] == '#': continue
        line = line.rstrip().split('\t')
        variant_id = '-'.join(line[:4])
        cadd_data[variant_id] = float(line[-1])
    return cadd_data

'''
result = {
    'v1':{
        'gnomad_af':
        'gnomad_ac':
        'gnomad_hom_af':
        'gnomad_hom_ac':
        'gnomad_hemi_ac':
    }
}
'''
def anno_gnomad(vars,path_to_gnomad):
    result = {}
    null = {
        'gnomad_af': None,
        'gnomad_ac': None,
        'gnomad_hom_af': None,
        'gnomad_hom_ac': None,
        'gnomad_hemi_ac': None,
    }
    for v in vars:
        
        if v.split('-')[0] not in VALID_CHROMOSOMES:
            result[v] = null
            continue
        covs = {
            'exome': gnomad_coverage(v,path_to_gnomad,mode='exome'),
            'genome':gnomad_coverage(v,path_to_gnomad,mode='genome'),
        }
        freqs = {
            'exome': gnomad_freqs(v,path_to_gnomad,mode='exome'),
            'genome':gnomad_freqs(v,path_to_gnomad,mode='genome'),
        }
        if not freqs['exome'] and not freqs['genome'] and not covs['exome'] and not covs['genome']:
            result[v] = null
            continue
        ac = hom_ac = af = hom_af = an = 0.
        hemi_ac = None
        for m in ['exome', 'genome']:
            if freqs[m]:
                ac += freqs[m]['AC']
                hom_ac += freqs[m]['Hom']
                an += freqs[m]['AN']
                if 'Hemi' in freqs[m]:
                    hemi_ac = hemi_ac + freqs[m]['Hemi'] if hemi_ac != None else freqs[m]['Hemi']
        if ac: af = ac / an
        if hom_ac: hom_af = hom_ac * 2 / an
        result[v] = {
            'gnomad_af':af,
            'gnomad_ac':ac,
            'gnomad_hom_af':hom_af,
            'gnomad_hom_ac':hom_ac,
            'gnomad_hemi_ac':hemi_ac,
            'gnomad_an':an,
        }
    return result

'''
class for single variant
'''
class Variant(object):
    def __init__(self, db_conn, variant_id=None, build='hg19'):
        '''
        attributes:
            variant_id
            db_conn
            build
        '''
        # initiate db
        _initiate_db(db_conn)
        self.variant_id = variant_id if variant_id else logging.warning('Need variant_id')
        self.cleaned_variant_id = find_leftmost_synonymous_variant(clean_variant(variant_id))
        self.build = build
        self.db_conn = db_conn
        # build not hg19? convert
        # always use self._v for downstream annotation
        self._v = self.cleaned_variant_id
        if build != 'hg19':
            self._v = _liftover(self._v,build,'hg19')

    @property
    def exac(self):
        # exac annotation for one variant.
        # check local database first. If not,
        #   use CommonFuncs to annotate exac, then store in database
        if getattr(self, '_exac', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM variants WHERE id=?',(self._v,))
            db_var = dict_factory(db_c,db_c.fetchone())
            if db_var == None or db_var['exac'] == None:
                # query web
                print('querying the web for exac')
                exac = anno_exac(self._v)
                # insert into database
                update_db(
                           self.db_conn,
                           'variants',
                           ['exac'],
                           {self._v:[json.dumps(exac)]}
                           )
            else:
                exac = json.loads(db_var['exac'])
            self._exac = exac
        return self._exac

    @property
    def kaviar_af(self):
        # kaviar annotation for one variant.
        # check local database first. If not,
        #   use CommonFuncs to annotate kaviar, then store in database
        # -1 means NA, None means not populated yet
        if getattr(self, '_kaviar_af', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM variants WHERE id=?',(self._v,))
            db_var = dict_factory(db_c,db_c.fetchone())
            if db_var == None or db_var['kaviar_af'] == None:
                # query web
                print('querying the web for kaviar')
                kaviar_af = anno_kaviar([self._v])[self._v]
                if kaviar_af == None: kaviar_af = -1
                # insert into database
                update_db(
                           self.db_conn,
                           'variants',
                           ['kaviar_af'],
                           {self._v:[kaviar_af]}
                           )
            else:
                kaviar_af = db_var['kaviar_af']
            self._kaviar_af = kaviar_af
        return self._kaviar_af

    @property
    def cadd_phred(self):
        # annotate cadd using database or a provided cadd_file
        if getattr(self, '_cadd_phred', None) is None:
            db_c = self.db_conn.cursor()
            db_c.execute('SELECT * FROM variants WHERE id=?',(self._v,))
            db_var = dict_factory(db_c,db_c.fetchone())
            cadd_phred = None
            if db_var == None or db_var['cadd_phred'] == None:
                if getattr(self,'cadd_file', None) is None:
                    # can't annotate
                    raise ValueError(self.variant_id+' cannot annotate since the cadd_phred for this variant is not in db, and cadd file is not provided. (can do this: object.cadd_file=file)')
                else:
                    # the file data are already cached in the object?
                    if getattr(self,'_cadd_data',None) is None:
                        self._cadd_data = _read_cadd(self.cadd_file)
                    cadd_phred = self._cadd_data.get(self._v, None)
                    if cadd_phred == None:
                        logging.warning(self.variant_id+' not found in cadd file provided')
                    else:
                        # write to db
                        print('write cadd_phred to database')
                        update_db(
                            self.db_conn,
                            'variants',
                            ['cadd_phred'],
                            {self._v:[cadd_phred]}
                            )
            else:
                cadd_phred = db_var['cadd_phred']
            self._cadd_phred = cadd_phred
        return self._cadd_phred

'''
class for an array of variants
'''
class Variants(object):
    def __init__(self, db_conn, vars=None, path_to_gnomad=None, build='hg19'):
        # initiate db
        _initiate_db(db_conn)
        self.variants = vars
        self.build = build
        self.db_conn = db_conn
        self.path_to_gnomad = path_to_gnomad
        self.cleaned_variants = {i:find_leftmost_synonymous_variant(clean_variant(i)) for i in vars}
        # as single Variant, use self._v for downstream annotation
        self._v = copy.copy(self.cleaned_variants)
        if build != 'hg19':
            self._v = {i:_liftover(self._v[i],build,'hg19') for i in list(self._v.values())}
        
    @property
    def exac(self):
        # Check local database first. If not,
        #   use CommonFuncs to annotate exac, then store in database
        if getattr(self, '_exac', None) is None:
            # check database
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'variants',list(self._v.values()))
            data = {}
            new_vars = {}
            exac = {}
            for i in result:
                temp = dict_factory(db_c,i)
                if temp['exac']:
                    data[temp['id']] = json.loads(temp['exac'])
            for k,v in self._v.items():
                if v in data and data[v] != None:
                    exac[k] = data[v]
                else:
                    # not in database, push to array for later query
                    new_vars[k] = v
            if new_vars:
                print('querying web')
                new_result = anno_exac_bulk(list(new_vars.values()))
                # update database
                update_db(
                           self.db_conn,
                           'variants',
                           ['exac'],
                           {k:[json.dumps(v)] for k,v in new_result.items()}
                           )
                # populate exac
                for k,v in new_vars.items():
                    exac[k] = new_result[v]
            self._exac = exac
        return self._exac

    @property
    def kaviar_af(self):
        # -1 means NA, None means not populated yet
        if getattr(self, '_kaviar_af', None) is None:
            # check database
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'variants',list(self._v.values()))
            data = {}
            new_vars = {}
            kaviar = {}
            for i in result:
                temp = dict_factory(db_c,i)
                data[temp['id']] = temp['kaviar_af']
            for k,v in self._v.items():
                if v in data and data[v] != None:
                    kaviar[k] = data[v]
                else:
                    # not in database, push to array for later query
                    new_vars[k] = v
            if new_vars:
                print('querying web')
                new_result = anno_kaviar(list(new_vars.values()))
                # change None to -1 in new_result
                for k in new_result:
                    if new_result[k] == None: new_result[k] = -1
                # update database
                update_db(
                           self.db_conn,
                           'variants',
                           ['kaviar_af'],
                           {k:[v] for k,v in new_result.items()}
                           )
                # populate kaviar_af
                for k,v in new_vars.items():
                    kaviar[k] = new_result[v]
            self._kaviar_af = kaviar
        return self._kaviar_af

    @property
    def cadd_phred(self):
        if getattr(self, '_cadd_phred', None) is None:
            # check database
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'variants',list(self._v.values()))
            data = {}
            new_vars = {}
            cadd_phred = {}
            for i in result:
                temp = dict_factory(db_c, i)
                data[temp['id']] = temp['cadd_phred']
            for k,v in self._v.items():
                if v in data and data[v] != None:
                    cadd_phred[k] = data[v]
                else:
                    new_vars[k] = v
            if new_vars:
                if getattr(self,'cadd_file', None) is None:
                    # can't annotate
                    logging.warning(', '.join(new_vars) + ' cannot annotate since the cadd_phred for this variant is not in db, and cadd file is not provided. (can do this: object.cadd_file=file)')
                    # write un-annotated cadd_phred as None
                    for k in new_vars:
                        cadd_phred[k] = None
                else:
                    self._cadd_data = _read_cadd(self.cadd_file)
                    new_result = {}
                    for k,v in new_vars.items():
                        if v in self._cadd_data:
                            new_result[v] = self._cadd_data[v]
                        cadd_phred[k] = self._cadd_data.get(v,None)
                    # update database
                    if new_result:
                        print('write cadd to database')
                        update_db(
                            self.db_conn,
                            'variants',
                            ['cadd_phred'],
                            {k:[v] for k,v in new_result.items()}
                            )
            self._cadd_phred = cadd_phred
        return self._cadd_phred

    @property
    def gnomad(self):
        # Check local database first. If not,
        #   use CommonFuncs to annotate gnomad, then store in database
        if getattr(self, '_gnomad', None) is None:
            if not self.path_to_gnomad:
                raise ValueError('Required to provide a path to gnomad for annotation')
            # check database
            db_c = self.db_conn.cursor()
            result = batch_query(db_c,'variants',list(self._v.values()))
            data = {}
            new_vars = {}
            gnomad = {}
            for i in result:
                temp = dict_factory(db_c,i)
                if temp['gnomad']:
                    data[temp['id']] = json.loads(temp['gnomad'])
            for k,v in self._v.items():
                if v in data and data[v] != None:
                    gnomad[k] = data[v]
                else:
                    # not in database, push to array for later query
                    new_vars[k] = v
            if new_vars:
                print('querying gnomad')
                # need to divide vars according to their chroms
                new_result = {}
                for chrom_vars in get_chrom_vars(new_vars.values()):
                    new_result.update(gnomad_utils.overall_freqs(chrom_vars, self.path_to_gnomad))
                # update database
                update_db(
                           self.db_conn,
                           'variants',
                           ['gnomad'],
                           {k:[json.dumps(v)] for k,v in new_result.items()}
                           )
                # populate exac
                for k,v in new_vars.items():
                    gnomad[k] = new_result.get(v,None)
            self._gnomad = gnomad
        return self._gnomad
