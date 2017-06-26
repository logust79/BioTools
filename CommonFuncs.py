#
#
'''
Some self use and common functions
'''
#
#
import time
import requests
from bs4 import BeautifulSoup
from collections import defaultdict
import json
import mygene
import re
import tabix
import os

'''
constants
'''
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X','Y']

'''
request ensembl for bases based on location
'''
def find_bases(chrom,start,end=None,build='hg19',strand=1):
    # translate build
    if build=='hg19':
        server = "http://grch37.rest.ensembl.org"
    elif build=='hg38':
        server = "http://rest.ensembl.org"
    end = end or start
    ext = '''/sequence/region/human/%(chrom)s:%(start)s..%(end)s:%(strand)s''' % locals()
    attempt = 5
    while attempt:
        try:
            r = requests.get(server+ext, headers={'Content-Type':'application/json' })
            time.sleep(0.05)
            break
        except requests.HTTPError:
            print('query ensembl connectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    if not r.ok:
        return r.raise_for_status()
    decoded = r.json()
    return str(decoded['seq'])

'''
    liftover between different human genome builds, say hg38 to hg19
'''
def liftover(v, frm, to):
    import pyliftover #pyliftover is slow!
    # note that pyliftover is 0 based
    # First frm-to pair may take time to download the data from UCSC
    # return a list of tuple.
    lo = pyliftover.LiftOver(frm, to)
    chrom,pos,ref,alt = v.split('-')
    results = lo.convert_coordinate('chr'+chrom, int(pos)-1)
    if not results:
        return []
    return ['-'.join([i[0][3:],str(i[1]+1),ref,alt]) for i in results]

'''
    this function cleans messy variant format
    var = '1-94512001-GTT-GAT'
    print clean_variant(var)
    prints
    1-94512002-T-A
'''
def clean_variant(v,build='hg19'):
    # sometimes variant has funny format, which has more - than expected, such as 1-117122294---TCT.
    #  use find_bases to fill in the gap
    if v.count('-') == 4:
        if v[-1] == '-':
            # deletion
            chrom,pos,ref,rubbish,rubbish = v.split('-')
            pos = int(pos)-1
            common_base = find_bases(chrom,pos,build=build)
            ref = common_base + ref
            alt = common_base
        else:
            # insertion
            chrom,pos,ref,rubbish,alt = v.split('-')
            pos = int(pos)-1
            common_base = find_bases(chrom,pos,build=build)
            ref = common_base
            alt = common_base + alt
    else:
        chrom,pos,ref,alt = v.split('-')
        pos = int(pos)
    if len(ref) < len(alt):
        ran = range(len(ref))
    else:
        ran = range(len(alt))
    # insert
    for e in ran:
        ref_e = len(ref) - e - 1
        alt_e = len(alt) - e - 1
        if ref[ref_e] != alt[alt_e]: break
    for b in ran:
        if ref[b] != alt[b] or len(ref[b:ref_e+1]) == 1 or len(alt[b:alt_e+1]) == 1:
            break
    return '-'.join([chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]])

'''
use harvard's rest service to query ExAC_freq of a variant
anno_exac('1-123-G-C')
'''

def anno_exac(v):
    rest_url = 'http://exac.hms.harvard.edu/rest'
    service = 'variant'
    attempt = 5
    while attempt:
        try:
            r = requests.get('/'.join([rest_url,service,v]))
            break
        except requests.ConnectionError:
            print('query exac connectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    exac_anno = r.json()
    return exac_anno
    # good to return the json object
'''
    if 'allele_freq' in exac_anno['variant']:
        return exac_anno['variant']['allele_freq']
    if exac_anno['base_coverage'] and exac_anno['base_coverage'][0]['has_coverage']:
        return 0
    return None
'''

def anno_exac_bulk(vars, chunk_size=100):
# chop into 100 chunks and send to exac annotation
    vars_array = _chop_array(vars, chunk_size)
    result = {}
    for vs in vars_array:
        print(vs)
        result.update(_anno_exac_bulk_100(vs))
    return result

def _anno_exac_bulk_100(vars):
    # limit to 100 variants
    rest_url = 'http://exac.hms.harvard.edu/rest/bulk/variant'
    attempt = 5
    while attempt:
        try:
            r = requests.post(rest_url, data=json.dumps(vars))
            break
        except requests.ConnectionError:
            print('query exac connectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    exac_anno = r.json()
    return exac_anno

def _chop_array(arr, size=100):
    for i in range(0, len(arr), size):
        yield arr[i:i + size]

'''
this function queries kaviar for allele frequency
it queries a chromosome at a time for multiple locations
vars = ['1-123-G-C','1-234-C-T','2-234-T-GA']
result = anno_kaviar(vars)
hg19
'''
def anno_kaviar(vars, chunk_size=100):
    # collapse on chroms
    chroms = {} # {1:[{pos:123,ref:A,alt:T},]}
    for v in vars:
        chrom,pos,ref,alt = v.split('-')
        chroms[chrom] = chroms.get(chrom,[])
        chroms[chrom].append({'pos':pos,'ref':ref,'alt':alt})
    # get kaviar result
    uri = 'http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl'
    args = {
        'frz' : 'hg19',
        'onebased': '1',
        'onebased_output': '1', # on the website the default is somehow 0
        'chr': '',
        'platform_specificity':'none',
        'format': 'text', # json is bugged.
        'pos':''
    }
    kaviar = []
    br=0
    for c in chroms:
        positions = list(set([i['pos'] for i in chroms[c]]))
        ind = 0
        while True:
            ind += chunk_size
            if ind > len(positions):
                position = ', '.join(positions[ind-chunk_size:len(positions)])
                if not position: break
                br = 1
            else:
                position = ', '.join(positions[ind-chunk_size:ind])
            print('process %s variants' % min(ind,len(positions)))
            args['pos'] = position
            args['chr'] = c
            print(args)
            r = requests.get(uri,params=args)
            # parse it
            soup = BeautifulSoup(r.content, 'html.parser')
            header = []
            for l in soup.body.pre.get_text().split('\n'):
                if not l:
                    continue
                if l[:3] == '#Ch':
                    # parse header
                    header = l[1:].split('\t')
                elif l[0] != '#':
                    # get data
                    l = l.split('\t')
                    pos = l[1]
                    if len(l) < len(header):
                        l.extend(['']*(len(header)-len(l)))
                    kaviar.append(dict([(header[i],l[i]) for i in range(len(header))]))
            #print pos
            if br: break
    result = {}
    for v in vars:
        chrom,pos,ref,alt = v.split('-')
        end = int(pos)+len(ref)-1
        match = [i for i in kaviar if i['Chrom'] == 'chr'+chrom and i['Position']==pos and i['End'] == str(end) and i['Variant']==alt]
        result[v] = float(match[0]['AF']) if match and match[0]['AF'] else None
    return result

def my_gene(gene_id):
    mg = mygene.MyGeneInfo()
    return mg.getgene(gene_id,fields='all')
def my_genes(gene_ids):
    mg = mygene.MyGeneInfo()
    return mg.getgenes(gene_ids,fields='all')
def my_genes_by_symbol(symbols,species=None):
    mg = mygene.MyGeneInfo()
    result = mg.querymany(symbols, scopes='symbol', species=species,fields='all')
    # which ones are not found
    not_found = [i['query'] for i in result if i.get('notfound',False) == True]
    # query again on alias
    result.extend(mg.querymany(not_found, scopes='alias', species=species,fields='all'))
    return result

def obo_parser(obofile):
    term_head = "[Term]"
    #Keep the desired object data here
    all_objects = {}

    def add_object(d):
        #Ignore obsolete objects
        if "is_obsolete" in d:
            return

        #Gather desired data into a single list,
        # and store it in the main all_objects dict
        key = d["id"][0]
        is_a = d["is_a"]
        alt_id = d["alt_id"]
        #Remove the next line if you want to keep the is_a description info
        is_a = [s.partition(' ! ')[0] for s in is_a]
        all_objects[key] = {
            'name':d["name"],
            'is_a':is_a,
            'alt_id':alt_id,
        }

    #A temporary dict to hold object data
    current = defaultdict(list)

    with open(obofile,'r') as f:
        #Skip header data
        for line in f:
            if line.rstrip() == term_head:
                break

        for line in f:
            line = line.rstrip()
            if not line:
                #ignore blank lines
                continue
            if line == term_head:
                #end of term
                add_object(current)
                current = defaultdict(list)
            else:
                #accumulate object data into current
                key, _, val = line.partition(": ")
                current[key].append(val)

    if current:
        add_object(current)
    # convert to dict
    for k,v in all_objects.items():
        all_objects[k] = {'name':v['name'],'is_a':v['is_a'],'alt_id':v['alt_id']}

    return all_objects

'''
check if ensembl id is active
'''
def check_ensemblId(ensemblId):
    url = 'http://rest.ensembl.org/lookup/id/'+ensemblId
    attempt = 5
    while attempt:
        try:
            r = requests.get(url, headers={ "Content-Type" : "application/json"})
            time.sleep(0.05)
            break
        except requests.HTTPError:
            print('query ensembl connectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    if not r.ok:
        #print(r.raise_for_status())
        return False
    decoded = r.json()
    if decoded.get("seq_region_name",None) in VALID_CHROMOSOMES:
        return True
    else:
        return False

'''
gnomad annotation. Currently there's no API service. Relying on local gnomad files
'''
#path = '/media/jing/SEAGATE/db/gnomAD'
#release = 'release-170228'

'''
coverage
'''
def gnomad_coverage(v,path_to_gnomad,mode='exome'):
    v = clean_variant(v)
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','coverage','exacv2.chr'+chrom+'.cov.txt.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','coverage','gnomad.chr'+chrom+'.cov.txt.gz')
    else:
        msg = "mode only accepts 'exome' or 'genome'"
        raise ValueError(msg)
    tb = tabix.open(file)
    r = tb.query(chrom, int(pos)-1, int(pos))
    r = list(r)
    if not r:
        # not covered
        return None
    r = r[0]
    return {a:b for a,b in zip(header,r)}

'''
freqs
'''
def gnomad_freqs(v,path_to_gnomad,mode='exome'):
    v = clean_variant(v)
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','id','ref','alt','quality','filter','info']
    chrom,pos,ref,alt = v.split('-')
    if mode == 'exome':
        file = os.path.join(path_to_gnomad,'exomes','vcf','gnomad.exomes.r2.0.1.sites.vcf.gz')
    elif mode == 'genome':
        file = os.path.join(path_to_gnomad,'genomes','vcf','gnomad.genomes.r2.0.1.sites.'+chrom+'.vcf.gz')
    tb = tabix.open(file)
    records = tb.query(chrom, int(pos)-1, int(pos))

    for r in records:
        if not r: return None
        data = {a:b for a,b in zip(header,r)}
        # find the variant
        g_alts = data['alt'].split(',')
        alt_ind = None
        for ind,this_alt in enumerate(g_alts):
            v_id = clean_variant('-'.join([data['chrom'],data['pos'],data['ref'],this_alt]))
            if v_id == v:
                alt_ind = ind

        if alt_ind == None:
            continue

        # parse info
        # no need for annotation?

        info = data['info'].split(';CSQ=A')[0] # 1 for annotation
        info = info.split(';')
        info_dict = {}
        for i in info:
            if not '=' in i: continue
            a,b = i.split('=')
            b = b.split(',')
            ind = min(alt_ind,len(b)-1)
            b = b[ind]
            # convert to number if possible
            try:
                if '.' in b:
                    b = float(b)
                else:
                    b = int(b)
            except ValueError:
                pass
            info_dict[a] = b
        info_dict['FILTER'] = info_dict['AS_FilterStatus']
        return info_dict
                

    return None
