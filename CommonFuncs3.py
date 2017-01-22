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
import pyliftover #pyliftover is slow!

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
        build = 'grch37'
    elif build=='hg38':
        build = 'grch38'
    end = end or start
    server = "http://%s.rest.ensembl.org" % build
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
            pos = int(pos) - 1
            common_base = find_bases(chrom,pos,build=build)
            ref = common_base
            alt = common_base + alt
    else:
        chrom,pos,ref,alt = v.split('-')
        pos = int(pos)
    if len(ref) < len(alt):
        ran = list(range(len(ref)))
    else:
        ran = list(range(len(alt)))
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

def anno_exac_bulk(vars):
# chop into 100 chunks and send to exac annotation
    vars_array = _chop_array(vars, size=100)
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
def anno_kaviar(vars):
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
            ind += 100
            print('process %s variants' % ind)
            if ind > len(positions):
                position = ', '.join(positions[ind-100:len(positions)])
                br = 1
            else:
                position = ', '.join(positions[ind-100:ind])
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

'''
mygenes use biomart to map gene ids, but sometimes it misses genes, such as entrez_id=6315
in this case, can use gene_card for mapping
gene_card accepts either entrez_id or symbol/aliases
sometimes it returns search result, if more than one matches
sometimes it returns the matching result
returns ensembl_id if has one, otherwise None
'''
def _find_ensembl_id_from_gene_cards(entrez_id,soup):
    # find the div
    the_div = [s for s in soup.findAll("div", {"class" : "gc-subsection"})
                        if '<h3>External' in str(s)][0]
    # then find the li
    lis = the_div.findAll("li")
    entrez_li = [i for i in lis if 'Entrez Gene:' in str(i)][0]
    this_entrez_id = entrez_li.a.get_text()
    
    if this_entrez_id != entrez_id:
        # not this one, return None
        return None
    # find ensembl_id
    ensembl_li = [i for i in lis if 'Ensembl:' in str(i)]
    if not ensembl_li:
        # it has no ensembl id, return 'NA'
        return 'NA'
    return ensembl_li[0].a.get_text()

def gene_cards(entrez_id):
    print('---querying gene_cards. might take some time---')
    # need selenium to scrape. might need the most recent version of geckodriver
    # https://github.com/mozilla/geckodriver/releases
    from selenium import webdriver
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
    from selenium.common.exceptions import TimeoutException
    # search the page
    base_url = 'http://www.genecards.org'
    url = base_url + ('/Search/Identifier?queryString=%s' % entrez_id)
    driver = webdriver.Firefox()
    driver.wait = WebDriverWait(driver, 5)
    driver.get(url)
    
    # is this a list of matches, or the actual gene page?
    # check title if it has 'GeneCards Search Results'
    raw = driver.page_source
    if re.search(r'GeneCards Search Results</title>',raw):
        # search result
        # parse it and get hrefs to match pages
        soup = BeautifulSoup(raw, 'html.parser')
        tds = soup.findAll("td", { "class" : "gc-gene-symbol" })
        hrefs = [base_url+t.a['href'] for t in tds]
        
        # go to each page and check out
        for h in hrefs:
            driver.get(h)
            raw = driver.page_source
            soup = BeautifulSoup(raw, 'html.parser')
            # find external ids section
            ensembl_id = _find_ensembl_id_from_gene_cards(entrez_id,soup)
            if not ensembl_id: continue
            driver.quit()
            return ensembl_id
    else:
        # match page
        soup = BeautifulSoup(raw, 'html.parser')
        # find external ids section
        ensembl_id = _find_ensembl_id_from_gene_cards(entrez_id,soup)
        driver.quit()
        return ensembl_id

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
        #Remove the next line if you want to keep the is_a description info
        is_a = [s.partition(' ! ')[0] for s in is_a]
        all_objects[key] = d["name"] + is_a

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
        all_objects[k] = {'name':v[0],'is_a':v[1:]}

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
        print(r.raise_for_status())
        return False
    decoded = r.json()
    if decoded.get("seq_region_name",None) in VALID_CHROMOSOMES:
        return True
    else:
        return False
