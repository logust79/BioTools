#
#
'''
Some self use and common functions
'''
#
#
import time
from requests.exceptions import HTTPError
import json
import os
import re
import requests
from collections import defaultdict

'''
constants
'''
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X','Y']

'''
divide variants according to their chrom and chunk size
'''
def get_chrom_vars(variants, chunk_size=1e5):
    result = {}
    for v in sorted(variants, key=lambda x: int(x.split('-')[1])):
        chrom,pos,_,_ = v.split('-')
        pos = int(pos)
        if chrom in result:
            if pos > result[chrom][-1]['start'] + chunk_size:
                result[chrom].append({'start':pos, 'variants':[v]})
            else:
                result[chrom][-1]['variants'].append(v)
        else:
            result[chrom] = [{'start':pos, 'variants':[v]}]
    for chrom, chunks in result.items():
        for chunk in chunks:
            yield chunk['variants']

'''
get chrom, start, stop for a group of variants
'''
def get_chrom_start_stop(vs):
    chrom = vs[0].split('-')[0]
    vs_arrays = [v.split('-') for v in vs]
    starts = [int(v[1]) for v in vs_arrays]
    stops = [starts[i] + len(vs_arrays[i][2]) for i in range(len(vs))]
    start = min(starts) - 1
    stop = max(stops) + 1
    return chrom,start,stop

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
            if r.ok:
                break
        except requests.HTTPError:
            print('query ensembl HTTPError, retry')
            attempt -= 1
            time.sleep(2)
        except requests.ConnectionError:
            print('query ensembl ConnectionError, retry')
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    '''
    if not r.ok:
        return r.raise_for_status()
    '''
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
def clean_variant(v,build='hg19',human_ref_pysam=None):
    # sometimes variant has funny format, which has more - than expected, such as 1-117122294---TCT.
    #  use find_bases to fill in the gap if human_ref_pysam is not provided
    if v.count('-') == 4:
        if v[-1] == '-':
            # deletion
            chrom,pos,ref,rubbish,rubbish = v.split('-')
            pos = int(pos)-1
            if human_ref_pysam:
                common_base = human_ref_pysam.fetch(chrom, pos-1, pos)
            else:
                common_base = find_bases(chrom,pos,build=build)
            ref = common_base + ref
            alt = common_base
        else:
            # insertion
            chrom,pos,ref,rubbish,alt = v.split('-')
            pos = int(pos) - 1
            if human_ref_pysam:
                common_base = human_ref_pysam.fetch(chrom, pos-1, pos)
            else:
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
    for b in ran:
        if ref[b] != alt[b] or len(ref[b:]) == 1 or len(alt[b:]) == 1:
            break
    if len(ref) == len(alt):
        for e in ran:
            ref_e = len(ref) - e - 1
            alt_e = len(alt) - e - 1
            if ref[ref_e] != alt[alt_e]: break
    else:
        ref_e = len(ref)
        alt_e = len(alt)
    return '-'.join([chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]])

def find_start_of_repeat(string, start, length):
    '''
    string: GCAGAGAGAG
    start: 5
    length: 2 #GA
    return 1 # the repeat starts from after 1:AG
    ===
    if no repeat, return start
    '''
    ind = start
    result = string[:start]+string[start+length:]
    while ind >= 0:
        ind -= 1
        if string[:ind] + string[ind+length:] != result:
            return ind
    return ind


def find_leftmost_synonymous_variant(variant, padding=200, build='hg19', human_ref_pysam=None):
    '''
    Only necessary for indel!
    find all synonymous variants given variant
    padding is how far you would like to search left and right of the change
    if human_ref_pysam is None, use find_base to query ensembl
    '''
    mode,pattern = None,None
    chrom, pos, ref, alt = variant.split('-')
    pos = int(pos)+1
    # removing commong base
    if ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    if len(ref) and not alt:
        pattern = ref
        mode = 'del'
    elif len(alt) and not ref:
        pattern = alt
        mode = 'in'
    else:
        # it's not indel, return
        return variant
    if human_ref_pysam:
        string = human_ref_pysam.fetch(chrom, pos-padding-1, pos+len(pattern)+padding-1)
    else:
        string = find_bases(chrom, pos-padding, pos+len(pattern)+padding, build=build)
    ind = find_start_of_repeat(string, padding, len(pattern))
    new_pos = pos - padding + ind
    missing_base = string[ind]
    pattern = string[ind+1:ind+len(pattern)+1]

    if mode == 'del':
        return '-'.join([chrom, str(new_pos), missing_base+pattern, missing_base])
    elif mode == 'in':
        return '-'.join([chrom, str(new_pos), missing_base, missing_base+pattern])
    else:
        msg = 'Cannot derive mode!'
        raise ValueError(msg)
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
    from bs4 import BeautifulSoup
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
    import mygene
    mg = mygene.MyGeneInfo()
    return mg.getgene(gene_id,fields='all')
def my_genes(gene_ids):
    import mygene
    mg = mygene.MyGeneInfo()
    return mg.getgenes(gene_ids,fields='all')
def my_genes_by_symbol(symbols,species=None):
    import mygene
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

def add_pop_freqs(infile, outfile, options):
    '''
    add pop freqs, such as gnomad, bravo and kaviar
    '''
    import pysam
    import gzip
    options['human_ref_pysam'] = pysam.FastaFile(options['human_ref'])
    fields = []
    if 'gnomad_path' in options['pop_freqs']:
        fields.extend(['gnomad_af', 'gnomad_hom_f'])
    if 'bravo_vcf' in options['pop_freqs']:
        fields.extend(['bravo_af','bravo_hom_f'])
    if 'kaviar_vcf' in options['pop_freqs']:
        fields.append('kaviar_af')
    line_cache = []
    variant_cache = {}
    header = []
    with gzip.open(infile, 'r') as inf, \
            gzip.open(outfile, 'w') as outf:
        for line in inf:
            if line.startswith('##'):
                outf.write(line)
            elif line.startswith('#'):
                # add an INFO line
                info_line = construct_pop_info(fields)
                outf.write(info_line)
                outf.write(line)
                header = line[1:].rstrip().split('\t')
            else:
                # write to cache
                line_cache.append(line)
                row = line.split('\t')
                for alt in row[4].split(','):
                    v_id = utils.clean_variant(
                            '-'.join([row[0],row[1],row[3],alt]),
                            human_ref_pysam = options['human_ref_pysam']
                    )
                    variant_cache[v_id]={}

                if len(line_cache) >= options['pop_freqs']['cache_size']:
                    # dump cache
                    pop_annotate(line_cache, variant_cache, header, fields, outf, options)
                    # empty caches
                    line_cache = []
                    variant_cache = {}

        # finally dump cache
        if line_cache:
            pop_annotate(line_cache, variant_cache, header, fields, outf, options)
            line_cache = []
            variant_cache = {}

def pop_annotate(line_cache, variant_cache, header, fields, outf, options):
    import gnomad_utils, bravo_utils, kaviar_utils
    # annotate
    # gnomad
    if 'gnomad_path' in options['pop_freqs']:
        gnomads = gnomad_utils.overall_freqs(
                list(variant_cache.keys()),
                options['pop_freqs']['gnomad_path']
        )
        for variant in variant_cache:
            af = gnomads[variant]['gnomad_af']
            if af is None:
                af = ''
            hom_f = gnomads[variant]['gnomad_hom_f']
            if hom_f is None:
                hom_f = ''
            variant_cache[variant]['gnomad_af'] = af

            variant_cache[variant]['gnomad_hom_f'] = hom_f

    # bravo
    if 'bravo_vcf' in options['pop_freqs']:
        bravos = bravo_utils.bravo(
                list(variant_cache.keys()),
                options['pop_freqs']['bravo_vcf']
        )
        for variant in variant_cache:
            if variant not in bravos:
                variant_cache[variant]['bravo_af'] = ''
                variant_cache[variant]['bravo_hom_f'] = ''
            else:
                variant_cache[variant]['bravo_af'] = \
                    bravos[variant]['af']
                variant_cache[variant]['bravo_hom_f'] = \
                    bravos[variant]['Hom']*2 / bravos[variant]['an']
    # kaviar
    if 'kaviar_vcf' in options['pop_freqs']:
        kaviars = kaviar_utils.kaviar(
                list(variant_cache.keys()),
                options['pop_freqs']['kaviar_vcf']
        )
        for variant in variant_cache:
            if variant not in kaviars:
                variant_cache[variant]['kaviar_af'] = ''
            else:
                variant_cache[variant]['kaviar_af'] = \
                    kaviars[variant]['af']

    for line in line_cache:
        record = dict(zip(header, line.rstrip().split('\t')))
        INFO = record['INFO']
        pop_info = []
        for alt in record['ALT'].split(','):
            v_id = clean_variant('-'.join([
                record['CHROM'],
                record['POS'],
                record['REF'],
                alt,
            ]), human_ref_pysam = options['human_ref_pysam'])
            pop_info.append('|'.join(
                [alt] + [str(variant_cache[v_id][f]) for f in fields]
            ))
        pop_info = 'POPF=' + ','.join(pop_info)
        new_INFO = ';'.join([INFO, pop_info])
        record['INFO'] = new_INFO
        new_line = '\t'.join([record[h] for h in header]) + '\n'
        outf.write(new_line)

def construct_pop_info(fields):
    info_line = (
            '##INFO=<ID=POPF,Number=.,Type=String,Description="'
            'Population frequency such as gnomAD and Bravo. '
            'Format: Allele'
    )
    # add fields
    info_line = '|'.join([info_line] + fields)
    # closing
    info_line += '">\n'
    return info_line

