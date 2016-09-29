#
#
'''
Some self use and common functions
'''
#
#
import requests
from bs4 import BeautifulSoup

'''
this function cleans messy variant format
    var = '1-94512001-GTT-GAT'
    print clean_variant(var)
prints
    1-94512002-T-A
'''
def clean_variant(v):
    chrom,pos,ref,alt = v.split('-')
    pos = int(pos)
    if len(ref) < len(alt):
        ran = range(len(ref))
    else:
        ran = range(len(alt))
    # insert
    for b in ran:
        if ref[b] != alt[b]:
            break
    for e in ran:
        ref_e = len(ref) - e - 1
        alt_e = len(alt) - e - 1
        if ref[ref_e] != alt[alt_e] or len(ref[b:ref_e+1]) == 1 or len(alt[b:alt_e+1]) == 1:
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
            print 'query exac connectionError, retry'
            attempt -= 1
            time.sleep(2)
    if r.status_code == 404: return None
    exac_anno = r.json()
    if 'allele_freq' in exac_anno['variant']:
        return exac_anno['variant']['allele_freq']
    if exac_anno['base_coverage'][0]['has_coverage']:
        return 0
    return None

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
    ind=0
    br=0
    for c in chroms:
        positions = list(set([i['pos'] for i in chroms[c]]))
        while True:
            ind += 100
            print 'process %s variants' % ind
            if ind > len(positions):
                position = ', '.join(positions[ind-100:len(positions)])
                br = 1
            else:
                position = ', '.join(positions[ind-100:ind])
            args['pos'] = position
            args['chr'] = c
            r = requests.get(uri,params=args)
            # parse it
            soup = BeautifulSoup(r.content, 'html.parser')
            header = []
            for l in soup.body.pre.string.split('\n'):
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
        result[v] = float(match[0]['AF']) if match else None
    return result
