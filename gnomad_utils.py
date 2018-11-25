'''
annotate gnomad population frequencies
!!!NOTE!!! that you need to normalise all variants before using this script, e.g.
maps = {variant:utils.clean_variant(variant) for variant in variants}
gnomads = gnomad_utils.overall_freqs(list(maps.values()),path_to_gnomad,group=True)
result = {variant:gnomads.get(maps[variant],None) for variant in maps}
'''
from __future__ import print_function, division
import sys
import pysam
import os
import CommonFuncs

class BadVariantException(Exception): pass
#path = '/cluster/project8/vyp/gnomad_data'
VALID_CHROMOSOMES = [str(i) for i in range(1,23)] + ['X']
POPS = ('AFR','NFE','AMR','EAS','ASJ','SAS','OTH','FIN')
'''
coverage, look at different parts of the ref
'''
def coverage(vs,path_to_gnomad,mode,chrom,start,stop):
    # pysam does not support header yet. hard code it
    header = ['chrom','pos','mean','median',1,5,10,15,20,25,30,50,100,]
    result = {v:None for v in vs}
    vcf = None
    if mode == 'exome':
        vcf = os.path.join(path_to_gnomad,'exomes','coverage','exacv2.chr'+chrom+'.cov.txt.gz')
    elif mode == 'genome':
        vcf = os.path.join(path_to_gnomad,'genomes','coverage','gnomad.chr'+chrom+'.cov.txt.gz')
    else:
        msg = "mode only accepts 'exome' or 'genome'"
        raise ValueError(msg)
    try:
        tb = pysam.TabixFile(vcf)
        rs = tb.fetch(chrom, start-1, stop+1)
    except OSError:
        return result
    except ValueError:
        return result

    cov_record = {}
    for r in rs:
        this = {a:b for a,b in zip(header,r.split('\t'))}
        cov_record[int(this['pos'])] = this
    for v in vs:
        _,pos,ref,_ = v.split('-')
        pos = int(pos)
        end = pos + len(ref)
        for i in range(pos,end):
            if i not in cov_record:
                result[v] = None
                break
            if result[v] is None:
                result[v] = {pos: cov_record[pos]}
            else:
                result[v][pos] = cov_record[pos]
            
    return result

'''
exome freqs
'''
def freqs(vs,path_to_gnomad,mode,chrom,start,stop):
    # pytabix does not support header yet. hard code it
    header = ['chrom','pos','id','ref','alt','quality','filter','info']
    if mode == 'exome':
        vcf = os.path.join(path_to_gnomad,'exomes','vcf','gnomad.exomes.r2.0.1.sites.vcf.gz')
    elif mode == 'genome':
        vcf = os.path.join(path_to_gnomad,'genomes','vcf','gnomad.genomes.r2.0.1.sites.'+chrom+'.vcf.gz')
    result = {v:{} for v in vs}
    try:
        tb = pysam.TabixFile(vcf)
        records = tb.fetch(chrom, start-1, stop+1)
    except OSError:
        return result
    except ValueError:
        return result
    #tb = tabix.open(vcf)
    #records = tb.query(chrom, int(pos)-1, int(pos))

    for r in records:
        if not r: return result
        data = {a:b for a,b in zip(header,r.split('\t'))}
        # find the variant
        g_alts = data['alt'].split(',')
        alt_ind = None
        v_ids = []
        for ind,this_alt in enumerate(g_alts):
            v_id = CommonFuncs.clean_variant('-'.join([data['chrom'],data['pos'],data['ref'],this_alt]))
            if v_id in vs:
                v_ids.append((v_id, ind))
        if not v_ids:
            continue

        # parse info
        # no need for annotation?

        info = data['info'].split(';CSQ=A')[0] # 1 for annotation
        info = info.split(';')
        for i in info:
            if not '=' in i: continue
            a,b = i.split('=')
            b = b.split(',')
            for v_ind in range(len(v_ids)):
                ind = min(v_ids[v_ind][1],len(b)-1)
                c = b[ind]
                # convert to number if possible
                try:
                    if '.' in c:
                        c = float(c)
                    else:
                        c = int(c)
                except ValueError:
                    pass
                result[v_ids[v_ind][0]][a] = c
        for v in v_ids:
            result[v[0]]['filter'] = data['filter']
                
    return result


'''
simple query on overall allele freq (or homozygote frequency). if covered, return at least 0. if not, return None. 
Assuming all variants are on the same chrom, and close together. If not, please
query a single variant (in a list) at a time
'''
def overall_freqs(vs,path_to_gnomad):
    # get start, stop and chrom if block
    chrom,start,stop = CommonFuncs.get_chrom_start_stop(vs)
    result = {}
    null = {
        'gnomad_af': None,
        'gnomad_ac': None,
        'gnomad_hom_f': None, # if X chrom, this takes into account of hemi
        'gnomad_hom_c': None,
        'gnomad_hemi_c': None,
        'filters':{'exome':None,'genome':None},
        'pop_filter':[],
        'most_freq_pops':[],
    }
    covs = {
        'exome': coverage(vs,path_to_gnomad,'exome',chrom,start,stop),
        'genome':coverage(vs,path_to_gnomad,'genome',chrom,start,stop),
    }
    fs = {
        'exome': freqs(vs,path_to_gnomad,'exome',chrom,start,stop),
        'genome':freqs(vs,path_to_gnomad,'genome',chrom,start,stop),
    }
    for v in vs:

        if v.split('-')[0] not in VALID_CHROMOSOMES:
            result[v] = null
            continue

        if not fs['exome'][v] and not fs['genome'][v] and not covs['exome'][v] and not covs['genome'][v]:
            result[v] = null
            continue
        ac = hom_c = af = hom_f = an = 0.
        hemi_c = None
        pop_filter = []
        filters = {'exome':None,'genome':None}
        # also check population frequencies to remove any variants 
        # with big af(>0.01)/hom_f(0) discrepancy, such as 1-144931607-C-T
        pops = {p:{'Hom':0,'Hemi':0,'AC':0,'AN':0} for p in POPS}
        for m in ['exome', 'genome']:
            if fs[m][v]:
                ac += fs[m][v]['AC']
                hom_c += fs[m][v]['Hom']
                an += fs[m][v]['AN']
                if 'Hemi' in fs[m][v]:
                    hemi_c = hemi_c + fs[m][v]['Hemi'] if hemi_c != None else fs[m][v]['Hemi']
                filters[m] = fs[m][v]['filter']
                for p in POPS:
                    for kk in pops[p]:
                        try:
                            this_c = fs[m][v].get('{}_{}'.format(kk,p), 0)
                            # sometimes on X this_c is a dot
                            if this_c == '.':
                                this_c = 0
                            pops[p][kk] += this_c
                        except TypeError:
                            print(v)
                            print(p,kk)
                            print(pops[p][kk])
                            print(fs[m][v].get('{}_{}'.format(kk,p), 0))
                            raise

        max_pop = ([], -1)
        for p in pops:
            if pops[p]['AC']:
                pop_af = pops[p]['AC'] / pops[p]['AN']
                if pop_af:
                    if pop_af > max_pop[1]:
                        max_pop = ([p], pop_af)
                    elif pop_af == max_pop[1]:
                        max_pop[0].append(p)
                if pop_af > 0.01 and pops[p]['Hemi'] == 0 and pops[p]['Hom'] == 0:
                    pop_filter.append(p)

        if ac: af = ac / an
        if hom_c:
            if hemi_c is not None and isinstance(hemi_c, int):
                hom_f = (hom_c * 2 + hemi_c) / an
            else:
                hom_f = hom_c * 2 / an

        result[v] = {
            'gnomad_af':af,
            'gnomad_ac':ac,
            'gnomad_hom_f':hom_f,
            'gnomad_hom_c':hom_c,
            'gnomad_hemi_c':hemi_c,
            'gnomad_an':an,
            'filters':filters,
            'pop_filter':pop_filter,
            'most_freq_pops':max_pop[0],
        }

    return result

if __name__ == '__main__':
    vs = ['M-150-T-C','1-12140-GCAT-C','1-12141-CAT-C','1-12143-T-C','1-40705-C-T','1-165389-CA-C','1-107189335-C-CTTTT']
    p2g = '/cluster/project8/vyp/gnomad_data'
    #print(freqs(vs[1:-1],p2g,'genome','1',12140,165390))
    print(overall_freqs(vs[1:-1],p2g))
    vs = ['1-55516888-G-GA','1-55516887-GG-GGA']
    print(overall_freqs(vs,p2g))
