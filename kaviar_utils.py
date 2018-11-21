'''
query kaviar to get af. no hom_f is available
'''
from __future__ import division, print_function
import gzip
import pysam
import CommonFuncs

def get_variants_from_tbx(fetch, variants, header):
    variants_set = set(variants)
    result = {}
    for line in fetch:
        row = line.split('\t')
        record = dict(zip(header, row))
        info = {}
        for field in record['INFO'].split(';'):
            field = field.split('=')
            info[field[0]] = field[1]
        for ind,alt in enumerate(record['ALT'].split(',')):
            v_id = CommonFuncs.clean_variant('-'.join([
                record['CHROM'],
                record['POS'],
                record['REF'],
                alt,
            ]))
            if v_id not in variants:
                continue
            # get all info in 'INFO'
            result[v_id] = result.get(v_id,{
                'filter': None,
                'af': None,
                'ac': 0,
                'an': 0
            })

            result[v_id] = {
                    'filter': None,
                    'ac': int(info['AC'].split(',')[ind]) + result[v_id]['ac'],
                    'an': int(info['AN']) + result[v_id]['an'],
            }
            result[v_id]['af'] = result[v_id]['ac'] / result[v_id]['an']
    return result

def kaviar(variants, kaviar_vcf, group=True):
    '''
    if group is True, should be just one chrom and all variants should
    be close together, so it will query in one go
    '''
    result = {}
    if len(variants) == 0:
        return result
    # get header first
    header = []
    with gzip.open(kaviar_vcf,'rt') as inf:
        for line in inf:
            if line.startswith('##'):
                continue
            row = line.rstrip().split()
            if not header:
                header = row
                # strip of # from #CHROM
                header[0] = header[0][1:]
                # all of these vcf files have just one sample.
                break

    tbx = pysam.TabixFile(kaviar_vcf)

    if group:
        # find chrom and range
        chrom = variants[0].split('-')[0]
        positions = [(int(v.split('-')[1]), v) for v in variants]
        positions = sorted(positions, key=lambda x: x[0])
        min_pos = positions[0][0]
        max_pos = positions[-1][0]
        try:
            fetch = tbx.fetch(chrom, min_pos-1, max_pos+1)
        except ValueError:
            result = {}
        else:
            result = get_variants_from_tbx(fetch, variants, header)
    else:
        # deal with variants one by one
        for variant in variants:
            chrom = variant.split('-')[0]
            pos = int(variant.split('-')[1])
            try:
                fetch = tbx.fetch(chrom, pos-1, pos+1)
            except ValueError:
                result = {}
            else:
                this = get_variants_from_tbx(
                        fetch,
                        [variant],
                        header)
                if variant in this:
                    result[variant] = this[variant]

    return result

if __name__ == '__main__':
    vs = ['1-10002-A-C','1-10003-A-T','1-884091-C-CACCCTGGTCCCCCTGGTCCCTTTGGCCCTGCACCTGGCTGG']
    vcf = '/cluster/project8/vyp/kaviar/Kaviar-160204-Public-hg19-trim.vcf.gz'
    result = kaviar(vs, vcf, group=True)
    print(result)
