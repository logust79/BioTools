'''
get het and hom freq from Bravo
https://bravo.sph.umich.edu/
Note that the hg19 version was liftovered from hg38 using Picard (by Jing Yu)
Also note that bravo file uses chr1, not 1 for chrom coding
'''
from __future__ import division,print_function
import pysam
import gzip

def get_chrom(variant_id):
    chrom = variant_id.split('-')[0]
    chrom = 'chr{}'.format(chrom)
    return chrom

def get_variants_from_tbx(fetch, variants, header):
    variants_set = set(variants)
    result = {}
    for line in fetch:
        row = line.split('\t')
        record = dict(zip(header, row))
        v_id = '-'.join([
            record['CHROM'].replace('chr',''),
            record['POS'],
            record['REF'],
            record['ALT']
        ])
        if v_id not in variants:
            continue
        # get all info in 'INFO'
        info = {}
        for field in record['INFO'].split(';'):
            field = field.split('=')
            info[field[0]] = field[1]
        result[v_id] = {
                'filter': record['FILTER'],
                'af': float(info['AF']),
                'ac': int(info['AC']),
                'an': int(info['AN']),
                'het': int(info['Het']),
                'Hom': int(info['Hom']), # easy enought, bravo use Hom for chrX
        }
    return result

def bravo(variants, bravo_vcf, group = True):
    '''
    if group is True, should be just one chrom and all variants should
    be close together, so it will query in one go
    '''
    result = {}
    if len(variants) == 0:
        return result
    # get header first
    header = []
    with gzip.open(bravo_vcf,'rt') as inf:
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

    tbx = pysam.TabixFile(bravo_vcf)

    if group:
        # find chrom and range
        chrom = get_chrom(variants[0])
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
            chrom = get_chrom(variant)
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
    vs = ['1-10051-A-AC','1-10132-T-C','1-123-T-C']
    vcf = '/cluster/project8/vyp/bravo-dbsnp-hg19.vcf.gz'
    result = bravo(vs, vcf, group=False)
    print(result)
