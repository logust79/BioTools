import unittest
from CommonFuncs import *
import os

class CommonFuncsTestCase(unittest.TestCase):
    def test_clean_variant(self):
        case = clean_variant('13-95363811---GCG')
        self.assertEqual(case,'13-95363811-C-CGCG')
        case = clean_variant('13-95363811---GCG', 'hg38')
        self.assertEqual(case,'13-95363811-C-CGCG')
        case = clean_variant('2-220400050-TGTGTGTGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTACGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTAC-GGTGTGTGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTACGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTAC')
        self.assertEqual(case, '2-220400050-T-G')
        case = clean_variant('13-95363809-GG--')
        self.assertEqual(case, '13-95363808-CGG-C')

    def test_find_bases(self):
        case = find_bases('13', 95363809)
        self.assertEqual(case,'G')
        case = find_bases('13', 95363809, 95363830)
        self.assertEqual(case, 'GGCGGCGGCGGCGGCGGCAGCG')
        case = find_bases('13', 95363809, 95363830, strand = -1)
        self.assertEqual(case, 'CGCTGCCGCCGCCGCCGCCGCC')
        case = find_bases('13', 95363809, 95363830, 'hg38')
        self.assertEqual(case, 'CACAGTTCAATGACAACAGGGA')

    def test_liftover(self):
        case = liftover('2-220400050-T-G','hg19','hg38')[0]
        self.assertEqual(case,'2-219535328-T-G')
        case = liftover('2-219535328-T-G','hg38','hg19')[0]
        self.assertEqual(case,'2-220400050-T-G')
    
    def test_anno_exac(self):
        case = anno_exac('2-220400050-T-G')['variant']['allele_freq']
        self.assertEqual(case, 2.4078979051288225e-05)
        case = anno_exac('5-178413161-C-G')
        self.assertTrue(case['any_covered'])
        self.assertFalse('allele_freq' in case['variant'])
        case = anno_exac('1-5994956-G-A')
        self.assertFalse(case['any_covered'])

    def test_anno_exac_bulk(self):
        vs = {
                '1-94496602-G-T': 0.02449307534291962,
                '1-94544234-T-C': 0.2552553542009885,
                '8-43051670-C-T': 0.0005546311702717693,
                }
        case1 = anno_exac_bulk(list(vs.keys()))
        case2 = anno_exac_bulk(list(vs.keys()), chunk_size=2)
        for k in vs:
            self.assertEqual(vs[k],case1[k]['variant']['allele_freq'])
            self.assertEqual(vs[k],case2[k]['variant']['allele_freq'])

    def test_anno_kaviar(self):
        vs = {
                '1-94496602-G-T': 0.02449,
                '1-94544234-T-C': 0.239336,
                '8-43051670-C-T': 0.000328,
                }
        case1 = anno_kaviar(list(vs.keys()))
        case2 = anno_kaviar(list(vs.keys()), chunk_size=2)
        for k in vs:
            self.assertEqual(vs[k],case1[k])
            self.assertEqual(vs[k],case2[k])

    def test_my_gene(self):
        gene_id = 'ENSG00000198691'
        case = my_gene(gene_id)
        self.assertEqual(case['symbol'], 'ABCA4')

    def test_my_genes(self):
        gene_ids = ['ENSG00000042781','ENSG00000198691']
        case = my_genes(gene_ids)
        self.assertEqual(case[0]['symbol'],'USH2A')
        self.assertEqual(case[1]['symbol'],'ABCA4')

    def test_my_genes_by_symbol(self):
        case = my_genes_by_symbol('ABCA4','human')
        self.assertEqual(case[0]['ensembl']['gene'],'ENSG00000198691')

    def test_obo_parser(self):
        obo_file = os.path.join('tests','data','hp.obo')
        case = obo_parser(obo_file)
        self.assertEqual(case['HP:0000556']['name'][0],'Retinal dystrophy')
    
    def test_check_ensemblId(self):
        good_gene = 'ENSG00000198691'
        case = check_ensemblId(good_gene)
        self.assertTrue(case)

if __name__ == '__main__':
    unittest.main()
