import unittest
from CommonFuncs import *

class CommonFuncsTestCase(unittest.TestCase):
    def test_clean_variant(self):
        case = clean_variant('13-95363811---GCG')
        self.assertEqual(case,'13-95363810-G-GGCG')
        case = clean_variant('13-95363811---GCG', 'hg38')
        self.assertEqual(case,'13-95363810-A-AGCG')
        case = clean_variant('2-220400050-TGTGTGTGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTACGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTAC-GGTGTGTGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTACGTGTGTGTGGGGGGTGGCTGTGTGACTCTGTGTGTAC')
        self.assertEqual(case, '2-220400050-T-G')
        case = clean_variant('13-95363809-GG--')
        self.assertEqual(case, '13-95363808-CGG-C')

    def test_find_baseis(self):
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

if __name__ == '__main__':
    unittest.main()
