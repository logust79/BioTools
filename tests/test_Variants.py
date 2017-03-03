import unittest
from Variants import *
import os
import sqlite3
import sys
sys.path.append('tests')
from helper import *

class VariantsTestCase(unittest.TestCase):
    def setUp(self):
        mkdir_p('db')
        self.db = sqlite3.connect(os.path.join('db','test.db'))
    def test_variant(self):
        case = Variant(self.db, '13-95363811---GCG')
        self.assertEqual(parse_exac(case.exac), 0)
        case = Variant(self.db, '2-220400050-T-G')
        self.assertEqual(parse_exac(case.exac),2.4078979051288225e-05)
    def test_variants(self):
        vs = {
                '1-94496602-G-T': {
                    'exac_af': 0.02449307534291962,
                    'kaviar_af': 0.02449,
                    'cadd_phred': 10.03,
                },
                '1-94544234-T-C': {
                    'exac_af': 0.2552553542009885,
                    'kaviar_af': 0.239336,
                    'cadd_phred': 1.131,
                },
                '8-43051670-C-T': {
                    'exac_af': 0.0005546311702717693,
                    'kaviar_af': 0.000328,
                    'cadd_phred': 4.887,
                }
        }
        case = Variants(self.db, vs.keys())
        case.cadd_file = os.path.join('tests','data','cadd.tsv')
        for k,v in vs.items():
            self.assertEqual(parse_exac(case.exac[k]), v['exac_af'])
            self.assertEqual(case.kaviar_af[k], v['kaviar_af'])
            self.assertEqual(case.cadd_phred[k], v['cadd_phred'])

if __name__ == '__main__':
    unittest.main()
