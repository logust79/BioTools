import unittest
from Genes import *
import os
import sqlite3
import sys
sys.path.append('tests')
from helper import *

class GenesTestCase(unittest.TestCase):
    def setUp(self):
        mkdir_p('db')
        self.db = sqlite3.connect(os.path.join('db','test.db'))

    def test_gene(self):
        gene_id = 'ENSG00000198691'
        case = Gene(self.db, gene_id)
        self.assertEqual(case.symbol, 'ABCA4')
        self.assertEqual("{0:.3f}".format(case.pLI), '0.000')
        self.assertEqual(case.entrez_id, '24')
        self.assertEqual(case.mis_z, -1.498)
        self.assertEqual(case.genomic_pos_hg19['start'], 94458393)
        self.assertEqual(case.genomic_pos['start'], 93992835)
        self.assertTrue('STGD' in case.alias)

if __name__ == '__main__':
    unittest.main()
