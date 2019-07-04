from helper import *
import unittest
from Genes import *
import os
import sqlite3
import sys
sys.path.append('tests')


class GenesTestCase(unittest.TestCase):
    def setUp(self):
        mkdir_p('db')
        self.db = sqlite3.connect(os.path.join('db', 'test.db'))

    def test_gene(self):
        gene_id = 'ENSG00000198691'
        case = Gene(self.db, gene_id)
        self.assertEqual(case.symbol, 'ABCA4')
        self.assertEqual("{0:.3f}".format(case.pLI), '0.000')
        self.assertEqual(case.entrez_id, '24')
        self.assertEqual("{0:.3f}".format(case.mis_z), '-1.498')
        self.assertEqual(case.genomic_pos_hg19['start'], 94458393)
        self.assertEqual(case.genomic_pos['start'], 93992834)
        self.assertTrue('STGD' in case.alias)

    def test_genes(self):
        gene_id = 'ENSG00000198691'
        gene_ids = ['ENSG00000042781', 'ENSG00000198691']
        case = Genes(self.db, gene_ids)
        self.assertEqual(case.symbol[gene_id], 'ABCA4')
        self.assertEqual("{0:.3f}".format(case.pLI[gene_id]), '0.000')
        self.assertEqual(case.entrez_id[gene_id], '24')
        self.assertEqual("{0:.3f}".format(case.mis_z[gene_id]), '-1.498')
        self.assertEqual(case.genomic_pos_hg19[gene_id]['start'], 94458393)
        self.assertEqual(case.genomic_pos[gene_id]['start'], 93992834)
        self.assertTrue('STGD' in case.alias[gene_id])

        self.assertEqual(case.entrezIds_to_ensemblIds(
            ['24', '7399'])['24'][0], 'ENSG00000198691')
        self.assertEqual(case.symbols_to_ensemblIds(
            ['ABCA4', 'USH2A'])['ABCA4'], 'ENSG00000198691')


if __name__ == '__main__':
    unittest.main()
