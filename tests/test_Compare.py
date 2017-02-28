import unittest
from helper import *
import os

class CompareTestCase(unittest.TestCase):
    def setUp(self):
        print("Testing %s" % self._testMethodName)
        self.datapath = os.path.join('tests','data')
    def test_equal(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new1.xlsx'),'Sheet1','key',['field1','field5'])
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(case['<>'], dict())

    def test_change(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new2.xlsx'),'Sheet1','key',['field1','field5'])
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(sorted(case['<>'].keys()), [2,3,4,5])
        self.assertEqual(case['<>'][4]['change'], ['field1'])

    def test_delete(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new3.xlsx'),'Sheet1','key',['field1','field5'])
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([3]))
        self.assertEqual(case['<>'], dict())

    def test_add(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new4.xlsx'),'Sheet1','key',['field1','field5'])
        self.assertEqual(case['+'], set([6]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(case['<>'], dict())

    def test_add_delete(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new5.xlsx'),'Sheet1','key',['field1','field5'])
        self.assertEqual(case['+'], set([6]))
        self.assertEqual(case['-'], set([3]))
        self.assertEqual(case['<>'], dict())

if __name__ == '__main__':
    unittest.main()
