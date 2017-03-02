import unittest
from helper import *
import os
import json
from numpy import NaN as nan

class CompareTestCase(unittest.TestCase):
    def setUp(self):
        print("Testing %s" % self._testMethodName)
        self.datapath = os.path.join('tests','data')
    def test_equal(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new1.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(case['<>'], dict())

    def test_changes(self):
        answer = json.dumps({"1": {"change": {"field5": {"to": 0.2, "from": 0.3}}}, "2": {"change": {"field5": {"to": 0.01, "from": NaN}}}, "4": {"change": {"field1": {"to": "d", "from": "a"}, "field5": {"to": 0.04, "from": 0.004}}}, "5": {"change": {"field1": {"to": "e", "from": "a"}, "field5": {"to": 0.0004, "from": 0.0003}}}})
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new2.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(json.dumps(case['<>']), answer)

    def test_change_methods(self):
        answer = json.dumps({"1": {"change": {"field4": {"to": "PASS", "from": "VQSR"}}}, "2": {"change": {"field5": {"to": 0.01, "from": NaN}}}, "4": {"change": {"field1": {"to": "d", "from": "a"}, "field5": {"to": 0.04, "from": 0.004}}}, "5": {"change": {"field1": {"to": "e", "from": "a"}, "field5": {"to": 0.0004, "from": 0.0003}}}})
        field5_closure = field5_cb_factory(0.01)
        case =  compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new2.xlsx'),'Sheet1','key',{'field1':None, 'field4':field4_cb, 'field5':field5_closure})
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(json.dumps(case['<>']), answer)

    def test_delete(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new3.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        self.assertEqual(case['+'], set([]))
        self.assertEqual(case['-'], set([3]))
        self.assertEqual(case['<>'], dict())

    def test_add(self):
        answer = json.dumps({"3": {"change": {"field5": {"to": 0.1, "from": NaN}}}})
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new4.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        self.assertEqual(case['+'], set([6]))
        self.assertEqual(case['-'], set([]))
        self.assertEqual(json.dumps(case['<>']), answer)

    def test_add_delete(self):
        case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_new5.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        self.assertEqual(case['+'], set([6]))
        self.assertEqual(case['-'], set([3]))
        self.assertEqual(case['<>'], dict())

    def test_duplicate_key_error(self):
        msg = "second array's key column's elements have duplicates"
        try:
            case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_bad1.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        except ValueError as e:
            self.assertEqual(e.message, msg)

    def test_field_error(self):
        msg = "field5"
        try:
            case = compare_excel(os.path.join(self.datapath,'compare_old.xlsx'),os.path.join(self.datapath,'compare_bad2.xlsx'),'Sheet1','key',{'field1':None,'field5':None})
        except KeyError as e:
            self.assertEqual(e.message, msg)

if __name__ == '__main__':
    unittest.main()
