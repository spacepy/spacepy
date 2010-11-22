#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

import unittest
import spacepy.toolbox as tb
import glob, os, datetime, numpy
from numpy import array

class PickleAssembleTests(unittest.TestCase):
    
    def setUp(self):
        super(PickleAssembleTests, self).setUp()
                
        D1 = {}
        D1['names'] = ['John', 'Joe', 'Joyce']
        D1['TAI'] = [1,2,3]
        D2 = D1.copy()
        D2['TAI'] = [4,5,6]
        D3 = D1.copy()
        D3['TAI'] = [7,8,9]
        self.D1 = D1
        self.D2 = D2
        self.D3 = D3
        self.all = {'names':['John', 'Joe', 'Joyce', 'John', 'Joe', 'Joyce', 'John', 'Joe', 'Joyce'],
                    'TAI':[1,2,3,4,5,6,7,8,9]}
        
        try:  # make sure test file is gone before test
            os.remove('test_pickle_1.pkl')
            os.remove('test_pickle_2.pkl')
            os.remove('test_pickle_3.pkl')
        except:
            pass
         

    def tearDown(self):
        super(PickleAssembleTests, self).tearDown()
        try:  # make sure test file is gone before test
            os.remove('test_pickle_1.pkl')
            os.remove('test_pickle_2.pkl')
            os.remove('test_pickle_3.pkl')
            os.remove('test_all.pkl')
        except:
            pass


    def testSaveLoadPickle(self):
        """savePickle should write a pickle to disk and loadPickle should load it"""
        
        tb.savepickle('test_pickle_1.pkl', self.D1)
        files = glob.glob('*.pkl')
        self.assertTrue('test_pickle_1.pkl' in files)
        DD = tb.loadpickle('test_pickle_1.pkl')
        self.assertEqual(self.D1, DD)      

    def test_assemble(self):
        tb.savepickle('test_pickle_1.pkl', self.D1)
        tb.savepickle('test_pickle_2.pkl', self.D2)
        tb.savepickle('test_pickle_3.pkl', self.D3)
        expected = self.all
        result = tb.assemble(['test_pickle_1.pkl', 'test_pickle_2.pkl', 'test_pickle_3.pkl'] , 'test_all.pkl', sortkey=None)
        for key in result:
            result[key] = result[key].tolist()
        
        self.assertEqual(expected, result)



class SimpleFunctionTests(unittest.TestCase):
    def setUp(self):
        super(SimpleFunctionTests, self).setUp()

    def tearDown(self):
        super(SimpleFunctionTests, self).tearDown()

    def testfeq_equal(self):
        """feq should return true when they are equal"""
        val1 = 1.1234
        val2 = 1.1235
        self.assertTrue(tb.feq(val1, val2, 0.0001))

    def testfeq_notequal(self):
        """feq should return false when they are not equal"""
        val1 = 1.1234
        val2 = 1.1235
        self.assertFalse(tb.feq(val1, val2, 0.000005))

    def test_medAbsDev(self):
        """medAbsDev should return a known range for given random input"""
        data = numpy.random.normal(0, 1, 100000)
        real_ans = 0.7
        ans = tb.medAbsDev(data)
        self.assertAlmostEqual(ans, real_ans, places=1)
        
    def test_binHisto(self):
        """binHisto should return know answer for known input"""
        input = range(0, 101)
        real_ans = (21.47300748096567, 5.0)
        ans = tb.binHisto(input)
        self.assertEqual(ans, real_ans)

    def test_logspace(self):
        """logspace should return know answer for known input"""
        real_ans = array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])
        ans = tb.logspace(1, 100, 5)
        for i, val in enumerate(real_ans):
            self.assertAlmostEqual(val, ans[i], places=4)

    def test_arraybin(self):
        """arraybin should return know answer for known input"""
        real_ans = [(array([0, 1, 2, 3, 4]),), (array([5, 6, 7, 8, 9]),)]
        ans = tb.arraybin(numpy.arange(10), [4.2])
        for i in range(2):
            self.assertTrue( (ans[i][0] == real_ans[i][0]).all()) 
        

    def test_pmm(self):
        """pmm should give known output for known input"""
        data = [[1,3,5,2,5,6,2], array([5,9,23,24,6]), [6,23,12,67.34] ]
        real_ans = [[[1,6]], [[5, 24]], [[6, 67.34]]] 
        for i, val in enumerate(real_ans):
            self.assertEqual(val, tb.pmm(data[i]))
        self.assertEqual([[1, 6], [5, 24], [6.0, 67.340000000000003]], tb.pmm(*data) )

    def test_listUniq(self):
        """listUniq should give known output for known input"""
        data = [[1,2,3], [2,3,1], [1,1,1], [1,2,3,1]]
        real_ans = [[1,2,3], [2,3,1], [1], [1,2,3]]
        for i, val in enumerate(real_ans):
            self.assertEqual(val, tb.listUniq(data[i]))


    def test_leap_year(self):
        """Leap_year should give known output for known input"""
        data = ( 1993 + array(range(10)), 1900, [1993 + val for val in range(10)] )
        real_ans = ( array([False, False, False,  True, False, False, False,  True, False, False], dtype=bool),
                     False,
                     [False, False, False,  True, False, False, False,  True, False, False] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), tb.leap_year(data[i]).tolist())
            else:
                self.assertEqual(val, tb.leap_year(data[i]))
        real_ans = ( array([365, 365, 365,  366, 365, 365, 365,  366, 365, 365]),
                     365,
                     [365, 365, 365,  366, 365, 365, 365,  366, 365, 365] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), tb.leap_year(data[i], True).tolist())
            else:
                self.assertEqual(val, tb.leap_year(data[i], True))

class tFunctionTests(unittest.TestCase):
    def setUp(self):
        super(tFunctionTests, self).setUp()
        tdelta1 = [datetime.timedelta(hours=val) for val in range(100)]
        tdelta2 = [datetime.timedelta(hours=val) for val in range(-20, 20)]

        dt1 = datetime.datetime(2000, 11, 12)
        dt_a = [dt1 + val for val in tdelta1]
        dt_b = [dt1 + val for val in tdelta2]
        self.dt_a = dt_a
        self.dt_b = dt_b
        

    def tearDown(self):
        super(tFunctionTests, self).tearDown()
    
    def test_tOverlap(self):
        """tOverlap should return a known value for known input"""
        real_ans = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
               [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39])
        ans = tb.tOverlap(self.dt_a, self.dt_b)
        self.assertEqual(real_ans, ans)

    def test_tCommon(self):
        """tCommon should return a known value for known input"""
        real_ans = (array([ True,  True,  True,  True,  True,  True,  True,  True,  True,
                            True,  True,  True,  True,  True,  True,  True,  True,  True,
                            True,  True, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False,
                            False, False, False, False, False, False, False, False, False, False], dtype=bool),
                    array([False, False, False, False, False, False, False, False, False,
                           False, False, False, False, False, False, False, False, False,
                           False, False,  True,  True,  True,  True,  True,  True,  True,
                           True,  True,  True,  True,  True,  True,  True,  True,  True,
                           True,  True,  True,  True], dtype=bool))
        ans = tb.tCommon(self.dt_a, self.dt_b)
        self.assertEqual(real_ans[0].tolist(), ans[0].tolist())
        self.assertEqual(real_ans[1].tolist(), ans[1].tolist())



if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    ## suite = unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)


    unittest.main()





