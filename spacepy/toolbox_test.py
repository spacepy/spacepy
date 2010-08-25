
import unittest
import toolbox as tb
import glob
import os
import datetime
from numpy import array
import numpy

class PickleTests(unittest.TestCase):

    def setUp(self):
        super(PickleTests, self).setUp()
        try:  # make sure test file is gone before test
            os.remove('test_pickle_1.pbin')
        except:
            pass
         

    def tearDown(self):
        super(PickleTests, self).tearDown()
        try:  # make sure test file is gone before test
            os.remove('test_pickle_1.pbin')
        except:
            pass


    def testSaveLoadPickle(self):
        """savePickle should write a pickle to disk and loadPickle should load it"""
        
        D = {}
        D['names'] = ['John', 'Joe', 'Joyce']
        D['TAI'] = [1,2,3]
        self.D = D
        tb.savepickle('test_pickle_1.pbin', D)
        files = glob.glob('*.pbin')
        self.assertTrue('test_pickle_1.pbin' in files)
        DD = tb.loadpickle('test_pickle_1.pbin')
        self.assertEqual(D, DD)


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
        self.assertAlmostEqual(ans, real_ans, 1)
        
    def test_binHisto(self):
        """binHisto should return know answer for known input"""
        input = range(100)
        real_ans = (21.328903431315652, 5.0)
        ans = tb.binHisto(input)
        self.assertEqual(ans, real_ans)

    def test_logspace(self):
        """logspace should return know answer for known input"""
        real_ans = array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])
        ans = tb.logspace(1, 100, 5)
        for i, val in enumerate(real_ans):
            self.assertAlmostEqual(val, ans[i], 4)

    def test_arraybin(self):
        """arraybin should return know answer for known input"""
        real_ans = [(array([0, 1, 2, 3, 4]),), (array([5, 6, 7, 8, 9]),)]
        ans = tb.arraybin(numpy.arange(10), [4.2])
        for i in range(2):
            self.assertTrue( (ans[i][0] == real_ans[i][0]).all()) 
        

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
	unittest.main()





