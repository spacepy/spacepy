# -*- coding: utf-8 -*-

import unittest
# import spacepy.toolbox as tb
import spacepy.time as t
import glob
import os
import datetime
from numpy import array
import numpy

        

class tFunctionTests(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass
    
    def test_doy2date(self):
        """doy2date should return a known value for known input"""
        inval = ( (2000, 1),
                  (2001, 34),
                  (2006, 34),
                  (2008, 60),
                  (2008, 366) )
        real_ans = ( (1, 1),
                     (2, 3),
                     (2, 3),
                     (2, 29),
                     (12, 31) )
        for i, val in enumerate(inval):
            ans = t.doy2date(*val)
            self.assertEqual(real_ans[i] , ans)

    def test_tickrange(self):
        """tickrange should return a known value for known input"""
        inval = ( ('2002-02-01T00:00:00', '2002-02-04T00:00:00', 1),
                  ('2002-02-01T00:00:00', '2002-02-04T00:00:00', 0.5) )
        real_ans = ( ['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00', 
                      '2002-02-04T00:00:00'],
                     ['2002-02-01T00:00:00', '2002-02-01T12:00:00', '2002-02-02T00:00:00', 
                      '2002-02-02T12:00:00', '2002-02-03T00:00:00', '2002-02-03T12:00:00',
                      '2002-02-04T00:00:00'] )

        for i, val in enumerate(inval):
            ans = t.tickrange(*val)
            self.assertEqual(real_ans[i], ans.ISO)
                             

    def test_sec2hms(self):
        """sec2hms should return a known value for known input"""
        inval = ( (30, False, False),
                  (86401, False, False),
                  (86401, False, True),
                  (30.3, True, False),
                  (3599, False, False) )
        real_ans = ( [0, 0, 30],
                     [24, 0, 1],
                     [0, 0, 1],
                     [0, 0, 30],
                     [0, 59, 59] )
        for i, val in enumerate(inval):
            ans = t.sec2hms(*val)
            self.assertEqual(real_ans[i], ans)

class classTests(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_Tickdelta(self):
        """Tickdelta should function"""
        tst = t.Tickdelta(hours=12)
        self.assertEqual(43200.0, tst.seconds)

    def test_subTicktock(self):
        """a ticktock minus a ticktock is a tickdelta"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        n2 = t.Ticktock('2002-02-01T00:00:00', 'ISO')
        self.assertTrue(isinstance(n2 - n1, list))
        self.assertTrue(isinstance((n2 - n1)[0], t.Tickdelta))
        self.assertAlmostEqual(28.47443287, (n1-n2)[0].days, places=8)

    def test_subTickdelta(self):
        """a ticktock minus a Tickdelta is a ticktock"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        de = t.Tickdelta(hours=12, seconds=2)
        self.assertEqual(t.Ticktock( '2002-02-28T23:23:09', 'ISO'), n1-de)

    def test_TickTock_with_xrange(self):
        t0 = 1663236947
        iter_ex = range(t0, t0+5000, 500)
        range_ex = list(range(t0, t0+5000, 500)) 
        self.assertEqual(t.Ticktock(iter_ex, 'TAI').TAI, t.Ticktock(range_ex, 'TAI').TAI)
        
    def test_append(self):
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-03', '2002-01-04'])
        expected = t.Ticktock(['2002-01-01', '2002-01-02', '2002-01-03', '2002-01-04'])
        actual_1 = t1.append(t2)
        actual_2 = t1.append(t2.convert('UTC'))
        numpy.testing.assert_equal(expected.RDT, actual_1.RDT)
        numpy.testing.assert_equal(expected.RDT, actual_2.RDT)
        
if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    ## suite = unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)


    unittest.main()





