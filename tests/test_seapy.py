#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for SeaPy

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import datetime as dt
import unittest
import warnings
import numpy as np
import numpy.testing as ntest

import spacepy_testing
try:
    import seapy
except ImportError:
    from spacepy import seapy

__all__ = ['SEATestsUniform', 'SEATests2dUniform', 'SEATestsUniWithBad', 'SeaClassExceptions']


class SEATestsUniform(unittest.TestCase):
    """Tests of the sea method using uniform input"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATestsUniform, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        time = list(range(200))
        self.epochs = [20,40,60,80,100,120,140,160,180]
        with warnings.catch_warnings():
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea(self.unidata, time, self.epochs, verbose=False)
            self.obj.sea()

    def testMeanUniform(self):
        """Check superposed means on uniform input"""
        ntest.assert_array_equal(self.obj.semean, \
              [self.testval]*(int(self.obj.window)*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input"""
        ntest.assert_array_equal(self.obj.semedian, \
              [self.testval]*(int(self.obj.window)*2+1))

    def testMeanMedEquality(self):
        """For uniform input mean and median are same"""
        ntest.assert_array_equal(self.obj.semedian, self.obj.semean)

    def testDatetimeEquiv(self):
        """Test of equivalence of serial and datetime handling"""
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs, \
                                window=window, delta=delta,verbose=False)
            compobj.sea()
        ntest.assert_array_equal(self.obj.semedian, compobj.semedian)
        ntest.assert_array_equal(self.obj.semean, compobj.semean)

    def testSeaLen(self):
        """len should return number of epochs"""
        self.assertEqual(len(self.obj), len(self.epochs))

    def testRandomEpochsNoArgs(self):
        """Random epochs should have requested number"""
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea = self.obj.random()
        self.assertEqual(len(newsea), len(self.epochs))

    def testRandomEpochsArgs(self):
        """Random epochs should have requested number"""
        n_req = 27
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea = self.obj.random(n_req)
        self.assertEqual(len(newsea), n_req)

    def testRandomType(self):
        """Random epochs should have requested number"""
        newsea = self.obj.random()
        self.assertEqual(type(newsea), type(self.obj))

    def testRandomBeforeSea(self):
        n_req = 27
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs, \
                                window=window, delta=delta,verbose=False)
            newsea = compobj.random(n_req)
        self.assertEqual(len(newsea), n_req)

    def testRandomBoundType(self):
        """Test bound_type attriubte is set correctly """
        newsea = self.obj.random()
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea.sea(mad=True, quartiles=False)
        self.assertEqual(newsea.bound_type, 'mad')
        newsea = self.obj.random()
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea.sea(ci=True, ci_quan='mean')
        self.assertEqual(newsea.bound_type, 'ci')

    def testSeaCIFunc(self):
        '''Use Sea object to test mean and function passing'''
        time = list(range(200))
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj1 = seapy.Sea(self.unidata, time, self.epochs, \
                                window=window, delta=delta,verbose=False)
            compobj1.sea(ci=95, ci_quan=np.mean)
            compobj2 = seapy.Sea(self.unidata, time, self.epochs, \
                                window=window, delta=delta,verbose=False)
            compobj2.sea(ci=95, ci_quan='mean')
        ntest.assert_allclose(compobj1.bound_low, compobj2.bound_low)

    def testSeaDict(self):
        '''Test seadict grouping'''
        namelist = ['O1','O2']
        sd = seapy.seadict([self.obj, self.obj], namelist)
        for key in sd.keys():
            self.assertTrue(key in namelist)

    def testSeaDictFail(self):
        '''Test seadict fails with bad inputs'''
        namelist = ['O1']
        namelist2 = ['O2']
        with self.assertRaises(ValueError):
            sd = seapy.seadict([self.obj, self.obj], namelist)
        with self.assertRaises(ValueError):
            sd = seapy.seadict({'O1': self.obj, 'O2': self.obj}, namelist2)

    def testSeaPlotShowFalse(self):
        '''Test that plot method (show=False) returns axes'''
        ax = self.obj.plot(show=False)
        self.assertTrue(hasattr(ax, '_axes_class'))


class SEATestsUniWithBad(unittest.TestCase):
    """Tests of sea method's handling of badvals"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATestsUniWithBad, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        #insert badvals
        for ind in range(30,180,16):
            self.unidata[ind] = -99
        time = list(range(200))
        self.epochs = [20,40,60,80,100,120,140,160,180]
        self.obj = seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        self.obj.sea(badval=-99)

    def testMeanUniform(self):
        """Check superposed means on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semean, \
              [self.testval]*(int(self.obj.window)*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semedian, \
              [self.testval]*(int(self.obj.window)*2+1))

    def testMeanMedEquality(self):
        """For uniform input mean and median are same with bad data"""
        ntest.assert_array_equal(self.obj.semedian, self.obj.semean)

    def testDatetimeEquiv(self):
        """Test of equivalence of serial and datetime handling"""
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs, \
                                window=window, delta=delta,verbose=False)
        compobj.sea(badval=-99)

        ntest.assert_array_equal(self.obj.semedian, compobj.semedian)
        ntest.assert_array_equal(self.obj.semean, compobj.semean)

class SeaClassExceptions(unittest.TestCase):
    """Tests of the exception handling in Sea class"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SeaClassExceptions, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        self.time = list(range(200))
        self.epochs = [20,40,60,80,100,120,140,160,180]

    def testRestoreEpochs(self):
        """Check that restoreepochs fails with no bad epochs"""
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea(self.unidata, self.time, \
                       self.epochs, verbose=False)
        re_fun = self.obj.restoreepochs
        self.assertRaises(AttributeError, re_fun)

class SEATests2dUniform(unittest.TestCase):
    """Tests of the sea method using uniform input"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATests2dUniform, self).setUp()

        self.testval = 5
        self.unidata = np.ones([200,200])
        self.unidata.fill(self.testval)
        time = list(range(200))
        self.epochs = [20,40,60,80,100,120,140,160,180]
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea2d(self.unidata, time, self.epochs, verbose=False)
            self.obj.sea()

    def testSea2dLen(self):
        '''Test that object length is as expected'''
        self.assertEqual(len(self.obj), len(self.epochs))

    def testSea2dMean(self):
        '''Mean of semean attr should be same as uniform test val'''
        self.assertEqual(np.mean(self.obj.semean), self.testval)

    def test2dRandomType(self):
        """Random object should have same type as parent"""
        newsea = self.obj.random()
        self.assertEqual(type(newsea), type(self.obj))

if __name__ == '__main__':
    unittest.main()
