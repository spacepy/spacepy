#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for SeaPy

Copyright ©2010 Los Alamos National Security, LLC.
"""


import unittest
import datetime as dt
import numpy.testing as ntest

try:
    import seapy
except ImportError:
    from spacepy import seapy

class SEATestsUniform(unittest.TestCase):
    """Tests of the sea method using uniform input"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATestsUniform, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        time = range(200)
        self.epochs = [20,40,60,80,100,120,140,160,180]
        self.obj = seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        self.obj.sea()

    def testMeanUniform(self):
        """Check superposed means on uniform input"""
        ntest.assert_array_equal(self.obj.semean, \
              [self.testval]*(self.obj.window*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input"""
        ntest.assert_array_equal(self.obj.semedian, \
              [self.testval]*(self.obj.window*2+1))

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
        compobj = seapy.Sea(self.unidata, time, epochs, \
                   window=window, delta=delta,verbose=False)
        compobj.sea()

        ntest.assert_array_equal(self.obj.semedian, compobj.semedian)
        ntest.assert_array_equal(self.obj.semean, compobj.semean)

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
        time = range(200)
        self.epochs = [20,40,60,80,100,120,140,160,180] 
        self.obj = seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        self.obj.sea(badval=-99)

    def testMeanUniform(self):
        """Check superposed means on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semean, \
              [self.testval]*(self.obj.window*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semedian, \
              [self.testval]*(self.obj.window*2+1))

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
        self.time = range(200)
        self.epochs = [20,40,60,80,100,120,140,160,180]

    def testRestoreEpochs(self):
        """Check that restoreepochs fails with no bad epochs"""
        self.obj = seapy.Sea(self.unidata, self.time, \
                   self.epochs, verbose=False)
        re_fun = self.obj.restoreepochs
        self.assertRaises(AttributeError, re_fun)


if __name__ == '__main__':
    unittest.main()
