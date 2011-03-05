#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Unit test suite for SeaPy"""

__version__ = "0.0"
__author__ = "Steve Morley <smorley@lanl.gov>"


import unittest
import numpy.testing as ntest

try:
    import seapy
except ImportError:
    from spacepy import seapy

class SEATestsUniform(unittest.TestCase):
    """Tests of the sea method"""

    def setUp(self):
        super(SEATestsUniform, self).setUp()
        
	self.testval = 5
        unidata = [self.testval]*200
        time = range(200)
        epochs = [20,40,60,80,100,120,140,160,180]
        self.obj = seapy.Sea(unidata, time, epochs, verbose=False)
        self.obj.sea()

    def testMeanUniform(self):
        """Check superposed means on uniform input"""
        ntest.assert_array_equal(self.obj.semean, [self.testval]*(self.obj.window*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input"""
        ntest.assert_array_equal(self.obj.semedian, [self.testval]*(self.obj.window*2+1))

    def testMeanMedEquality(self):
        """For uniform input mean and median are same"""
        ntest.assert_array_equal(self.obj.semedian, self.obj.semean)



if __name__ == '__main__':
    unittest.main()
