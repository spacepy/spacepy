# -*- coding: utf-8 -*-
"""
Test suite for dates

As this is an exention module should have tests in C but I am too lazy

Copyright ©2011 Los Alamos National Security, LLC.
"""
import unittest
import datetime

from matplotlib.dates import date2num as date2num_mpl
import numpy

import spacepy.time as spt

class datesTests(unittest.TestCase):
    def setUp(self):
        super(datesTests, self).setUp()
        self.dt = datetime.datetime(2000, 12, 13, 04, 54, 34)

    def tearDown(self):
        super(datesTests, self).tearDown()

    def test_input(self):
        """There is some input checking"""
        self.assertRaises(TypeError, spt.date2num, self.dt, self.dt)
        self.assertRaises(ValueError, spt.date2num, 'bad in')
        self.assertRaises(ValueError, spt.date2num, 8675309)
        self.assertRaises(ValueError, spt.date2num, [self.dt, 1])        
    
    def test_mpl_same(self):
        """the C version should give the same result as the matplotlib version"""
        self.assertAlmostEqual(spt.date2num(self.dt), date2num_mpl(self.dt))

    def test_array(self):
        """output is an array if the input is an iterable"""
        numpy.testing.assert_allclose(spt.date2num([self.dt]*10), date2num_mpl([self.dt]*10))


if __name__ == "__main__":
    unittest.main()
