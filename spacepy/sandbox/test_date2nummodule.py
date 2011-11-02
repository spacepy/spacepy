# -*- coding: utf-8 -*-
"""
Test suite for date2num

As this is an exention module should have tests in C but I am too lazy

Copyright ©2011 Los Alamos National Security, LLC.
"""
import unittest
import datetime

from matplotlib.dates import date2num as date2num_mpl
import numpy

import date2num

class date2numTests(unittest.TestCase):
    def setUp(self):
        super(date2numTests, self).setUp()
        self.dt = datetime.datetime(2000, 12, 13, 04, 54, 34)

    def tearDown(self):
        super(date2numTests, self).tearDown()

    def test_input(self):
        """There is some input checking"""
        self.assertRaises(TypeError, date2num.date2num, self.dt, self.dt)
        self.assertRaises(ValueError, date2num.date2num, 'bad in')
        self.assertRaises(ValueError, date2num.date2num, 8675309)
    
    def test_mpl_same(self):
        """the C version should give the same result as the matplotlib version"""
        self.assertAlmostEqual(date2num.date2num(self.dt), date2num_mpl(self.dt))

    def test_array(self):
        """output is an array if the input is an iterable"""
        numpy.testing.assert_allclose(date2num.date2num([self.dt]*10), date2num_mpl([self.dt]*10))


if __name__ == "__main__":
    unittest.main()
