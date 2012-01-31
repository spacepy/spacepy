# -*- coding: utf-8 -*-
"""
Test suite for dates

As this is an exention module should have tests in C but I am too lazy

Copyright 2011 Los Alamos National Security, LLC.
"""
import unittest
import datetime

from matplotlib.dates import date2num as date2num_mpl
from matplotlib.dates import num2date as num2date_mpl
import numpy

import spacepy.time as spt

class date2numTests(unittest.TestCase):
    def setUp(self):
        super(date2numTests, self).setUp()
        self.dt = datetime.datetime(2000, 12, 13, 04, 54, 34)

    def tearDown(self):
        super(date2numTests, self).tearDown()

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

    def test_2darray(self):
        """2d arrys should come back n the same shape"""
        numpy.testing.assert_allclose(numpy.asarray(spt.date2num([self.dt]*10)).reshape((5,2)),
                                      date2num_mpl([self.dt]*10).reshape((5,2)))


class num2dateTests(unittest.TestCase):
    def setUp(self):
        super(num2dateTests, self).setUp()
        self.num = 734443.7143487677
        self.dt = datetime.datetime(2011, 11, 2, 17, 8, 39, 733525)

    def tearDown(self):
        super(num2dateTests, self).tearDown()
        
    def test_input(self):
        """There is some input checking"""
        self.assertRaises(TypeError, spt.num2date, 8675309, 8675309)
        self.assertRaises(ValueError, spt.num2date, 'bad in')
        self.assertRaises(ValueError, spt.num2date, self.dt)
        self.assertRaises(ValueError, spt.num2date, [8675309, 'bad']) 

    def test_mpl_same(self):
        """the C version should give the same result as the matplotlib version"""
        self.assertEqual(spt.num2date(self.num),
                         num2date_mpl(self.num).replace(tzinfo=None))

    def test_array(self):
        """output is an array if the input is an iterable"""
        numpy.testing.assert_array_equal(spt.num2date([self.num]*10), [val.replace(tzinfo=None) for val in num2date_mpl([self.num]*10)])
        numpy.testing.assert_array_equal(spt.num2date(numpy.asarray([self.num]*10)), [val.replace(tzinfo=None) for val in num2date_mpl([self.num]*10)])

    def test_2darray(self):
        """there was a segfault on 2d input, don't come back'"""
        numpy.testing.assert_array_equal(spt.num2date(numpy.asarray([self.num]*10).reshape((5,2))), numpy.asarray([val.replace(tzinfo=None) for val in num2date_mpl([self.num]*10)]).reshape((5,2)))

    def test_intypes(self):
        """Can input float, long, or int"""
        self.assertEqual(spt.num2date(int(self.num)), datetime.datetime(2011, 11, 2, 0, 0))
        self.assertEqual(spt.num2date(long(self.num)), datetime.datetime(2011, 11, 2, 0, 0))
        self.assertEqual(spt.num2date(float(long(self.num))), datetime.datetime(2011, 11, 2, 0, 0))
        self.assertEqual(spt.num2date([int(self.num)][0]), datetime.datetime(2011, 11, 2, 0, 0))
        self.assertEqual(spt.num2date([long(self.num)][0]), datetime.datetime(2011, 11, 2, 0, 0))
        self.assertEqual(spt.num2date([float(long(self.num))][0]), datetime.datetime(2011, 11, 2, 0, 0))


if __name__ == "__main__":
    unittest.main()
