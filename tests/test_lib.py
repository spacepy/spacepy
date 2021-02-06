#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for libspacepy

Copyright 2014 Los Alamos National Security, LLC.
"""

import ctypes
import unittest
import warnings

import numpy
import numpy.testing

import spacepy_testing
import spacepy.lib


class SpacepyLibFindTests(unittest.TestCase):
    """Tests that libspacepy was built/exists"""

    def setUp(self):
        warnings.simplefilter('always')

    def testExists(self):
        """Make sure we're finding libspacepy"""
        self.assertTrue(spacepy.lib.have_libspacepy)

@unittest.skipUnless(spacepy.lib.have_libspacepy, 'libspacepy not found')
class SpacepyLibTests(unittest.TestCase):
    """Tests for libspacepy functions

    These are not deep functional tests. They're a chance to exercise the C
    code directly to make sure the library is loading, memory is safe, etc.
    The Python code which uses the C has much more robust test of output.
    """

    def setUp(self):
        warnings.simplefilter('always')

    def testHypotTB(self):
        """Test hypoteneuse"""
        self.assertEqual(5., spacepy.lib.hypot_tb(numpy.array([3., 4.]), 2))

    def testAACI(self):
        """Test confidence interval of association analysis"""
        dtype = numpy.dtype('uint' + str(ctypes.sizeof(ctypes.c_long) * 8))
        n_lags = 5
        n_p1 = 10
        n_surr = 8
        n_assoc = numpy.zeros(shape=(n_lags, n_p1), dtype=dtype)
        #Just constructing something that could be reasonable...
        #First five processes: single association where the lag index
        #is the p1 index
        for i in range(5):
            n_assoc[i, i] = 1
        #Next 5: one assoc where lag is p1-1, two where lag is p1
        for i in range(5, 10):
            n_assoc[i-5, i] = 2
            if i > 5:
                n_assoc[i-6, i] = 1
        n_assoc = n_assoc.reshape((-1,))
        surr_assoc_total = numpy.zeros(shape=(n_lags * n_surr,),
                                       dtype=dtype)
        #Force the seed
        seeds = numpy.zeros(shape=(n_lags,), dtype=dtype)
        spacepy.lib.aa_ci(n_assoc.ctypes.data_as(spacepy.lib.ulptr),
                          surr_assoc_total.ctypes.data_as(spacepy.lib.ulptr),
                          ctypes.c_ulong(n_lags), ctypes.c_ulong(n_p1),
                          ctypes.c_ulong(n_surr),
                          seeds.ctypes.data_as(spacepy.lib.ulptr),
                          ctypes.c_int(0))
        #This is just a regression test, really...
        expected = numpy.array([
            5, 2, 6, 1, 2, 6, 1, 4, 1, 9, 0, 3, 2, 3, 6, 2, 3, 9, 1, 4, 6, 1, 6,
            0, 4, 6, 6, 6, 3, 3, 5, 3, 3, 0, 5, 0, 5, 5, 4, 6
        ], dtype=dtype)
        numpy.testing.assert_array_equal(surr_assoc_total, expected)

    def testAssoc(self):
        """Test main association analysis"""
        np_long = numpy.dtype('int' + str(ctypes.sizeof(ctypes.c_long) * 8))
        np_double = numpy.dtype('float' +
                                str(ctypes.sizeof(ctypes.c_double) * 8))
        p2 = numpy.require(numpy.array([1, 2, 20, 40], dtype=np_double),
                           requirements='CAW')
        p1 = numpy.require(numpy.array([1, 5, 15], dtype=np_double),
                           requirements='CAW')
        lags = numpy.require(numpy.array([0, 2], dtype=np_double),
                             requirements='CAW')
        n_assoc = numpy.empty(shape=(len(lags), len(p1)), dtype=np_long)
        spacepy.lib.assoc(p2.ctypes.data_as(spacepy.lib.dptr),
                          p1.ctypes.data_as(spacepy.lib.dptr),
                          lags.ctypes.data_as(spacepy.lib.dptr),
                          n_assoc.ctypes.data_as(spacepy.lib.lptr),
                          1., len(p2), len(p1), len(lags))
        expected = numpy.array([[2, 0, 0],
                                [1, 0, 0]], dtype=np_long)
        numpy.testing.assert_array_equal(n_assoc, expected)

    def testEuler(self):
        """Test Euler tracing"""
        gridx = numpy.arange(10, dtype=ctypes.c_double)
        gridy = numpy.arange(10, dtype=ctypes.c_double)
        maxstep = 100
        outx = numpy.empty(shape=(maxstep,), dtype=ctypes.c_double)
        outy = numpy.empty(shape=(maxstep,), dtype=ctypes.c_double)
        fieldx = numpy.empty((10, 10), dtype=ctypes.c_double)
        fieldy = numpy.empty((10, 10), dtype=ctypes.c_double)
        #Constant velocity just out in x
        fieldx[...] = 1
        fieldy[...] = 0
        count = spacepy.lib.cEuler(10, 10, maxstep, 1, 3, 3,
                                   gridx, gridy, fieldx, fieldy, outx, outy)
        self.assertEqual(7, count)
        numpy.testing.assert_array_equal(
            outx[:count], [3, 4, 5, 6, 7, 8, 9])
        numpy.testing.assert_array_equal(
            outy[:count], [3, 3, 3, 3, 3, 3, 3])


if __name__ == '__main__':
    unittest.main()
