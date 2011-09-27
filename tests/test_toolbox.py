#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for toolbox

Copyright ©2010 Los Alamos National Security, LLC.
"""

import datetime
import glob
import itertools
import math
import os
import random
import unittest

import numpy
from numpy import array
from scipy import inf
import spacepy.toolbox as tb
import matplotlib.pyplot as plt
import spacepy.time as st
from matplotlib.text import Text

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

        # make sure test file is gone before test
        self.tearDown()

    def tearDown(self):
        super(PickleAssembleTests, self).tearDown()
        try:  # make sure test file is gone before test
            os.remove('test_pickle_1.pkl')
        except:
            pass
        try:
            os.remove('test_pickle_2.pkl')
        except:
            pass
        try:
            os.remove('test_pickle_3.pkl')
        except:
            pass
        try:
            os.remove('test_all.pkl')
        except:
            pass
        try:
            os.remove('test_pickle_1.pkl.gz')
        except:
            pass
        try:
            os.remove('test_pickle_2.pkl.gz')
        except:
            pass
        try:
            os.remove('test_pickle_3.pkl.gz')
        except:
            pass
        try:
            os.remove('test_all.pkl.gz')
        except:
            pass

    def testSaveLoadPickle(self):
        """savePickle should write a pickle to disk and loadPickle should load it"""
        tb.savepickle('test_pickle_1.pkl', self.D1)
        files = glob.glob('*.pkl')
        self.assertTrue('test_pickle_1.pkl' in files)
        DD = tb.loadpickle('test_pickle_1.pkl')
        self.assertEqual(self.D1, DD)

    def testSaveLoadPickleCompress(self):
        """savePickle should write a pickle to disk and loadPickle should load it (compressed)"""
        tb.savepickle('test_pickle_1.pkl', self.D1, compress=True)
        files = glob.glob('*.pkl.gz')
        self.assertTrue('test_pickle_1.pkl.gz' in files)
        DD = tb.loadpickle('test_pickle_1.pkl')
        self.assertEqual(self.D1, DD)
        # this doesn't but should also work
        # DD = tb.loadpickle('test_pickle_1.pkl.gz')
        # self.assertEqual(self.D1, DD)

    def test_assemble(self):
        tb.savepickle('test_pickle_1.pkl', self.D1)
        tb.savepickle('test_pickle_2.pkl', self.D2)
        tb.savepickle('test_pickle_3.pkl', self.D3)
        expected = self.all
        result = tb.assemble('test_pickle_[1-3].pkl', 'test_all.pkl', sortkey=None, verbose=False)
        for key in result:
            result[key] = result[key].tolist()
        self.assertEqual(expected, result)


class SimpleFunctionTests(unittest.TestCase):
    def test_mlt2rad(self):
        """mlt2rad should have known output for known input"""
        self.assertAlmostEqual(-2.8797932657906435, tb.mlt2rad(1))
        self.assertAlmostEqual(3.6651914291880918, tb.mlt2rad(26))
        val = tb.mlt2rad([1, 2])
        ans = [-2.8797932657906435, -2.6179938779914944]
        numpy.testing.assert_almost_equal(val, ans)
        self.assertAlmostEqual(0.26179938779914941, tb.mlt2rad(1, midnight=True))
        self.assertAlmostEqual(6.8067840827778854, tb.mlt2rad(26, midnight=True))
        val = tb.mlt2rad([1, 2], midnight=True)
        ans = [0.26179938779914941, 0.52359877559829882]
        numpy.testing.assert_almost_equal(val, ans)

    def test_interpol(self):
        """interpol should give known results"""
        ans = array([ 0.5,  1.5,  2.5,  3.5,  4.5])
        x = numpy.arange(10)
        y = numpy.arange(10)
        numpy.testing.assert_equal(ans, tb.interpol(numpy.arange(5)+0.5, x, y))
        # now test with baddata
        ans = numpy.ma.masked_array([0.5, 1.5, 2.5, 3.5, 4.5],
            mask = [False,  True,  True, False, False], fill_value = 1e+20)
        numpy.testing.assert_equal(ans, tb.interpol(numpy.arange(5)+0.5, x, y, baddata=2))
        # test wrap hour
        y = range(24)*2
        x = range(len(y))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 23.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 23.5], x, y, wrap='hour'))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 1.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 1.5], x, y)) # as a regression don't need wrap
        # test wrap lon
        y = range(360)*2
        x = range(len(y))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 359.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 359.5], x, y, wrap='lon'))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 10.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 370.5], x, y)) # as a regression don't need wrap
        # test wrap arb
        y = range(14)*2
        x = range(len(y))
        real_ans = [1.5, 10.5, 13.5]
        numpy.testing.assert_almost_equal(real_ans, tb.interpol([1.5, 10.5, 13.5], x, y, wrap=14).compressed())
        real_ans = [1.5, 10.5, 1.5]
        numpy.testing.assert_almost_equal(real_ans, tb.interpol([1.5, 10.5, 15.5], x, y)) # as a regression don't need wrap

    def test_normalize(self):
        """normalize should give known results"""
        self.assertEqual([0.0, 0.5, 1.0], tb.normalize([1,2,3]))
        numpy.testing.assert_equal(array([0.0, 0.5, 1.0]), tb.normalize(array([1,2,3])))

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
        numpy.testing.assert_almost_equal(tb.binHisto([100]*10), (3.3333333333333335, 3.0))
        numpy.testing.assert_almost_equal(tb.binHisto([100]), (1.0, 1.0))

    def test_logspace(self):
        """logspace should return know answer for known input"""
        real_ans = array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])
        numpy.testing.assert_almost_equal(real_ans, tb.logspace(1, 100, 5) , 4)
        t1 = datetime.datetime(2000, 1, 1)
        t2 = datetime.datetime(2000, 1, 2)
        real_ans = [datetime.datetime(1999, 12, 31, 23, 59, 59, 999989),
            datetime.datetime(2000, 1, 1, 2, 39, 59, 994207),
            datetime.datetime(2000, 1, 1, 5, 19, 59, 989762),
            datetime.datetime(2000, 1, 1, 7, 59, 59, 986897),
            datetime.datetime(2000, 1, 1, 10, 39, 59, 985369),
            datetime.datetime(2000, 1, 1, 13, 19, 59, 985431),
            datetime.datetime(2000, 1, 1, 15, 59, 59, 986830),
            datetime.datetime(2000, 1, 1, 18, 39, 59, 989808),
            datetime.datetime(2000, 1, 1, 21, 19, 59, 994124),
            datetime.datetime(2000, 1, 2, 0, 0, 0, 30)]
        ans = tb.logspace(t1, t2, 10)
        try:
            from matplotlib.dates import date2num
        except ImportError:
            pass
        else:
            numpy.testing.assert_almost_equal(date2num(real_ans), date2num(ans) , 4)

    def test_linspace(self):
        """linspace should return know answer for known input"""
        # should exactly match here since it is the same
        numpy.testing.assert_equal(tb.linspace(10, 100, 200), numpy.linspace(10,100, 200))
        t1 = datetime.datetime(2000, 1, 1)
        t2 = datetime.datetime(2000, 1, 10)
        real_ans = [datetime.datetime(2000, 1, 1, 0, 0),
             datetime.datetime(2000, 1, 3, 6, 0),
             datetime.datetime(2000, 1, 5, 12, 0),
             datetime.datetime(2000, 1, 7, 18, 0),
             datetime.datetime(2000, 1, 10, 0, 0)]
        ans = tb.linspace(t1, t2, 5)
        try:
            from matplotlib.dates import date2num
        except ImportError:
            pass
        else:
            numpy.testing.assert_almost_equal(date2num(real_ans), date2num(ans) , 4)

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
        leaps = [1600, 1604, 1608, 1612, 1616, 1620, 1624, 1628, 1632, 1636,
            1640, 1644, 1648, 1652, 1656, 1660, 1664, 1668, 1672, 1676,
            1680, 1684, 1688, 1692, 1696, 1704, 1708, 1712, 1716, 1720,
            1724, 1728, 1732, 1736, 1740, 1744, 1748, 1752, 1756, 1760,
            1764, 1768, 1772, 1776, 1780, 1784, 1788, 1792, 1796, 1804,
            1808, 1812, 1816, 1820, 1824, 1828, 1832, 1836, 1840, 1844,
            1848, 1852, 1856, 1860, 1864, 1868, 1872, 1876, 1880, 1884,
            1888, 1892, 1896, 1904, 1908, 1912, 1916, 1920, 1924, 1928,
            1932, 1936, 1940, 1944, 1948, 1952, 1956, 1960, 1964, 1968,
            1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008,
            2012, 2016, 2020, 2024, 2028, 2032, 2036, 2040]
        for i in range(1600, 2041):
            if i in leaps:
                self.assertTrue(tb.leap_year(i))
            else:
                self.assertFalse(tb.leap_year(i))

        data = ( 1993 + array(range(10)), 1900, [1993 + val for val in range(10)] )
        real_ans = ( array([365, 365, 365,  366, 365, 365, 365,  366, 365, 365]),
                     365,
                     [365, 365, 365,  366, 365, 365, 365,  366, 365, 365] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), tb.leap_year(data[i], True))
            else:
                self.assertEqual(val, tb.leap_year(data[i], True))
        real_ans = ( array([False, False, False,  True, False, False, False, True, False, False]),
                     False,
                     [False, False, False, True, False, False, False,  True, False, False] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), tb.leap_year(data[i], False))
            else:
                self.assertEqual(val, tb.leap_year(data[i], False))

    def test_rad2mlt(self):
        """rad2mlt should give known answers (regression)"""
        real_ans = array([ -6.        ,  -4.10526316,  -2.21052632,  -0.31578947,
         1.57894737,   3.47368421,   5.36842105,   7.26315789,
         9.15789474,  11.05263158,  12.94736842,  14.84210526,
        16.73684211,  18.63157895,  20.52631579,  22.42105263,
        24.31578947,  26.21052632,  28.10526316,  30.        ])
        numpy.testing.assert_almost_equal(real_ans, tb.rad2mlt(numpy.linspace(-numpy.pi*1.5, numpy.pi*1.5, 20)))
        real_ans = array([  6.        ,   7.89473684,   9.78947368,  11.68421053,
        13.57894737,  15.47368421,  17.36842105,  19.26315789,
        21.15789474,  23.05263158,  24.94736842,  26.84210526,
        28.73684211,  30.63157895,  32.52631579,  34.42105263,
        36.31578947,  38.21052632,  40.10526316,  42.        ])
        numpy.testing.assert_almost_equal(real_ans, tb.rad2mlt(numpy.linspace(-numpy.pi*1.5, numpy.pi*1.5, 20), midnight=True))


    def testIntSolve(self):
        """Find function input to give a desired integral value"""
        inputs = [[lambda x: x**2, 1, 0, 1000],
                  [lambda x: x / 2, 4, 0, 100],
                  [lambda x: math.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                          (5 * math.sqrt(2 * math.pi)), 0.6, -inf, inf],
                  ]
        outputs = [3.0 ** (1.0 / 3),
                   4,
                   1.266735515678999
                   ]
        for (input, output) in zip(inputs, outputs):
            self.assertAlmostEqual(output,
                                   tb.intsolve(*input),
                                   places=6)

    def testDistToList(self):
        """Convert probability distribution to list of values"""
        inputs = [[lambda x: x, 10, 0, 10],
                  [lambda x: math.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                   (5 * math.sqrt(2 * math.pi)), 20, -inf, inf],
                  ]
        outputs = [[2.2360679774998005,
                    3.8729833462074126,
                    5.0,
                    5.91607978309959,
                    6.708203932499401,
                    7.416198487095684,
                    8.06225774829855,
                    8.66025403784434,
                    9.219544457292841,
                    9.746794344808976],
                   [-9.799819946289062, -7.197657585144043,
                    -5.751746901776642, -4.672946454957128,
                    -3.777075132355094, -2.9888006299734116,
                    -2.268810950219631, -1.5931968204677105,
                    -0.9455921351909637, -0.31353388726711273,
                    0.3135339021682739, 0.9455921053886414,
                    1.5931968688964844, 2.268810987472534,
                    2.988800525665283, 3.7770752906799316,
                    4.672946453094482, 5.751747131347656,
                    7.197657346725464, 9.799819946289062]
                   ]
        for (input, output) in zip(inputs, outputs):
            numpy.testing.assert_almost_equal(output, tb.dist_to_list(*input))
        ans = [0.22360679689512963,
                0.3872983331400981,
                0.49999999808849793,
                0.5916079760482376,
                0.6708203906853853,
                0.7416198458743466,
                0.8062257717476484,
                0.8660254004736174,
                0.9219544422046511,
                0.9746794307546907]
        # test the "if max is None:" part
        numpy.testing.assert_almost_equal(ans, tb.dist_to_list(lambda x: x, 10, 0))

    def testBinCenterToEdges(self):
        """Convert a set of bin centers to bin edges"""
        inputs = [[1, 2, 3, 4, 5],
                  [12.4, 77, 100],
                  ]
        outputs = [[0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
                   [-19.9, 44.7, 88.5, 111.5],
                   ]
        for (input, output) in itertools.izip(inputs, outputs):
            self.assertEqual(output, tb.bin_center_to_edges(input))

    def test_hypot(self):
        """hypot should have known output"""
        invals = [ [3, 4], range(3, 6), range(3,10), [-1,2,3] ]
        ans = [ 5, 7.0710678118654755, 16.73320053068151, 3.7416573867739413 ]
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(*tst))
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(tst))
        self.assertEqual(5.0, tb.hypot(5.0))

    def testThreadJob(self):
        """Multithread the square of an array"""
        numpy.random.seed(0)
        a = numpy.random.randint(0, 100, [1000000])
        b = numpy.empty([1000000], dtype='int64')
        expected = a ** 2
        def targ(ina, outa, start, count):
            outa[start:start + count] = ina[start:start + count] ** 2
        tb.thread_job(len(a), 0, targ, a, b)
        self.assertEqual(list(expected), list(b))

    def testThreadMap(self):
        """Multithread summing a bunch of arrays"""
        numpy.random.seed(0)
        inputs = [numpy.random.randint(0, 100, [100000]) for i in range(100)]
        totals = tb.thread_map(numpy.sum, inputs)
        expected = [numpy.sum(i) for i in inputs]
        self.assertEqual(expected, totals)

    def test_dictree(self):
        """dictree has known output (None)"""
        a = {'a':1, 'b':2, 'c':{'aa':11, 'bb':22}}
        self.assertEqual(tb.dictree(a), None)
        self.assertEqual(tb.dictree(a, attrs=True), None)
        self.assertEqual(tb.dictree(a, verbose=True), None)
        self.assertRaises(TypeError, tb.dictree, 'bad')

    def test_geomspace(self):
        """geomspace shopuld give known output"""
        ans = [1, 10, 100, 1000]
        numpy.testing.assert_array_equal(tb.geomspace(1, 10, 1000), ans)
        ans = [1, 10.0, 100.0]
        numpy.testing.assert_array_equal(tb.geomspace(1, stop = 100, num=3), ans)
        ans = [1, 10, 100]
        numpy.testing.assert_array_equal(tb.geomspace(1, ratio = 10, num=3), ans)
        # there was a roiunding issue that this test catches
        ans = [1, 3.1622776601683795, 10.000000000000002]
        numpy.testing.assert_allclose(tb.geomspace(1, stop = 10, num=3), ans)

class tFunctionTests(unittest.TestCase):
    def setUp(self):
        super(tFunctionTests, self).setUp()
        dt1 = datetime.datetime(2000, 11, 12)
        self.dt_a = [dt1 + datetime.timedelta(hours=val)
                     for val in range(100)]
        self.dt_b = [dt1 + datetime.timedelta(hours=val)
                     for val in range(-20, 20)]
        self.dt_b2 = [dt1 + datetime.timedelta(hours=val)
                     for val in range(-20, -2)]

    def tearDown(self):
        super(tFunctionTests, self).tearDown()

    def test_tOverlap(self):
        """tOverlap should return a known value for known input"""
        real_ans = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                     11, 12, 13, 14, 15, 16, 17, 18, 19],
                    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                     30, 31, 32, 33, 34, 35, 36, 37, 38, 39])
        ans = tb.tOverlap(self.dt_a, self.dt_b)
        self.assertEqual(real_ans, ans)
        self.assertEqual( (None, None), tb.tOverlap(self.dt_a, self.dt_b2) )

    def test_tOverlap_random(self):
        """Shuffle input before calling tOverlap"""
        real_ans = ([1, 5, 6, 10, 15, 16, 18, 24, 29, 30, 43,
                     46, 47, 51, 53, 55, 56, 64, 67, 74],
                    [1, 2, 6, 7, 10, 12, 13, 14, 15, 17, 18,
                     19, 24, 27, 28, 30, 32, 35, 37, 38])
        random.seed(0)
        random.shuffle(self.dt_a)
        random.shuffle(self.dt_b)
        ans = tb.tOverlap(self.dt_a, self.dt_b)
        numpy.testing.assert_array_equal(real_ans, ans)

    def test_tOverlapHalf(self):
        """Get overlap of only one list"""
        real_ans = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                    30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
        ans = tb.tOverlapHalf(self.dt_a, self.dt_b)
        self.assertEqual(real_ans, ans)

    def test_tOverlapHalf_random(self):
        """Shuffle input before calling tOverlapHalf"""
        real_ans = [1, 2, 6, 7, 10, 12, 13, 14, 15, 17, 18,
                     19, 24, 27, 28, 30, 32, 35, 37, 38]
        random.seed(0)
        random.shuffle(self.dt_a)
        random.shuffle(self.dt_b)
        ans = tb.tOverlapHalf(self.dt_a, self.dt_b)
        self.assertEqual(real_ans, ans)

    def test_tOverlapSorted(self):
        """Exploit the sorting for a fast tOverlap"""
        real_ans = ([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                     11, 12, 13, 14, 15, 16, 17, 18, 19],
                    [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                     30, 31, 32, 33, 34, 35, 36, 37, 38, 39])
        ans = tb.tOverlap(self.dt_a, self.dt_b, presort=True)
        numpy.testing.assert_array_equal(real_ans, ans)

    def test_tOverlapHalfSorted(self):
        """Get overlap of only one list, exploiting the sort"""
        real_ans = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                    30, 31, 32, 33, 34, 35, 36, 37, 38, 39]
        ans = tb.tOverlapHalf(self.dt_a, self.dt_b, presort=True)
        numpy.testing.assert_array_equal(real_ans, ans)

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
        real_ans2 = ([datetime.datetime(2000, 11, 12, 0, 0),
            datetime.datetime(2000, 11, 12, 1, 0, ),
            datetime.datetime(2000, 11, 12, 2, 0, ),
            datetime.datetime(2000, 11, 12, 3, 0, ),
            datetime.datetime(2000, 11, 12, 4, 0, ),
            datetime.datetime(2000, 11, 12, 5, 0, ),
            datetime.datetime(2000, 11, 12, 6, 0, ),
            datetime.datetime(2000, 11, 12, 7, 0, ),
            datetime.datetime(2000, 11, 12, 8, 0, ),
            datetime.datetime(2000, 11, 12, 9, 0, ),
            datetime.datetime(2000, 11, 12, 10, 0, ),
            datetime.datetime(2000, 11, 12, 11, 0, ),
            datetime.datetime(2000, 11, 12, 12, 0, ),
            datetime.datetime(2000, 11, 12, 13, 0, ),
            datetime.datetime(2000, 11, 12, 14, 0, ),
            datetime.datetime(2000, 11, 12, 15, 0, ),
            datetime.datetime(2000, 11, 12, 16, 0, ),
            datetime.datetime(2000, 11, 12, 17, 0, ),
            datetime.datetime(2000, 11, 12, 18, 0, ),
            datetime.datetime(2000, 11, 12, 19, 0, )],
           [datetime.datetime(2000, 11, 12, 0, 0, ),
            datetime.datetime(2000, 11, 12, 1, 0, ),
            datetime.datetime(2000, 11, 12, 2, 0, ),
            datetime.datetime(2000, 11, 12, 3, 0, ),
            datetime.datetime(2000, 11, 12, 4, 0, ),
            datetime.datetime(2000, 11, 12, 5, 0, ),
            datetime.datetime(2000, 11, 12, 6, 0, ),
            datetime.datetime(2000, 11, 12, 7, 0, ),
            datetime.datetime(2000, 11, 12, 8, 0, ),
            datetime.datetime(2000, 11, 12, 9, 0, ),
            datetime.datetime(2000, 11, 12, 10, 0, ),
            datetime.datetime(2000, 11, 12, 11, 0, ),
            datetime.datetime(2000, 11, 12, 12, 0, ),
            datetime.datetime(2000, 11, 12, 13, 0, ),
            datetime.datetime(2000, 11, 12, 14, 0, ),
            datetime.datetime(2000, 11, 12, 15, 0, ),
            datetime.datetime(2000, 11, 12, 16, 0, ),
            datetime.datetime(2000, 11, 12, 17, 0, ),
            datetime.datetime(2000, 11, 12, 18, 0, ),
            datetime.datetime(2000, 11, 12, 19, 0, )])
        ans = tb.tCommon(self.dt_a, self.dt_b, mask_only=False)
        self.assertEqual(real_ans2[0], ans[0])
        self.assertEqual(real_ans2[1], ans[1])
        # test ts1 being an array
        ans = tb.tCommon(array(self.dt_a), self.dt_b, mask_only=False)
        numpy.testing.assert_equal(real_ans2[0], ans[0])
        numpy.testing.assert_equal(real_ans2[1], ans[1])
        # test ts2 being an array
        ans = tb.tCommon(self.dt_a, array(self.dt_b), mask_only=False)
        numpy.testing.assert_equal(real_ans2[0], ans[0])
        numpy.testing.assert_equal(real_ans2[1], ans[1])

    def test_smartTimeTicks(self):
        """smartTimeTicks should give known output (regression)"""
        # hits all the different cases
        # else
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 10), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%d %b', fmt.fmt)
        # elif nHours < 4:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 1), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nHours < 24:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 13), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nHours < 12:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 11), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # if nHours < 1:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 0, 30), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nHours < 48:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 2, 0, 30), 20)
        Mtick, mtick, fmt = tb.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)

    def test_windowMean(self):
        """windowMean should give known results (regression)"""
        wsize = datetime.timedelta(days=1)
        olap = datetime.timedelta(hours=12)
        data = [10, 20]*50
        time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
        outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
        od_ans = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]
        ot_ans = [datetime.datetime(2001, 1, 1, 12, 0),
                    datetime.datetime(2001, 1, 2, 0, 0),
                    datetime.datetime(2001, 1, 2, 12, 0),
                    datetime.datetime(2001, 1, 3, 0, 0),
                    datetime.datetime(2001, 1, 3, 12, 0),
                    datetime.datetime(2001, 1, 4, 0, 0),
                    datetime.datetime(2001, 1, 4, 12, 0)]
        numpy.testing.assert_allclose(od_ans, outdata)
        self.assertEqual(ot_ans, outtime)
        outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap)
        od_ans = [14.8, 14.8, 14.8, 14.8, 14.8, 14.8, 14.8]
        ot_ans = [datetime.datetime(2001, 1, 1, 12, 30),
                    datetime.datetime(2001, 1, 2, 0, 30),
                    datetime.datetime(2001, 1, 2, 12, 30),
                    datetime.datetime(2001, 1, 3, 0, 30),
                    datetime.datetime(2001, 1, 3, 12, 30),
                    datetime.datetime(2001, 1, 4, 0, 30),
                    datetime.datetime(2001, 1, 4, 12, 30)]
        numpy.testing.assert_allclose(od_ans, outdata)
        self.assertEqual(ot_ans, outtime)

        time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
        time[50:] = [val + datetime.timedelta(days=2) for val in time[50:]]
        outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
        od_ans = [ 15.,  15.,  15.,  15.,  15.,  numpy.nan,  numpy.nan,  15.,  15.,  15.,  15.]
        ot_ans = [datetime.datetime(2001, 1, 1, 12, 0),
                  datetime.datetime(2001, 1, 2, 0, 0),
                  datetime.datetime(2001, 1, 2, 12, 0),
                  datetime.datetime(2001, 1, 3, 0, 0),
                  datetime.datetime(2001, 1, 3, 12, 0),
                  datetime.datetime(2001, 1, 4, 0, 0),
                  datetime.datetime(2001, 1, 4, 12, 0),
                  datetime.datetime(2001, 1, 5, 0, 0),
                  datetime.datetime(2001, 1, 5, 12, 0),
                  datetime.datetime(2001, 1, 6, 0, 0),
                  datetime.datetime(2001, 1, 6, 12, 0)]
        numpy.testing.assert_allclose(od_ans, outdata)
        self.assertEqual(ot_ans, outtime)

        # now test the pointwise
        outdata, outtime = tb.windowMean(data, winsize=24, overlap=12)
        od_ans = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]
        ot_ans = [12.0, 24.0, 36.0, 48.0, 60.0, 72.0, 84.0]
        numpy.testing.assert_allclose(ot_ans, outtime)
        numpy.testing.assert_allclose(od_ans, outdata)
        # winsize tests
        outdata, outtime = tb.windowMean(data, winsize=24.6, overlap=12)
        od_ans, ot_ans = tb.windowMean(data, winsize=24.6, overlap=12)
        numpy.testing.assert_allclose(ot_ans, outtime)
        numpy.testing.assert_allclose(od_ans, outdata)
        outdata, outtime = tb.windowMean(data, winsize=0.4)
        od_ans, ot_ans = tb.windowMean(data, winsize=1.0)
        numpy.testing.assert_allclose(ot_ans, outtime)
        numpy.testing.assert_allclose(od_ans, outdata)
        outdata, outtime = tb.windowMean(data, winsize=1.0, overlap=2)
        od_ans, ot_ans = tb.windowMean(data, winsize=1.0, overlap=0)
        numpy.testing.assert_allclose(ot_ans, outtime)
        numpy.testing.assert_allclose(od_ans, outdata)

    def test_windowMeanInputs(self):
        """windowMean does some input checking (regression)"""
        wsize = datetime.timedelta(days=1)
        olap = datetime.timedelta(hours=12)
        data = [10, 20]*50
        time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
        self.assertRaises(ValueError, tb.windowMean, data[1:], time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
        self.assertRaises(TypeError, tb.windowMean, data, time, winsize='bad', overlap=olap, st_time=datetime.datetime(2001,1,1))
        self.assertRaises(TypeError, tb.windowMean, data, time, winsize=wsize, overlap='bad', st_time=datetime.datetime(2001,1,1))
        self.assertRaises(TypeError, tb.windowMean, data, time, overlap=olap, st_time=datetime.datetime(2001,1,1))
        olap = datetime.timedelta(days=2)
        self.assertRaises(ValueError, tb.windowMean, data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
        time = range(len(time))
        self.assertRaises(TypeError, tb.windowMean, data, time, overlap=olap, st_time=datetime.datetime(2001,1,1))


class ArrayBinTests(unittest.TestCase):
    """Tests for arraybin function"""

    def testNormalList(self):
        """Non-pathological cases, lists as input"""
        inputs = [[range(10), [4.2]],
                  [[5, 6, 7, 8, 9, 10, 11, 12], [7, 11]],
                  [[5, 6, 7, 8, 9, 10, 11, 12], [7, 11, 11.5, 13]],
                  ]
        outputs = [[[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]],
                   [[0, 1], [2, 3, 4, 5], [6, 7]],
                   [[0, 1], [2, 3, 4, 5], [6], [7], []],
                    ]
        for (input, output) in zip(inputs, outputs):
            self.assertEqual(output,
                             tb.arraybin(*input))

class PlottingTests(unittest.TestCase):
    """Tests for plotting functionality"""

    def test_applySmartTimeTicks(self):
        """applySmartTimeTicks should have known behaviour"""
        ticks = st.tickrange('2002-02-01T00:00:00', '2002-02-10T00:00:00', deltadays = 1)
        y = range(len(ticks))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        line = ax.plot(ticks.UTC, y)
        tb.applySmartTimeTicks(ax, ticks.UTC)
        plt.draw()
        # should not have moved the ticks
        real_ans = numpy.array([ 730882.,  730883.,  730884.,  730885.,  730886.,  730887.,
        730888.,  730889.,  730890.,  730891.])
        numpy.testing.assert_allclose(real_ans, ax.get_xticks())
        # should have named them 01 Feb, 02 Feb etc
        real_ans = ["Text(0,0,u'01 Feb')", "Text(0,0,u'02 Feb')", "Text(0,0,u'03 Feb')", "Text(0,0,u'04 Feb')",
                    "Text(0,0,u'05 Feb')", "Text(0,0,u'06 Feb')", "Text(0,0,u'07 Feb')", "Text(0,0,u'08 Feb')",
                    "Text(0,0,u'09 Feb')", "Text(0,0,u'10 Feb')"]
        ans = [str(ax.xaxis.get_majorticklabels()[i]) for i in range(len(ax.xaxis.get_majorticklabels()))]
        numpy.testing.assert_array_equal(real_ans, ans)
        plt.close()


if __name__ == "__main__":
    unittest.main()
