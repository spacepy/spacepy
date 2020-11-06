#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for toolbox

Copyright 2010-2014 Los Alamos National Security, LLC.
"""

import time, datetime
import glob, os, sys
import shutil
import random
import re
import tempfile
try:
    import StringIO
except:
    import io as StringIO
import unittest
import warnings
from contextlib import contextmanager
try:
    import __builtin__ as builtins #python 2.x
except ImportError:
    import builtins #python 3.x

import numpy
from numpy import array
from scipy import inf
from scipy.stats import poisson

import spacepy
import spacepy.toolbox as tb
import spacepy.lib

#Py3k compatibility renamings
try:
    xrange
except NameError:
    xrange = range

__all__ = ['PickleAssembleTests', 'SimpleFunctionTests', 'TBTimeFunctionTests',
           'ArrayBinTests']

@contextmanager
def mockRawInput3(mock):
    original_raw_input = builtins.input
    builtins.input = lambda: mock
    yield
    builtins.input = original_raw_input

@contextmanager
def mockRawInput2(mock):
    original_raw_input = builtins.raw_input
    builtins.raw_input = lambda: mock
    yield
    builtins.raw_input = original_raw_input

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
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        super(PickleAssembleTests, self).tearDown()
        try:
            shutil.rmtree(self.tempdir)
        except OSError:
            pass

    def testSaveLoadPickle(self):
        """savePickle should write a pickle to disk and loadPickle should load it"""
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_1.pkl'), self.D1)
        files = glob.glob(os.path.join(self.tempdir, '*.pkl'))
        self.assertTrue(os.path.join(self.tempdir, 'test_pickle_1.pkl') in files)
        DD = tb.loadpickle(os.path.join(self.tempdir, 'test_pickle_1.pkl'))
        self.assertEqual(self.D1, DD)

    def testSaveLoadPickleCompress(self):
        """savePickle should write a pickle to disk and loadPickle should load it (compressed)"""
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_1.pkl'), self.D1, compress=True)
        files = os.listdir(self.tempdir)
        self.assertTrue('test_pickle_1.pkl.gz' in files)
        self.assertFalse('test_pickle_1.pkl' in files)
        DD = tb.loadpickle(os.path.join(self.tempdir,'test_pickle_1.pkl'))
        self.assertEqual(self.D1, DD)
        DD = tb.loadpickle(os.path.join(self.tempdir,'test_pickle_1.pkl.gz'))
        self.assertEqual(self.D1, DD)
        # Save without specifying compression, make sure saves compressed
        # (because compressed file already exists)
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_1.pkl'), self.D1)
        files = os.listdir(self.tempdir)
        self.assertTrue('test_pickle_1.pkl.gz' in files)
        self.assertFalse('test_pickle_1.pkl' in files)

    def test_assemble(self):
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_1.pkl'), self.D1)
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_2.pkl'), self.D2)
        tb.savepickle(os.path.join(self.tempdir, 'test_pickle_3.pkl'), self.D3)
        expected = self.all
        result = tb.assemble(os.path.join(self.tempdir, 'test_pickle_[1-3].pkl'), os.path.join(self.tempdir, 'test_all.pkl'), sortkey=None, verbose=False)
        for key in result:
            result[key] = result[key].tolist()
        self.assertEqual(expected, result)


class SimpleFunctionTests(unittest.TestCase):

    def test_humansort(self):
        """human_sort should give known answers"""
        dat = ['1.10', '1.2', '1.3', '1.20']
        self.assertEqual(
            ['1.2', '1.3', '1.10', '1.20'],
            tb.human_sort(dat))
        dat = ['r1.txt', 'r10.txt', 'r2.txt']
        dat.sort() # Standard Python sort
        self.assertEqual(['r1.txt', 'r10.txt', 'r2.txt'], dat)
        self.assertEqual(['r1.txt', 'r2.txt', 'r10.txt'],
                         tb.human_sort(dat))
        dat = [5, 1, 3, -1]
        self.assertEqual([-1, 1, 3, 5], tb.human_sort(dat))

    def test_quaternionDeprecation(self):
        """Make sure deprecated quaternion functions work"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            tst = tb.quaternionNormalize([0.707, 0, 0.707, 0.2])
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual(
            'moved to spacepy.coordinates',
            str(w[0].message))
        numpy.testing.assert_array_almost_equal(
            [0.693,  0.,  0.693,  0.196], tst, decimal=2)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            tst = tb.quaternionRotateVector([0.7071, 0, 0, 0.7071],
                                            [0, 1, 0])
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual(
            'moved to spacepy.coordinates',
            str(w[0].message))
        numpy.testing.assert_array_almost_equal(
            [0, 0, 1], tst, decimal=5)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            tst = tb.quaternionMultiply([1., 0, 0, 0],
                                        [0., 0, 0, 1], scalarPos='first')
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual(
            'moved to spacepy.coordinates',
            str(w[0].message))
        numpy.testing.assert_array_equal([0, 0, 0, 1], tst)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            tst = tb.quaternionConjugate([.707, 0, .707, 0.2])
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual(
            'moved to spacepy.coordinates',
            str(w[0].message))
        numpy.testing.assert_array_equal([-.707, 0, -.707, 0.2], tst)

    def test_indsFromXrange(self):
        """indsFromXrange should have known result"""
        foo = xrange(23, 39)
        self.assertEqual([23, 39], tb.indsFromXrange(foo))
        foo = xrange(5)
        self.assertEqual([0, 5], tb.indsFromXrange(foo))

    def test_indsFromXrange_zerolen(self):
        """indsFromXrange with zero length range"""
        foo = xrange(20, 20)  # empty
        self.assertEqual([20, 20], tb.indsFromXrange(foo))

    def test_interweave(self):
        """interweave should have known result"""
        a = numpy.arange(5)
        b = numpy.arange(5, 10)
        numpy.testing.assert_equal(numpy.vstack((a,b)).reshape((-1,),order='F'),
                                   tb.interweave(a, b))
        numpy.testing.assert_equal(array([0, 5, 1, 6, 2, 7, 3, 8, 4, 9]),
                                   tb.interweave(a, b))

    def test_getNamedPath(self):
        """getNamedPath should have known result"""
        curloc = os.path.dirname(os.path.abspath(__file__))
        tmpdir = os.path.join(curloc, 'tmp', 'test1', 'test2')
        os.makedirs(tmpdir)
        ans = ['tests', 'tmp', 'test1']
        os.chdir(tmpdir)
        res = tb.getNamedPath('test1').split(os.path.sep)[-3:]
        self.assertEqual(ans[0], res[0][0:len(ans[0])])
        self.assertEqual(ans[1], res[1][0:len(ans[0])])
        self.assertEqual(res[0],
                         tb.getNamedPath(res[0]).split(os.path.sep)[-1])
        os.chdir(curloc)
        os.removedirs(tmpdir)

    def test_getNamedPath_badInput(self):
        """getNamedPath should return None if directory does not exist"""
        import string, random
        len_fn = 16 #should be exceedingly unlikely to exist...
        dname = ''.join(random.choice(string.ascii_uppercase + 
                        string.digits) for _ in range(len_fn))
        res = tb.getNamedPath(dname)
        self.assertTrue(res is None)

    def test_progressbar(self):
        """progressbar shouldhave a known output"""
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        self.assertEqual(tb.progressbar(0, 1, 100), None)
        self.assertEqual(tb.progressbar(100, 1, 100), None)
        result = output.getvalue()
        output.close()
        self.assertEqual(result, "\rDownload Progress ...0%"
                         "\rDownload Progress ...100%\n")
        sys.stdout = realstdout

    def test_progressbar_bigblock(self):
        """progressbar should not go over 100% with big blocks"""
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        try:
            for i in range(4):
                self.assertEqual(tb.progressbar(i, 100, 205), None)
            result = output.getvalue()
        finally:
            sys.stdout = realstdout
            output.close()
        self.assertEqual(
            result,
            "\rDownload Progress ...0%"
            "\rDownload Progress ...49%"
            "\rDownload Progress ...98%"
            "\rDownload Progress ...100%"
            "\n"
        )

    def test_query_yes_no(self):
        '''query_yes_no should return known strings for known input'''
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        majorVersion = sys.version_info[0]
        if  majorVersion>2:
            mockRaw = mockRawInput3
        else:
            mockRaw = mockRawInput2
        with mockRaw('y'):
            self.assertEqual(tb.query_yes_no('yes?'), 'yes')
        with mockRaw('n'):
            self.assertEqual(tb.query_yes_no('no?'), 'no')
        with mockRaw(''):
            self.assertEqual(tb.query_yes_no('no?', default='no'), 'no')
        output.close()
        sys.stdout = realstdout

    def test_query_yes_no_badDefault(self):
        '''query_yes_no should return error for bad args'''
        self.assertRaises(ValueError, tb.query_yes_no, '', default='bad')

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

    def test_interpol_baddata(self):
        """interpol should give known results in presence of fill values"""
        ans = array([ 0.5,  1.5,  2.5,  3.5,  4.5])
        x = numpy.arange(10)
        y = numpy.arange(10)
        numpy.testing.assert_equal(ans, tb.interpol(numpy.arange(5)+0.5, x, y))
        # now test with baddata
        ans = numpy.ma.masked_array([0.5, 1.5, 2.5, 3.5, 4.5],
            mask = [False,  True,  True, False, False], fill_value = 1e+20)
        numpy.testing.assert_equal(ans, tb.interpol(numpy.arange(5)+0.5, x, y, baddata=2))
        #test with baddata at end of array
        ans = array([1.0, 9.0])
        numpy.testing.assert_equal(ans, tb.interpol([-1,12], x, y, baddata=0))

    def test_interpol_keywords(self):
        """Interpol should give known results with hour and lon keyword wrapping"""
        # test wrap hour
        y = list(range(24))*2
        x = list(range(len(y)))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 23.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 23.5], x, y, wrap='hour'))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 1.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 1.5], x, y)) # as a regression don't need wrap
        # test wrap lon
        y = list(range(360))*2
        x = list(range(len(y)))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 359.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 359.5], x, y, wrap='lon'))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 10.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 370.5], x, y)) # as a regression don't need wrap

    def test_interpol_arb(self):
        """Interpol should give known results for arbitrary float/int wrapping"""
        # test wrap arb
        y = list(range(14))*2
        x = list(range(len(y)))
        real_ans = [1.5, 10.5, 13.5]
        numpy.testing.assert_almost_equal(real_ans, tb.interpol([1.5, 10.5, 13.5], x, y, wrap=14).compressed())
        real_ans = [1.5, 10.5, 1.5]
        numpy.testing.assert_almost_equal(real_ans, tb.interpol([1.5, 10.5, 15.5], x, y)) # as a regression don't need wrap
        # and test wrap with a float
        real_ans = [1.5, 10.5, 13.5]
        numpy.testing.assert_almost_equal(real_ans, tb.interpol([1.5, 10.5, 13.5], x, y, wrap=14.0).compressed())
        # test wrap lon using 360.0 as input
        y = list(range(360))*2
        x = list(range(len(y)))
        real_ans = numpy.ma.masked_array([1.5, 10.5, 359.5],
            mask = False, fill_value = 1e+20)
        numpy.testing.assert_equal(real_ans, tb.interpol([1.5, 10.5, 359.5], x, y, wrap=360.0))

    def test_normalize(self):
        """normalize should give known results, default range"""
        numpy.testing.assert_array_almost_equal(array([0.0, 0.5, 1.0]), tb.normalize(array([1,2,3])))

    def test_normalize_small(self):
        """normalize should give known results, smaller range"""
        numpy.testing.assert_array_almost_equal(array([0.1 , 0.45, 0.8 ]), tb.normalize(array([1,2,3]), low=0.1, high=0.8))

    def test_normalize_large(self):
        """normalize should give known results, larger range"""
        numpy.testing.assert_array_almost_equal(array([1.0, 6.0, 11.0]), tb.normalize(array([1,2,3]), low=1, high=11))

    def test_normalize_nan(self):
        """normalize should give known results, larger range with nan"""
        numpy.testing.assert_array_almost_equal(array([1.0, 6.0, 11.0,numpy.nan]), tb.normalize(array([1,2,3,numpy.nan]),
                                                                                   low=1, high=11))
        numpy.testing.assert_array_almost_equal(array([1.0, 6.0, numpy.nan, 11.0]), tb.normalize(array([1,2,numpy.nan, 3,]),
                                                                                    low=1, high=11))

    def testfeq_equal(self):
        """feq should return true when they are equal"""
        val1 = 1.1234
        val2 = 1.1235
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            self.assertTrue(tb.feq(val1, val2, 0.0001))
            numpy.testing.assert_array_equal(
                [False, True, False, False],
                tb.feq([1., 2., 3., 4.],
                       [1.25, 2.05, 2.2, 500.1],
                       0.1)
            )
        self.assertEqual(2, len(w))
        for this_w in w:
            self.assertEqual(DeprecationWarning, this_w.category)
            self.assertEqual(DeprecationWarning, this_w.category)
            self.assertEqual('use numpy.isclose', str(this_w.message))

    def testfeq_notequal(self):
        """feq should return false when they are not equal"""
        val1 = 1.1234
        val2 = 1.1235
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always', category=DeprecationWarning)
            self.assertFalse(tb.feq(val1, val2, 0.000005))
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual('use numpy.isclose', str(w[0].message))

    def test_medAbsDev(self):
        """medAbsDev should return a known range for given random input"""
        data = numpy.random.normal(0, 1, 100000)
        real_ans = 0.7
        ans = tb.medAbsDev(data)
        self.assertAlmostEqual(ans, real_ans, places=1)

    def test_medAbsDev_scale(self):
        """medAbsDev should return a known range for given random input"""
        data = numpy.random.normal(0, 1, 100000)
        real_ans = 0.7*1.4826
        ans = tb.medAbsDev(data, scale=True)
        self.assertAlmostEqual(ans, real_ans, places=1)

    def test_binHisto(self):
        """binHisto should return know answer for known input"""
        input = list(range(0, 101))
        real_ans = (21.47300748096567, 5.0)
        ans = tb.binHisto(input)
        self.assertEqual(ans, real_ans)
        numpy.testing.assert_almost_equal(tb.binHisto([100]*10), (3.3333333333333335, 3.0))
        numpy.testing.assert_almost_equal(tb.binHisto([100]), (1.0, 1.0))
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        numpy.testing.assert_almost_equal(tb.binHisto([100], verbose=True), (1.0, 1.0))
        result = output.getvalue()
        output.close()
        self.assertEqual(result, "Used sqrt rule\n")
        sys.stdout = realstdout
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        numpy.testing.assert_almost_equal(tb.binHisto([90, 100]*10, verbose=True), (7.3680629972807736, 1.0))
        result = output.getvalue()
        output.close()
        self.assertEqual(result, "Used F-D rule\n")
        sys.stdout = realstdout

    def testBootHisto(self):
        """Bootstrap histogram known output for known input"""
        numpy.random.seed(28420)
        data = numpy.random.randn(1000)
        bin_edges, ci_low, ci_high, sample = spacepy.toolbox.bootHisto(
            data, n=1000, seed=28420)
        numpy.testing.assert_allclose(
            [-3.17371804, -2.40693385, -1.64014966, -0.87336547,
             -0.10658127, 0.66020292,  1.42698711,  2.1937713 ,
             2.9605555 ,  3.72733969, 4.49412388],
            bin_edges, rtol=1e-5)
        numpy.testing.assert_equal(
            [  7,  45, 149, 283, 272, 162,  71,   9,   0,   2], sample)
        #This is a very coarse chunk-check to allow for variations in
        #RNGs; by using the same seed it should be deterministic on a
        #given RNG. Values from 100000 iterations
        numpy.testing.assert_allclose(
            [3.,  35., 131., 260., 249., 143.,  58.,   4.,   0.,   0.],
            ci_low, atol=2, rtol=1e-2)
        numpy.testing.assert_allclose(
            [12.,  56., 168., 307., 295., 181.,  85.,  14.,   0.,   5.],
            ci_high, atol=2, rtol=1e-2)

    def test_logspace(self):
        """logspace should return know answer for known input"""
        try:
            from matplotlib.dates import date2num
        except ImportError: # just pass if matplotlib is not installed
            return
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
        numpy.testing.assert_almost_equal(date2num(real_ans), date2num(ans) , 4)

    def test_linspace(self):
        """linspace should return know answer for known input"""
        try:
            from matplotlib.dates import date2num
        except ImportError:  # just pass if matplotlib is not installed
            return
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
        numpy.testing.assert_almost_equal(date2num(real_ans), date2num(ans) , 4)

    def test_linspace_bug(self):
        """This catches a linspace datetime bug with 0-d arrays (regression)"""
        try:
            from matplotlib.dates import date2num
        except ImportError:  # just pass if matplotlib is not installed
            return
        t1 = numpy.array(datetime.datetime(2000, 1, 1))
        t2 = numpy.array(datetime.datetime(2000, 1, 10))
        real_ans = [datetime.datetime(2000, 1, 1, 0, 0),
             datetime.datetime(2000, 1, 3, 6, 0),
             datetime.datetime(2000, 1, 5, 12, 0),
             datetime.datetime(2000, 1, 7, 18, 0),
             datetime.datetime(2000, 1, 10, 0, 0)]
        ans = tb.linspace(t1, t2, 5)
        numpy.testing.assert_almost_equal(date2num(real_ans), date2num(ans) , 4)

    def test_pmm(self):
        """pmm should give known output for known input"""
        data = [[1,3,5,2,5,6,2], array([5,9,23,24,6]), [6,23,12,67.34] ]
        real_ans = [[[1,6]], [[5, 24]], [[6, 67.34]]]
        for i, val in enumerate(real_ans):
            self.assertEqual(val, tb.pmm(data[i]))
        self.assertEqual([[1, 6], [5, 24], [6.0, 67.340000000000003]], tb.pmm(*data) )

    def test_pmm_object(self):
        """Test pmm handling of object arrays"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            data = [array([5,9,23,24,6]).astype(object),
                    [datetime.datetime(2000, 3, 1, 0, 1), datetime.datetime(2000, 2, 28), datetime.datetime(2000, 3, 1)],
                    numpy.array(['foo', 'bar', 'baz'], dtype=object),
                    ]
            real_ans = [[[5, 24]],
                        [[datetime.datetime(2000, 2, 28), datetime.datetime(2000, 3, 1, 0, 1)]],
                        [['bar', 'foo']],
            ]
            for i, val in enumerate(real_ans):
                self.assertEqual(val, tb.pmm(data[i]))

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
                  [lambda x: numpy.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                          (5 * numpy.sqrt(2 * numpy.pi)), 0.6, -inf, inf],
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
                  [lambda x: numpy.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                   (5 * numpy.sqrt(2 * numpy.pi)), 20],
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

    def testBinCenterToEdges(self):
        """Convert a set of bin centers to bin edges"""
        inputs = [[1, 2, 3, 4, 5],
                  [12.4, 77, 100],]
        outputs = [[0.5, 1.5, 2.5, 3.5, 4.5, 5.5],
                   [-19.9, 44.7, 88.5, 111.5],]
        for i, val in enumerate(inputs):
            numpy.testing.assert_almost_equal(outputs[i], tb.bin_center_to_edges(val))

    def testBinEdgesToCenterToEdges_datetimeroundtrip(self):
        """Convert a set of datetime bin edges to centers and back to bin edges"""
        inputs = [datetime.datetime(2012,9,3,n) for n in range(10)]
        computedOut = tb.bin_edges_to_center(inputs)
        roundtripResults = tb.bin_center_to_edges(computedOut)
        self.assertTrue((inputs == roundtripResults).all())

    def testBinEdgesToCenter(self):
        """Convert a set of bin edges to bin centers"""
        inputs = [[1, 2, 3, 4, 5],
                  [1,2,3,7,10,20],
                  ]
        outputs = [[1.5, 2.5, 3.5, 4.5],
                   [1.5, 2.5, 5, 8.5, 15],]
        for i, val in enumerate(inputs):
            numpy.testing.assert_almost_equal(outputs[i], tb.bin_edges_to_center(val))

    def test_hypot(self):
        """hypot should have known output, call Python or C (as default)"""
        invals = [ [3, 4], list(range(3, 6)), list(range(3,10)), [-1,2,3] ]
        ans = [ 5, 7.0710678118654755, 16.73320053068151, 3.7416573867739413 ]
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(*tst))
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(tst))
        self.assertEqual(5.0, tb.hypot(5.0))
        # Same as first set of inputs but explicitly np arrays
        invals = [numpy.asarray([3, 4]), numpy.array([3., 4.]),
                  numpy.arange(3, 6), numpy.arange(3, 10),
                  numpy.asarray([-1, 2, 3])]
        ans = [5, 5.,
               7.0710678118654755, 16.73320053068151,
               3.7416573867739413]
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(*tst))
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(tst))
        # Few tests of array subclasses
        invals = [spacepy.dmarray([3, 4]), spacepy.dmarray([3., 4.]),
                  spacepy.dmarray([-1, 2, 3])]
        ans = [5, 5.,
               3.7416573867739413]
        for i, tst in enumerate(invals):
            self.assertAlmostEqual(ans[i], tb.hypot(*tst))

    @unittest.skipUnless(spacepy.lib.have_libspacepy, 'libspacepy not found')
    def test_hypot_python(self):
        """Explicitly call Python version of hypot if had C before"""
        try:
            spacepy.lib.have_libspacepy = False
            self.test_hypot() # Repeat tests without the C code
        finally:
            spacepy.lib.have_libspacepy = True

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
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        self.assertEqual(tb.dictree(a), None)
        self.assertEqual(tb.dictree(a, attrs=True), None)
        self.assertEqual(tb.dictree(a, verbose=True), None)
        sys.stdout = realstdout
        result = output.getvalue()
        output.close()
        self.assertRaises(TypeError, tb.dictree, 'bad')
        expected = """+
|____a
|____b
|____c
     |____aa
     |____bb
+
|____a
|____b
|____c
     |____aa
     |____bb
+
|____a (int)
|____b (int)
|____c (dict [2])
     |____aa (int)
     |____bb (int)
"""
        self.assertEqual(expected, result)

    def test_geomspace(self):
        """geomspace should give known output"""
        ans = [1, 10, 100, 1000]
        numpy.testing.assert_array_equal(tb.geomspace(1, 10, 1000), ans)
        ans = [1, 10.0, 100.0]
        numpy.testing.assert_almost_equal(tb.geomspace(1, stop = 100, num=3), ans)
        ans = [1, 10, 100]
        numpy.testing.assert_array_equal(tb.geomspace(1, ratio = 10, num=3), ans)
        # there was a rounding issue that this test catches
        ans = [1, 3.1622776601683795, 10.000000000000002]
        numpy.testing.assert_almost_equal(tb.geomspace(1, stop = 10, num=3), ans)

    def test_isview(self):
        """isview should have known output"""
        a = numpy.arange(100)
        b = a[0:10]
        self.assertTrue(tb.isview(b))
        self.assertFalse(tb.isview(a))
        self.assertFalse(tb.isview([1, 2, 3]))
        numpy.testing.assert_array_equal(tb.isview(a, [1,2,3]), [False, False]) # a bit of a pathological case
        numpy.testing.assert_array_equal(tb.isview(b, a), [True, True])
        numpy.testing.assert_array_equal(tb.isview(b, a[2:]), [True, False])
        numpy.testing.assert_array_equal(tb.isview(a, a), [False, False])
        numpy.testing.assert_array_equal(tb.isview(a[:], a), [True, True])
        numpy.testing.assert_array_equal(tb.isview([1,2,3], 4), [False, False])

    def test_do_with_timeout(self):
        """Check for timeout"""
        def testfunc(x):
            time.sleep(1)
            return x + 1
        self.assertEqual(6, tb.do_with_timeout(2.0, testfunc, 5))
        self.assertRaises(tb.TimeoutError, tb.do_with_timeout,
                          0.5, testfunc, 5)

    def test_do_with_timeout_exception(self):
        """Check for timeout"""
        def testfunc(x):
            foo = ['hi', 'there']
            return foo[x]
        self.assertRaises(IndexError, tb.do_with_timeout,
                          0.5, testfunc, 5)

    def test_timeout_check_call(self):
        """Make sure check_call replacement handles timeout"""
        #This test definitely doesn't work on Windows, and the function
        #itself probably doesn't, either
        if sys.platform != 'win32':
            self.assertEqual(0, tb.timeout_check_call(10.0, 'sleep 2',
                                                      shell=True))
            self.assertRaises(tb.TimeoutError, tb.timeout_check_call,
                              1.0, 'sleep 5', shell=True)

    def test_unique_columns(self):
        """Make sure unique_columns gives the expected answers"""
        a = numpy.array([[1, 1, 1, 0, 0, 0],
                            [0, 1, 1, 1, 0, 0],
                            [0, 1, 1, 1, 0, 0],
                            [1, 1, 1, 0, 0, 0],
                            [1, 1, 1, 1, 1, 0]])
        ans1 = numpy.array([[0, 1, 1, 1, 0, 0],
                            [1, 1, 1, 0, 0, 0],
                            [1, 1, 1, 1, 1, 0]])
        ans0 = numpy.array([[0, 0, 0, 0, 0],
                            [0, 0, 0, 0, 1],
                            [0, 1, 1, 0, 1],
                            [1, 0, 0, 1, 1],
                            [1, 1, 1, 1, 1]])
        numpy.testing.assert_array_equal(ans1, tb.unique_columns(a, axis=1))
        numpy.testing.assert_array_equal(ans0, tb.unique_columns(a, axis=0))

    def test_poisson_fit(self):
        """Make sure that we get the right Poisson fit answer"""
        numpy.random.seed(8675309)
        ans = 20
        data = poisson.rvs(ans, size=1000)
        res = tb.poisson_fit(data)
        self.assertEqual(ans, numpy.round(res.x))

    def test_assemble_qindenton_daily(self):
        """Assemble OMNI data structure from Qin-Denton daily files"""
        dailydir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data', 'qindenton_daily')
        omnidata = tb._assemble_qindenton_daily(dailydir)
        self.assertEqual(
            ['ByIMF', 'Bz1', 'Bz2', 'Bz3', 'Bz4', 'Bz5', 'Bz6', 'BzIMF',
             'DOY', 'Dst', 'G1', 'G2', 'G3', 'Kp', 'Pdyn', 'Qbits',
             'RDT', 'UTC', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'akp3',
             'dens', 'velo'],
            sorted(omnidata.keys()))
        numpy.testing.assert_array_equal(
            [14,  16,  18,  17,  17,  12,  10,  11,  11,   4,   0,  -8, -13,
             -13, -12, -11,  -9,  -5,   2,   6,   9,   9,  10,   6,  -2,  -8,
             -11, -11, -12, -20, -27, -29, -26, -22, -17, -22, -25, -29, -33,
             -36, -41, -47, -47, -39, -35, -41, -42, -45],
            omnidata['Dst'])
        numpy.testing.assert_array_equal(
            2, omnidata['Qbits']['W6'])


class TBTimeFunctionTests(unittest.TestCase):
    def setUp(self):
        super(TBTimeFunctionTests, self).setUp()
        dt1 = datetime.datetime(2000, 11, 12)
        self.dt_a = [dt1 + datetime.timedelta(hours=val)
                     for val in range(100)]
        self.dt_b = [dt1 + datetime.timedelta(hours=val)
                     for val in range(-20, 20)]
        self.dt_b2 = [dt1 + datetime.timedelta(hours=val)
                     for val in range(-20, -2)]

    def tearDown(self):
        super(TBTimeFunctionTests, self).tearDown()

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
        random.shuffle(self.dt_a, lambda:round(random.random(), 9))
        random.shuffle(self.dt_b, lambda:round(random.random(), 9))
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
        #Cut down the python3 rng to the same precision as python2
        random.shuffle(self.dt_a, lambda:round(random.random(), 9))
        random.shuffle(self.dt_b, lambda:round(random.random(), 9))
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

    def test_eventTimer(self):
        """eventTimer should behave in a known way"""
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        t1 = time.time()
        time.sleep(0.25)
        t2 = tb.eventTimer('', t1)
        sys.stdout = realstdout
        result = output.getvalue()
        output.close()
        #There may be some overhead that pushes this past 0.25, but should
        #never be less
        self.assertTrue(0.25 <= t2 - t1 < 0.28)
        m = re.match(r"^\('(\d\.\d\d)', ''\)\n$", result)
        self.assertTrue(m)
        self.assertTrue(0.25 <= float(m.group(1)) < 0.28)

    def test_windowMean_outputTimes(self):
        '''windowMean should return a known set of output times for a given set of input times and windows'''
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(hours=1)
            olap = datetime.timedelta(0)
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(minutes=n*15) for n in range(48)]
            #last time is datetime.datetime(2001, 1, 1, 11, 45), so last hourly bin should
            #cover 1100 to 1200 and be centered at 1130
            data= [10]*48
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
            od_ans = [10]*12
            ot_ans = [datetime.datetime(2001, 1, 1, 0, 30),
                      datetime.datetime(2001, 1, 1, 1, 30),
                      datetime.datetime(2001, 1, 1, 2, 30),
                      datetime.datetime(2001, 1, 1, 3, 30),
                      datetime.datetime(2001, 1, 1, 4, 30),
                      datetime.datetime(2001, 1, 1, 5, 30),
                      datetime.datetime(2001, 1, 1, 6, 30),
                      datetime.datetime(2001, 1, 1, 7, 30),
                      datetime.datetime(2001, 1, 1, 8, 30),
                      datetime.datetime(2001, 1, 1, 9, 30),
                      datetime.datetime(2001, 1, 1, 10, 30),
                      datetime.datetime(2001, 1, 1, 11, 30)]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)

    def test_windowMean1(self):
        """windowMean should give known results 1(regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
            od_ans = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]
            ot_ans = [datetime.datetime(2001, 1, 1, 12, 0),
                      datetime.datetime(2001, 1, 2, 0, 0),
                      datetime.datetime(2001, 1, 2, 12, 0),
                      datetime.datetime(2001, 1, 3, 0, 0),
                      datetime.datetime(2001, 1, 3, 12, 0),
                      datetime.datetime(2001, 1, 4, 0, 0),
                      datetime.datetime(2001, 1, 4, 12, 0),
                      datetime.datetime(2001, 1, 5, 0, 0),
                      datetime.datetime(2001, 1, 5, 12, 0),
                      ]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)

    def test_windowMean2(self):
        """windowMean should give known results 2(regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap)
            od_ans = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]
            ot_ans = [datetime.datetime(2001, 1, 1, 12, 30),
                      datetime.datetime(2001, 1, 2, 0, 30),
                      datetime.datetime(2001, 1, 2, 12, 30),
                      datetime.datetime(2001, 1, 3, 0, 30),
                      datetime.datetime(2001, 1, 3, 12, 30),
                      datetime.datetime(2001, 1, 4, 0, 30),
                      datetime.datetime(2001, 1, 4, 12, 30),
                      datetime.datetime(2001, 1, 5, 0, 30),
                      datetime.datetime(2001, 1, 5, 12, 30)]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)

    def test_windowMean3(self):
        """windowMean should give known results 3(regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            time[50:] = [val + datetime.timedelta(days=2) for val in time[50:]]
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
            od_ans = [ 15.,  15.,  15.,  15.,  15.,  numpy.nan,  numpy.nan,  15.,  15.,  15.,  15., 15., 15.]
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
                      datetime.datetime(2001, 1, 6, 12, 0),
                      datetime.datetime(2001, 1, 7, 0, 0),
                      datetime.datetime(2001, 1, 7, 12, 0),
                      ]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)

    def test_windowMean4(self):
        """windowMean should give known results 4(regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            time[50:] = [val + datetime.timedelta(days=2) for val in time[50:]]
            # now test the pointwise
            outdata, outtime = tb.windowMean(data, winsize=24, overlap=12)
            od_ans = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]
            ot_ans = [12.0, 24.0, 36.0, 48.0, 60.0, 72.0, 84.0]
            numpy.testing.assert_almost_equal(ot_ans, outtime)
            numpy.testing.assert_almost_equal(od_ans, outdata)

    def test_windowMean5(self):
        """windowMean should give known results 5(regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            time[50:] = [val + datetime.timedelta(days=2) for val in time[50:]]
            # winsize tests
            outdata, outtime = tb.windowMean(data, winsize=24.6, overlap=12)
            od_ans, ot_ans = tb.windowMean(data, winsize=24.6, overlap=12)
            numpy.testing.assert_almost_equal(ot_ans, outtime)
            numpy.testing.assert_almost_equal(od_ans, outdata)
            outdata, outtime = tb.windowMean(data, winsize=0.4)
            od_ans, ot_ans = tb.windowMean(data, winsize=1.0)
            numpy.testing.assert_almost_equal(ot_ans, outtime)
            numpy.testing.assert_almost_equal(od_ans, outdata)
            outdata, outtime = tb.windowMean(data, winsize=1.0, overlap=2)
            od_ans, ot_ans = tb.windowMean(data, winsize=1.0, overlap=0)
            numpy.testing.assert_almost_equal(ot_ans, outtime)
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(8, len(w))

    def test_windowMean_op(self):
        """windowMean should give known results (regression)"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(hours=12)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1), op=len)
            od_ans = [24, 24, 24, 24, 24, 24, 24, 16, 4]
            ot_ans = [datetime.datetime(2001, 1, 1, 12, 0),
                      datetime.datetime(2001, 1, 2, 0, 0),
                      datetime.datetime(2001, 1, 2, 12, 0),
                      datetime.datetime(2001, 1, 3, 0, 0),
                      datetime.datetime(2001, 1, 3, 12, 0),
                      datetime.datetime(2001, 1, 4, 0, 0),
                      datetime.datetime(2001, 1, 4, 12, 0),
                      datetime.datetime(2001, 1, 5, 0, 0),
                      datetime.datetime(2001, 1, 5, 12, 0)]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)

    def test_windowMean_op2(self):
        """windowMean should give expected sums with len as passed function"""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            wsize = datetime.timedelta(days=1)
            olap = datetime.timedelta(0)
            data = [10, 20]*50
            time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
            #four points lie on the 5th, so the sum for that day is 4
            outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1), op=len)
            od_ans = [24, 24, 24, 24, 4]
            ot_ans = [datetime.datetime(2001, 1, 1, 12, 0),
                      datetime.datetime(2001, 1, 2, 12, 0),
                      datetime.datetime(2001, 1, 3, 12, 0),
                      datetime.datetime(2001, 1, 4, 12, 0),
                      datetime.datetime(2001, 1, 5, 12, 0)]
            numpy.testing.assert_almost_equal(od_ans, outdata)
            self.assertEqual(ot_ans, outtime)
            
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
        time = list(range(len(time)))
        self.assertRaises(TypeError, tb.windowMean, data, time, overlap=olap, st_time=datetime.datetime(2001,1,1))

    def test_windowMean_2d(self):
        """windowMean flattens and extends the time var over all the other dims"""
        wsize = 2
        olap = 0
        data = numpy.arange(20).reshape(10,2)
        out = tb.windowMean(data, winsize=wsize, overlap=olap)
        ansd = [  1.5,   5.5,   9.5,  13.5]
        numpy.testing.assert_almost_equal(ansd, out[0])
        anst = [ 1.,  3.,  5.,  7.]
        numpy.testing.assert_almost_equal(anst, out[1])


class ArrayBinTests(unittest.TestCase):
    """Tests for arraybin function"""

    def testNormalList(self):
        """Non-pathological cases, lists as input"""
        inputs = [[list(range(10)), [4.2]],
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



if __name__ == "__main__":
    unittest.main()
