#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's data manager

Copyright 2015-2020 contributors
"""

import datetime
import itertools
import ntpath
import os
import os.path
import posixpath
import shutil
import sys
import tempfile
import unittest
import warnings

import numpy.random
import numpy.testing

import spacepy.datamanager


__all__ = ["RePathTests", "DataManagerFunctionTests",
           "DataManagerBinningTests"]


class DataManagerClassTests(unittest.TestCase):
    def test_files_matching(self):
        """Files matching a format"""
        pth = os.path.dirname(os.path.abspath(__file__))
        dirlist = [os.path.join(pth, 'data', 'datamanager_test', '1'),
                   os.path.join(pth, 'data', 'datamanager_test', '2'),
                   ]
        cases = [
            {'files1': ['rbspa_ect-hope-sci-L2_20150409_v4.0.0.cdf',
                        'rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf'],
             'files2': ['rbspa_ect-hope-sci-L2_20150409_v4.1.0.cdf'],
             'regex': r'rbspa_ect-hope-sci-L2_%Y%m%d_v(\d\.){3}cdf',
             'dt': None,
             'descend': False,
         },
            {'files1': ['rbspa_ect-hope-sci-L2_20150409_v4.0.0.cdf',
                        'rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf',
                        'foo/rbspa_ect-hope-sci-L2_20150409_v5.0.0.cdf',
                        '2015/rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf',
                        ],
             'files2': ['rbspa_ect-hope-sci-L2_20150409_v4.1.0.cdf',
                        '2015/rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf'],
             'regex': r'rbspa_ect-hope-sci-L2_%Y%m%d_v(\d\.){3}cdf',
             'dt': None,
             'descend': True,
         },
            {'files1': ['2015/rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf'],
             'files2': ['2015/rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf'],
             'regex': r'%Y/rbspa_ect-hope-sci-L2_%Y%m%d_v(\d\.){3}cdf',
             'dt': None,
             'descend': False,
         },
            ]
        for c in cases:
            dm = spacepy.datamanager.DataManager(dirlist, c['regex'],
                                                 c['descend'])
            expected = [os.path.join(dirlist[0], f) for f in c['files1']]
            expected.extend([os.path.join(dirlist[1], f) for f in c['files2']])
            #We're not going to be crazy and put slashes in the tests....
            expected = [e.replace(posixpath.sep, os.path.sep)
                        for e in expected]
            output = list(dm.files_matching(c['dt']))
            numpy.testing.assert_array_equal(
                sorted(expected), sorted(output), 'Case:' + str(c))


class RePathTests(unittest.TestCase):
    def test_path_split(self):
        """Verify the simple splitting works"""
        path = "foo/bar/baz"
        self.assertEqual(["foo", "bar", "baz"],
                         spacepy.datamanager.RePath.path_split(path))
        path = "/foo/rbspa_ect-hope-sci-L2_20150409_v5.0.0.cdf"
        self.assertEqual(
            ["/", "foo", "rbspa_ect-hope-sci-L2_20150409_v5.0.0.cdf"],
            spacepy.datamanager.RePath.path_split(path))

    def test_path_split_win(self):
        """Simple splitting on Windows paths"""
        if sys.platform != 'win32':
            oldpath = os.path
            os.path = ntpath #Pretend we're on Windows
        try:
            path = "foo\\bar\\baz"
            self.assertEqual(
                ["foo", "bar", "baz"],
                spacepy.datamanager.RePath.path_split(path, native=True))
            path = "C:\\foo\\rbspa_ect-hope-sci-L2_20150409_v5.0.0.cdf"
            self.assertEqual(
                ["C:\\", "foo", "rbspa_ect-hope-sci-L2_20150409_v5.0.0.cdf"],
                spacepy.datamanager.RePath.path_split(path, native=True))
        finally:
            if sys.platform != 'win32':
                os.path = oldpath

    def test_path_slice(self):
        """Verify slicing on a path"""
        path = "foo/bar/baz"
        self.assertEqual("bar",
                         spacepy.datamanager.RePath.path_slice(path, 1))
        self.assertEqual("foo/baz",
                         spacepy.datamanager.RePath.path_slice(
                             path, 0, step=2))
        self.assertEqual("bar/baz",
                         spacepy.datamanager.RePath.path_slice(
                             path, 1, 3))

    def test_path_slice_win(self):
        """Verify slicing on a path, Windows"""
        if sys.platform != 'win32':
            oldpath = os.path
            os.path = ntpath #Pretend we're on Windows
        try:
            path = "foo\\bar\\baz"
            self.assertEqual(
                "bar",
                spacepy.datamanager.RePath.path_slice(path, 1, native=True))
            self.assertEqual(
                "foo\\baz",
                spacepy.datamanager.RePath.path_slice(
                    path, 0, step=2, native=True))
            self.assertEqual(
                "bar\\baz",
                spacepy.datamanager.RePath.path_slice(path, 1, 3, native=True))
        finally:
            if sys.platform != 'win32':
                os.path = oldpath

    def test_path_match(self):
        """Verify matching a path regex"""
        r = spacepy.datamanager.RePath(r'directory/%y/file%Y%m%d_v\d\.cdf')
        self.assertTrue(r.match('directory/99/file19990502_v0.cdf',
                                datetime.datetime(1999, 5, 2)))
        self.assertTrue(r.match('directory/99/file19990502_v0.cdf'))
        self.assertFalse(r.match('directory/99/file19990502_v0.cdf',
                                 datetime.datetime(1999, 5, 3)))

    def test_path_match_start(self):
        """Verify matching beginning of a path regex"""
        r = spacepy.datamanager.RePath(r'directory/%y/file%Y%m%d_v\d\.cdf')
        self.assertTrue(r.match('directory/99',
                                datetime.datetime(1999, 5, 2), 'start'))
        self.assertTrue(r.match('directory/99', where='start'))
        self.assertFalse(r.match('directory/99',
                                 datetime.datetime(2000, 5, 2), 'start'))

    def test_path_match_end(self):
        """Verify matching end of a path regex"""
        r = spacepy.datamanager.RePath(r'%y/file%Y%m%d_v\d\.cdf')
        self.assertTrue(r.match('directory/99/file19990502_v0.cdf',
                                datetime.datetime(1999, 5, 2), 'end'))
        self.assertTrue(r.match('directory/99/file19990502_v0.cdf',
                                where='end'))
        self.assertFalse(r.match('directory/99/file19990502_v0.cdf',
                                 datetime.datetime(2000, 5, 2), 'end'))
        r = spacepy.datamanager.RePath(
            r'rbspa_ect-hope-sci-L2_%Y%m%d_v(\d\.){3}cdf')
        self.assertTrue(r.match(
            '2015/rbspa_ect-hope-sci-L2_20150410_v4.0.0.cdf',
            where='end'))

    def test_path_match_toplevel(self):
        """Single toplevel dir"""
        #This is the reduced case of Case 3 in test_files_matching,
        #which fails on Windows
        r = spacepy.datamanager.RePath(
            r'%Y/rbspa_ect-hope-sci-L2_%Y%m%d_v(\d\.){3}cdf')
        self.assertTrue(r.match('2015', None, 'start'))


class DataManagerFunctionTests(unittest.TestCase):
    def test_insert_fill(self):
        """Verify insert_fill with relative tolerances"""
        e = [1, 2, 3, 5, 6, 7]
        d = [1, 2, 3, 3, 2, 1]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5, 6, 7]
        expected_df = [1, 2, 3, -1, 3, 2, 1]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1, 2], [3, 4], [4, 3], [2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2], [3, 4], [-1, -1], [4, 3], [2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2, -1, 3, 4], [4, 3, -1, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1., 2, 3, 4], [4, 3, 2, 1]] #has to be float for nan output
        ef, df = spacepy.datamanager.insert_fill(e, d)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2, numpy.nan, 3, 4], [4, 3, numpy.nan, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 3, 4]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4]
        expected_df = [[1, 2, 3, 4], [4, 3, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, 1, i, 0) for i in (1, 2, 4, 5)]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, 1, i, 0)
                       for i in (1, 2, 3, 4, 5)]
        expected_df = [[1, 2, -1, 3, 4], [4, 3, -1, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, 1, 0, 0, 0, i)
             for i in itertools.chain(range(1, 4), range(5, 8), range(9, 12))]
        d = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, 1, 0, 0, 0, i)
                       for i in range(1, 12)]
        expected_df = [1, 2, 3, -1, 5, 6, 7, -1, 9, 10, 11]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, i, 0, 0, 0, 0)
             for i in itertools.chain(range(1, 4), range(5, 8), range(9, 12))]
        d = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, i, 0, 0, 0, 0)
                       for i in range(1, 12)]
        expected_df = [1, 2, 3, -1, 5, 6, 7, -1, 9, 10, 11]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

    def test_insert_fill_abs(self):
        """Verify insert_fill with absolute tolerances"""
        e = [1, 2, 3, 5, 6, 7]
        d = [1, 2, 3, 3, 2, 1]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1, absolute=1.2)
        expected_ef = [1, 2, 3, 4, 5, 6, 7]
        expected_df = [1, 2, 3, -1, 3, 2, 1]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, 1, 0, i, 0, 0)
             for i in itertools.chain(range(1, 4), range(5, 8), range(9, 12))]
        d = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        ef, df = spacepy.datamanager.insert_fill(
            e, d, -1, absolute=datetime.timedelta(seconds=90.0))
        expected_ef = [datetime.datetime(2000, 1, 1, 0, i, 0, 0)
                       for i in range(1, 12)]
        expected_df = [1, 2, 3, -1, 5, 6, 7, -1, 9, 10, 11]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

    def test_apply_idx(self):
        """Verify apply_idx"""
        numpy.random.seed(0)
        indata = numpy.random.randint(0, 10000, (50, 3, 12, 8))
        #indata will by sorted in dimension 2 (12) by orderby's
        #sort order for the same dimension 0 (50)
        orderby = numpy.random.randint(0, 10000, (50, 12))
        idx = numpy.argsort(orderby, axis=1)
        outdata = spacepy.datamanager.apply_index(indata, idx)
        for i in range(50):
            numpy.testing.assert_array_equal(numpy.sort(orderby[i, :]),
                                             orderby[i, idx[i]])
        for j in range(3):
            for l in range(8):
                for i in range(50):
                    numpy.testing.assert_array_equal(indata[i, j, idx[i], l],
                                                     outdata[i, j, :, l])

    def test_values_to_steps(self):
        inval = [[1, 3, 5, 7, 6, 4, 2, 0],
                 [0, 1, 2, 3, 3, 2, 1, 0],
                 [50, 50, 100, 100, 100, 100, 50, 50],
             ]
        outval = spacepy.datamanager.values_to_steps(inval, axis=-1)
        expected = [[1, 3, 5, 7, 6, 4, 2, 0],
                    [0, 1, 2, 3, 3, 2, 1, 0],
                    [0, 0, 1, 1, 1, 1, 0, 0],
                ]
        numpy.testing.assert_array_equal(expected, outval)

    def test_flatten_idx(self):
        #Really easy cases: just do zero for everything, make sure the offset
        #is okay (before testing the stride)
        inval = numpy.zeros((5, 4), dtype=numpy.int64)
        outval = spacepy.datamanager.flatten_idx(inval, axis=1).reshape(5, 4)
        for i in range(5):
            self.assertTrue((outval[i, :] == i * 4).all())

        inval = numpy.zeros((5, 4), dtype=numpy.int64)
        outval = spacepy.datamanager.flatten_idx(inval, axis=0)
        checkval = numpy.empty((5, 4))
        for j in range(4):
            checkval[0, j] = j
            checkval[1:, j] = -1
        self.assertTrue((checkval.flatten()[outval].reshape(5, 4) ==
                [0, 1, 2, 3]).all())

        #Easy cases
        inval = [[0, 1, 2], [0, 1, 2]]
        outval = spacepy.datamanager.flatten_idx(inval, axis=1)
        self.assertTrue((outval == [0, 1, 2, 3, 4, 5]).all())

        inval = [[0, 0, 0], [1, 1, 1]]
        outval = spacepy.datamanager.flatten_idx(inval, axis=0)
        numpy.testing.assert_array_equal([0, 1, 2, 3, 4, 5], outval)

        #Big and fancy and difficult
        numpy.random.seed(0)
        inval = numpy.random.randint(10000, size=(5, 6, 7))

        idx = numpy.argsort(inval, axis=0)
        sorted = inval.ravel()[spacepy.datamanager.flatten_idx(idx, axis=0)
        ].reshape(inval.shape)
        self.assertEqual(inval.shape, sorted.shape)
        for j in range(6):
            for k in range(7):
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[:, j, k]),
                    inval[idx[:, j, k], j, k])
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[:, j, k]),
                    sorted[:, j, k])

        idx = numpy.argsort(inval, axis=1)
        sorted = inval.ravel()[spacepy.datamanager.flatten_idx(idx, axis=1)
        ].reshape(inval.shape)
        assert(sorted.shape == inval.shape)
        for i in range(5):
            for k in range(7):
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[i, :, k]),
                    inval[i, idx[i, :, k], k])
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[i, :, k]),
                    sorted[i, :, k])

        idx = numpy.argsort(inval, axis=2)
        sorted = inval.ravel()[spacepy.datamanager.flatten_idx(idx, axis=2)
        ].reshape(inval.shape)
        assert(sorted.shape == inval.shape)
        for i in range(5):
            for j in range(6):
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[i, j, :]),
                    inval[i, j, idx[i, j, :]])
                numpy.testing.assert_array_equal(
                    numpy.sort(inval[i, j, :]),
                    sorted[i, j, :])
        #Add a test that assignment works?

    def test_axis_index(self):
        shape = (6, 3, 4, 12, 9)
        axis = 2
        output = spacepy.datamanager.axis_index(shape, axis)
        expected = numpy.arange(shape[axis])
        assert(output.shape == shape)
        for i in range(6):
            for j in range(3):
                for l in range(12):
                    for m in range(9):
                        numpy.testing.assert_array_equal(expected,
                                                         output[i, j, :, l, m])

    def test_rev_index(self):
        numpy.random.seed(0)
        inval = numpy.random.randint(10000, size=(5, 6, 7))
        for axis in (0, 1, 2):
            idx = numpy.argsort(inval, axis=axis)
            idx_rev = spacepy.datamanager.rev_index(idx, axis=axis)
    #apply_index doesn't allow a choice of axis...maybe later
    #        assert((inval == spacepy.datamanager.apply_index(
    #            spacepy.datamanager.apply_index(inval, idx, axis=axis),
    #            idx_rev, axis=axis)).all())
            self.assertTrue(
                (inval.ravel()[spacepy.datamanager.flatten_idx(idx, axis)]
                 [spacepy.datamanager.flatten_idx(idx_rev, axis)]
                 .reshape(inval.shape) == inval).all())


class DataManagerBinningTests(unittest.TestCase):
    """Test of binning and helper functions"""

    def testFindShape(self):
        """Test reshaping of smaller to larger"""
        # Large, small, reshaped
        cases = [
            [(5, 6, 10), (5, 6, 10), (5, 6, 10)],
            [(5, 6, 10), (5, 1, 10), (5, 1, 10)],
            [(5, 6, 10), (5, 10), (5, 1, 10)],
            [(5, 6, 10), (5,), (5, 1, 1)],
            [(5, 6, 10), (5, 6), (5, 6, 1)],
            [(5, 6, 10), (6, 10), (1, 6, 10)],
            ]
        for large, small, reshaped in cases:
            newshape = spacepy.datamanager._find_shape(large, small)
            self.assertEqual(reshaped, newshape)
            if small == reshaped: # No actual reshaping
                continue
            numpy.random.seed(0xdeadbeef)
            # Make sure reshaping doesn't affect indexing once the
            # new dims are removed.
            oldarray = numpy.random.randint(100, size=small)
            newarray = numpy.reshape(oldarray, reshaped)
            idx = tuple([0 if i == 1 else slice(None) for i in newshape])
            numpy.testing.assert_array_equal(
                oldarray,
                newarray[idx])

    def testRebinSimple(self):
        """Test rebinning of an array, simple case"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 6, 10)
        bindata = numpy.random.rand(10) * 10

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    expected[i, j, binno] = numpy.mean(indata[i, j, :][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, bintype='unc')
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    if count:
                        total = numpy.sum(indata[i, j, :][idx] ** 2)
                        expected[i, j, binno] = numpy.sqrt(total) / count
                    else:
                        expected[i, j, binno] = numpy.nan
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, bintype='count')
        expected = numpy.empty(dtype=numpy.int64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    expected[i, j, binno] = count
        numpy.testing.assert_array_equal(
            expected, rebinned)

    def testRebinClipped(self):
        """Test rebinning of an array with out-of-range data"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 6, 10)
        bindata = numpy.random.rand(10) * 10
        bindata[0] = -99 # Force out-of-range bin

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    expected[i, j, binno] = numpy.mean(indata[i, j, :][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins, clip=True)
        # The out-of-range values, when clipped, go to lowest bin...
        # so clip in advance and recalc our expected for that bin.
        bindata[0] = 0
        binno = 0
        idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
        for i in range(100):
            for j in range(6):
                expected[i, j, binno] = numpy.mean(indata[i, j, :][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

    def testRebinSimpleWeighted(self):
        """Test rebinning of an array with weights specified"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 6, 10)
        bindata = numpy.random.rand(10) * 10
        weights = numpy.arange(.1, 1.1, 0.1)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, weights=weights)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    total = numpy.sum(indata[i, j, :][idx] * weights[idx])
                    count = numpy.sum(weights[idx])
                    expected[i, j, binno] = total / count if count \
                                            else numpy.nan
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)
        # Multiplying all the weights by a factor shouldn't change answer
        newrebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, weights=weights * 10)
        numpy.testing.assert_allclose(
            rebinned, newrebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins,
                                weights=weights, bintype='unc')
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    total = numpy.sum((indata[i, j, :][idx]
                                                  * weights[idx]) ** 2)
                    count = numpy.sum(weights[idx])
                    expected[i, j, binno] = numpy.sqrt(total) / count if count \
                                            else numpy.nan
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)
        # Multiplying all the weights by a factor shouldn't change answer
        newrebinned = spacepy.datamanager.rebin(indata, bindata, bins,
                                   weights=weights * 10, bintype='unc')
        numpy.testing.assert_allclose(
            rebinned, newrebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins,
                                weights=weights, bintype='count')
        expected = numpy.empty(dtype=numpy.float, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    count = numpy.sum(weights[idx])
                    expected[i, j, binno] = count
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)
        # Multiplying all the weights by 10 should increase counts by 10
        newrebinned = spacepy.datamanager.rebin(indata, bindata, bins,
                                   weights=weights * 10, bintype='count')
        numpy.testing.assert_allclose(
            rebinned * 10, newrebinned, atol=1e-20)

    def testRebinSimpleNaN(self):
        """Test rebinning of an array, NaN input"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 6, 10)
        bindata = numpy.random.rand(10) * 10
        indata[0, 0, 0] = numpy.nan
        indata[5, ...] = numpy.nan

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    with warnings.catch_warnings(record=True):
                        expected[i, j, binno] \
                            = numpy.nanmean(indata[i, j, :][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)
        self.assertTrue(numpy.isnan(rebinned[5, ...]).all())

    def testRebinSimpleDiffAxis(self):
        """Test rebinning of an array, simple case, different axis"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 10, 6)
        bindata = numpy.random.rand(10) * 10

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins, axis=1)
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    expected[i, binno, j] = numpy.mean(indata[i, :, j][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, axis=1, bintype='unc')
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    if count:
                        total = numpy.sum(indata[i, :, j][idx] ** 2)
                        expected[i, binno, j] = numpy.sqrt(total) / count
                    else:
                        expected[i, binno, j] = numpy.nan
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, axis=1, bintype='count')
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.int64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            idx = (bins[binno] <= bindata) & (bindata < bins[binno + 1])
            for i in range(100):
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    expected[i, binno, j] = count
        numpy.testing.assert_array_equal(
            expected, rebinned)

    def testRebin2DDiffAxis(self):
        """Test rebinning of an array, 2D bin data, axis specified"""
        numpy.random.seed(0xdeadbeef)
        bins = numpy.arange(0, 11, 2)
        indata = numpy.random.rand(100, 10, 6)
        bindata = numpy.reshape(numpy.random.rand(100 * 10) * 10,
                                (100, 10))

        rebinned = spacepy.datamanager.rebin(indata, bindata, bins, axis=1)
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            for i in range(100):
                idx = (bins[binno] <= bindata[i, :])\
                       & (bindata[i, :] < bins[binno + 1])
                for j in range(6):
                    with warnings.catch_warnings(record=True):
                        expected[i, binno, j] = numpy.mean(indata[i, :, j][idx])
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, axis=1, bintype='unc')
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            for i in range(100):
                idx = (bins[binno] <= bindata[i, :])\
                       & (bindata[i, :] < bins[binno + 1])
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    if count:
                        total = numpy.sum(indata[i, :, j][idx] ** 2)
                        expected[i, binno, j] = numpy.sqrt(total) / count
                    else:
                        expected[i, binno, j] = numpy.nan
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)

        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, axis=1, bintype='count')
        self.assertEqual((100, 5, 6), rebinned.shape)
        expected = numpy.empty(dtype=numpy.int64, shape=(100, 5, 6))
        for binno in range(len(bins) - 1):
            for i in range(100):
                idx = (bins[binno] <= bindata[i, :])\
                      & (bindata[i, :] < bins[binno + 1])
                for j in range(6):
                    count = numpy.sum(numpy.require(idx, dtype=numpy.int64))
                    expected[i, binno, j] = count
        numpy.testing.assert_array_equal(
            expected, rebinned)

    def testRebinBindataDelta(self):
        """Rebin data where the binning data have deltas"""
        bins = numpy.arange(6)
        indata = numpy.array([1, 2, 3, 4, 5])
        bindata = numpy.array([1, 1.5, 1, 3, 4])
        # This means input ranges
        # A 0.5 to 1.5: 1
        # B 1.0 to 2.0: 2
        # C 0.5 to 1.5: 3
        # D 2.5 to 3.5: 4
        # E 3.5 to 4.5: 5
        #Outputs are 0-1, 1-2, 2-3, 3-4, 4-5
        expected = numpy.array([
            2, 2, 4, 4.5, 5])
        rebinned = spacepy.datamanager.rebin(indata, bindata, bins,
                                             bindatadelta=0.5)
        numpy.testing.assert_array_equal(
            expected, rebinned)
        # Try this specifying an array of bins
        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, bindatadelta=numpy.repeat([0.5], 5))
        numpy.testing.assert_array_equal(
            expected, rebinned)

    def testRebinBindataDelta2D(self):
        """Rebin data where the binning data have deltas, multi-D"""
        bins = numpy.arange(0, 9, 4)
        indata = numpy.array([[40., 56], [38, 93], [51, 60], [91, 42]])
        bindata = numpy.array([[0.5, 1], [2.5, 3.5], [5, 6], [6, 6.5]])
        bindeltas = numpy.array([[0.5, 1], [1.5, 1.5], [1, 1], [1, 0.5]])
        expected = numpy.array([[2, 5. / 3], [2, 7. / 3]])
        counts = spacepy.datamanager.rebin(
            indata, bindata, bins, bindatadelta=bindeltas, axis=0,
            bintype='count')
        numpy.testing.assert_allclose(expected, counts)
        expected = numpy.array([[39, 70.8], [71, 57]])
        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, bindatadelta=bindeltas, axis=0)
        numpy.testing.assert_allclose(
            expected, rebinned)
        # Treat these values as uncertainties
        expected = numpy.array([[27.586228448267445, 50.12783657809302],
                                [52.15841255253078, 34.084229401257566]])
        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, bindatadelta=bindeltas, axis=0,
            bintype='unc')
        numpy.testing.assert_allclose(
            expected, rebinned)

    def testRebinBindataAxis0Irregular(self):
        """Rebin data on axis 0 with irregular sampling"""
        bins = numpy.arange(0, 9, 4)
        indata = numpy.array([[40., 56], [38, 93], [51, 60], [91, 42]])
        bindata = numpy.array([[0.5, 1], [2.5, 3.5], [4, 6], [6, 6.5]])
        expected = numpy.array([[39, 74.5], [71, 51]])
        rebinned = spacepy.datamanager.rebin(
            indata, bindata, bins, axis=0)
        numpy.testing.assert_array_equal(
            expected, rebinned)

    def testRebinListInput(self):
        """Test rebinning of an array, list instead of array input"""
        bins = [0, 4, 8]
        indata = [[40, 38, 51, 91], [56, 93, 60, 42]]
        bindata = [[0.5, 2.5, 4, 6], [1, 3.5, 6, 6.5]]
        expected = numpy.empty(dtype=numpy.float64, shape=(100, 6, 5))
        expected = [[39, 71], [74.5, 51]]
        rebinned = spacepy.datamanager.rebin(indata, bindata, bins)
        numpy.testing.assert_allclose(
            expected, rebinned, atol=1e-20)


if __name__ == "__main__":
    unittest.main()
