#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for pycdf

Copyright 2010-2014 Los Alamos National Security, LLC.
"""

try:
    from collections.abc import Callable
except ImportError:
    from collections import Callable
import ctypes
import datetime
import gc
import hashlib
import operator
import os, os.path
try:
    import cPickle as pickle
except:
    import pickle
import re
import shutil
import sys
import tempfile
import unittest
import warnings

import matplotlib.dates
import numpy
import numpy.testing
import spacepy_testing
from spacepy import datamodel
import spacepy.pycdf as cdf
import spacepy.pycdf.const as const

__all__ = ['NoCDF', 'MakeCDF', 'CDFTestsBase', 'CDFTests', 'ColCDFTests',
           'OpenCDF', 'ReadCDF', 'ReadColCDF', 'ChangeCDFBase', 'ChangeCDF',
           'ChangezVar', 'ChangeAttr', 'ChangeColCDF']


class est_tz(datetime.tzinfo):
    """Eastern Standard timezone (no daylight time"""

    def utcoffset(self, dt):
        """Offset from UTC"""
        return datetime.timedelta(hours=-5)

    def dst(self, dt):
        """Minute offset for DST"""
        return datetime.timedelta(0)

    def tzname(self, dt):
        """Name of this time zone"""
        return 'EST'


class NoCDF(unittest.TestCase):
    """Tests that do not involve a CDF file"""
    def testErrorMessage(self):
        """Displays correct messages for exceptions"""
        exceptdic = { cdf.const.CDF_OK:
                      'CDF_OK: Function completed successfully.',
                      cdf.const.ATTR_EXISTS:
                      'ATTR_EXISTS: Named attribute already exists.',
                      cdf.const.CDF_CREATE_ERROR:
                      'CDF_CREATE_ERROR: Creation failed - error from file system.',
                      }
        for status, message in list(exceptdic.items()):
            try:
                raise cdf.CDFError(status)
            except cdf.CDFError:
                (type, val, traceback) = sys.exc_info()
                self.assertEqual(val.__str__(), message)
            else:
                self.assertTrue(False, 'Should have raised a CDFError: ' + message)

    def testHypersliceReorder(self):
        """Reorders sequences to switch array majority"""
        input = [[1, 2, 3, 4, 5], [3, -5, 6, 12], ]
        output = [[1, 5, 4, 3, 2], [3, 12, 6, -5], ]
        for (inp, outp) in zip(input, output):
            self.assertEqual(cdf._Hyperslice.reorder(inp).tolist(),
                             outp)

    def testHypersliceconvert(self):
        """Converts start/stop/step to CDF intervals"""
        input = [[None, None, None, 5],
                 [1, 4, None, 5],
                 [-5, -1, 1, 5],
                 [-1, -5, 1, 5],
                 [-1, -5, -1, 5],
                 [-1, -6, -1, 5],
                 [-1, None, -1, 5],
                 [-1, -20, -1, 5],
                 [-4, 0, -6, 10],
                 [-10, 10, 4, 10],
                 [-10, -6, 9, 10],
                 [-6, -9, -7, 10],
                 [-4, -9, -2, 10],
                 [-2, -1, -2, 10],
                 [-3, 4, -1, 10],
                 [10, -17, 10, 20],
                 [-6, -15, -10, 20],
                 ]
        output = [[0, 5, 1, False],
                  [1, 3, 1, False],
                  [0, 4, 1, False],
                  [0, 0, 1, False],
                  [1, 4, 1, True],
                  [0, 5, 1, True],
                  [0, 5, 1, True],
                  [0, 5, 1, True],
                  [6, 1, 6, True],
                  [0, 3, 4, False],
                  [0, 1, 9, False],
                  [4, 1, 7, True],
                  [2, 3, 2, True],
                  [10, 0, 2, True],
                  [5, 3, 1, True],
                  [10, 0, 10, False],
                  [14, 1, 10, True],
                  ]
        for (inp, outp) in zip(input, output):
            result = cdf._Hyperslice.convert_range(*inp)
            self.assertEqual(tuple(outp), result,
                             str(tuple(outp)) + ' != ' + str(result) +
                             ' for input ' + str(inp))

    def testHypersliceDimensions(self):
        """Find dimensions of an array"""
        data = [[[2, 3], [4, 5], [6, 7]],
                [[8, 9], [0, 1], [2, 3]],
                [[4, 5], [6, 7], [8, 9]],
                [[0, 1], [2, 3], [4, 5]],
                ]
        self.assertEqual(cdf._Hyperslice.dimensions(data),
                         (4, 3, 2))

        data = [[[2, 3], [4, 5], [6, 7]],
                [[8, 9], [0, 1], [2, 3]],
                [[4, 5], [6, 7],],
                [[0, 1], [2, 3], [4, 5]],
                ]
        messages = ('Data must be well-formed, regular array of number, '
                    'string, or datetime',
                    'setting an array element with a sequence.',
                )
        try:
            cdf._Hyperslice.dimensions(data)
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertTrue(str(v) in messages)
        else:
            self.fail('Should raise ValueError: ' + messages[0])

        self.assertEqual(cdf._Hyperslice.dimensions('hi'),
                         ())

    def testExpandEllipsis(self):
        """Tests of hyperslice's expand ellipsis method"""
        q = numpy.array([1, 2, 3]) #comparison fails if not the SAME array
        inputs = [((Ellipsis, 0), 2),
                  ((0,), 1),
                  (Ellipsis, 3),
                  (([1, 2, 3], 4, Ellipsis, 0), 5),
                  ((q, 4, Ellipsis, 0), 5),
                  ]
        expected = [(slice(None), 0),
                    (0,),
                    (slice(None, None, None),
                     slice(None, None, None),
                     slice(None, None, None)),
                    ([1, 2, 3], 4, slice(None), slice(None), 0),
                    (q, 4, slice(None), slice(None), 0),
                    ]
        for i, e in zip(inputs, expected):
            output = cdf._Hyperslice.expand_ellipsis(*i)
            #Later versions of numpy compare poorly with mixed-type arrays,
            #so have to loop over (which makes figuring out which one failed
            #pretty difficult!)
            self.assertEqual(len(e), len(output))
            for idx in range(len(e)):
                if isinstance(output[idx], numpy.ndarray):
                    numpy.testing.assert_array_equal(output[idx], e[idx])
                else:
                    self.assertEqual(e[idx], output[idx])

    def testExpandEllipsisError(self):
        """Test hyperslice expand ellipsis with too many indices"""
        self.assertRaises(IndexError,
                          cdf._Hyperslice.expand_ellipsis,
                          (Ellipsis, 0, 0, 0), 2)
        self.assertRaises(IndexError,
                          cdf._Hyperslice.expand_ellipsis,
                          (Ellipsis, Ellipsis,), 2)

    def testTT2000ToDatetime(self):
        if not cdf.lib.supports_int8:
            self.assertRaises(NotImplementedError, cdf.lib.tt2000_to_datetime,
                              1)
            return
        epochs = [284040066184000000,
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(dt, cdf.lib.tt2000_to_datetime(epoch))
        result = cdf.lib.v_tt2000_to_datetime(numpy.array(epochs))
        expected = numpy.array(dts)
        numpy.testing.assert_array_equal(expected, result)

    def testEpoch16ToDatetime(self):
        epochs = [[63397987199.0, 999999999999.0],
                  [-1.0, -1.0],
                  [0.0, 0.0],
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               datetime.datetime(9999, 12, 13, 23, 59, 59, 999999),
               datetime.datetime(9999, 12, 13, 23, 59, 59, 999999),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(dt, cdf.lib.epoch16_to_datetime(*epoch))
        result = cdf.lib.v_epoch16_to_datetime(numpy.array(epochs))
        expected = numpy.array(dts)
        numpy.testing.assert_array_equal(expected, result)

    def testEpochToDatetime(self):
        epochs = [63397987200000.0,
                  -1.0,
                  0.0,
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               datetime.datetime(9999, 12, 31, 23, 59, 59, 999000),
               datetime.datetime(9999, 12, 13, 23, 59, 59, 999000),
               ]
        #3.7.1 updated the "invalid" value to 12/31, but only in one
        #case; rest remain 12/13.
        if cdf.lib.version[0:3] < (3, 7, 1):
            dts[1] = datetime.datetime(
                9999, 12, 13, 23, 59, 59, 999000)
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(dt, cdf.lib.epoch_to_datetime(epoch))
        result = cdf.lib.v_epoch_to_datetime(numpy.array(epochs))
        expected = numpy.array(dts)
        numpy.testing.assert_array_equal(expected, result)

    def testDatetimeToTT2000(self):
        if not cdf.lib.supports_int8:
            self.assertRaises(NotImplementedError, cdf.lib.datetime_to_tt2000,
                              datetime.datetime(2009, 1, 1))
            return
        epochs = [284040066184000000,
                  284040066184000000,
                  -9223372036854775808]
        dts = [datetime.datetime(2009, 1, 1),
               datetime.datetime(2008, 12, 31, 19, tzinfo=est_tz()),
               datetime.datetime(9999, 12, 31, 23, 59, 59, 999999)
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(epoch, cdf.lib.datetime_to_tt2000(dt))
        result = cdf.lib.v_datetime_to_tt2000(numpy.array(dts))
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)

    def testDatetimeToEpoch16(self):
        epochs = [(63397987200.0, 0.0),
                  (63397987200.0, 0.0),
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               datetime.datetime(2008, 12, 31, 19, tzinfo=est_tz()),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(epoch, cdf.lib.datetime_to_epoch16(dt))
        result = cdf.lib.v_datetime_to_epoch16(numpy.array(dts))
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)
        result = cdf.lib.v_datetime_to_epoch16(dts[0])
        numpy.testing.assert_array_equal(epochs[0], result)

    def testDatetimeToEpoch(self):
        epochs = [63397987200000.0,
                  63397987200000.0,
                  63397987200001.0,
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               datetime.datetime(2008, 12, 31, 19, tzinfo=est_tz()),
               datetime.datetime(2009, 1, 1, 0, 0, 0, 501),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(epoch, cdf.lib.datetime_to_epoch(dt))
        result = cdf.lib.v_datetime_to_epoch(numpy.array(dts))
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)

    def testEpochToTT2000(self):
        """Epoch to TT2000"""
        if not cdf.lib.supports_int8:
            self.assertRaises(NotImplementedError, cdf.lib.epoch_to_tt2000,
                              63397987200000.0)
            return
        epochs = [63397987200000.0,
                  63397987200001.0,
                  ]
        tt2000s = [284040066184000000,
                   284040066185000000,
                   ]
        for (epoch, tt2000) in zip(epochs, tt2000s):
            self.assertEqual(
                tt2000, cdf.lib.epoch_to_tt2000(epoch))
        result = cdf.lib.v_epoch_to_tt2000(epochs)
        expected = numpy.array(tt2000s)
        numpy.testing.assert_array_equal(expected, result)

    def testTT2000ToEpoch(self):
        """TT2000 to Epoch"""
        if not cdf.lib.supports_int8:
            self.assertRaises(NotImplementedError, cdf.lib.tt2000_to_epoch,
                              284040066184000000)
            return
        epochs = [63397987200000.0,
                  63397987200001.0,
                  ]
        tt2000s = [284040066184000000,
                   284040066185000000,
                   ]
        for (epoch, tt2000) in zip(epochs, tt2000s):
            self.assertEqual(
                epoch, cdf.lib.tt2000_to_epoch(tt2000))
        result = cdf.lib.v_tt2000_to_epoch(tt2000s)
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)

    def testEpoch16ToTT2000(self):
        epochs = [[63397987199.0, 999999999999.0],
                  [63400665600.0, 100000000.0],
                  ]
        tt2000s = [284040065183999999,
                   286718466184100000,
                   ]
        for (epoch, tt2000) in zip(epochs, tt2000s):
            self.assertEqual(tt2000, cdf.lib.epoch16_to_tt2000(*epoch))
        result = cdf.lib.v_epoch16_to_tt2000(numpy.array(epochs))
        expected = numpy.array(tt2000s)
        numpy.testing.assert_array_equal(expected, result)

    def testTT2000ToEpoch16(self):
        epochs = [[63366364799.0, 999999999000.0],
                  [63400665600.0, 100000000.0],
                  ]
        tt2000s = [252417665183999999,
                   286718466184100000,
                   ]
        result = cdf.lib.v_tt2000_to_epoch16(numpy.array(tt2000s))
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)
        for (epoch, tt2000) in zip(epochs, tt2000s):
            self.assertEqual(epoch, list(cdf.lib.tt2000_to_epoch16(tt2000)))

    def testTT2000ToEpoch16Stress(self):
        """Test TT2000toEpoch16 repeatedly to make sure it doesn't crash"""
        epochs = [[63366364799.0, 999999999000.0],
                  [63400665600.0, 100000000.0],
                  ]
        tt2000s = [252417665183999999,
                   286718466184100000,
                   ]
        #this is basically 2x what reliably fails in the Win32/py2.7 case,
        #where this problem was found
        for i in range(10):
            for (epoch, tt2000) in zip(epochs, tt2000s):
                self.assertEqual(epoch, list(cdf.lib.tt2000_to_epoch16(tt2000)))

    def testTT2000ToEpoch16Filled(self):
        """Test conversion of TT2000 fill value to epoch 16"""
        self.assertEqual(
            (-1.e31, -1.e31),
            cdf.lib.tt2000_to_epoch16(const.FILLED_TT2000_VALUE))

    def testEpochtoEpoch16(self):
        """Convert an Epoch to Epoch16"""
        epochs = [63397987200000.0,
                  63397987200001.0,
                  ]
        epoch16s = [(63397987200.0, 0.0),
                   (63397987200.0, 1000000000.0),
                   ]
        for (epoch, epoch16) in zip(epochs, epoch16s):
            numpy.testing.assert_array_equal(
                epoch16, cdf.lib.epoch_to_epoch16(epoch))
        result = cdf.lib.epoch_to_epoch16(epochs)
        expected = numpy.array(epoch16s)
        numpy.testing.assert_array_equal(expected, result)

    def testEpochtoNum(self):
        """Convert Epoch to matplotlib number"""
        dts = [datetime.datetime(2012, 1, 3, 23, 12),
               datetime.datetime(1857, 7, 23, 12, 1),
               datetime.datetime(1934, 2, 7, 23, 5),
               ]
        epochs = cdf.lib.v_datetime_to_epoch(dts)
        expected = matplotlib.dates.date2num(dts)
        result = cdf.lib.epoch_to_num(epochs)
        numpy.testing.assert_array_equal(expected, result)

    def testEpoch16toEpoch(self):
        """Convert an Epoch16 to Epoch"""
        epochs = [63397987200000.0,
                  63397987200001.0,
                  ]
        epoch16s = [(63397987200.0, 0.0),
                   (63397987200.0, 1000000000.0),
                   ]
        for (epoch, epoch16) in zip(epochs, epoch16s):
            self.assertEqual(epoch, cdf.lib.epoch16_to_epoch(epoch16))
        result = cdf.lib.epoch16_to_epoch(epoch16s)
        expected = numpy.array(epochs)
        numpy.testing.assert_array_equal(expected, result)

    def testDatetimeEpoch16RT(self):
        """Roundtrip datetimes to epoch16s and back"""
        dts = [datetime.datetime(2008, 12, 15, 3, 12, 5, 1000),
               datetime.datetime(1821, 1, 30, 2, 31, 5, 23000),
               datetime.datetime(2050, 6, 5, 15, 0, 5, 0),
               ]
        for dt in dts:
            self.assertEqual(dt, cdf.lib.epoch16_to_datetime(
                *cdf.lib.datetime_to_epoch16(dt)))

    def testDatetimeEpochRT(self):
        """Roundtrip datetimes to epochs and back"""
        if not cdf.lib.supports_int8:
            return
        dts = [datetime.datetime(2008, 12, 15, 3, 12, 5, 1000),
               datetime.datetime(1821, 1, 30, 2, 31, 5, 23000),
               datetime.datetime(2050, 6, 5, 15, 0, 5, 0),
               ]
        for dt in dts:
            self.assertEqual(dt, cdf.lib.epoch_to_datetime(
                cdf.lib.datetime_to_epoch(dt)))

    def testDatetimeTT2000RT(self):
        """Roundtrip datetimes to TT2000 and back"""
        if not cdf.lib.supports_int8:
            return
        dts = [datetime.datetime(2008, 12, 15, 3, 12, 5, 1000),
               datetime.datetime(1821, 1, 30, 2, 31, 5, 23000),
               datetime.datetime(2050, 6, 5, 15, 0, 5, 0),
               ]
        #CDF library bug: some dates don't work on 32-bit
        if sys.maxsize <= 2 ** 31:
            dts[1] = datetime.datetime(1945, 1, 30, 2, 31, 5, 23000)
        for dt in dts:
            self.assertEqual(dt, cdf.lib.tt2000_to_datetime(
                cdf.lib.datetime_to_tt2000(dt)))

    def testIgnoreErrors(self):
        """Call the library and ignore particular error"""
        nbytes = ctypes.c_long(0)
        status = cdf.lib.call(cdf.const.GET_, cdf.const.DATATYPE_SIZE_,
                              ctypes.c_long(100), ctypes.byref(nbytes),
                              ignore=(cdf.const.BAD_DATA_TYPE,))
        self.assertEqual(cdf.const.BAD_DATA_TYPE, status)

    def testVersion(self):
        """Check library's version"""
        self.assertTrue(cdf.lib.version[0] in (2, 3))
        self.assertTrue(cdf.lib.version[1] in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
        self.assertTrue(cdf.lib.version[2] in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9))
        self.assertTrue(re.match(r'( |a|\d)?', str(cdf.lib.version[3])))
        if cdf.lib.version == (3, 3, 0, ' '):
            self.assertTrue(cdf.lib._del_middle_rec_bug)
        elif cdf.lib.version == (3, 3, 1, ' '):
            self.assertTrue(cdf.lib._del_middle_rec_bug)
        elif cdf.lib.version == (3, 4, 0, '0'):
            self.assertTrue(cdf.lib._del_middle_rec_bug)
        elif cdf.lib.version == (3, 4, 1, '0'):
            self.assertFalse(cdf.lib._del_middle_rec_bug)

    def testTypeGuessing(self):
        """Guess CDF types based on input data"""
        self.longMessage = True #doesn't help on 2.6, but that's dying
        samples = [[1, 2, 3, 4],
                   [[1.2, 1.3, 1.4], [2.2, 2.3, 2.4]],
                   ['hello', 'there', 'everybody'],
                   datetime.datetime(2009, 1, 1),
                   datetime.datetime(2009, 1, 1, 12, 15, 12, 1000),
                   datetime.datetime(2009, 1, 1, 12, 15, 12, 1),
                   [1.0],
                   0.0,
                   numpy.array([1, 2, 3], dtype=numpy.int32),
                   numpy.array([1, 2, 4], dtype=numpy.float64),
                   numpy.array([1, 2, 5], dtype=numpy.int64),
                   2 ** 62,
                   -1.0,
                   numpy.array([1, 2, 6], dtype='<u2'),
                   numpy.array([1, 2, 7], dtype='>u2'),
                   numpy.int64(-1 * 2 ** 63),
                   numpy.int32(-1 * 2 ** 31),
                   -1 * 2 ** 31,
                   numpy.array([5, 6, 7], dtype=numpy.uint8),
                   [4611686018427387904],
                   numpy.array([1], dtype=object),
                   ]
        type8 = [((4,), [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4, const.CDF_INT8,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((2, 3), [const.CDF_FLOAT, const.CDF_REAL4,
                           const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_CHAR, const.CDF_UCHAR], 9),
                 ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH,
                       const.CDF_EPOCH16], 1),
                 ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH,
                       const.CDF_EPOCH16], 1),
                 ((), [const.CDF_TIME_TT2000, const.CDF_EPOCH16,
                       const.CDF_EPOCH], 1),
                 ((1,), [const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_INT4], 1),
                 ((3,), [const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_INT8], 1),
                 ((), [const.CDF_INT8, const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_UINT2], 1),
                 ((3,), [const.CDF_UINT2], 1),
                 ((), [const.CDF_INT8], 1),
                 ((), [const.CDF_INT4], 1),
                 ((), [const.CDF_INT4, const.CDF_INT8,
                       const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_UINT1, const.CDF_UCHAR], 1),
                 ((1,), [const.CDF_INT8,
                       const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((1,), [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4, const.CDF_INT8,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ]
        types = [((4,), [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((2, 3), [const.CDF_FLOAT, const.CDF_REAL4,
                           const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_CHAR, const.CDF_UCHAR], 9),
                 ((), [const.CDF_EPOCH, const.CDF_EPOCH16], 1),
                 ((), [const.CDF_EPOCH, const.CDF_EPOCH16], 1),
                 ((), [const.CDF_EPOCH16, const.CDF_EPOCH], 1),
                 ((1,), [const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_INT4], 1),
                 ((3,), [const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_UINT2], 1),
                 ((3,), [const.CDF_UINT2], 1),
                 ((), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((), [const.CDF_INT4], 1),
                 ((), [const.CDF_INT4, const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((3,), [const.CDF_UINT1, const.CDF_UCHAR], 1),
                 ((1,), [const.CDF_FLOAT, const.CDF_REAL4,
                       const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ((1,), [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ]
        self.assertRaises(ValueError, cdf._Hyperslice.types, [object()])
        if cdf.lib.supports_int8: #explicitly test backward-compatible
            cdf.lib.supports_int8 = False
            try:
                for (s, t) in zip(samples, types):
                    t = (t[0], [i.value for i in t[1]], t[2])
                    self.assertEqual(t, cdf._Hyperslice.types(s),
                                     msg='Input ' + str(s))
            finally: #don't leave the library hashed if test fails
                cdf.lib.supports_int8 = True
            types = type8
        for (s, t) in zip(samples, types):
            t = (t[0], [i.value for i in t[1]], t[2])
            self.assertEqual(t, cdf._Hyperslice.types(s),
                             msg='Input ' + str(s))

    @unittest.skipIf(cdf.lib.version[0] < 3,
                     "Not supported with CDF library < 3")
    def testMinMaxTT2000(self):
        """Get min/max values for TT2000 types"""
        minval, maxval = cdf.lib.get_minmax(const.CDF_TIME_TT2000)
        #Make sure the minimum isn't just plain invalid
        self.assertTrue(minval < datetime.datetime(9999, 1, 1))

    def testMinMaxFloat(self):
        """Get min/max values for a float"""
        with spacepy_testing.assertDoesntWarn(
                self, 'always',
                r'Conversion of the second argument of issubdtype',
                FutureWarning, r'spacepy'):
            minval, maxval = cdf.lib.get_minmax(const.CDF_FLOAT)
        self.assertAlmostEqual(-3.4028234663853e+38, minval, places=-30)
        self.assertAlmostEqual(3.4028234663853e+38, maxval, places=-30)

    def testMinMaxInt(self):
        """Get min/max values for an integer"""
        minval, maxval = cdf.lib.get_minmax(const.CDF_INT1)
        self.assertEqual(-128, minval)
        self.assertEqual(127, maxval)

    def testConcatCDF(self):
        """Read from two sequential CDFs"""
        td = tempfile.mkdtemp()
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        warnings.filterwarnings(
            'ignore', r'^No type specified for time input.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(os.path.join(td, 'one.cdf'), create=True) as cdffile:
                cdffile.attrs['gattrone'] = 1
                cdffile.attrs['gattrtwo'] = 2
                cdffile.attrs['gattrthree'] = 3
                cdffile['var1'] = numpy.array([1, 2, 3], dtype=numpy.float32)
                cdffile['var2'] = numpy.array([11, 12, 13], dtype=numpy.float32)
                cdffile['var3'] = numpy.array([21, 22, 23], dtype=numpy.float32)
                cdffile['Epoch'] = numpy.array([datetime.datetime(2010, 1, i)
                                            for i in range(1, 4)])
                cdffile['var2'].attrs['foo'] = 'file1'
                cdffile.new(
                    'var4', data=numpy.array([99, 100], dtype=numpy.float32),
                    recVary=False)
            with cdf.CDF(os.path.join(td, 'two.cdf'),
                                   create=True) as cdffile:
                cdffile.attrs['gattrone'] = 1
                cdffile.attrs['gattrtwo'] = 1
                cdffile.attrs['gattrfour'] = 4
                cdffile['var1'] = numpy.array([4, 5, 6], dtype=numpy.float32)
                cdffile['var2'] = numpy.array([14, 15, 16], dtype=numpy.float32)
                cdffile['var3'] = numpy.array([24, 25, 26], dtype=numpy.float32)
                cdffile['var2'].attrs['foo'] = 'file2'
                cdffile['Epoch'] = numpy.array([datetime.datetime(2010, 1, i)
                                            for i in range(4, 7)])
                cdffile.new('var4', data=numpy.array(
                    [101, 102], dtype=numpy.float32), recVary=False)
            with cdf.CDF(os.path.join(td, 'one.cdf')) as cdf1:
                with cdf.CDF(os.path.join(td, 'two.cdf')) as cdf2:
                    data = cdf.concatCDF(
                        [cdf1, cdf2], ['var1', 'var2', 'var4', 'Epoch'],
                        raw=True)
        finally:
            del warnings.filters[0:2]
            shutil.rmtree(td)
        self.assertEqual(
            ['gattrone', 'gattrthree', 'gattrtwo'],
            sorted(data.attrs.keys()))
        self.assertEqual(['Epoch', 'var1', 'var2', 'var4'], sorted(data.keys()))
        self.assertEqual(['foo'], list(data['var2'].attrs.keys()))
        self.assertEqual(b'file1', data['var2'].attrs['foo']) #raw variable
        numpy.testing.assert_array_equal(
            data['var1'][...],
            numpy.array([1, 2, 3, 4, 5, 6], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            data['var2'][...],
            numpy.array([11, 12, 13, 14, 15, 16], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            data['var4'][...],
            numpy.array([99, 100], dtype=numpy.float32))
        numpy.testing.assert_array_equal(
            data['Epoch'][...],
            cdf.lib.v_datetime_to_tt2000([datetime.datetime(2010, 1, i)
                                                   for i in range(1, 7)]))


class MakeCDF(unittest.TestCase):
    def setUp(self):
        self.testdir = tempfile.mkdtemp()
        self.testfspec = os.path.join(self.testdir, 'foo.cdf')
        self.testmaster = os.path.join(spacepy_testing.testsdir,
                                       'po_l1_cam_testc.cdf')

    def tearDown(self):
        shutil.rmtree(self.testdir)

    def testOpenCDFNew(self):
        """Create a new CDF"""

        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            newcdf = cdf.CDF(self.testfspec, '')
        finally:
            del warnings.filters[0]
        self.assertTrue(os.path.isfile(self.testfspec))
        self.assertFalse(newcdf.readonly())
        newcdf.close()
        os.remove(self.testfspec)

    def testCreateCDFKeyword(self):
        """Create a CDF specifying the create keyword"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            newcdf = cdf.CDF(self.testfspec, create=True)
        finally:
            del warnings.filters[0]
        self.assertTrue(os.path.isfile(self.testfspec))
        self.assertFalse(newcdf.readonly())
        newcdf.close()
        os.remove(self.testfspec)

    def testOpenCDFNonexistent(self):
        """Open a CDF which doesn't exist"""

        self.assertRaises(cdf.CDFError, cdf.CDF, self.testfspec)

    def testOpenCDFNoMaster(self):
        """Open a CDF from a master CDF which doesn't exist"""

        self.assertRaises(IOError, cdf.CDF, self.testfspec, 'nonexist.cdf')

    def testCDFNewMajority(self):
        """Creates a new CDF and changes majority"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            newcdf = cdf.CDF(self.testfspec, '')
        finally:
            del warnings.filters[0]
        newcdf.col_major(True)
        self.assertTrue(newcdf.col_major())
        newcdf.col_major(False)
        self.assertFalse(newcdf.col_major())
        newcdf.close()
        os.remove(self.testfspec)

    def testCreateCDFFromMaster(self):
        """Create a CDF from a master"""
        newcdf = cdf.CDF(self.testfspec, self.testmaster)
        self.assertTrue('ATC' in newcdf)
        self.assertFalse(newcdf.readonly())
        newcdf.close()
        os.remove(self.testfspec)

    def testCreateCDFBackward(self):
        """Try a backward-compatible CDF"""
        cdf.lib.set_backward(True)
        newcdf = cdf.CDF(self.testfspec, '')
        (ver, rel, inc) = newcdf.version()
        backward = newcdf.backward
        newcdf.close()
        os.remove(self.testfspec)
        self.assertEqual(2, ver)
        self.assertTrue(backward)

        cdf.lib.set_backward(False)
        newcdf = cdf.CDF(self.testfspec, '')
        (ver, rel, inc) = newcdf.version()
        backward = newcdf.backward
        newcdf.close()
        os.remove(self.testfspec)
        self.assertEqual(3, ver)
        self.assertFalse(backward)

    def testNewEPOCHAssign(self):
        """Create a new epoch variable by assigning to a CDF element"""
        cdf.lib.set_backward(True)
        newcdf = cdf.CDF(self.testfspec, '')
        data = [datetime.datetime(2000, 1, 1, 0, 0, 0, 999999),
                datetime.datetime(2001, 1, 1, 0, 0, 0, 999999)]
        newcdf['newzVar'] = data
        newtype = newcdf['newzVar'].type()
        newdata = newcdf['newzVar'][...]
        newcdf.close()
        os.remove(self.testfspec)
        self.assertEqual(const.CDF_EPOCH.value, newtype)
        numpy.testing.assert_array_equal(
            [datetime.datetime(2000, 1, 1, 0, 0, 1),
             datetime.datetime(2001, 1, 1, 0, 0, 1)],
            newdata)
        cdf.lib.set_backward(False)  # Revert to default

    def testCreateCDFLeak(self):
        """Make a CDF that doesn't get collected"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            newcdf = cdf.CDF(self.testfspec, '')
        finally:
            del warnings.filters[0]
        newcdf.close()
        gc.collect()
        old_garblen = len(gc.garbage)
        del newcdf
        os.remove(self.testfspec)
        gc.collect()
        new_garblen = len(gc.garbage)
        self.assertEqual(old_garblen, new_garblen)

    def testCreateCDFFromSpaceData(self):
        """Make a CDF from a Spacedata"""
        sd = datamodel.SpaceData(
            {
            'Epoch': datamodel.dmarray([datetime.datetime(2011, 1, 1),
                                        datetime.datetime(2011, 1, 2)],
                                       attrs={'min':
                                              datetime.datetime(2011, 1, 1)}),
            'flux': datamodel.dmarray([5.0, 6.0], dtype=numpy.float64,
                                      attrs={'type': 'data'}),
            },
            attrs={'project': 'junk'}
            )
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        warnings.filterwarnings(
            'ignore', r'^No type specified for time input.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            cdf.CDF.from_data(self.testfspec, sd)
        finally:
            del warnings.filters[0:2]
        with cdf.CDF(self.testfspec) as cdffile:
            self.assertEqual(['project'], list(cdffile.attrs.keys()))
            self.assertEqual(['min'], list(cdffile['Epoch'].attrs.keys()))
            self.assertEqual(['type'], list(cdffile['flux'].attrs.keys()))
            numpy.testing.assert_array_equal([datetime.datetime(2011, 1, 1)],
                                             cdffile['Epoch'].attrs['min'])
            self.assertEqual('data', cdffile['flux'].attrs['type'])
            numpy.testing.assert_array_equal(
                [datetime.datetime(2011, 1, 1), datetime.datetime(2011, 1, 2)],
                cdffile['Epoch'][...])
            numpy.testing.assert_array_equal(
                [5.0, 6.0], cdffile['flux'][...])
            self.assertEqual(cdffile['flux'].dtype, numpy.float64)

    def testEPOCH16inBackward(self):
        """Create backward-compatible CDF with EPOCH16"""
        msg = 'Cannot use EPOCH16, INT8, or TIME_TT2000 ' \
            'in backward-compatible CDF'
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            newcdf = cdf.CDF(self.testfspec, '')
        finally:
            del warnings.filters[0]
        try:
            newcdf.new('foo', type=const.CDF_EPOCH16)
        except ValueError:
            self.assertEqual(msg, str(sys.exc_info()[1]))
        else:
            self.fail('Should have raised ValueError: ' + msg)
        newcdf.close()
        os.remove(self.testfspec)

    def testInt64inBackward(self):
        """Create backward-compatible CDF with INT8"""
        if not cdf.lib.supports_int8:
            return
        msg = 'Data requires EPOCH16, INT8, or TIME_TT2000; ' \
            'incompatible with backward-compatible CDF'
        cdf.lib.set_backward(True)
        newcdf = cdf.CDF(self.testfspec, '')
        try:
            newcdf.new('foo', data=numpy.array([1,2,3], dtype=numpy.int64))
        except ValueError:
            self.assertEqual(msg, str(sys.exc_info()[1]))
        else:
            self.fail('Should have raised ValueError: ' + msg)
        newcdf.close()
        os.remove(self.testfspec)
        cdf.lib.set_backward(False)  # Revert to default

    def testEPOCH16AttrinBackward(self):
        """Create backward-compatible CDF with EPOCH16 attribute"""
        cdf.lib.set_backward(True)
        newcdf = cdf.CDF(self.testfspec, '')
        try:
            newcdf.attrs['foo'] = datetime.datetime(
                9999, 12, 31, 23, 59, 59, 999999)
            self.assertEqual(cdf.const.CDF_EPOCH.value,
                             newcdf.attrs['foo'].type(0))
            self.assertEqual(
                datetime.datetime(9999, 12, 31, 23, 59, 59, 999000),
                newcdf.attrs['foo'][0])
            newcdf.attrs.new('bar', datetime.datetime(
                9999, 12, 31, 23, 59, 59, 999999))
            self.assertEqual(cdf.const.CDF_EPOCH.value,
                             newcdf.attrs['bar'].type(0))
            self.assertEqual(
                datetime.datetime(9999, 12, 31, 23, 59, 59, 999000),
                newcdf.attrs['bar'][0])
        finally:
            cdf.lib.set_backward(True)  # Revert to default
            newcdf.close()
            os.remove(self.testfspec)

    def testEntryType(self):
        """Entry type should match variable type in some cases"""
        #This is very hard to reproduce, thus creating a new CDF just for it
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, '') as f:
                f.new('one', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('two', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('three', data=numpy.array([1, 2, 3], dtype=numpy.uint8))
                self.assertEqual(const.CDF_UINT1.value, f['three'].type())
                for k in f:
                    f[k].attrs['foo'] = 5
                self.assertNotEqual(const.CDF_FLOAT.value,
                                    f['three'].attrs.type('foo'))
                self.assertEqual(const.CDF_UINT1.value,
                                 f['three'].attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEntryType2(self):
        """Entry type should match variable if no One True entry type"""
        #This is very hard to reproduce, thus creating a new CDF just for it
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, '') as f:
                f.new('one', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('two', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('three', data=numpy.array([1, 2, 3], dtype=numpy.uint8))
                self.assertEqual(const.CDF_UINT1.value, f['three'].type())
                f['one'].attrs.new('foo', 5, type=const.CDF_INT2)
                f['two'].attrs.new('foo', 5, type=const.CDF_INT4)
                f['three'].attrs['foo'] = 5
                self.assertNotEqual(const.CDF_FLOAT.value,
                                    f['three'].attrs.type('foo'))
                self.assertEqual(const.CDF_UINT1.value,
                                 f['three'].attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEntryType3(self):
        """Entry type should not match variable type in other cases"""
        #Another hard to reproduce
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, '') as f:
                f.new('one', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('three', data=numpy.array([1, 2, 3], dtype=numpy.uint8))
                f['one'].attrs.new('foo', data=5, type=const.CDF_INT2)
                f['three'].attrs['foo'] = 5
                self.assertEqual(const.CDF_INT2.value,
                                 f['three'].attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEntryType3WithNew(self):
        """Entry type should not match variable type in other cases"""
        #Another hard to reproduce
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, '') as f:
                f.new('one', data=numpy.array([1, 2, 3], dtype=numpy.float32))
                f.new('three', data=numpy.array([1, 2, 3], dtype=numpy.uint8))
                f['one'].attrs.new('foo', data=5, type=const.CDF_INT2)
                f['three'].attrs.new('foo', 5)
                self.assertEqual(const.CDF_INT2.value,
                                 f['three'].attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEntryType4(self):
        """Another case where Entry type should match variable type"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, create=True) as f:
                v = f.new('newvar', data=[1, 2, 3])
                v.attrs['foo'] = 5
                self.assertEqual(v.type(), v.attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEntryType4MultiElements(self):
        """Entry with multiple elements, type should match variable type"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, create=True) as f:
                v = f.new('newvar', data=[1, 2, 3])
                v.attrs['foo'] = [5, 3]
                self.assertEqual(v.type(), v.attrs.type('foo'))
        finally:
            del warnings.filters[0]

    def testEmptyNRV(self):
        """Read an empty NRV variable, should be empty"""
        #This is strictly a READ test, but creating a new CDF and
        #new variable is the easiest way to get to it
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(self.testfspec, '') as f:
                v = f.new('nrv_test', recVary=False, dims=[5, 3],
                          type=const.CDF_INT1)
                hslice = cdf._Hyperslice(v, (0, 0))
                self.assertEqual(3, hslice.dims)
                #This is 1 for both NRV and RV but it still raises index error,
                #in the actual __getitem__
                numpy.testing.assert_array_equal(hslice.counts, [1, 1, 1])
                numpy.testing.assert_array_equal(hslice.degen, [True, True, True])
                numpy.testing.assert_array_equal(hslice.dimsizes, [0, 5, 3])
                self.assertRaises(IndexError, operator.getitem, v, 0)

                hslice = cdf._Hyperslice(v, Ellipsis)
                self.assertEqual(3, hslice.dims)
                self.assertEqual((0, 0, 1, False),
                                 hslice.convert_range(None, None, None, 0))
                #For RV, this is zero, since it's a slice.
                #For NRV, this is a 1, since there's an implicit 0,
                #at the front.
                numpy.testing.assert_array_equal(hslice.counts, [1, 5, 3])
                numpy.testing.assert_array_equal(hslice.degen, [True, False, False])
                numpy.testing.assert_array_equal(hslice.dimsizes, [0, 5, 3])
                data = v[...]
                self.assertEqual((0, 0), data.shape)

                #One more test: NRV scalar with no records
                v = f.new('nrv_scalar', recVary=False, dims=[],
                          type=const.CDF_INT1)
                hslice = cdf._Hyperslice(v, Ellipsis)
                data = v[...]
                #TODO: This is an awful special case, but it's impossible to
                #have a SCALAR with no value! i.e. you cannot be both
                #zero-dimensional and empty
                self.assertEqual((0,), data.shape)
        finally:
            del warnings.filters[0]

    def testNoSetBackward(self):
        """Warn if create a CDF without explicitly setting backward/not"""
        # Awkward, but need to make sure the default state at load "knows" that
        # set_backward has not been called.
        cdf.lib = cdf.Library(libpath=cdf.lib, library=cdf._library)
        self.assertFalse(cdf.lib._explicit_backward)
        with spacepy_testing.assertWarns(
                self, 'always',
                r'spacepy\.pycdf\.lib\.set_backward not called\; making'
                r' v3-compatible CDF\.$',
                DeprecationWarning, r'spacepy\.pycdf$'):
            cdf.CDF(self.testfspec, create=True).close()
        with cdf.CDF(self.testfspec) as f:
            ver, rel, inc = f.version()
        self.assertEqual(3, ver) # Still the default

    def testSetBackward(self):
        """But no warn if explicit set"""
        # Awkward, but need to make sure the default state at load "knows" that
        # set_backward has not been called.
        cdf.lib = cdf.Library(libpath=cdf.lib, library=cdf._library)
        self.assertFalse(cdf.lib._explicit_backward)
        cdf.lib.set_backward(True)
        with spacepy_testing.assertDoesntWarn(
                self, 'always', category=DeprecationWarning,
                module=r'spacepy\.pycdf$'):
            cdf.CDF(self.testfspec, create=True).close()
        with cdf.CDF(self.testfspec) as f:
            ver, rel, inc = f.version()
        self.assertEqual(2, ver)
        # Revert to the default
        cdf.lib.set_backward(False)

    def testSetBackwardFalse(self):
        """But no warn if explicit set"""
        cdf.lib = cdf.Library(libpath=cdf.lib, library=cdf._library)
        self.assertFalse(cdf.lib._explicit_backward)
        cdf.lib.set_backward(False)
        with spacepy_testing.assertDoesntWarn(
                self, 'always', category=DeprecationWarning,
                module=r'spacepy\.pycdf$'):
            cdf.CDF(self.testfspec, create=True).close()
        with cdf.CDF(self.testfspec) as f:
            ver, rel, inc = f.version()
        self.assertEqual(3, ver)


class CDFTestsBase(unittest.TestCase):
    """Base class for tests involving existing CDF, column or row major"""
    def __init__(self, *args, **kwargs):
        self.testfile = os.path.join(tempfile.gettempdir(), self.testbase)
        assert(self.calcDigest(self.testmaster) == self.expected_digest)
        super(CDFTestsBase, self).__init__(*args, **kwargs)

    @staticmethod
    def calcDigest(file):
        m = hashlib.md5()
        with open(file, 'rb') as f:
            m.update(f.read())
        return m.hexdigest()


class CDFTests(CDFTestsBase):
    """Tests that involve an existing CDF, read or write"""
    testmaster = os.path.join(spacepy_testing.testsdir, 'po_l1_cam_test.cdf')
    testbase = 'test.cdf'
    expected_digest = '94515e62d38a31ad02f6d435274cbfe7'


class ColCDFTests(CDFTestsBase):
    """Tests that involve an existing column-major CDF, read or write"""
    testmaster = os.path.join(spacepy_testing.testsdir, 'po_l1_cam_testc.cdf')
    testbase = 'testc.cdf'
    expected_digest = '7728439e20bece4c0962a125373345bf'


class OpenCDF(CDFTests):
    """Tests that open a CDF"""
    def setUp(self):
        shutil.copy(self.testmaster, self.testfile)

    def tearDown(self):
        os.remove(self.testfile)

    def testopenUnicode(self):
        """Opens a CDF providing a Unicode name"""
        try:
            cdffile = cdf.CDF(unicode(self.testfile))
        except NameError: #Py3k, all strings are unicode
            cdffile = cdf.CDF(self.testfile)
        cdffile.close()
        del cdffile

    def testcreateMaster(self):
        """Creates a new CDF from a master"""
        testfspec = 'foo.cdf'
        new = cdf.CDF(testfspec, self.testfile)
        new.close()
        self.assertTrue(os.path.isfile(testfspec))
        self.assertEqual(self.calcDigest(testfspec), self.calcDigest(self.testfile))
        os.remove(testfspec)

    def testcreateMasterExisting(self):
        """Creates a new CDF from a master, on top of an existing"""
        testfspec = 'foo.cdf'
        open(testfspec, 'w').close()
        errstr = 'CDF_EXISTS: The CDF named already exists.'
        try:
            new = cdf.CDF(testfspec, self.testfile)
        except cdf.CDFError:
            self.assertEqual(sys.exc_info()[1].__str__(),
                             errstr)
        else:
            self.fail('Should have raised CDFError: ' +
                      errstr)
        os.remove(testfspec)

    def testContextManager(self):
        expected = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
                    'RateScalerNames', 'SectorRateScalerNames',
                    'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
                    'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
                    'MajorNumbers', 'MeanCharge', 'Epoch', 'Epoch2D',
                    'String1D', 'StringUnpadded']
        with cdf.CDF(self.testfile) as f:
            names = list(f.keys())
        self.assertEqual(expected, names)
        self.assertRaises(cdf.CDFError, f.close)

    def testOpenCDFLeak(self):
        """Open a CDF that doesn't get collected"""
        cdffile = cdf.CDF(self.testfile)
        cdffile.close()
        gc.collect()
        old_garblen = len(gc.garbage)
        del cdffile
        gc.collect()
        new_garblen = len(gc.garbage)
        self.assertEqual(old_garblen, new_garblen)


class ReadCDF(CDFTests):
    """Tests that read an existing CDF, but do not modify it."""
    longMessage = True
    testbase = 'test_ro.cdf'
    varnames = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
               'RateScalerNames', 'SectorRateScalerNames',
               'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
               'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
               'MajorNumbers', 'MeanCharge', 'Epoch', 'Epoch2D',
               'String1D', 'StringUnpadded']

    def __init__(self, *args, **kwargs):
        super(ReadCDF, self).__init__(*args, **kwargs)
        #Unittest docs say 'the order in which the various test cases will be
        #run is determined by sorting the test function names with the built-in
        #cmp() function'
        testnames = [name for name in dir(self)
                     if name[0:4] == 'test' and
                     isinstance(getattr(self,name), Callable)]
        self.last_test = max(testnames)

    def setUp(self):
        super(ReadCDF, self).setUp()
        if not os.path.exists(self.testfile):
            shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)

    def tearDown(self):
        self.cdf.close()
        del self.cdf
        if self._testMethodName == self.last_test:
            os.remove(self.testfile)
        super(ReadCDF, self).tearDown()

    def testGetATC(self):
        """Get ATC zVar using subscripting"""
        atc = self.cdf['ATC']
        self.assertEqual(type(atc), cdf.Var)

    def testGetATCByNum(self):
        """Get ATC zVar using subscripting by variable number"""
        atc = self.cdf[0]
        self.assertEqual(type(atc), cdf.Var)
        self.assertEqual(atc.name(), 'ATC')

    def testGetAllzVars(self):
        """Check getting a list of zVars"""
        expectedNames = self.varnames
        names = [zVar.name() for zVar in self.cdf.values()]
        self.assertEqual(names, expectedNames)

    def testGetAllVarNames(self):
        """Getting a list of zVar names"""
        expectedNames = self.varnames
        names = list(self.cdf.keys())
        self.assertEqual(expectedNames, names)

    def testGetVarNum(self):
        self.assertEqual(0, self.cdf['ATC']._num())

    def testCDFIterator(self):
        expected = self.varnames
        self.assertEqual(expected, [i for i in self.cdf])
        a = self.cdf.__iter__()
        a.send(None)
        self.assertEqual('SectorNumbers', a.send('SpinNumbers'))
        try:
            res = a.next()
        except AttributeError:
            res = next(a)
        self.assertEqual('RateScalerNames', res)

    def testRecCount(self):
        """Get number of records in a zVariable"""
        self.assertEqual(len(self.cdf['ATC']), 100)
        self.assertEqual(len(self.cdf['MeanCharge']), 100)
        self.assertEqual(len(self.cdf['SpinNumbers']), 1)

    def testMajority(self):
        """Get majority of the CDF"""
        self.assertFalse(self.cdf.col_major())

    def testgetndims(self):
        """Get number of dimensions in zVar"""
        expected = {'ATC': 0, 'PhysRecNo': 0, 'SpinNumbers': 1,
                    'SectorNumbers': 1, 'RateScalerNames': 1,
                    'SectorRateScalerNames': 1,
                    'SectorRateScalersCounts': 3, 'SectorRateScalersCountsSigma': 3,
                    'SpinRateScalersCounts': 2, 'SpinRateScalersCountsSigma': 2}
        for i in expected:
            self.assertEqual(self.cdf[i]._n_dims(), expected[i])

    def testgetdimsizes(self):
        """Get size of dimensions in zVar"""
        expected = {'ATC': [], 'PhysRecNo': [], 'SpinNumbers': [18],
                    'SectorNumbers': [32], 'RateScalerNames': [16],
                    'SectorRateScalerNames': [9],
                    'SectorRateScalersCounts': [18, 32, 9],
                    'SectorRateScalersCountsSigma': [18, 32, 9],
                    'SpinRateScalersCounts': [18, 16],
                    'SpinRateScalersCountsSigma': [18, 16]}
        for i in expected:
            self.assertEqual(self.cdf[i]._dim_sizes(), expected[i])

    def testShape(self):
        """Get numpy-like shape (n_recs plus dimensions) in zvar"""
        expected = {'ATC': (100,), 'PhysRecNo': (100,), 'SpinNumbers': (18,),
                    'SectorNumbers': (32,), 'RateScalerNames': (16,),
                    'SectorRateScalerNames': (9,),
                    'SectorRateScalersCounts': (100, 18, 32, 9),
                    'SectorRateScalersCountsSigma': (100, 18, 32, 9),
                    'SpinRateScalersCounts': (100, 18, 16),
                    'SpinRateScalersCountsSigma': (100, 18, 16)}
        for i in expected:
            self.assertEqual(self.cdf[i].shape, expected[i])
            
    def testgetrecvary(self):
        """Get record variance of zVar"""
        expected = {'ATC': True, 'PhysRecNo': True, 'SpinNumbers': False,
                    'SectorNumbers': False, 'RateScalerNames': False,
                    'SectorRateScalerNames': False,
                    'SectorRateScalersCounts': True,
                    'SectorRateScalersCountsSigma': True,
                    'SpinRateScalersCounts': True,
                    'SpinRateScalersCountsSigma': True}
        for i in expected:
            self.assertEqual(self.cdf[i].rv(), expected[i])

    def testHyperslices(self):
        slices = {'ATC': 1,
                  'PhysRecNo': slice(10, 2, -2),
                  'SpinNumbers': slice(2, None, 2),
                  'SectorRateScalersCounts': (slice(3, 6, None),
                                              slice(None, None, None),
                                              slice(None, None, None)),
                  'SpinRateScalersCounts': (Ellipsis, slice(-1, None, -1)),
                  'MeanCharge': (0, -1)
                  } #Slice objects indexed by variable
        #Expected results [dims, dimsizes, starts, counts, intervals, degen, rev]
        #indexed by variable
        expected = {'ATC': [1, [100], [1], [1], [1], [True], [False]],
                    'PhysRecNo': [1, [100], [4], [4], [2], [False], [True]],
                    'SpinNumbers': [2, [1, 18], [0, 2], [1, 8], [1, 2],
                                    [True, False], [False, False]],
                    'SectorRateScalersCounts': [4, [100, 18, 32, 9],
                                                [0, 3, 0, 0], [100, 3, 32, 9],
                                                [1, 1, 1, 1],
                                                [False, False, False, False],
                                                [False, False, False, False]],
                    'SpinRateScalersCounts': [3, [100, 18, 16],
                                              [0, 0, 0], [100, 18, 16],
                                              [1, 1, 1], [False, False, False],
                                              [False, False, True]],
                    'MeanCharge': [2, [100, 16], [0, 15], [1, 1], [1, 1],
                                   [True, True], [False, False]]
                    }
        for i in expected:
            zvar = self.cdf[i]
            sliced = cdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts.tolist(), sliced.intervals,
                      sliced.degen.tolist(), sliced.rev.tolist())
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)
        self.assertRaises(IndexError, cdf._Hyperslice,
                          self.cdf['ATC'], (1, 2))
        self.assertRaises(IndexError, cdf._Hyperslice,
                          self.cdf['ATC'], 800)
        self.assertRaises(IndexError, cdf._Hyperslice,
                          self.cdf['ATC'], -1000)

    def testHyperslices2(self):
        """Additional checks: converting python slices to CDF counts, etc."""
        slices = {'ATC': Ellipsis,
                  } #Slice objects indexed by variable
        #Expected results [dims, dimsizes, starts, counts, intervals, degen, rev]
        #indexed by variable
        expected = {'ATC': [1, [100], [0], [100], [1], [False], [False]],
                    }
        for i in expected:
            zvar = self.cdf[i]
            sliced = cdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts, sliced.intervals, sliced.degen,
                      sliced.rev)
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)

    def testHypersliceExpand(self):
        """Expand a slice to store the data passed in"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._Hyperslice(zvar, slice(0, None, 1))
        self.assertEqual(100, sliced.counts[0])
        sliced.expand(list(range(110)))
        self.assertEqual(110, sliced.counts[0])
        sliced = cdf._Hyperslice(zvar, slice(0, 100, 2))
        sliced.expand(list(range(110)))
        self.assertEqual(50, sliced.counts[0])

    def testHypersliceExpectedDims(self):
        """Find dimensions expected by a slice"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._Hyperslice(zvar, slice(0, None, 1))
        self.assertEqual([100], sliced.expected_dims())
        sliced.expand(list(range(110)))
        self.assertEqual([110], sliced.expected_dims())
        sliced = cdf._Hyperslice(zvar, slice(0, 100, 2))
        sliced.expand(list(range(110)))
        self.assertEqual([50], sliced.expected_dims())

        zvar = self.cdf['SpinRateScalersCounts']
        sliced = cdf._Hyperslice(zvar, (slice(None, None, None),
                                               slice(None, None, 2),
                                               slice(0, None, 3)))
        self.assertEqual([100, 9, 6], sliced.expected_dims())

        zvar = self.cdf['SpinNumbers']
        sliced = cdf._Hyperslice(zvar, 2)
        self.assertEqual([1, 18], sliced.dimsizes)

    def testCDFTypes(self):
        """Look up variable type from the CDF"""
        expected = {'ATC': cdf.const.CDF_EPOCH16,
                    'PhysRecNo': cdf.const.CDF_INT4,
                    'SpinNumbers': cdf.const.CDF_CHAR,
                    'MeanCharge': cdf.const.CDF_FLOAT,
                    'Epoch': cdf.const.CDF_EPOCH,
                    }
        for i in expected:
            self.assertEqual(expected[i].value,
                             self.cdf[i].type())

    def testNPTypesInternal(self):
        """Look up numpy type to match internal representation of variable"""
        expected = {'ATC': numpy.dtype((numpy.float64, 2)),
                    'PhysRecNo': numpy.int32,
                    'SpinNumbers': numpy.dtype('S2'),
                    'MeanCharge': numpy.float32,
                    'Epoch': numpy.float64,
                    }
        for i in expected:
            self.assertEqual(expected[i],
                             self.cdf[i]._np_type())
            self.assertEqual(expected[i],
                             self.cdf.raw_var(i).dtype)

    def testNPTypes(self):
        """Look up numpy types to match variable"""
        expected = {'ATC': numpy.dtype('O'),
                    'PhysRecNo': numpy.int32,
                    'SpinNumbers': (numpy.dtype('S2') if str is bytes
                                    else numpy.dtype('U2')),
                    'MeanCharge': numpy.float32,
                    'Epoch': numpy.dtype('O'),
                    }
        for i in expected:
            self.assertEqual(expected[i],
                             self.cdf[i].dtype)

    def testSubscriptVariable(self):
        """Refer to an array by subscript"""
        numpy.testing.assert_array_equal([3, 25, 47],
                                         self.cdf['PhysRecNo'][0:5:2])
        numpy.testing.assert_array_equal([1094, 1083, 1072, 1061],
                                         self.cdf['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         self.cdf['SpinRateScalersCounts'][41, 2, 15])

    def testIncompleteSubscript(self):
        """Get data from a variable with a less-than-complete specification"""
        chargedata = self.cdf['MeanCharge'][0] #Should be the first record
        self.assertEqual(len(chargedata), 16)
        SpinRateScalersCounts = self.cdf['SpinRateScalersCounts'][...]
        self.assertEqual(100, len(SpinRateScalersCounts))

    def testEmptyResults(self):
        """Request an empty slice from a variable"""
        data = self.cdf['SectorRateScalersCounts'][1:1]
        self.assertEqual((0,18, 32, 9), data.shape)
        self.assertEqual(data.dtype, numpy.float32)

    def testReadEpochs(self):
        """Read an Epoch16 value"""
        expected = datetime.datetime(1998, 1, 15, 0, 0, 5, 334662)
        self.assertEqual(expected,
                         self.cdf['ATC'][0])
        expected = [datetime.datetime(1998, 1, 15, 0, 6, 48, 231),
                    datetime.datetime(1998, 1, 15, 0, 8, 30, 157015),
                    datetime.datetime(1998, 1, 15, 0, 10, 12, 313815),
                    datetime.datetime(1998, 1, 15, 0, 11, 54, 507400)
                    ]
        numpy.testing.assert_array_equal(
            expected, self.cdf['ATC'][4:8])

    def testReadEpoch8(self):
        """Read an Epoch value"""
        expected = datetime.datetime(1998, 1, 15, 0, 0, 0, 0)
        self.assertEqual(expected,
                         self.cdf['Epoch'][0])
        expected = [datetime.datetime(1998, 1, 15, 0, 4, 0, 0),
                    datetime.datetime(1998, 1, 15, 0, 5, 0, 0),
                    datetime.datetime(1998, 1, 15, 0, 6, 0, 0),
                    datetime.datetime(1998, 1, 15, 0, 7, 0, 0),
                    ]
        numpy.testing.assert_array_equal(
            expected, self.cdf['Epoch'][4:8])

    def testRead2DEpoch(self):
        """Read an Epoch16 variable with nonzero dimension"""
        expected = [[datetime.datetime(2000, 1, 1),
                     datetime.datetime(2000, 1, 1, 1)],
                    [datetime.datetime(2000, 1, 2),
                     datetime.datetime(2000, 1, 2, 1)],
                    [datetime.datetime(2000, 1, 3),
                     datetime.datetime(2000, 1, 3, 1)],
                    ]
        numpy.testing.assert_array_equal(expected, self.cdf['Epoch2D'][...])

    def testRead1DString(self):
        """Read a string with nonzero dimension"""
        expected = [['A', 'B', 'C'], ['D', 'E', 'F']]
        numpy.testing.assert_array_equal(expected, self.cdf['String1D'][...])

    def testnElems(self):
        """Read number of elements in a string variable"""
        self.assertEqual(2, self.cdf['SpinNumbers'].nelems())
        self.assertEqual(2, self.cdf['SectorNumbers'].nelems())

    def testSubscriptString(self):
        """Refer to a string array by subscript"""
        numpy.testing.assert_array_equal(
            ['0', '1', '2', '3', '4', '5', '6', '7',
             '8', '9', '10', '11', '12', '13', '14', '15',
             '16', '17'],
            numpy.char.rstrip(self.cdf['SpinNumbers'][:]))

    def testSubscriptIrregString(self):
        """Refer to a variable-length string array by subscript"""
        expected = ['H+', 'He+', 'He++', 'O<=+2', 'O>=+3', 'CN<=+2',
                    'H0', 'He0', 'CNO0', 'CN>=+3', 'Ne-Si', 'S-Ni',
                    '3He', 'D', 'Molecules', 'Others']
        out = self.cdf['RateScalerNames'][:]
        numpy.testing.assert_array_equal(expected, numpy.char.rstrip(out))

    def testGetAllNRV(self):
        """Get an entire non record varying variable"""
        numpy.testing.assert_array_equal(
            ['0', '1', '2', '3', '4', '5', '6', '7',
             '8', '9', '10', '11', '12', '13', '14', '15',
             '16', '17'],
            numpy.char.rstrip(self.cdf['SpinNumbers'][...]))

    def testGetsingleNRV(self):
        """Get single element of non record varying variable"""
        self.assertEqual('0 ',
                         self.cdf['SpinNumbers'][0])

    def testGetSingleNotArray(self):
        """Get a single element of a variable, not array type"""
        res = self.cdf['PhysRecNo'][0]
        self.assertEqual(3, res)
        self.assertFalse(isinstance(res, numpy.ndarray))

    def testcharType(self):
        """Get a CDF_CHAR variable and make sure it's a string"""
        self.assertTrue(isinstance(self.cdf['SpinNumbers'][0], str))

    def testGetVarUnicode(self):
        name = 'ATC'
        try:
            name = unicode(name)
        except NameError: #Py3k, all strings are unicode
            pass
        self.assertEqual(cdf.Var(self.cdf, name).name(), 'ATC')

    def testGetAllData(self):
        data = self.cdf.copy()
        expected = sorted(self.varnames)
        self.assertEqual(expected,
                         sorted([i for i in data]))

    def testCDFGetItem(self):
        """Look up a variable in CDF as a dict key"""
        result = self.cdf['ATC']
        self.assertEqual('ATC', result.name())
        self.assertRaises(KeyError, self.cdf.__getitem__, 'noexist')

    def testCDFlen(self):
        """length of CDF (number of zVars)"""
        result = len(self.cdf)
        self.assertEqual(len(self.varnames), result)

    def testReadonlyDefault(self):
        """CDF should be opened RO by default"""
        message = 'READ_ONLY_MODE: CDF is in read-only mode.'
        try:
            self.cdf['PhysRecNo']._delete()
        except cdf.CDFError:
            (type, val, traceback) = sys.exc_info()
            self.assertEqual(str(val), message)
        else:
            self.fail('Should have raised CDFError: '+ message)

    def testgAttrContains(self):
        """Check for attribute in global attrlist"""
        self.assertTrue('Mission_group' in self.cdf.attrs)
        self.assertFalse('DEPEND_0' in self.cdf.attrs)
        self.assertFalse('notanattratall' in self.cdf.attrs)

    def testzAttrContains(self):
        """Check for attribute in variable attrlist"""
        attrs = self.cdf['PhysRecNo'].attrs
        self.assertFalse('Mission_group' in attrs)
        self.assertTrue('DEPEND_0' in attrs)
        self.assertFalse('notanattratall' in attrs)
        self.assertFalse('LABL_PTR_1' in attrs)

    def testzEntryType(self):
        """Get the type of a zEntry"""
        names = ['DEPEND_0', 'VALIDMAX', ]
        numbers = [1, 0, ]
        types = [cdf.const.CDF_CHAR, cdf.const.CDF_EPOCH16, ]
        for (name, number, cdf_type) in zip(names, numbers, types):
            attribute = cdf.zAttr(self.cdf, name)
            actual_type = attribute.type(number)
            self.assertEqual(actual_type, cdf_type.value,
                             'zAttr ' + name + ' zEntry ' + str(number) +
                             ' ' + str(cdf_type.value) + ' != ' +
                             str(actual_type))
        self.assertEqual(cdf.const.CDF_CHAR.value,
                         self.cdf['PhysRecNo'].attrs.type('DEPEND_0'))

    def testgEntryType(self):
        """Get the type of a gEntry"""
        names = ['PI_name', 'Project', ]
        numbers = [0, 0, ]
        types = [cdf.const.CDF_CHAR, cdf.const.CDF_CHAR, ]
        for (name, number, cdf_type) in zip(names, numbers, types):
            attribute = cdf.gAttr(self.cdf, name)
            actual_type = attribute.type(number)
            self.assertEqual(actual_type, cdf_type.value,
                             'gAttr ' + name + ' gEntry ' + str(number) +
                             ' ' + str(cdf_type.value) + ' != ' +
                             str(actual_type))

    def testgEntryTypeBadNumber(self):
        """Get the type of nonexistent gEntry"""
        attr = cdf.gAttr(self.cdf, 'PI_name')
        self.assertRaises(IndexError, attr.type, 9999)

    def testzEntryNelems(self):
        """Get number of elements of a zEntry"""
        names = ['DEPEND_0', 'VALIDMAX', ]
        numbers = [1, 0, ]
        nelems = [3, 1, ]
        for (name, number, nelem) in zip(names, numbers, nelems):
            attribute = cdf.zAttr(self.cdf, name)
            actual_number = attribute._entry_len(number)
            self.assertEqual(actual_number, nelem,
                             'zAttr ' + name + ' zEntry ' + str(number) +
                             ' ' + str(nelem) + ' != ' + str(actual_number))

    def testgEntryNelems(self):
        """Get number of elements of a gEntry"""
        names = ['PI_name', 'Project', ]
        numbers = [0, 0, ]
        nelems = [8, 44, ]
        for (name, number, nelem) in zip(names, numbers, nelems):
            attribute = cdf.gAttr(self.cdf, name)
            actual_number = attribute._entry_len(number)
            self.assertEqual(actual_number, nelem,
                             'gAttr ' + name + ' gEntry ' + str(number) +
                             ' ' + str(nelem) + ' != ' + str(actual_number))

    def testzAttrLen(self):
        """Get number of zEntries for a zAttr"""
        names = ['DEPEND_0', 'VALIDMAX', ]
        lengths = [6, 8, ]
        for (name, length) in zip(names, lengths):
            attribute = cdf.zAttr(self.cdf, name)
            actual_length = len(attribute)
            self.assertEqual(actual_length, length,
                             'zAttr ' + name +
                             ' ' + str(length) + ' != ' + str(actual_length))

    def testgAttrLen(self):
        """Get number of gEntries for a gAttr"""
        names = ['PI_name', 'Project', ]
        lengths = [1, 1, ]
        for (name, length) in zip(names, lengths):
            attribute = cdf.gAttr(self.cdf, name)
            actual_length = len(attribute)
            self.assertEqual(actual_length, length,
                             'gAttr ' + name +
                             ' ' + str(length) + ' != ' + str(actual_length))

    def testGetAllgEntries(self):
        """Get all gEntries in a gAttr"""
        self.assertEqual(['T. Fritz'], self.cdf.attrs['PI_name'][:])
        self.assertEqual(['T. Fritz'], self.cdf.attrs['PI_name'][...])

    def testzEntryValue(self):
        """Get the value of a zEntry"""
        names = ['DEPEND_0', 'VALIDMAX', ]
        numbers = [1, 0]
        values = ['ATC', datetime.datetime(2009, 1, 1)]
        for (name, number, value) in zip(names, numbers, values):
            attribute = cdf.zAttr(self.cdf, name)
            entry = attribute._get_entry(number)
            self.assertEqual(value, entry)

    def testzEntryEpoch(self):
        """Get the value of an Epoch zEntry"""
        expected = datetime.datetime(2008, 12, 31, 23, 59, 59, 999000)
        actual = self.cdf['Epoch'].attrs['VALIDMAX']
        self.assertEqual(expected, actual)

    def testgEntryValue(self):
        """Get the value of a gEntry"""
        names = ['Project', 'TEXT', ]
        numbers = [0, 0, ]
        values = ['ISTP>International Solar-Terrestrial Physics',
                  'Polar CAMMICE Level One intermediate files', ]
        for (name, number, value) in zip(names, numbers, values):
            attribute = cdf.gAttr(self.cdf, name)
            entry = attribute._get_entry(number)
            self.assertEqual(value, entry)

    def testzAttrSlice(self):
        """Slice a zAttribute"""
        entries = cdf.zAttr(self.cdf, 'DEPEND_0')[6:10:2]
        values = [entry for entry in entries]
        self.assertEqual(['ATC', 'ATC'], values)

    def testgAttrSlice(self):
        """Slice a gAttribute"""
        entry = cdf.gAttr(self.cdf, 'Instrument_type')[0]
        value = entry
        self.assertEqual('Particles (space)', value)

    def testzAttrMaxIdx(self):
        """Find max index of a zAttr"""
        self.assertEqual(11,
            cdf.zAttr(self.cdf, 'DEPEND_0').max_idx())

    def testgAttrMaxIdx(self):
        """Find max index of a gAttr"""
        self.assertEqual(0,
            cdf.gAttr(self.cdf, 'Mission_group').max_idx())
        self.assertEqual(-1,
            cdf.gAttr(self.cdf, 'HTTP_LINK').max_idx())

    def testzEntryExists(self):
        """Checks for existence of a zEntry"""
        attribute = cdf.zAttr(self.cdf, 'DEPEND_0')
        self.assertFalse(attribute.has_entry(0))
        self.assertTrue(attribute.has_entry(6))

    def testGetBadzEntry(self):
        message = "'foobar: NO_SUCH_ATTR: Named attribute not found in this CDF.'"
        try:
            attrib = self.cdf['ATC'].attrs['foobar']
        except KeyError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should raise KeyError: ' + message)

        message = "'DEPEND_0: no such attribute for variable ATC'"
        try:
            attrib = self.cdf['ATC'].attrs['DEPEND_0']
        except KeyError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should raise KeyError: ' + message)

    def testgEntryExists(self):
        """Checks for existence of a gEntry"""
        attribute = cdf.gAttr(self.cdf, 'TEXT')
        self.assertTrue(attribute.has_entry(0))
        self.assertFalse(attribute.has_entry(6))

    def testzAttrIterator(self):
        """Iterate through all zEntries of a zAttr"""
        expected = ['ATC'] * 6
        attrib = cdf.zAttr(self.cdf, 'DEPEND_0')
        self.assertEqual(expected, [i for i in attrib])
        a = attrib.__iter__()
        a.send(None)
        res = a.send(1)
        self.assertEqual('ATC', res)
        try:
            res = a.next()
        except AttributeError:
            res = next(a)
        self.assertEqual('ATC', res)

    def testzAttrRevIterator(self):
        """Iterate backwards through all zEntries of a zAttr"""
        expected = ['ATC'] * 6
        attrib = cdf.zAttr(self.cdf, 'DEPEND_0')
        output = [entry for entry in reversed(attrib)]
        self.assertEqual(expected, output)

    def testgAttrIterator(self):
        """Iterate through all gEntries of a gAttr"""
        expected = ['ISTP>International Solar-Terrestrial Physics']
        attrib = cdf.gAttr(self.cdf, 'Project')
        self.assertEqual(expected, [i for i in attrib])

    def testzAttribList(self):
        """Get a zAttrib from the list on the CDF"""
        attrlist = cdf.zAttrList(self.cdf['PhysRecNo'])
        self.assertEqual(attrlist['DEPEND_0'], 'ATC')
        self.assertRaises(KeyError, attrlist.__getitem__, 'Data_type')

    def testgAttribList(self):
        """Get a gAttrib from the list on the CDF"""
        attrlist = cdf.gAttrList(self.cdf)
        self.assertEqual(attrlist['Data_type'][0],
                         'H0>High Time Resolution')
        self.assertRaises(KeyError, attrlist.__getitem__, 'DEPEND_0')

    def testAttribScope(self):
        """Get the variable/global scope of an attribute"""
        self.assertFalse(cdf.gAttr(self.cdf, 'CATDESC').global_scope())
        self.assertFalse(cdf.gAttr(self.cdf, 'DEPEND_0').global_scope())
        self.assertTrue(cdf.gAttr(self.cdf, 'Data_type').global_scope())

    def testAttribNumber(self):
        """Get the number of an attribute"""
        self.assertEqual(25, cdf.zAttr(self.cdf, 'CATDESC').number())
        self.assertEqual(26, cdf.zAttr(self.cdf, 'DEPEND_0').number())
        self.assertEqual(3, cdf.gAttr(self.cdf, 'Data_type').number())

    def testzAttribListLen(self):
        """Number of zAttrib for a zVar"""
        self.assertEqual(10, len(cdf.zAttrList(self.cdf['ATC'])))
        self.assertEqual(8, len(cdf.zAttrList(self.cdf['PhysRecNo'])))

    def testgAttribListLen(self):
        """Number of gAttrib in a CDF"""
        self.assertEqual(25, len(cdf.gAttrList(self.cdf)))

    def testzAttribsonVar(self):
        """Check zAttribs as an attribute of Var"""
        self.assertEqual(10, len(self.cdf['ATC'].attrs))
        self.assertEqual(8, len(self.cdf['PhysRecNo'].attrs))

    def testgAttribsonCDF(self):
        """Check gAttribs as an attribute of CDF"""
        self.assertEqual(25, len(self.cdf.attrs))

    def testzAttribListIt(self):
        """Iterate over keys in a zAttrList"""
        attrlist = cdf.zAttrList(self.cdf['PhysRecNo'])
        self.assertEqual(['CATDESC', 'DEPEND_0', 'FIELDNAM', 'FILLVAL',
                          'FORMAT', 'VALIDMIN', 'VALIDMAX', 'VAR_TYPE'],
                         list(attrlist))

    def testgAttribListIt(self):
        """Iterate over keys in a gAttrList"""
        attrlist = cdf.gAttrList(self.cdf)
        self.assertEqual(['Project', 'Source_name', 'Discipline',
                          'Data_type', 'Descriptor',
                          'File_naming_convention', 'Data_version',
                          'PI_name', 'PI_affiliation', 'TEXT',
                          'Instrument_type', 'Mission_group',
                          'Logical_source',
                          'Logical_file_id', 'Logical_source_description',
                          'Time_resolution', 'Rules_of_use', 'Generated_by',
                          'Generation_date', 'Acknowledgement', 'MODS',
                          'ADID_ref', 'LINK_TEXT', 'LINK_TITLE',
                          'HTTP_LINK',
                          ],
                         list(attrlist))

    def testgAttrsViaMeta(self):
        """Get gAttrs from the meta property"""
        self.assertEqual(
            'ISTP>International Solar-Terrestrial Physics',
            self.cdf.meta['Project'][0])
        self.assertTrue(self.cdf.attrs is self.cdf.meta)

    def testzAttrsViaMeta(self):
        """Get zAttrs from the meta property"""
        self.assertEqual(
            'ATC',
            self.cdf['ATC'].meta['FIELDNAM'])
        v = self.cdf['ATC']
        self.assertTrue(v.attrs is v.meta)

    def testzAttribListCopy(self):
        """Make a copy of a zAttr list"""
        attrs = self.cdf['PhysRecNo'].attrs
        attrcopy = attrs.copy()
        self.assertEqual(attrs, attrcopy)
        self.assertFalse(attrs is attrcopy)

    def testgAttribListCopy(self):
        """Copy a gAttr list"""
        attrs = self.cdf.attrs
        attrcopy = attrs.copy()
        for key in attrs:
            self.assertEqual(attrs[key][:], attrcopy[key])
            self.assertFalse(attrs[key] is attrcopy[key])

    def testgAttribListSame(self):
        """Are two instances of attributes from a CDF the same?"""
        attrs = self.cdf.attrs
        self.assertTrue(attrs is self.cdf.attrs)

    def testzAttribListSame(self):
        """Are two instances of attributes from a zVar the same?"""
        zv = self.cdf['PhysRecNo']
        attrs = zv.attrs
        self.assertTrue(attrs is zv.attrs)

    def testzVarCopy(self):
        """Make a copy of an entire zVar"""
        zvar = self.cdf['PhysRecNo']
        zvarcopy = zvar.copy()
        self.assertFalse(zvar is zvarcopy)
        for i in range(len(zvar)):
            self.assertEqual(zvar[i], zvarcopy[i])
        for i in zvarcopy.attrs:
            self.assertEqual(zvar.attrs[i], zvarcopy.attrs[i])
        numpy.testing.assert_array_equal(zvar[...], zvarcopy[...])

    def testVarCopyNRV(self):
        """VarCopy of NRV makes an NRV variable"""
        varcopy = self.cdf['RateScalerNames'].copy()
        testdir = tempfile.mkdtemp()
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(os.path.join(testdir, 'temp.cdf'), create=True) as f:
                f['newvar'] = varcopy
                self.assertFalse(f['newvar'].rv())
        finally:
            del warnings.filters[0]
            shutil.rmtree(testdir)

    @unittest.skipIf(cdf.lib.version[0] < 3,
                     "Not supported with CDF library < 3")
    def testVarCopyCDFType(self):
        """Assigning from VarCopy preserves CDF type"""
        varcopy = self.cdf['ATC'].copy()
        #Remove sub-second resolution so a normal Epoch would do
        for i in range(len(varcopy)):
            varcopy[i] = varcopy[i].replace(microsecond=0)
        testdir = tempfile.mkdtemp()
        try:
            with cdf.CDF(os.path.join(testdir, 'temp.cdf'), create=True) as f:
                f['newvar'] = varcopy
                self.assertEqual(self.cdf['ATC'].type(),
                                 f['newvar'].type())
        finally:
            shutil.rmtree(testdir)

    def testVarCopyMungeCDFType(self):
        """Change CDF type of VarCopy before assignment"""
        varcopy = self.cdf['MeanCharge'].copy()
        varcopy.set('type', const.CDF_DOUBLE)
        testdir = tempfile.mkdtemp()
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF(os.path.join(testdir, 'temp.cdf'), create=True) as f:
                f.new('newvar', data=varcopy)
                self.assertEqual(const.CDF_DOUBLE.value,
                                 f['newvar'].type())
                f['newvar2'] = varcopy
                self.assertEqual(const.CDF_DOUBLE.value,
                                 f['newvar2'].type())
        finally:
            del warnings.filters[0]
            shutil.rmtree(testdir)

    def testVarCopyBadAssign(self):
        """Assign to invalid CDF metadata"""
        varcopy = self.cdf['MeanCharge'].copy()
        self.assertRaises(KeyError, varcopy.set, 'foo', 'bar')

    def testVarCopyCompressThunk(self):
        """Check that compress overloading works"""
        varcopy = self.cdf['MeanCharge'].copy()
        actual = varcopy.compress()
        self.assertEqual(0, actual[0].value)
        self.assertEqual(0, actual[1])
        foo = varcopy.compress([False] * 17 + [True], axis=0)
        numpy.testing.assert_allclose(foo, self.cdf['MeanCharge'][17:18, :])

    def testVarCopyPickle(self):
        """Pickling a VarCopy preserves information"""
        varcopy = self.cdf['RateScalerNames'].copy()
        testdir = tempfile.mkdtemp()
        testfile = os.path.join(testdir, 'foo.pkl')
        try:
            with open(testfile, 'wb') as f:
                pickle.dump(varcopy, f)
            with open(testfile, 'rb') as f:
                varcopy2 = pickle.load(f)
        finally:
            shutil.rmtree(testdir)
        self.assertFalse(varcopy2.rv())

    def testCDFCopy(self):
        """Make a copy of an entire CDF"""
        cdfcopy = self.cdf.copy()
        self.assertFalse(cdfcopy is self.cdf)
        for key in self.cdf:
            orig = self.cdf[key][...]
            cp = cdfcopy[key]
            self.assertEqual(orig.shape, cp.shape)
            numpy.testing.assert_array_equal(orig.flatten(), cp.flatten())
            self.assertFalse(self.cdf[key] is cdfcopy[key])
        for key in self.cdf.attrs:
            self.assertEqual(self.cdf.attrs[key][:], cdfcopy.attrs[key])
            self.assertNotEqual(self.cdf.attrs[key], cdfcopy.attrs[key])

    def testSliceCDFCopy(self):
        """Slice a copy of a CDF"""
        cdfcopy = self.cdf.copy()
        numpy.testing.assert_array_equal([3, 25, 47],
                                         cdfcopy['PhysRecNo'][0:5:2])
        numpy.testing.assert_array_equal([1094, 1083, 1072, 1061],
                                         cdfcopy['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         cdfcopy['SpinRateScalersCounts'][41, 2, 15])

    def testVarString(self):
        """Convert a variable to a string representation"""
        expected = {'StringUnpadded': 'CDF_CHAR*6 [12]', 'String1D': 'CDF_CHAR*1 [2, 3]', 'SectorRateScalerNames': 'CDF_CHAR*9 [9] NRV', 'PhysRecNo': 'CDF_INT4 [100]', 'RateScalerNames': 'CDF_CHAR*9 [16] NRV', 'SpinRateScalersCountsSigma': 'CDF_FLOAT [100, 18, 16]', 'SectorRateScalersCountsSigma': 'CDF_FLOAT [100, 18, 32, 9]', 'SpinRateScalersCounts': 'CDF_FLOAT [100, 18, 16]', 'SpinNumbers': 'CDF_CHAR*2 [18] NRV', 'Epoch': 'CDF_EPOCH [11]', 'SectorRateScalersCounts': 'CDF_FLOAT [100, 18, 32, 9]', 'MeanCharge': 'CDF_FLOAT [100, 16]', 'SectorNumbers': 'CDF_CHAR*2 [32] NRV', 'MajorNumbers': 'CDF_CHAR*2 [11] NRV', 'Epoch2D': 'CDF_EPOCH16 [3, 2]', 'ATC': 'CDF_EPOCH16 [100]'}
        actual = dict([(varname, str(zVar))
                       for (varname, zVar) in self.cdf.items()])
        self.assertEqual(expected, actual)

    def testCDFString(self):
        """Convert a CDF to a string representation"""
        expected = 'ATC: CDF_EPOCH16 [100]\nEpoch: CDF_EPOCH [11]\nEpoch2D: CDF_EPOCH16 [3, 2]\nMajorNumbers: CDF_CHAR*2 [11] NRV\nMeanCharge: CDF_FLOAT [100, 16]\nPhysRecNo: CDF_INT4 [100]\nRateScalerNames: CDF_CHAR*9 [16] NRV\nSectorNumbers: CDF_CHAR*2 [32] NRV\nSectorRateScalerNames: CDF_CHAR*9 [9] NRV\nSectorRateScalersCounts: CDF_FLOAT [100, 18, 32, 9]\nSectorRateScalersCountsSigma: CDF_FLOAT [100, 18, 32, 9]\nSpinNumbers: CDF_CHAR*2 [18] NRV\nSpinRateScalersCounts: CDF_FLOAT [100, 18, 16]\nSpinRateScalersCountsSigma: CDF_FLOAT [100, 18, 16]\nString1D: CDF_CHAR*1 [2, 3]\nStringUnpadded: CDF_CHAR*6 [12]'
        actual = str(self.cdf)
        self.assertEqual(expected, actual)

    def testgAttrListString(self):
        """Convert a list of gattributes to a string"""
        expected = 'ADID_ref: NSSD0241 [CDF_CHAR]\nAcknowledgement: \nData_type: H0>High Time Resolution [CDF_CHAR]\nData_version: 1 [CDF_CHAR]\nDescriptor: CAM>Charge and Mass Magnetospheric Ion Composition Experiment [CDF_CHAR]\nDiscipline: Space Physics>Magnetospheric Science [CDF_CHAR]\nFile_naming_convention: source_datatype_descriptor [CDF_CHAR]\nGenerated_by: BU Energetic Particle Group [CDF_CHAR]\nGeneration_date: 20100625 [CDF_CHAR]\nHTTP_LINK: \nInstrument_type: Particles (space) [CDF_CHAR]\nLINK_TEXT: \nLINK_TITLE: \nLogical_file_id: polar_h0_cam_00000000_v01 [CDF_CHAR]\nLogical_source: polar_h0_cam [CDF_CHAR]\nLogical_source_description: PO_L1_CAM [CDF_CHAR]\nMODS: \nMission_group: Polar [CDF_CHAR]\nPI_affiliation: Boston University [CDF_CHAR]\nPI_name: T. Fritz [CDF_CHAR]\nProject: ISTP>International Solar-Terrestrial Physics [CDF_CHAR]\nRules_of_use: \nSource_name: POLAR>POLAR PLASMA LABORATORY [CDF_CHAR]\nTEXT: Polar CAMMICE Level One intermediate files [CDF_CHAR]\n      another entry to simply pad it out [CDF_CHAR]\nTime_resolution: millisecond [CDF_CHAR]'
        self.assertEqual(expected, str(self.cdf.attrs))
        self.assertEqual('<gAttrList:\n' + expected + '\n>',
                         repr(self.cdf.attrs))

    def testzAttrListString(self):
        """Convert a list of zAttributes to a string"""
        expected = {
            'ATC': 'CATDESC: Absolute Time Code [CDF_CHAR]\nFIELDNAM: ATC [CDF_CHAR]\nFILLVAL: 9999-12-31 23:59:59.999999 [CDF_EPOCH16]\nLABLAXIS: UT [CDF_CHAR]\nMONOTON: INCREASE [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 2009-01-01 00:00:00 [CDF_EPOCH16]\nVALIDMIN: 1996-01-01 00:00:00 [CDF_EPOCH16]\nVAR_NOTES: Time when data in this master started accumulating. [CDF_CHAR]\nVAR_TYPE: support_data [CDF_CHAR]',
            'Epoch': 'CATDESC: Standard CDF Epoch time (8 byte) [CDF_CHAR]\nFIELDNAM: UTC [CDF_CHAR]\nFILLVAL: 9999-12-31 23:59:59.999000 [CDF_EPOCH]\nMONOTON: INCREASE [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 2008-12-31 23:59:59.999000 [CDF_EPOCH]\nVALIDMIN: 1996-01-01 00:00:00 [CDF_EPOCH]\nVAR_TYPE: support_data [CDF_CHAR]',
            'Epoch2D': '',
            'MajorNumbers': 'CATDESC: major frame number within the TM Master [CDF_CHAR]\nFIELDNAM: Major Frame Number [CDF_CHAR]\nFORMAT: A2 [CDF_CHAR]\nVAR_TYPE: metadata [CDF_CHAR]',
            'MeanCharge': 'CATDESC: Mean charge state [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nFIELDNAM: avg charge [CDF_CHAR]\nFILLVAL: -1e+31 [CDF_FLOAT]\nFORMAT: F3.1 [CDF_CHAR]\nLABL_PTR_1: RateScalerNames [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nUNITS: e [CDF_CHAR]\nVALIDMAX: 8.0 [CDF_FLOAT]\nVALIDMIN: 1.0 [CDF_FLOAT]\nVAR_NOTES: Mean charge state in each rate scaler. For the ENTIRE master period (i.e. summed over all energies), based on COUNTS. [CDF_CHAR]\nVAR_TYPE: support_data [CDF_CHAR]',
            'PhysRecNo': 'CATDESC: LZ record number for first major in this master [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nFIELDNAM: physical record [CDF_CHAR]\nFILLVAL: -2147483648 [CDF_INT4]\nFORMAT: I8 [CDF_CHAR]\nVALIDMAX: 20000 [CDF_INT4]\nVALIDMIN: 0 [CDF_INT4]\nVAR_TYPE: metadata [CDF_CHAR]',
            'RateScalerNames': 'CATDESC: Species found in each rate scaler [CDF_CHAR]\nFIELDNAM: Rate Scaler Names [CDF_CHAR]\nFORMAT: A10 [CDF_CHAR]\nVAR_NOTES: From J. Fennell revision 1997/02/28 [CDF_CHAR]\nVAR_TYPE: metadata [CDF_CHAR]',
            'SectorNumbers': 'CATDESC: Data accumulation sector number within the spin [CDF_CHAR]\nFIELDNAM: Sector Number [CDF_CHAR]\nFORMAT: A3 [CDF_CHAR]\nVAR_TYPE: metadata [CDF_CHAR]',
            'SectorRateScalersCounts': 'CATDESC: Counts in the per-sector rate scalers [CDF_CHAR]\nDELTA_MINUS_VAR: SectorRateScalersCountsSigma [CDF_CHAR]\nDELTA_PLUS_VAR: SectorRateScalersCountsSigma [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nDEPEND_1: SpinNumbers [CDF_CHAR]\nDEPEND_2: SectorNumbers [CDF_CHAR]\nDEPEND_3: SectorRateScalerNames [CDF_CHAR]\nDISPLAY_TYPE: time_series [CDF_CHAR]\nFIELDNAM: Sector rate scaler counts [CDF_CHAR]\nFILLVAL: -1e+31 [CDF_FLOAT]\nFORMAT: E6.2 [CDF_CHAR]\nLABL_PTR_1: SpinNumbers [CDF_CHAR]\nLABL_PTR_2: SectorNumbers [CDF_CHAR]\nLABL_PTR_3: SectorRateScalerNames [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 1e+06 [CDF_FLOAT]\nVALIDMIN: 0.0 [CDF_FLOAT]\nVAR_NOTES: Total counts accumulated over one sector (divide by SectorLength for rate, subtracting 58ms for sector 0). [CDF_CHAR]\nVAR_TYPE: data [CDF_CHAR]',
            'SectorRateScalersCountsSigma': 'CATDESC: Uncertainty in counts in the per-sector rate scalers. [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nFIELDNAM: Sector rate scaler uncertainty [CDF_CHAR]\nFILLVAL: -1e+31 [CDF_FLOAT]\nFORMAT: E12.2 [CDF_CHAR]\nLABL_PTR_1: SpinNumbers [CDF_CHAR]\nLABL_PTR_2: SectorNumbers [CDF_CHAR]\nLABL_PTR_3: SectorRateScalerNames [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 1e+12 [CDF_FLOAT]\nVALIDMIN: 0.0 [CDF_FLOAT]\nVAR_NOTES: Combines uncertainty from RS compression and Poisson stats. Total counts accumulated over one sector (divide by SectorLength for rate. Subtract 58ms for sector 0) [CDF_CHAR]\nVAR_TYPE: support_data [CDF_CHAR]',
            'SectorRateScalerNames': 'CATDESC: Species found in each per-sector rate scaler [CDF_CHAR]\nFIELDNAM: Sector Rate Scaler Names [CDF_CHAR]\nFORMAT: A10 [CDF_CHAR]\nVAR_NOTES: From J. Fennell revision 1997/02/28 [CDF_CHAR]\nVAR_TYPE: metadata [CDF_CHAR]',
            'SpinNumbers': 'CATDESC: Spin number within the TM Master [CDF_CHAR]\nFIELDNAM: Spin Number [CDF_CHAR]\nFORMAT: A3 [CDF_CHAR]\nVAR_TYPE: metadata [CDF_CHAR]',
            'SpinRateScalersCounts': 'CATDESC: Counts in the per-spin rate scalers [CDF_CHAR]\nDELTA_MINUS_VAR: SpinRateScalersCountsSigma [CDF_CHAR]\nDELTA_PLUS_VAR: SpinRateScalersCountsSigma [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nDEPEND_1: SpinNumbers [CDF_CHAR]\nDEPEND_2: RateScalerNames [CDF_CHAR]\nDISPLAY_TYPE: time_series [CDF_CHAR]\nFIELDNAM: Spin rate scaler number counts [CDF_CHAR]\nFILLVAL: -1e+31 [CDF_FLOAT]\nFORMAT: E6.2 [CDF_CHAR]\nLABL_PTR_1: SpinNumbers [CDF_CHAR]\nLABL_PTR_2: RateScalerNames [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 1e+06 [CDF_FLOAT]\nVALIDMIN: 0.0 [CDF_FLOAT]\nVAR_NOTES: Total counts accumulated over one spin (divide by SectorLength*32-58ms for rate). [CDF_CHAR]\nVAR_TYPE: data [CDF_CHAR]',
            'SpinRateScalersCountsSigma': 'CATDESC: Uncertainty in counts in the per-spin rate scalers. [CDF_CHAR]\nDEPEND_0: ATC [CDF_CHAR]\nFIELDNAM: Spin rate scaler uncertainty [CDF_CHAR]\nFILLVAL: -1e+31 [CDF_FLOAT]\nFORMAT: E6.2 [CDF_CHAR]\nLABL_PTR_1: SpinNumbers [CDF_CHAR]\nLABL_PTR_2: RateScalerNames [CDF_CHAR]\nSCALETYP: linear [CDF_CHAR]\nVALIDMAX: 1e+06 [CDF_FLOAT]\nVALIDMIN: 0.0 [CDF_FLOAT]\nVAR_NOTES: Combines uncertainty from RS compression and Poisson stats. Total counts accumulated over one spin (divide by SectorLength*32-58ms for rate). [CDF_CHAR]\nVAR_TYPE: support_data [CDF_CHAR]',
            'String1D': '',
            'StringUnpadded': '',
            }
        actual = dict([(varname, str(zVar.attrs))
                        for (varname, zVar) in self.cdf.items()])
        #Py3k and 2k display the floats differently,
        #as do numpy and Python
        ignorelist = ('SectorRateScalersCountsSigma',
                      'SpinRateScalersCountsSigma',
                      'SpinRateScalersCounts',
                      'SectorRateScalersCounts',
                      'MeanCharge',
                      )
        for k in ignorelist:
            del expected[k]
            del actual[k]
        self.assertEqual(expected, actual)
        for idx in expected:
            expected[idx] = '<zAttrList:\n' + expected[idx] + '\n>'
        actual = dict([(varname, repr(zVar.attrs))
                        for (varname, zVar) in self.cdf.items()])
        for k in ignorelist:
            del actual[k]
        self.assertEqual(expected, actual)

    def testReadClosedCDF(self):
        """Read a CDF that has been closed"""
        self.cdf.close()
        try:
            keylist = list(self.cdf.keys())
        except cdf.CDFError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(const.BAD_CDF_ID, v.status)
        else:
            self.fail('Should raise CDFError: BAD_CDF_ID')
        finally:
            self.cdf = cdf.CDF(self.testfile) #keep tearDown from failing

    def testStrClosedCDF(self):
        """String representation of CDF that has been closed"""
        self.cdf.close()
        cdflabel = str(self.cdf)
        try:
            self.assertEqual('Closed CDF', cdflabel[0:10])
            self.assertEqual(self.testbase, cdflabel[-len(self.testbase):])
        finally:
            self.cdf = cdf.CDF(self.testfile) #keep tearDown from failing

    def testStrClosedVar(self):
        """String representation of zVar in CDF that has been closed"""
        zVar = self.cdf['ATC']
        self.cdf.close()
        varlabel = str(zVar)
        try:
            self.assertEqual('zVar "ATC" in closed CDF ', varlabel[0:25])
            self.assertEqual(self.testbase, varlabel[-len(self.testbase):])
        finally:
            self.cdf = cdf.CDF(self.testfile) #keep tearDown from failing

    def testStrClosedAttribute(self):
        """String representation of attribute in CDF that has been closed"""
        attrib = self.cdf.attrs['Project']
        self.cdf.close()
        attrlabel = str(attrib)
        try:
            self.assertEqual('Attribute "Project" in closed CDF ', attrlabel[0:34])
            self.assertEqual(self.testbase, attrlabel[-len(self.testbase):])
        finally:
            self.cdf = cdf.CDF(self.testfile) #keep tearDown from failing

    def testStrClosedAttributeList(self):
        """String representation of attribute list in closed CDF"""
        al = self.cdf.attrs
        self.cdf.close()
        allabel = str(al)
        try:
            self.assertEqual('Attribute list in closed CDF ', allabel[0:29])
            self.assertEqual(self.testbase, allabel[-len(self.testbase):])
        finally:
            self.cdf = cdf.CDF(self.testfile) #keep tearDown from failing

    def testFancyIndexing(self):
        """See if fancy indexing works on a copy"""
        copy = self.cdf['PhysRecNo'].copy()
        numpy.testing.assert_array_equal(
            [3, 58],
            copy[numpy.array([True, False, False, False, False, True]
                             + [False] * 94)])

    def testReadEpoch16Raw(self):
        """Read an Epoch16 value, raw mode"""
        expected = numpy.array([63052041605.0, 334662000000.0],
                               dtype=numpy.float64)
        numpy.testing.assert_array_equal(
            expected, self.cdf.raw_var('ATC')[0])
        expected = numpy.array([
            [63052042008.0, 231000000.0],
            [63052042110.0, 157015000000.0],
            [63052042212.0, 313815000000.0],
            [63052042314.0, 507400000000.0]
            ], dtype=numpy.float64)
        numpy.testing.assert_array_equal(
            expected, self.cdf.raw_var('ATC')[4:8])

    def testReadEpoch8Raw(self):
        """Read an Epoch value, raw mode"""
        expected = 63052041600000.0
        self.assertEqual(expected,
                         self.cdf.raw_var('Epoch')[0])
        expected = numpy.array([63052041840000.0,
                                63052041900000.0,
                                63052041960000.0,
                                63052042020000.0,
                                ], dtype=numpy.float64)
        numpy.testing.assert_array_equal(
            expected, self.cdf.raw_var('Epoch')[4:8])

    def testReadEpoch8AttrRaw(self):
        """Read an Epoch attribute, raw mode"""
        expected = 62987673600000.0
        self.assertEqual(expected,
                         self.cdf.raw_var('Epoch').attrs['VALIDMIN'])

    def testReadCharRaw(self):
        """Read a string, raw mode"""
        #Verify that we're getting bytes, not unicode
        self.assertEqual('S',
                         self.cdf.raw_var('RateScalerNames')[...].dtype.char)

    def testReadCharConverted(self):
        """Read a string, not raw mode"""
        #verify getting unicode on py3k
        self.assertEqual('U' if str != bytes else 'S',
                         self.cdf['RateScalerNames'][...].dtype.char)

    def testDtypeChar(self):
        """Check dtype of a char"""
        if str is bytes: #py2k
            self.assertEqual('|S6', str(self.cdf['StringUnpadded'].dtype))
        else: #py3k
            self.assertEqual('|S6',
                             str(self.cdf.raw_var('StringUnpadded').dtype))
            self.assertEqual('U6', str(self.cdf['StringUnpadded'].dtype)[1:])

    def testReadCharUnpadded(self):
        """Read a string NOT padded out to full length"""
        self.assertEqual(
            '6', str(self.cdf.raw_var('StringUnpadded').dtype)[-1])
        self.assertEqual(
            '6', str(self.cdf.raw_var('StringUnpadded')[...].dtype)[-1])
        self.assertEqual('6', str(self.cdf['StringUnpadded'].dtype)[-1])
        self.assertEqual('6', str(self.cdf['StringUnpadded'][...].dtype)[-1])

    def testCachezVarNum(self):
        """Test basic reads of the zVar number cache"""
        self.assertFalse(b'ATC' in self.cdf._var_nums)
        self.assertEqual(0, self.cdf.var_num(b'ATC'))
        self.assertTrue(b'ATC' in self.cdf._var_nums)
        self.assertEqual(0, self.cdf._var_nums[b'ATC'])
        self.assertRaises(cdf.CDFError, self.cdf.var_num, b'foobar')

    def testDeleteCachezVarNum(self):
        """Test deleting something from the zvar cache"""
        #This is a bit backwards, but easiest way to fill the cache
        #without assuming it's being used
        for i in range(len(self.cdf)):
            n = self.cdf[i].name()
            if str is not bytes:
                n = n.encode('ascii')
            _ = self.cdf.var_num(n)
        #This is variable #8
        self.cdf.clear_from_cache(b'SpinRateScalersCounts')
        self.assertFalse(b'SpinRateScalersCounts' in self.cdf._var_nums)
        self.assertFalse(b'SpinRateScalersCountSigma' in self.cdf._var_nums) #9
        self.assertFalse(b'MajorNumbers' in self.cdf._var_nums) #10
        #7
        self.assertTrue(b'SectorRateScalersCountsSigma' in self.cdf._var_nums)
        self.assertEqual(10, self.cdf.var_num(b'MajorNumbers')) #Repopulated

    def testCachezAttrNum(self):
        """Test basic reads of the attr number cache"""
        self.assertFalse(b'Project' in self.cdf._attr_info)
        self.assertEqual((0, True), self.cdf.attr_num(b'Project'))
        self.assertTrue(b'Project' in self.cdf._attr_info)
        self.assertEqual((0, True), self.cdf._attr_info[b'Project'])
        self.assertRaises(cdf.CDFError, self.cdf.attr_num, b'foobar')
        #zAttr?
        self.assertEqual((41, False), self.cdf.attr_num(b'VALIDMIN'))

    def testDeleteCacheAttrNum(self):
        """Test deleting something from the attr cache"""
        #This is a bit backwards, but easiest way to fill the cache
        _ = [self.cdf.attr_num(self.cdf.attrs[i]._name)
             for i in range(len(self.cdf.attrs))]
        #This is attr #8
        self.cdf.clear_attr_from_cache(b'PI_affiliation')
        self.assertFalse(b'PI_affiliation' in self.cdf._attr_info)
        self.assertFalse(b'TEXT' in self.cdf._attr_info) #9
        self.assertFalse(b'Instrument_type' in self.cdf._attr_info) #10
        self.assertTrue(b'PI_name' in self.cdf._attr_info) #7
        #Repopulated
        self.assertEqual((10, True), self.cdf.attr_num(b'Instrument_type'))

    def testReadUnsetPad(self):
        """Test getting pad value on variable where it isn't set."""
        self.assertIs(self.cdf['Epoch'].pad(), None)

    def testPadNRV(self):
        """Read pad value for NRV"""
        # This wasn't explicitly set but that's what it gives
        self.assertEqual(self.cdf['SectorNumbers'].pad(), ' P')

    def testPrepare(self):
        """Test data conversion to numpy arrays"""
        # Data here are only prepared, CDF itself does not change
        # Each test is variable name, input data, expected prepared data,
        # expected dtype
        tests = [('ATC',
                  [datetime.datetime(1, 1, 1)],
                  numpy.array([[31622400.0, 0]]), numpy.float64),
                 ('MeanCharge',
                  [1, 2], numpy.array([1., 2.]), numpy.float32)
                 ]
        for vname, inputs, expected, dtype in tests:
            actual = self.cdf[vname]._prepare(inputs)
            self.assertIs(type(actual), numpy.ndarray, vname)
            self.assertEqual(actual.dtype, dtype, vname)
            numpy.testing.assert_array_equal(actual, expected)


class ReadColCDF(ColCDFTests):
    """Tests that read a column-major CDF, but do not modify it."""
    testbase = 'testc_ro.cdf'

    def __init__(self, *args, **kwargs):
        super(ReadColCDF, self).__init__(*args, **kwargs)
        #Unittest docs say 'the order in which the various test cases will be
        #run is determined by sorting the test function names with the built-in
        #cmp() function'
        testnames = [name for name in dir(self)
                     if name[0:4] == 'test' and
                     isinstance(getattr(self,name), Callable)]
        self.last_test = max(testnames)

    def setUp(self):
        super(ReadColCDF, self).setUp()
        if not os.path.exists(self.testfile):
            shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)

    def tearDown(self):
        self.cdf.close()
        del self.cdf
        if self._testMethodName == self.last_test:
            os.remove(self.testfile)
        super(ReadColCDF, self).tearDown()

    def testCMajority(self):
        """Get majority of the CDF"""
        self.assertTrue(self.cdf.col_major())

    def testCgetndims(self):
        """Get number of dimensions in zVar"""
        expected = {'ATC': 0, 'PhysRecNo': 0, 'SpinNumbers': 1,
                    'SectorNumbers': 1, 'RateScalerNames': 1,
                    'SectorRateScalerNames': 1,
                    'SectorRateScalersCounts': 3, 'SectorRateScalersCountsSigma': 3,
                    'SpinRateScalersCounts': 2, 'SpinRateScalersCountsSigma': 2}
        for i in expected:
            self.assertEqual(self.cdf[i]._n_dims(), expected[i])

    def testCgetdimsizes(self):
        """Get size of dimensions in zVar"""
        expected = {'ATC': [], 'PhysRecNo': [], 'SpinNumbers': [18],
                    'SectorNumbers': [32], 'RateScalerNames': [16],
                    'SectorRateScalerNames': [9],
                    'SectorRateScalersCounts': [18, 32, 9],
                    'SectorRateScalersCountsSigma': [18, 32, 9],
                    'SpinRateScalersCounts': [18, 16],
                    'SpinRateScalersCountsSigma': [18, 16]}
        for i in expected:
            self.assertEqual(self.cdf[i]._dim_sizes(), expected[i])

    def testColHyperslices(self):
        slices = {'ATC': 1,
                  'PhysRecNo': slice(10, 2, -2),
                  'SpinNumbers': slice(2, None, 2),
                  'SectorRateScalersCounts': (slice(3, 6, None),
                                              slice(None, None, None),
                                              slice(None, None, None)),
                  'SpinRateScalersCounts': (Ellipsis, slice(-1, None, -1)),
                  } #Slice objects indexed by variable
        #Expected results [dims, dimsizes, starts, counts, intervals, degen, rev]
        #indexed by variable
        expected = {'ATC': [1, [747], [1], [1], [1], [True], [False]],
                    'PhysRecNo': [1, [100], [4], [4], [2], [False], [True]],
                    'SpinNumbers': [2, [1, 18], [0, 2], [1, 8], [1, 2],
                                    [True, False], [False, False]],
                    'SectorRateScalersCounts': [4, [100, 18, 32, 9],
                                                [0, 3, 0, 0], [100, 3, 32, 9],
                                                [1, 1, 1, 1],
                                                [False, False, False, False],
                                                [False, False, False, False]],
                    'SpinRateScalersCounts': [3, [100, 18, 16],
                                              [0, 0, 0], [100, 18, 16],
                                              [1, 1, 1], [False, False, False],
                                              [False, False, True]],
                    }
        for i in expected:
            zvar = self.cdf[i]
            sliced = cdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts.tolist(), sliced.intervals,
                      sliced.degen.tolist(), sliced.rev.tolist())
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)

    def testColSubscriptVariable(self):
        """Refer to an column-major array by subscript"""
        #NB: Should be in SAME order as row-major,
        #since converted in convert_array
        numpy.testing.assert_array_equal([3, 25, 47],
                                         self.cdf['PhysRecNo'][0:5:2])
        numpy.testing.assert_array_equal([1094, 1083, 1072, 1061],
                                         self.cdf['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         self.cdf['SpinRateScalersCounts'][41, 2, 15])

    def testColSubscriptString(self):
        """Refer to a string array by subscript"""
        numpy.testing.assert_array_equal(
            ['0', '1', '2', '3', '4', '5', '6', '7',
             '8', '9', '10', '11', '12', '13', '14', '15',
             '16', '17'],
            numpy.char.rstrip(self.cdf['SpinNumbers'][:]))

    def testColSubscriptIrregString(self):
        """Refer to a variable-length string array by subscript"""
        expected = ['H+', 'He+', 'He++', 'O<=+2', 'O>=+3', 'CN<=+2',
                    'H0', 'He0', 'CNO0', 'CN>=+3', 'Ne-Si', 'S-Ni',
                    '3He', 'D', 'Molecules', 'Others']
        out = self.cdf['RateScalerNames'][:]
        numpy.testing.assert_array_equal(expected, numpy.char.rstrip(out))

    def testColReadEpochs(self):
        """Read an Epoch16 value"""
        expected = datetime.datetime(1998, 1, 15, 0, 0, 5, 334662)
        self.assertEqual(expected,
                         self.cdf['ATC'][0])
        expected = [datetime.datetime(1998, 1, 15, 0, 6, 48, 231),
                    datetime.datetime(1998, 1, 15, 0, 8, 30, 157015),
                    datetime.datetime(1998, 1, 15, 0, 10, 12, 313815),
                    datetime.datetime(1998, 1, 15, 0, 11, 54, 507400)
                    ]
        numpy.testing.assert_array_equal(expected,
                                         self.cdf['ATC'][4:8])

    def testContains(self):
        """See if variable exists in CDF"""
        self.assertTrue('ATC' in self.cdf)
        self.assertFalse('notthere' in self.cdf)

    def testgetdimsizescol(self):
        """Get size of dimensions in zVar, column-major"""
        expected = {'ATC': [], 'PhysRecNo': [], 'SpinNumbers': [18],
                    'SectorNumbers': [32], 'RateScalerNames': [16],
                    'SectorRateScalerNames': [9],
                    'SectorRateScalersCounts': [18, 32, 9],
                    'SectorRateScalersCountsSigma': [18, 32, 9],
                    'SpinRateScalersCounts': [18, 16],
                    'SpinRateScalersCountsSigma': [18, 16]}
        for i in expected:
            self.assertEqual(self.cdf[i]._dim_sizes(), expected[i])


class ChangeCDFBase(CDFTests):
    """Base for tests that modify an existing CDF"""
    def setUp(self):
        super(ChangeCDFBase, self).setUp()
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)
        self.cdf.readonly(False)

    def tearDown(self):
        warnings.filterwarnings(
            'ignore', message='^DID_NOT_COMPRESS.*$', category=cdf.CDFWarning,
            module='^spacepy.pycdf')
        self.cdf.close()
        del self.cdf
        del warnings.filters[0]
        os.remove(self.testfile)
        super(ChangeCDFBase, self).tearDown()

    
class ChangeCDF(ChangeCDFBase):
    """Tests that modify an existing CDF, not otherwise specified"""
    def testDeletezVar(self):
        """Delete a zVar"""
        self.cdf['PhysRecNo']._delete()
        self.assertRaises(KeyError, self.cdf.__getitem__, 'PhysRecNo')
        del self.cdf['ATC']
        self.assertFalse('ATC' in self.cdf)

    def testCreateScalarRV(self):
        """Create an RV with scalar data, check error message"""
        msg = 'Record-varying data cannot be scalar. ' \
              'Specify NRV with CDF.new() or put data in array.'
        try:
            self.cdf['value'] = 10
        except ValueError:
            actual = str(sys.exc_info()[1])
            self.assertEqual(msg, actual)
        else:
            self.fail('Should have raised ValueError ' + msg)

    def testSaveCDF(self):
        """Save the CDF and make sure it's different"""
        self.cdf['PhysRecNo']._delete()
        self.cdf.save()
        self.assertNotEqual(self.calcDigest(self.testfile), self.expected_digest)
        self.assertTrue(self.cdf._handle)
        self.cdf['ATC']

    def testReadonlySettable(self):
        """Readonly mode should prevent changes"""
        self.cdf.readonly(True)
        self.assertTrue(self.cdf.readonly())
        message = 'READ_ONLY_MODE: CDF is in read-only mode.'
        try:
            self.cdf['PhysRecNo']._delete()
        except cdf.CDFError:
            (type, val, traceback) = sys.exc_info()
            self.assertEqual(str(val), message)
        else:
            self.fail('Should have raised CDFError: '+ message)

    def testReadonlyDisable(self):
        """Turn off readonly and try to change"""
        self.cdf.readonly(True)
        self.assertTrue(self.cdf.readonly())
        self.cdf.readonly(False)
        try:
            self.cdf['PhysRecNo']._delete()
        except:
            (type, val, traceback) = sys.exc_info()
            self.fail('Raised exception ' + str(val))

    def testRenameVar(self):
        """Rename a variable"""
        zvar = self.cdf['PhysRecNo']
        zvardata = zvar[...]
        zvar.rename('foobar')
        numpy.testing.assert_array_equal(
            zvardata, self.cdf['foobar'][...])
        numpy.testing.assert_array_equal(
            zvardata, zvar[...])
        try:
            zvar = self.cdf['PhysRecNo']
        except KeyError:
            pass
        else:
            self.fail('Should have raised KeyError')

        try:
            zvar.rename('a' * 300)
        except cdf.CDFError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(v.status, cdf.const.BAD_VAR_NAME)
        else:
            self.fail('Should have raised CDFError')

    #This is renamed via a *different* reference, so the original
    #reference still has the old name, and referenced by name.
    @unittest.expectedFailure
    def testRenameSameVar(self):
        """Rename a variable while keeping a reference to it"""
        zvar = self.cdf['PhysRecNo']
        self.cdf['PhysRecNo'].rename('foobar')
        numpy.testing.assert_array_equal(
            zvar[...], self.cdf['foobar'][...])

    #A variable is deleted but we still have a reference to it
    #I'm not sure what would be expected behavior here, but
    #there's a NO_SUCH_VAR error raised.
    @unittest.expectedFailure
    def testDeleteAndAccess(self):
        """Delete a variable and then try to access it again"""
        zvar = self.cdf['PhysRecNo']
        number = zvar._num()
        del self.cdf['PhysRecNo']
        foo = zvar[...]

    #A completely different variable is made with the same name,
    #so the old reference will point to the new variable
    #(referenced by name)
    @unittest.expectedFailure
    def testDeleteSameName(self):
        """Delete a variable and re-make with same name"""
        zvar = self.cdf['PhysRecNo']
        number = zvar._num()
        del self.cdf['PhysRecNo']
        self.cdf['PhysRecNo'] = [1, 2, 3, 4]
        #Verify we have different variable numbers between two things
        #with the same name
        self.assertNotEqual(number, self.cdf['PhysRecNo']._num())
        #And since they're referring to different variables, the
        #contents should be different.
        #https://stackoverflow.com/questions/38506044/numpy-testing-assert-array-not-equal
        self.assertRaises(AssertionError, numpy.testing.assert_array_equal,
                          zvar[...], self.cdf['PhysRecNo'][...])

    def testDeleteWithRef(self):
        """Delete a variable while have a reference to another one"""
        zvar = self.cdf['PhysRecNo']
        zvardata = zvar[...]
        del self.cdf['ATC']
        numpy.testing.assert_array_equal(
            zvardata, zvar[...])

    def testDeleteAttrWithRef(self):
        """Delete a gAttr while having reference to another one"""
        #This is #19
        ack = self.cdf.attrs['Acknowledgement']
        ackdata = ack[...]
        #This is #3
        del self.cdf.attrs['Data_type']
        numpy.testing.assert_array_equal(ackdata, ack[...])
        
    def testNewVar(self):
        """Create a new variable"""
        self.cdf.new('newzVar', [[1, 2, 3], [4, 5, 6]],
                     const.CDF_INT4)
        self.assertTrue('newzVar' in self.cdf)
        zvar = self.cdf['newzVar']
        numpy.testing.assert_array_equal(
            [[1, 2, 3], [4, 5, 6]], zvar[...])
        self.assertEqual(2, len(zvar))
        self.assertEqual([3], zvar._dim_sizes())

    def testNewVarCompress(self):
        """Create a new, compressed variable"""
        data = numpy.arange(1000) #big enough to actually compress
        self.cdf.new('newzVar', data,
                     const.CDF_INT4,
                     compress=const.GZIP_COMPRESSION, compress_param=9)
        self.assertTrue('newzVar' in self.cdf)
        zvar = self.cdf['newzVar']
        numpy.testing.assert_array_equal(
            data, zvar[...])
        self.assertEqual(1000, len(zvar))
        (comptype, compparam) = zvar.compress()
        self.assertEqual(const.GZIP_COMPRESSION, comptype)
        self.assertEqual(9, compparam)

    def testNewVarFromdmarray(self):
        """Create a new variable with data in dmarray"""
        indata = datamodel.dmarray([1,2,3], dtype=numpy.int8,
                                   attrs={'name': 'bob'})
        self.cdf.new('newzVar', indata)
        numpy.testing.assert_array_equal(
            indata[...], self.cdf['newzVar'][...])
        self.assertEqual('bob', self.cdf['newzVar'].attrs['name'])
        self.assertEqual(numpy.int8, self.cdf['newzVar'].dtype)

    def testNewVarFromVarCompress(self):
        """Create a new variable from a variable, check compression"""
        self.cdf.new('newzvar1', compress=const.GZIP_COMPRESSION,
                     compress_param=8, data=numpy.arange(1000))
        zvar = self.cdf.new('newzvar2', data=self.cdf['newzvar1'])
        comptype, compparam = zvar.compress()
        self.assertEqual(const.GZIP_COMPRESSION.value,
                         comptype.value)
        self.assertEqual(8, compparam)

    def testNewVarFromVarCompress(self):
        """Create a new variable from a variable, change compression"""
        zvar = self.cdf.new('newzvar1', compress=const.GZIP_COMPRESSION,
                            data=self.cdf['SpinRateScalersCounts'])
        comptype, compparam = zvar.compress()
        self.assertEqual(const.GZIP_COMPRESSION.value,
                         comptype.value)
        self.assertEqual(5, compparam)

    def testNewVarFromdmarrayAssign(self):
        """Create a new variable by assigning from dmarray"""
        indata = datamodel.dmarray([1,2,3], dtype=numpy.int8,
                                   attrs={'name': 'bob'})
        self.cdf['newzVar'] = indata
        numpy.testing.assert_array_equal(
            indata[...], self.cdf['newzVar'][...])
        self.assertEqual('bob', self.cdf['newzVar'].attrs['name'])
        self.assertEqual(numpy.int8, self.cdf['newzVar'].dtype)

    def testNewVarAssign(self):
        """Create a new variable by assigning to CDF element"""
        self.cdf['newzVar'] = [[1, 2, 3], [4, 5, 6]]
        self.assertTrue('newzVar' in self.cdf)
        zvar = self.cdf['newzVar']
        numpy.testing.assert_array_equal(
            [[1, 2, 3], [4, 5, 6]], zvar[...])
        self.assertEqual(2, len(zvar))
        self.assertEqual([3], zvar._dim_sizes())

    def testNewVarApostrophe(self):
        """Assign to new variable with apostrophe in name"""
        self.cdf['new\'zVar'] = [[1, 2, 3], [4, 5, 6]]
        self.assertTrue('new\'zVar' in self.cdf)
        zvar = self.cdf['new\'zVar']
        numpy.testing.assert_array_equal(
            [[1, 2, 3], [4, 5, 6]], zvar[...])
        self.assertEqual(2, len(zvar))
        self.assertEqual([3], zvar._dim_sizes())

    def testNewVarUnicode(self):
        """Create a new variable with a unicode name"""
        if str is bytes:
            self.cdf['newzVar'.decode()] = [[1, 2, 3], [4, 5, 6]]
        else:
            self.cdf['newzVar'] = [[1, 2, 3], [4, 5, 6]]
        numpy.testing.assert_array_equal(
            [[1, 2, 3], [4, 5, 6]], self.cdf['newzVar'][...])

    def testNewVarTime(self):
        with spacepy_testing.assertWarns(
                self, 'always',
                r'No type specified for time input\; assuming'
                r' CDF_TIME_TT2000\.$',
                DeprecationWarning, r'spacepy\.pycdf$'):
            self.cdf['newzVar'] = [datetime.datetime(2010, 1, 1)]
        # Most of the type-guessing testing is in NoCDF, but this is here
        # because the warning of the default changing is associated with
        # creating a zVar.
        expected = cdf.const.CDF_TIME_TT2000.value if cdf.lib.supports_int8 \
                   else cdf.const.CDF_EPOCH.value
        self.assertEqual(expected, self.cdf['newzVar'].type())

    def testBadDataSize(self):
        """Attempt to assign data of the wrong size to a zVar"""
        try:
            self.cdf['MeanCharge'] = [1.0, 2.0, 3.0]
        except ValueError:
            pass
        else:
            self.fail('Should have raised ValueError')

    def testChangeVarType(self):
        """Change the type of a variable"""
        self.cdf['new'] = [-1, -2, -3]
        self.cdf['new'].type(const.CDF_UINT1)
        numpy.testing.assert_array_equal(
            [255, 254, 253], self.cdf['new'][...])

    def testNewVarNoData(self):
        """Create a new variable without providing any data"""
        self.assertRaises(ValueError, self.cdf.new, 'newvar')
        self.cdf.new('newvar', None, const.CDF_INT4)
        self.assertEqual([], self.cdf['newvar']._dim_sizes())

        self.cdf.new('newvar2', None, const.CDF_CHAR, dims=[])
        self.assertEqual(1, self.cdf['newvar2'].nelems())

    def testNewVarEmptyArray(self):
        """Create a variable with an empty numpy array"""
        self.cdf['newvar1'] = numpy.array([], dtype=numpy.int32)
        self.assertEqual(const.CDF_INT4.value, self.cdf['newvar1'].type())
        #This should fail, can't guess a CDF type with no data
        try:
            self.cdf['newvar2'] = numpy.array([], dtype=object)
        except ValueError:
            (t, val, traceback) = sys.exc_info()
            self.assertEqual('Cannot determine CDF type of empty object array.',
                             str(val))
        else:
            self.fail('Should have raised ValueError')
        self.assertFalse('newvar2' in self.cdf)

    def testNewVarDatetimeArray(self):
        """Create a variable with a datetime numpy array"""
        warnings.filterwarnings(
            'ignore', r'^No type specified for time input.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            self.cdf['newvar'] = numpy.array([datetime.datetime(2010, 1, 1)])
        finally:
            del warnings.filters[0]
        self.assertEqual((const.CDF_TIME_TT2000 if cdf.lib.supports_int8
                          else const.CDF_EPOCH).value,
                         self.cdf['newvar'].type())

    def testNewVarNRV(self):
        """Create a new non-record-varying variable"""
        self.cdf.new('newvar2', [1, 2, 3], recVary=False)
        self.assertFalse(self.cdf['newvar2'].rv())
        self.assertEqual([3], self.cdf['newvar2']._dim_sizes())
        numpy.testing.assert_array_equal(
            [1, 2, 3], self.cdf['newvar2'][...])

    def testChangeRV(self):
        """Change record variance"""
        zVar = self.cdf.new('newvar', dims=[], type=const.CDF_INT4)
        self.assertTrue(zVar.rv())
        zVar.rv(False)
        self.assertFalse(zVar.rv())
        zVar.rv(True)
        self.assertTrue(zVar.rv())

    def testChangeSparseRecordsPrev(self):
        """Change sparse records mode to PREV"""
        zVar = self.cdf.new('newvarSR', dims=[], type=const.CDF_INT4)
        self.assertEqual(zVar.sparse(), const.NO_SPARSERECORDS)
        zVar.sparse(const.PREV_SPARSERECORDS)
        self.assertEqual(zVar.sparse(), const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[1], 1);
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[2], 1);
    
    def testSparseMultidimPrev(self):
        """test multi-dimensional sparse record variables, PREV"""
        zVar = self.cdf.new('sr', dims=[2], type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = [1, 2];
        zVar[3] = [3, 4];
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            numpy.testing.assert_array_equal(
                [[1, 2], [1, 2], [1, 2], [3, 4]],
                zVar[...])

    def testSparseRecordsReadAll(self):
        """Read all records from a sparse variable"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        self.assertEqual(4, len(zVar))
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            numpy.testing.assert_array_equal(zVar[...],
                                             [1, 1, 1, 2]);

    def testSparseReadOffEnd(self):
        """Read past the end of a sparse variable"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[...] = [1, 2, 3]
        hs = cdf._Hyperslice(zVar, slice(None, 4, None))
        # Make sure the slice sizing is proper before read
        self.assertEqual([0], hs.starts)
        self.assertEqual([4], hs.counts)
        with spacepy_testing.assertWarns(
                self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                r'spacepy\.pycdf$'):
            numpy.testing.assert_array_equal(zVar[:4],
                                             [1, 2, 3, 3]);

    def testChangeSparseRecordsPad(self):
        """Change sparse records mode to PAD"""
        zVar = self.cdf.new('newvarSRPad', dims=[], type=const.CDF_INT4)
        zVar.sparse(const.PAD_SPARSERECORDS)
        self.assertEqual(zVar.sparse(), const.PAD_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        pad = zVar.pad()
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[1], pad);
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[2], pad);
        pad = zVar.pad(123)
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[1], pad);
        with spacepy_testing.assertWarns(*aw_args):
            self.assertEqual(zVar[2], pad);

    def testSparseMultidimPad(self):
        """test multi-dimensional sparse record variables, PAD"""
        zVar = self.cdf.new('sr', dims=[2], type=const.CDF_INT4)
        zVar.sparse(const.PAD_SPARSERECORDS)
        zVar[0] = [1, 2];
        zVar[3] = [3, 4];
        pad = zVar.pad(20)
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            numpy.testing.assert_array_equal(
                [[1, 2], [20, 20], [20, 20], [3, 4]],
                zVar[...])

    def testSparsePrevInsert(self):
        """Insert records in sparse records previous mode"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        with self.assertRaises(NotImplementedError) as cm:
            zVar.insert(2, 99)
        self.assertEqual('Sparse records do not support insertion.',
                         str(cm.exception))
        # Following is test for if this did work
        # Arguments to use for assertWarngs, since used a lot....
        #aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
        #           r'spacepy\.pycdf$')
        #with spacepy_testing.assertWarns(*aw_args):
        #    numpy.testing.assert_array_equal(zVar[0:5],
        #                                     [1, 1, 99, 99, 2])

    def testSparsePrevSingle(self):
        """Delete one record in sparse records previous mode"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        zVar[5] = 3;
        del zVar[3]
        # Arguments to use for assertWarngs, since used a lot....
        aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
                   r'spacepy\.pycdf$')
        with spacepy_testing.assertWarns(*aw_args):
            numpy.testing.assert_array_equal(zVar[0:6],
                                             [1, 1, 1, 1, 1, 3])

    def testSparsePrevDelete(self):
        """Delete records in sparse records previous mode"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        zVar[5] = 3;
        with self.assertRaises(NotImplementedError) as cm:
            del zVar[3:5]
        self.assertEqual('Sparse records do not support multi-record delete.',
                         str(cm.exception))
        # Following is test for if this did work
        # Arguments to use for assertWarngs, since used a lot....
        #aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
        #           r'spacepy\.pycdf$')
        #with spacepy_testing.assertWarns(*aw_args):
        #    numpy.testing.assert_array_equal(zVar[0:4],
        #                                     [1, 1, 1, 3])

    def testSparsePrevTruncate(self):
        """Truncate records in sparse records previous mode"""
        zVar = self.cdf.new('newvarSR', type=const.CDF_INT4)
        zVar.sparse(const.PREV_SPARSERECORDS)
        zVar[0] = 1;
        zVar[3] = 2;
        zVar[4] = 3;
        with self.assertRaises(NotImplementedError) as cm:
            zVar[0:] = [100, 101, 102]
        self.assertEqual('Sparse records do not support truncation on write.',
                         str(cm.exception))
        # Following is test for if this did work
        # Arguments to use for assertWarngs, since used a lot....
        #aw_args = (self, 'always', r'VIRTUAL_RECORD_DATA', cdf.CDFWarning,
        #           r'spacepy\.pycdf$')
        #with spacepy_testing.assertWarns(*aw_args):
        #    numpy.testing.assert_array_equal(zVar[0:4],
        #                                     [1000, 101, 102, 3])

    def testSparseOnCreate(self):
        """Specify sparseness when creating a variable"""
        zVar = self.cdf.new('newvar', data=[1, 2, 3, 4],
                            sparse=const.PAD_SPARSERECORDS, pad=99)
        self.assertEqual(const.PAD_SPARSERECORDS, zVar.sparse())
        self.assertEqual(99, zVar.pad())
        with spacepy_testing.assertWarns(self, 'always', r'VIRTUAL_RECORD_DATA',
                                         cdf.CDFWarning, r'spacepy\.pycdf$'):
            numpy.testing.assert_array_equal(zVar[0:5], [1, 2, 3, 4, 99])

    def testSparseCopy(self):
        """Make sure sparseness carries through to VarCopy"""
        zVar = self.cdf.new('newvar', data=[1, 2, 3, 4],
                            sparse=const.PAD_SPARSERECORDS, pad=99)
        cp = zVar.copy()
        self.assertEqual(const.PAD_SPARSERECORDS, cp.sparse())
        self.assertEqual(99, cp.pad())
        numpy.testing.assert_array_equal(zVar[...], [1, 2, 3, 4])

    def testNRVWritePad(self):
        """Write pad value for NRV"""
        v = self.cdf['SectorNumbers']
        v.pad('Q')
        self.assertEqual('Q', v.pad())

    def testNRVSparse(self):
        """Make NRV sparse variable"""
        v = self.cdf.new('newvar', recVary=False, type=const.CDF_INT1)
        with self.assertRaises(cdf.CDFError) as cm:
            v.sparse(const.PAD_SPARSERECORDS)
        self.assertEqual("CANNOT_SPARSERECORDS: Sparse records can't be"
                         " set/modified for the variable.", str(cm.exception))
        self.assertEqual(const.CANNOT_SPARSERECORDS, cm.exception.status)

    def testSparseNonConst(self):
        """Set sparseness without using the const module"""
        v = self.cdf.new('newvar', type=const.CDF_INT1)
        v.sparse(ctypes.c_long(0))
        self.assertEqual(0, v.sparse().value)
        v.sparse(0)
        self.assertEqual(0, v.sparse().value)

    def testChecksum(self):
        """Change checksumming on the CDF"""
        self.cdf.checksum(True)
        self.assertTrue(self.cdf.checksum())
        self.cdf.checksum(False)
        self.assertFalse(self.cdf.checksum())

    def testCompress(self):
        """Change compression on the CDF"""
        self.cdf.compress(const.GZIP_COMPRESSION)
        (comptype, parm) = self.cdf.compress()
        self.assertEqual(const.GZIP_COMPRESSION, comptype),
        self.assertEqual(5, parm)
        self.cdf.compress(const.NO_COMPRESSION)
        (comptype, parm) = self.cdf.compress()
        self.assertEqual(const.NO_COMPRESSION, comptype),
        self.assertEqual(0, parm)

    def testVarCompress(self):
        """Change compression on a variable"""
        zvar = self.cdf.new('newvar', type=const.CDF_INT4, dims=[])
        zvar.compress(const.GZIP_COMPRESSION)
        (comptype, parm) = zvar.compress()
        self.assertEqual(const.GZIP_COMPRESSION, comptype),
        self.assertEqual(5, parm)
        zvar.compress(const.NO_COMPRESSION)
        (comptype, parm) = zvar.compress()
        self.assertEqual(const.NO_COMPRESSION, comptype),
        self.assertEqual(0, parm)

    def testWarnings(self):
        """Bizarre way to force a warning"""
        attrnum = ctypes.c_long(0)
        msg = 'this is a very long string intended to get up to ' \
        '257 characters or so because the maximum length ' \
        'of an attribute name is 256 characters and ' \
        'attribute name truncated is just about the ONLY ' \
        'warning I can figure out how to raise in the CDF ' \
        'library and this is really a serious pain in just ' \
        'about every portion of the anatomy.'
        if not str is bytes:
            msg = msg.encode('ascii')
        with spacepy_testing.assertWarns(self, 'always', r'ATTR_NAME_TRUNC',
                                         cdf.CDFWarning, r'spacepy\.pycdf$'):
            self.cdf._call(cdf.const.CREATE_, cdf.const.ATTR_, msg,
                           cdf.const.GLOBAL_SCOPE, ctypes.byref(attrnum))

    def testAssignEmptyList(self):
        """Assign an empty list to a variable"""
        self.cdf['ATC'] = []
        self.assertEqual(0, len(self.cdf['ATC']))

    def testReadEmptyList(self):
        """Read from an empty variable"""
        self.cdf['ATC'] = []
        data = self.cdf['ATC'][...]
        self.assertEqual((0,), data.shape)
        self.assertEqual(object, data.dtype)

    def testCopyVariable(self):
        """Copy one variable to another"""
        varlist = list(self.cdf.keys())
        for name in varlist:
            oldvar = self.cdf[name]
            self.cdf[name + '_2'] = oldvar
            newvar = self.cdf[name + '_2']
            msg = 'Variable ' + name + ' failed.'
            self.assertEqual(oldvar._n_dims(), newvar._n_dims(), msg)
            self.assertEqual(oldvar._dim_sizes(), newvar._dim_sizes(), msg)
            self.assertEqual(oldvar.type(), newvar.type(), msg)
            self.assertEqual(oldvar.nelems(), newvar.nelems(), msg)
            self.assertEqual(oldvar.compress(), newvar.compress(), msg)
            self.assertEqual(oldvar.rv(), newvar.rv(), msg)
            self.assertEqual(oldvar.dv(), newvar.dv(), msg)
            numpy.testing.assert_array_equal(
                oldvar[...], newvar[...], msg)
            oldlist = oldvar.attrs
            newlist = newvar.attrs
            for attrname in oldlist:
                self.assertTrue(attrname in newlist)
                self.assertEqual(oldlist[attrname], newlist[attrname])
                self.assertEqual(oldlist.type(attrname),
                                 newlist.type(attrname))

    def testCloneVariable(self):
        """Clone a variable's type, dims, etc. to another"""
        varlist = list(self.cdf.keys())
        for name in varlist:
            oldvar = self.cdf[name]
            self.cdf.clone(oldvar, name + '_2', False)
            newvar = self.cdf[name + '_2']
            msg = 'Variable ' + name + ' failed.'
            self.assertEqual(oldvar._n_dims(), newvar._n_dims(), msg)
            self.assertEqual(oldvar._dim_sizes(), newvar._dim_sizes(), msg)
            self.assertEqual(oldvar.type(), newvar.type(), msg)
            self.assertEqual(oldvar.nelems(), newvar.nelems(), msg)
            self.assertEqual(oldvar.compress(), newvar.compress(), msg)
            self.assertEqual(oldvar.rv(), newvar.rv(), msg)
            self.assertEqual(oldvar.dv(), newvar.dv(), msg)
            if newvar.rv():
                self.assertEqual(0, newvar[...].size, msg)
            oldlist = oldvar.attrs
            newlist = newvar.attrs
            for attrname in oldlist:
                self.assertTrue(
                    attrname in newlist,
                    'Attribute {0} not found in copy of {1}'.format(
                    attrname, name))
                self.assertEqual(oldlist[attrname], newlist[attrname])
                self.assertEqual(oldlist.type(attrname),
                                 newlist.type(attrname))

    def testCloneVarCompressed(self):
        """Clone a variable and keep its compression"""
        self.cdf.new('IAmCompressed', data=numpy.array([1, 2, 3, 4, 5]),
                     compress=const.GZIP_COMPRESSION, compress_param=7)
        self.cdf.clone(self.cdf['IAmCompressed'], 'IAmCompressed_2')
        output = self.cdf['IAmCompressed_2'].compress()
        self.assertEqual((const.GZIP_COMPRESSION, 7), output)

    def testDimVariance(self):
        """Check and change dimension variance of a variable"""
        self.assertEqual([True],
            self.cdf['SpinNumbers'].dv())
        self.assertEqual([True, True, True],
            self.cdf['SectorRateScalersCounts'].dv())
        self.cdf.new('foobar', type=const.CDF_INT1,
                     dims=[2, 3], dimVarys=[True, False])
        self.assertEqual([True, False],
                         self.cdf['foobar'].dv())
        self.cdf['foobar'].dv([False, True])
        self.assertEqual([False, True],
                         self.cdf['foobar'].dv())

    def testAssignEpoch16Entry(self):
        """Assign to an Epoch16 entry"""
        warnings.filterwarnings(
            'ignore', r'^Assuming CDF_.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            self.cdf['ATC'].attrs['FILLVAL'] = datetime.datetime(2010,1,1)
        finally:
            del warnings.filters[0]
        self.assertEqual(datetime.datetime(2010,1,1),
                         self.cdf['ATC'].attrs['FILLVAL'])

    def testVarTrailingSpaces(self):
        """Cut trailing spaces from names of vars"""
        self.cdf['foobar  '] = [1, 2, 3]
        namelist = list(self.cdf.keys())
        self.assertTrue('foobar' in namelist)
        self.assertFalse('foobar  ' in namelist)

    def testAttrTrailingSpaces(self):
        """Cut trailing spaces from names of attributes"""
        self.cdf.attrs['hi '] = 'hello'
        namelist = list(self.cdf.attrs.keys())
        self.assertTrue('hi' in namelist)
        self.assertFalse('hi ' in namelist)

    def testInt8(self):
        """Create a new INT8 zVar"""
        self.cdf['foobar'] = numpy.array([1, 2, 3], dtype=numpy.int64)
        if cdf.lib.supports_int8:
            self.assertEqual(self.cdf['foobar'].type(), const.CDF_INT8.value)
        else:
            self.assertEqual(self.cdf['foobar'].type(), const.CDF_BYTE.value)

    def testTT2000New(self):
        """Create a new TT2000 zVar"""
        expected = [datetime.datetime(2010, 1, 1) + 
                    datetime.timedelta(days=i) for i in range(5)]
        if cdf.lib.supports_int8:
            self.cdf.new('foobar', data=expected, type=const.CDF_TIME_TT2000)
            self.assertEqual(self.cdf['foobar'].type(),
                             const.CDF_TIME_TT2000.value)
            numpy.testing.assert_array_equal(expected,
                                             self.cdf['foobar'][...])
        else:
            message = 'INT8 and TIME_TT2000 require CDF library 3.4.0'
            try:
                self.cdf.new('foobar', data=expected,
                             type=const.CDF_TIME_TT2000)
            except ValueError:
                self.assertEqual(message, str(sys.exc_info()[1]))
            else:
                self.fail('Should have raised ValueError: ' + message)

    def testUnicodeString(self):
        """Write Unicode to a string variable"""
        if str != bytes: #py3k:
            data = [ 'hi', 'there']
        else:
            data = ['hi'.decode(), 'there'.decode()]
        self.cdf['teststr'] = data
        expected = data
        out = self.cdf['teststr'][0:2]
        numpy.testing.assert_array_equal(expected, numpy.char.rstrip(out))

    def testFloatEpoch(self):
        """Write floats to an Epoch variable"""
        self.cdf.new('epochtest', type=const.CDF_EPOCH)
        data = numpy.array([62987673600000.0,
                            62987760000000.0], dtype=numpy.float64)
        self.cdf['epochtest'][:] = data
        numpy.testing.assert_array_equal(
            numpy.array([datetime.datetime(1996, 1, 1),
                         datetime.datetime(1996, 1, 2)]),
            self.cdf['epochtest'][:])

    def testEpochForceRaw(self):
        """Try to write datetime to a forced-raw Epoch"""
        self.cdf.new('epochtest', type=const.CDF_EPOCH)
        try:
            self.cdf.raw_var('epochtest')[:] = [
                datetime.datetime(1996, 1, 1),
                datetime.datetime(1996, 1, 2)]
        except (TypeError, ValueError):
            pass
        else:
            self.fail('Should have raised TypeError or ValueError')

    def testFloatEpoch16(self):
        """Write floats to an Epoch16 variable"""
        self.cdf.new('epochtest', type=const.CDF_EPOCH16)
        data = numpy.array([[62987673600.0, 0.0],
                            [62987760000.0, 0.0]], dtype=numpy.float64)
        self.cdf['epochtest'][:] = data
        numpy.testing.assert_array_equal(
            numpy.array([datetime.datetime(1996, 1, 1),
                         datetime.datetime(1996, 1, 2)]),
            self.cdf['epochtest'][:])

    def testAttrsRawEpoch16(self):
        """Assign float attribute to Epoch16"""
        self.cdf.new('epochtest', type=const.CDF_EPOCH16)
        data = numpy.array([[62987673600.0, 0.0],
                            [62987760000.0, 0.0]], dtype=numpy.float64)
        self.cdf['epochtest'][:] = data
        self.cdf.raw_var('epochtest').attrs.new(
            'FILLVAL', numpy.array([-1e31, -1e31]), const.CDF_EPOCH16)
        numpy.testing.assert_array_equal(
            numpy.array([-1e31, -1e31], dtype=numpy.float64),
            self.cdf.raw_var('epochtest').attrs['FILLVAL'])

    def testInt8TT2000(self):
        """Write integers to a TT2000 variable"""
        if not cdf.lib.supports_int8:
            return
        self.cdf.new('epochtest', type=const.CDF_TIME_TT2000)
        data = numpy.array([-126273537816000000,
                            -126187137816000000], dtype=numpy.int64)
        self.cdf['epochtest'][:] = data
        numpy.testing.assert_array_equal(
            numpy.array([datetime.datetime(1996, 1, 1),
                         datetime.datetime(1996, 1, 2)]),
            self.cdf['epochtest'][:])

    def testSetAttrListTypes(self):
        """Make sure that assigning an entire attrlist gets types right"""
        self.cdf.new('foo2', data=[0, 1], type=const.CDF_DOUBLE)
        a = {'FILLVAL': -1e31}
        self.cdf['foo2'].attrs = a
        self.assertEqual(const.CDF_DOUBLE.value,
            self.cdf['foo2'].attrs.type('FILLVAL'))

    def testCreateVarFromStringArray(self):
        """make a zvar from a numpy string array and get the size right"""
        inarray = numpy.array(['hi', 'there'], dtype='|S6')
        self.cdf['string6'] = inarray
        self.assertEqual(6, self.cdf['string6'].nelems())
        expected = numpy.require(inarray, dtype=str)
        outarray = numpy.char.rstrip(self.cdf['string6'][...])
        numpy.testing.assert_array_equal(expected, outarray)

    def testCreateVarFromUnicodeArray(self):
        """make a zvar from numpy string array in unicode"""
        if str is bytes: #Py2k, don't expect unicode handling of char
            return
        inarray = numpy.array(['hi', 'there'], dtype='U6')
        self.cdf['string62'] = inarray
        self.assertEqual(6, self.cdf['string62'].nelems())
        out = self.cdf['string62'][...]
        numpy.testing.assert_array_equal(inarray, numpy.char.rstrip(out))

    def testAppendgEntry(self):
        """Append to a gAttr"""
        self.cdf.attrs['foo'] = ['foo', 'bar', 'baz']
        self.cdf.attrs['foo'].append('qux')
        self.assertEqual(['foo', 'bar', 'baz', 'qux'],
                         self.cdf.attrs['foo'][:])
        self.cdf.attrs['foo'].new(5, type=const.CDF_INT2, number=5)
        self.cdf.attrs['foo'].append('quxx')
        self.assertEqual('quxx', self.cdf.attrs['foo'][6])

    def testInsertgEntry(self):
        """Insert into a gAttr"""
        self.cdf.attrs.new('foo')
        a = self.cdf.attrs['foo']
        a.new(0, const.CDF_INT2, number=0)
        a.new(1, const.CDF_INT4, number=1)
        a.new(2, const.CDF_UINT2, number=2)
        a.new(5, const.CDF_UINT4, number=5)
        a.insert(2, 10)
        self.assertEqual(const.CDF_INT2.value, a.type(0))
        self.assertEqual(0, a[0])
        self.assertEqual(const.CDF_INT4.value, a.type(1))
        self.assertEqual(1, a[1])
        self.assertEqual(10, a[2])
        self.assertEqual(const.CDF_UINT2.value, a.type(2))
        self.assertEqual(2, a[3])
        self.assertEqual(const.CDF_UINT4.value, a.type(6))
        self.assertEqual(5, a[6])
        self.assertFalse(a.has_entry(5))

    def testIterateEmpty(self):
        """Make an empty variable and iterate over it"""
        self.cdf.new('Empty', type=const.CDF_INT1)
        self.assertEqual([], list(self.cdf['Empty']))
        self.assertRaises(StopIteration,
                          next,
                          iter(self.cdf['Empty']))

    def testReadAllEmpty(self):
        """Make an empty variable and read all data"""
        self.cdf.new('Empty', type=const.CDF_INT1)
        self.assertEqual((0,), self.cdf['Empty'][...].shape)

    def testReadSomeEmpty(self):
        """Make an empty variable and read off end"""
        self.cdf.new('Empty', type=const.CDF_INT1)
        self.assertEqual((0,), self.cdf['Empty'][0:2].shape)

    def testLengthEmpty(self):
        """Make an empty variable and check its length"""
        self.cdf.new('Empty', type=const.CDF_INT1)
        self.assertEqual(0, len(self.cdf['Empty']))

    def testIndexEmpty(self):
        """Make an empty variable and index it"""
        self.cdf.new('Empty', type=const.CDF_INT1)
        try:
            self.cdf['Empty'][1]
        except IndexError:
            pass
        else:
            self.fail('Should have raised IndexError on 1')
        try:
            self.cdf['Empty'][0]
        except IndexError:
            pass
        else:
            self.fail('Should have raised IndexError on 0')

    def testIndexOneValue(self):
        """Put a single value into a var and index past it"""
        #This is basically to test that indexing empty at [0] is special case
        self.cdf.new('Single', data=[1], type=const.CDF_INT1)
        try:
            self.cdf['Single'][1]
        except IndexError:
            pass
        else:
            self.fail('Should have raised IndexError on 1')


class ChangezVar(ChangeCDFBase):
    """Tests that modify a zVar"""

    def testWriteSubscripted(self):
        """Write data to a slice of a zVar"""
        expected = ['0', '1', '99', '3', '98', '5', '97', '7',
                    '8', '9']
        self.cdf['SpinNumbers'][2:7:2] = ['99', '98', '97']
        numpy.testing.assert_array_equal(
            expected, numpy.char.rstrip(self.cdf['SpinNumbers'][0:10]))

        expected = self.cdf['SectorRateScalersCounts'][...]
        expected[4][5][5][8:3:-1] = [101.0, 102.0, 103.0, 104.0, 105.0]
        self.cdf['SectorRateScalersCounts'][4, 5, 5, 8:3:-1] = \
            [101.0, 102.0, 103.0, 104.0, 105.0]
        numpy.testing.assert_array_equal(
            expected, self.cdf['SectorRateScalersCounts'][...])

        self.cdf['PhysRecNo'] = [1, 2, 3]
        numpy.testing.assert_array_equal(
            [1, 2, 3], self.cdf['PhysRecNo'][...])

    def testWriteOffEnd(self):
        """Extend variable by writing off the end"""
        additional = [2000 + i for i in range(20)]
        expected = self.cdf['PhysRecNo'][0:95].tolist() + additional
        self.cdf['PhysRecNo'][95:] = additional
        self.assertEqual(115, len(self.cdf['PhysRecNo']))
        numpy.testing.assert_array_equal(expected, self.cdf['PhysRecNo'][:])

    def testWriteExtend(self):
        """Extend variable with explicit extend call"""
        additional = [2000 + i for i in range(20)]
        oldlen = len(self.cdf['PhysRecNo'])
        expected = self.cdf['PhysRecNo'][...].tolist() + additional
        self.cdf['PhysRecNo'].extend(additional)
        self.assertEqual(oldlen + 20, len(self.cdf['PhysRecNo']))
        numpy.testing.assert_array_equal(expected, self.cdf['PhysRecNo'][:])

    def testInsertRecord(self):
        """Insert a record into the middle of a variable"""
        PhysRecNoData = self.cdf['PhysRecNo'][...].tolist()
        PhysRecNoData[5:6] = [-1, -2, -3, -4]
        self.cdf['PhysRecNo'][5:6] = [-1, -2, -3, -4]
        self.assertEqual(103, len(self.cdf['PhysRecNo']))
        numpy.testing.assert_array_equal(
            PhysRecNoData, self.cdf['PhysRecNo'][...])

    def testzVarInsert(self):
        """Insert a record into a zVariable"""
        before = self.cdf['ATC'][:].tolist()
        self.cdf['ATC'].insert(100, datetime.datetime(2010, 12, 31))
        before.insert(100, datetime.datetime(2010, 12, 31))
        numpy.testing.assert_array_equal(before, self.cdf['ATC'][:])
        before = self.cdf['MeanCharge'][:].tolist()
        self.cdf['MeanCharge'].insert(20, [99] * 16)
        before.insert(20, [99] * 16)
        numpy.testing.assert_array_equal(before, self.cdf['MeanCharge'][:])

    def testWriteAndTruncate(self):
        """Write with insufficient data to fill all existing records"""
        expected = [-1 * i for i in range(20)]
        self.cdf['PhysRecNo'][:] = expected
        numpy.testing.assert_array_equal(
            expected, self.cdf['PhysRecNo'][:])

    def testWriteWrongSizeData(self):
        """Write with data sized or shaped differently from expected"""
        #This is a bit odd but numpy "shape" is sometimes Long...
        message = ('attempt to assign data of dimensions {} '
                  'to slice of dimensions (5,)').format(
                      str(numpy.array([1, 2, 3]).shape))
        try:
            self.cdf['SpinNumbers'][0:5] = [b'99', b'98', b'97']
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should have raised ValueError: ' + message)

        message = ('attempt to assign data of dimensions {} '
                  'to slice of dimensions (3, 6, 2)').format(
                      str(numpy.empty((2, 6, 2)).shape))

        try:
            self.cdf['SpinRateScalersCounts'][0:3, 12:, 0:4:2] = \
                [[[0, 1], [2, 3], [4, 5], [6, 7], [8, 9], [0, 1]],
                 [[2, 3], [4, 5], [6, 7], [8, 9], [0, 1], [2, 3]],
                 ]
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should have raised ValueError: ' + message)

    def testDeleteRecord(self):
        """Delete records from a variable"""
        oldlen = len(self.cdf['PhysRecNo'])
        PhysRecCopy = self.cdf['PhysRecNo'].copy()
        del self.cdf['PhysRecNo'][5]
        PhysRecCopy = numpy.append(PhysRecCopy[:5], PhysRecCopy[6:], 0)
        self.assertEqual(oldlen - 1, len(self.cdf['PhysRecNo']))
        numpy.testing.assert_array_equal(
            PhysRecCopy[0:15], self.cdf['PhysRecNo'][0:15])

        oldlen = len(self.cdf['ATC'])
        ATCCopy = self.cdf['ATC'].copy()
        del self.cdf['ATC'][0::2]
        self.assertEqual(int(oldlen / 2), len(self.cdf['ATC']))
        numpy.testing.assert_array_equal(ATCCopy[1::2], self.cdf['ATC'][...])

        message = 'Cannot delete records from non-record-varying variable.'
        try:
            del self.cdf['SpinNumbers'][0]
        except:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(t, TypeError)
            self.assertEqual(str(v), message)
        else:
            self.fail('Should have raised TypeError: ' + message)

        oldlen = len(self.cdf['SectorRateScalersCounts'])
        SectorRateScalersCountsCopy = \
                                    self.cdf['SectorRateScalersCounts'].copy()
        SectorRateScalersCountsCopy = numpy.delete(
            SectorRateScalersCountsCopy,
            list(range(*slice(-1,-5,-1).indices(
                        len(SectorRateScalersCountsCopy)))),
            0)
        del self.cdf['SectorRateScalersCounts'][-1:-5:-1]
        self.assertEqual(oldlen - 4, len(self.cdf['SectorRateScalersCounts']))
        self.assertEqual(
            SectorRateScalersCountsCopy[...].shape,
            self.cdf['SectorRateScalersCounts'][...].shape)
        numpy.testing.assert_array_equal(
            SectorRateScalersCountsCopy[...].flatten(),
            self.cdf['SectorRateScalersCounts'][...].flatten())

        oldlen = len(self.cdf['SectorRateScalersCounts'])
        del self.cdf['SectorRateScalersCounts'][-1:-5]
        self.assertEqual(oldlen, len(self.cdf['SectorRateScalersCounts']))

        message = 'Can only delete entire records.'
        try:
            del self.cdf['SpinRateScalersCounts'][0, 12, 0:5]
        except:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(t, TypeError)
            self.assertEqual(str(v), message)
        else:
            self.fail('Should have raised TypeError: ' + message)

class ChangeAttr(ChangeCDFBase):
    """Tests that modify Attributes and Attribute lists"""

    def testChangezEntry(self):
        """Write new or changed zEntry"""
        zvar = self.cdf['PhysRecNo']
        zvar.attrs['DEPEND_0'] = 'foobar'
        self.assertEqual('foobar', zvar.attrs['DEPEND_0'])
        self.assertEqual(const.CDF_CHAR.value,
                         cdf.zAttr(self.cdf,
                                   'DEPEND_0').type(zvar._num()))

        zvar.attrs['FILLVAL'] = [0, 1]
        numpy.testing.assert_array_equal([0,1], zvar.attrs['FILLVAL'])
        self.assertEqual(const.CDF_INT4.value,
                         cdf.zAttr(self.cdf,
                                   'FILLVAL').type(zvar._num()))

        message = 'Entry strings must be scalar.'
        try:
            zvar.attrs['CATDESC'] = ['hi', 'there']
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should have raised ValueError: ' + message)

        message = 'Entries must be scalar or 1D.'
        try:
            zvar.attrs['FILLVAL'] = [[1, 2], [3, 4]]
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should have raised ValueError: ' + message)

    def testNewzAttr(self):
        """Write a zEntry for a zAttribute that doesn't exist"""
        zvar = self.cdf['PhysRecNo']
        zvar.attrs['NEW_ATTRIBUTE'] = 1
        self.assertTrue('NEW_ATTRIBUTE' in zvar.attrs)
        self.assertEqual(1, zvar.attrs['NEW_ATTRIBUTE'])
        self.assertEqual(const.CDF_INT4.value,
                         cdf.zAttr(self.cdf,
                                   'NEW_ATTRIBUTE').type(zvar._num()))

        zvar.attrs['NEW_ATTRIBUTE2'] = [1, 2]
        numpy.testing.assert_array_equal([1, 2], zvar.attrs['NEW_ATTRIBUTE2'])
        self.assertEqual(const.CDF_INT4.value,
                         cdf.zAttr(self.cdf,
                                   'NEW_ATTRIBUTE2').type(zvar._num()))

        zvar = self.cdf['SpinNumbers']
        zvar.attrs['NEW_ATTRIBUTE3'] = 1
        self.assertEqual(1, zvar.attrs['NEW_ATTRIBUTE3'])
        self.assertEqual(const.CDF_BYTE.value,
                         cdf.zAttr(self.cdf,
                                   'NEW_ATTRIBUTE3').type(zvar._num()))

    def testDelzAttr(self):
        """Delete a zEntry"""
        del self.cdf['PhysRecNo'].attrs['DEPEND_0']
        self.assertFalse('DEPEND_0' in self.cdf['PhysRecNo'].attrs)
        #Make sure attribute still exists
        attrib = cdf.zAttr(self.cdf, 'DEPEND_0')

        del self.cdf['SectorRateScalersCounts'].attrs['DEPEND_3']
        self.assertFalse('DEPEND_3' in
                         self.cdf['SectorRateScalersCounts'].attrs)
        try:
            attrib = cdf.zAttr(self.cdf, 'DEPEND_3')
        except cdf.CDFError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(const.NO_SUCH_ATTR, v.status)
        else:
            self.fail('Should have raised CDFError')

    def testChangegAttr(self):
        """Change an existing gEntry"""
        self.cdf.attrs['Project'][0] = 'not much'
        self.assertEqual('not much',
                         self.cdf.attrs['Project'][0])

        warnings.filterwarnings(
            'ignore', r'^Assuming CDF_.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            self.cdf.attrs['Source_name'][0] = datetime.datetime(2009, 1, 1)
        finally:
            del warnings.filters[0]
        self.assertEqual([datetime.datetime(2009, 1, 1)],
                         self.cdf.attrs['Source_name'][:])

        self.cdf.attrs['Data_type'] = 'stuff'
        self.assertEqual('stuff',
                         self.cdf.attrs['Data_type'][0])
        self.cdf.attrs['Data_type'] = ['stuff', 'more stuff']
        self.assertEqual(['stuff', 'more stuff'],
                         self.cdf.attrs['Data_type'][:])

    def testNewgAttr(self):
        """Create a new gAttr by adding a gEntry"""
        self.cdf.attrs['new_attr'] = 1.5
        self.assertEqual([1.5],
                         self.cdf.attrs['new_attr'][:])

        self.cdf.attrs['new_attr2'] = []
        self.assertTrue('new_attr2' in self.cdf.attrs)
        self.cdf.attrs['new_attr2'][0:6:2] = [1, 2, 3]
        self.assertEqual([1, None, 2, None, 3],
                         self.cdf.attrs['new_attr2'][:])

        self.cdf.attrs['new_attr3'] = ['hello', 'there']
        self.assertEqual(['hello', 'there'],
                         self.cdf.attrs['new_attr3'][:])

    def testDelgAttr(self):
        """Delete a gEntry"""
        del self.cdf.attrs['TEXT'][0]
        self.assertTrue('TEXT' in self.cdf.attrs)

        del self.cdf.attrs['Project'][0]
        self.assertTrue('Project' in self.cdf.attrs)

        del self.cdf.attrs['PI_name']
        self.assertFalse('PI_name' in self.cdf.attrs)

    def testRenamegAttr(self):
        """Rename a gAttribute"""
        textcopy = self.cdf.attrs['TEXT'][:]
        self.cdf.attrs['TEXT'].rename('notTEXT')
        self.assertTrue('notTEXT' in self.cdf.attrs)
        self.assertFalse('TEXT' in self.cdf.attrs)
        self.assertEqual(textcopy, self.cdf.attrs['notTEXT'][:])

    def testRenamezAttr(self):
        """Rename a zAttribute"""
        prn_attrs = self.cdf['PhysRecNo'].attrs
        prn_depend = prn_attrs['DEPEND_0']
        mc_attrs = self.cdf['MeanCharge'].attrs
        mc_depend = mc_attrs['DEPEND_0']
        prn_attrs.rename('DEPEND_0', 'notDEPEND_0')
        self.assertTrue('notDEPEND_0' in prn_attrs)
        self.assertTrue('notDEPEND_0' in mc_attrs)
        self.assertFalse('DEPEND_0' in prn_attrs)
        self.assertFalse('DEPEND_0' in mc_attrs)
        self.assertEqual(prn_depend, prn_attrs['notDEPEND_0'])
        self.assertEqual(mc_depend, mc_attrs['notDEPEND_0'])

    def testChangegEntryType(self):
        """Change the type of a gEntry"""
        attrs = self.cdf.attrs
        attrs['new_attr'] = []
        attrs['new_attr'][0] = [ord('a'), ord('b'), ord('c')]
        attrs['new_attr'].type(0, const.CDF_CHAR)
        self.assertEqual(attrs['new_attr'][0], 'abc')
        try:
            attrs['new_attr'].type(0, const.CDF_INT2)
        except cdf.CDFError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(v.status, const.CANNOT_CHANGE)
        else:
            self.fail('Should have raised CDFError')

    def testChangezEntryType(self):
        """Change the type of a zEntry"""
        attrs = self.cdf['ATC'].attrs
        attrs['new_attr'] = [ord('a'), ord('b'), ord('c')]
        attrs.type('new_attr', const.CDF_CHAR)
        self.assertEqual(attrs['new_attr'], 'abc')
        self.assertEqual(const.CDF_CHAR.value,
                         attrs.type('new_attr'))

    def testzEntryTypeGuessing(self):
        """Guess the type of a zEntry"""
        v = self.cdf.new('newvar', data=[1, 2, 3])
        v.attrs['foo'] = 5
        self.assertEqual(v.type(), v.attrs.type('foo'))

    def testzEntryTypeGuessingMultiElements(self):
        """Guess the type of a zEntry with mulitple elements"""
        v = self.cdf.new('newvar', data=[1, 2, 3])
        v.attrs['foo'] = [5, 3]
        self.assertEqual(v.type(), v.attrs.type('foo'))

    def testgAttrNewEntry(self):
        """Create a new gEntry using Attr.new()"""
        attr = self.cdf.attrs['Project']
        #no type or number
        attr.new([0, 1, 2, 3])
        self.assertEqual(2, len(attr))
        numpy.testing.assert_array_equal([0, 1, 2, 3], attr[1])
        self.assertEqual(const.CDF_BYTE.value, attr.type(1))
        #explicit number
        attr.new('hello there', number=10)
        self.assertEqual(3, len(attr))
        self.assertEqual(10, attr.max_idx())
        self.assertEqual('hello there', attr[10])
        self.assertEqual(const.CDF_CHAR.value, attr.type(10))
        #explicit type and number
        attr.new(10, const.CDF_INT4, 15)
        self.assertEqual(4, len(attr))
        self.assertEqual(15, attr.max_idx())
        self.assertEqual(10, attr[15])
        self.assertEqual(const.CDF_INT4.value, attr.type(15))
        #explicit type
        attr.new([10, 11, 12, 13], const.CDF_REAL8)
        self.assertEqual(5, len(attr))
        numpy.testing.assert_array_equal([10.0, 11.0, 12.0, 13.0], attr[2])
        self.assertEqual(const.CDF_REAL8.value, attr.type(2))

    def testgAttrListNew(self):
        """Create a new gAttr and/or gEntry using gAttrList.new"""
        attrs = self.cdf.attrs
        attrs.new('new')
        self.assertTrue('new' in attrs)
        attrs.new('new2', [1, 2, 3])
        self.assertTrue('new2' in attrs)
        numpy.testing.assert_array_equal([1, 2, 3], attrs['new2'][0])
        attrs.new('new3', [1, 2, 3], const.CDF_INT4)
        self.assertTrue('new3' in attrs)
        numpy.testing.assert_array_equal([1, 2, 3], attrs['new3'][0])
        self.assertEqual(const.CDF_INT4.value, attrs['new3'].type(0))

    def testzAttrListNew(self):
        """Create a new zEntry using zAttrList.new"""
        attrs = self.cdf['ATC'].attrs
        attrs.new('new2', [1, 2, 3])
        self.assertTrue('new2' in attrs)
        numpy.testing.assert_array_equal([1, 2, 3], attrs['new2'])
        attrs.new('new3', [1, 2, 3], const.CDF_INT4)
        self.assertTrue('new3' in attrs)
        numpy.testing.assert_array_equal([1, 2, 3], attrs['new3'])
        self.assertEqual(const.CDF_INT4.value, attrs.type('new3'))

    def testAttrsFromDict(self):
        """Dump a bunch of attrs on a variable from a dict, using clone"""
        indict = { 'CATDESC': numpy.array([1, 2, 3], dtype=numpy.int32),
                   'b': 'hello',
                   }
        attrlist = self.cdf['ATC'].attrs
        attrlist.clone(indict)
        self.assertEqual(['CATDESC', 'b'], sorted(attrlist.keys()))
        numpy.testing.assert_array_equal(indict['CATDESC'],
                                         attrlist['CATDESC'])
        self.assertEqual('hello', attrlist['b'])
        types = {'CATDESC': const.CDF_INT4.value,
                 'b': const.CDF_CHAR.value,
                 }
        for k in types:
            self.assertEqual(types[k], attrlist.type(k))

    def testgAttrsAssign(self):
        """Assign to the attrs attribute of CDF"""
        self.cdf.attrs = {'foobar': ['global']}
        self.cdf.close()
        self.cdf = cdf.CDF(self.testfile) #reopen
        self.assertFalse(isinstance(self.cdf.attrs, dict))
        self.assertEqual(self.cdf.attrs['foobar'][0], 'global')
        self.assertEqual(len(self.cdf.attrs['foobar']), 1)
        self.assertFalse('TEXT' in self.cdf.attrs)

    def testzAttrsAssign(self):
        """Assign to the attrs attribute of variable"""
        self.cdf['ATC'].attrs = {'foobar': ['var']}
        self.cdf.close()
        self.cdf = cdf.CDF(self.testfile) #reopen
        self.assertFalse(isinstance(self.cdf['ATC'].attrs, dict))
        self.assertEqual(self.cdf['ATC'].attrs['foobar'], 'var')
        self.assertFalse('CATDESC' in self.cdf['ATC'].attrs)

    def testzAttrsAssignTimeType(self):
        """Assign a time type to a zAttr"""
        self.cdf['ATC'].attrs['testtime'] = datetime.datetime(2010, 1, 1)
        expected = cdf.const.CDF_EPOCH16.value # Matches var
        self.assertEqual(expected, self.cdf['ATC'].attrs.type('testtime'))
        with spacepy_testing.assertWarns(
                self, 'always',
                r'Assuming CDF_TIME_TT2000 for time input\.$',
                DeprecationWarning, r'spacepy\.pycdf$'):
            self.cdf['SectorRateScalersCounts'].attrs['testtime'] \
                = datetime.datetime(2010, 1, 1)
        # Assigned to attribute of non-time variable
        expected = cdf.const.CDF_TIME_TT2000.value if cdf.lib.supports_int8 \
                   else cdf.const.CDF_EPOCH.value
        self.assertEqual(
            expected,
            self.cdf['SectorRateScalersCounts'].attrs.type('testtime'))

    def testgAttrsAssignTimeType(self):
        """Assign a time type to a gAttr"""
        with spacepy_testing.assertWarns(
                self, 'always',
                r'Assuming CDF_TIME_TT2000 for time input\.$',
                DeprecationWarning, r'spacepy\.pycdf$'):
            self.cdf.attrs['testtime'] = datetime.datetime(2010, 1, 1)
        expected = cdf.const.CDF_TIME_TT2000.value if cdf.lib.supports_int8 \
                   else cdf.const.CDF_EPOCH.value
        self.assertEqual(expected, self.cdf.attrs['testtime'].type(0))

    def testzAttrsDelete(self):
        """Try to delete attrs attribute of variable, CDF"""
        try:
            del self.cdf['ATC'].attrs
        except AttributeError:
            pass
        else:
            self.fail('AttributeError not raised.')
        try:
            del self.cdf.attrs
        except AttributeError:
            pass
        else:
            self.fail('AttributeError not raised.')

    def testCopyAttr(self):
        """Assign a gAttribute to another"""
        self.cdf.attrs['new_attr'] = self.cdf.attrs['TEXT']
        old_attr = self.cdf.attrs['TEXT']
        new_attr = self.cdf.attrs['new_attr']
        for i in range(self.cdf.attrs['TEXT'].max_idx()):
            self.assertEqual(old_attr.has_entry(i),
                             new_attr.has_entry(i))
            if old_attr.has_entry(i):
                self.assertEqual(old_attr[i], new_attr[i])
                self.assertEqual(old_attr.type(i),
                                 new_attr.type(i))

    def testCloneAttrList(self):
        """Copy an entire attribute list from one CDF to another"""
        warnings.filterwarnings(
            'ignore', r'^spacepy\.pycdf\.lib\.set_backward not called.*',
            DeprecationWarning, r'^spacepy\.pycdf$')
        try:
            with cdf.CDF('attrcopy.cdf', '') as newcdf:
                newcdf.attrs['deleteme'] = ['hello']
                newcdf.attrs.clone(self.cdf.attrs)
                for attrname in self.cdf.attrs:
                    self.assertTrue(attrname in newcdf.attrs)
                    old_attr = self.cdf.attrs[attrname]
                    new_attr = newcdf.attrs[attrname]
                    self.assertEqual(old_attr.max_idx(),
                                     new_attr.max_idx())
                    for i in range(old_attr.max_idx()):
                        self.assertEqual(old_attr.has_entry(i),
                                         new_attr.has_entry(i))
                        if old_attr.has_entry(i):
                            self.assertEqual(old_attr[i], new_attr[i])
                            self.assertEqual(old_attr.type(i),
                                             new_attr.type(i))
                for attrname in newcdf.attrs:
                    self.assertTrue(attrname in self.cdf.attrs)
        finally:
            del warnings.filters[0]
            os.remove('attrcopy.cdf')

    def testClonezAttrList(self):
        """Copy entire attribute list from one zVar to another"""
        oldlist = self.cdf['ATC'].attrs
        newlist = self.cdf['PhysRecNo'].attrs
        newlist.clone(oldlist)
        for attrname in oldlist:
            self.assertTrue(attrname in newlist)
            self.assertEqual(oldlist[attrname], newlist[attrname])
            self.assertEqual(oldlist.type(attrname),
                             newlist.type(attrname))
        oldlist = self.cdf['Epoch'].attrs
        newlist = self.cdf['MeanCharge'].attrs
        newlist.clone(oldlist)
        for attrname in oldlist:
            self.assertTrue(attrname in newlist,
                            'Attribute {0} not found in copy.'.format(attrname)
                            )
            self.assertEqual(oldlist[attrname], newlist[attrname])
            self.assertEqual(oldlist.type(attrname),
                             newlist.type(attrname))

    def testUnicodeAttribute(self):
        """Assign unicode string to attributes"""
        self.cdf['ATC'].attrs['foo'] = u'C'
        self.assertEqual('C', self.cdf['ATC'].attrs['foo'])
        try:
            self.cdf['ATC'].attrs['foo2'] = u'\xb0C'
        except UnicodeEncodeError:
            pass
        else:
            self.fail('Should have raised UnicodeEncodeError')
        #This basically fails same way as numpy.array(u'\xb0C', dtype='|S2')
#        self.assertEqual(b'C', self.cdf['ATC'].attrs['foo2'])


class ChangeColCDF(ColCDFTests):
    """Tests that modify an existing colum-major CDF"""
    def __init__(self, *args):
        super(ChangeColCDF, self).__init__(*args)

    def setUp(self):
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)
        self.cdf.readonly(False)

    def tearDown(self):
        self.cdf.close()
        del self.cdf
        os.remove(self.testfile)

    def testWriteColSubscripted(self):
        """Write data to a slice of a zVar"""
        expected = ['0', '1', '99', '3', '98', '5', '97', '7',
                    '8', '9']
        self.cdf['SpinNumbers'][2:7:2] = ['99', '98', '97']
        numpy.testing.assert_array_equal(
            expected, numpy.char.rstrip(self.cdf['SpinNumbers'][0:10]))

        expected = self.cdf['SectorRateScalersCounts'][...]
        expected[4][5][5][8:3:-1] = [101.0, 102.0, 103.0, 104.0, 105.0]
        self.cdf['SectorRateScalersCounts'][4, 5, 5, 8:3:-1] = \
            [101.0, 102.0, 103.0, 104.0, 105.0]
        numpy.testing.assert_array_equal(
            expected, self.cdf['SectorRateScalersCounts'][...])


if __name__ == '__main__':
    unittest.main()
