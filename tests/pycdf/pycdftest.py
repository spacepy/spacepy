#!/usr/bin/env python

"""Unit test suite for pycdf"""

import ctypes
import datetime
import hashlib
import os, os.path
import shutil
import sys
import unittest
import time

try:
    type(callable)
except NameError:
    import collections
    def callable(obj):
        return isinstance(obj, collections.Callable)

import pycdf as cdf
import spacepy.pycdf.const as const


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
            self.assertEqual(cdf._pycdf._Hyperslice.reorder(inp), outp)

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
            result = cdf._pycdf._Hyperslice.convert_range(*inp)
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
        self.assertEqual(cdf._pycdf._Hyperslice.dimensions(data),
                         [4, 3, 2])

        data = [[[2, 3], [4, 5], [6, 7]],
                [[8, 9], [0, 1], [2, 3]],
                [[4, 5], [6, 7],],
                [[0, 1], [2, 3], [4, 5]],
                ]
        message = 'Data irregular in dimension 1'
        try:
            cdf._pycdf._Hyperslice.dimensions(data)
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should raise ValueError: ' + message)

        self.assertEqual(cdf._pycdf._Hyperslice.dimensions('hi'),
                         [])

    def testFlipMajority(self):
        """Changes the majority of an array"""
        #Code to generate this 5x5x5:
        #[[[random.randint(0,100) for i in range(5)]
        #  for j in range(5)] for k in range(5)]
        three_d = [[[90, 90, 96, 71, 90], [29, 18, 90, 78, 51],
                   [14, 29, 41, 25, 50], [73, 59, 83, 92, 24],
                   [10, 1, 4, 61, 54]],
                  [[40, 8, 0, 28, 47], [3, 98, 28, 9, 38],
                   [34, 95, 7, 87, 9], [11, 73, 71, 54, 69],
                   [42, 75, 82, 16, 73]],
                  [[88, 40, 5, 69, 41], [35, 15, 32, 68, 8],
                   [68, 74, 6, 30, 9], [86, 48, 52, 49, 100],
                   [8, 35, 26, 16, 61]],
                  [[49, 81, 57, 37, 98], [54, 64, 28, 21, 17],
                   [73, 100, 90, 8, 25], [40, 75, 52, 41, 40],
                   [42, 72, 55, 16, 39]],
                  [[24, 38, 26, 85, 25], [5, 98, 63, 29, 33],
                   [91, 100, 17, 85, 9], [59, 50, 50, 41, 82],
                   [21, 45, 65, 51, 90]]]
        flipped = cdf._pycdf._Hyperslice.flip_majority(three_d)
        for i in range(5):
            for j in range(5):
                for k in range(5):
                    self.assertEqual(three_d[i][j][k],
                                     flipped[k][j][i],
                                     'Original index ' +
                                     str(i) + ', ' +
                                     str(j) + ', ' +
                                     str(k) + ' mismatch ' +
                                     str(three_d[i][j][k]) + ' != ' +
                                     str(flipped[k][j][i]))

        #[[[[random.randint(0,100) for i in range(5)]
        #  for j in range(4)] for k in range(3)] for l in range(2)]
        four_d = [[[[14, 84, 79, 74, 45], [39, 47, 93, 32, 59],
                    [15, 47, 1, 84, 44], [13, 43, 13, 88, 3]],
                   [[65, 75, 36, 90, 93], [64, 36, 59, 39, 42],
                    [59, 85, 21, 88, 61], [64, 29, 62, 33, 35]],
                   [[46, 69, 3, 50, 44], [86, 15, 32, 17, 51],
                    [79, 20, 29, 10, 55], [29, 10, 79, 7, 58]]],
                  [[[20, 76, 81, 40, 85], [44, 56, 5, 83, 32],
                    [34, 88, 23, 57, 74], [24, 55, 83, 39, 60]],
                   [[79, 56, 5, 98, 29], [28, 50, 77, 33, 45],
                    [38, 82, 82, 28, 97], [42, 14, 56, 48, 38]],
                   [[58, 27, 38, 43, 25], [72, 91, 85, 44, 43],
                    [17, 57, 91, 19, 35], [98, 62, 61, 14, 60]]]]
        flipped = cdf._pycdf._Hyperslice.flip_majority(four_d)
        for i in range(2):
            for j in range(3):
                for k in range(4):
                    for l in range(5):
                        self.assertEqual(four_d[i][j][k][l],
                                         flipped[l][k][j][i],
                                         'Original index ' +
                                         str(i) + ', ' +
                                         str(j) + ', ' +
                                         str(k) + ', ' +
                                         str(l) + ' mismatch ' +
                                         str(four_d[i][j][k][l]) + ' != ' +
                                         str(flipped[l][k][j][i]))

        zero_d = 1
        flipped = cdf._pycdf._Hyperslice.flip_majority(zero_d)
        self.assertEqual(zero_d, flipped)

        one_d = [1, 2, 3, 4]
        flipped = cdf._pycdf._Hyperslice.flip_majority(one_d)
        self.assertEqual(one_d, flipped)

        two_d = [[6, 7, 48, 81], [61, 67, 90, 99], [71, 96, 58, 85],
                 [35, 31, 71, 73], [77, 41, 71, 92], [74, 89, 94, 64],
                 [64, 30, 66, 94]]
        flipped = cdf._pycdf._Hyperslice.flip_majority(two_d)
        for i in range(7):
            for j in range(4):
                self.assertEqual(two_d[i][j],
                                 flipped[j][i],
                                 'Original index ' +
                                 str(i) + ', ' +
                                 str(j) + ' mismatch ' +
                                 str(two_d[i][j]) + ' != ' +
                                 str(flipped[j][i]))

    def testEpoch16ToDatetime(self):
        epochs = [[63397987199.0, 999999999999.0],
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(dt, cdf.lib.epoch16_to_datetime(epoch))

    def testEpochToDatetime(self):
        epochs = [63397987200000.0,
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(dt, cdf.lib.epoch_to_datetime(epoch))

    def testDatetimeToEpoch16(self):
        epochs = [[63397987200.0, 0.0],
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(epoch, cdf.lib.datetime_to_epoch16(dt))

    def testDatetimeToEpoch(self):
        epochs = [63397987200000.0,
                  ]
        dts = [datetime.datetime(2009, 1, 1),
               ]
        for (epoch, dt) in zip(epochs, dts):
            self.assertEqual(epoch, cdf.lib.datetime_to_epoch(dt))

    def testDatetimeEpoch16RT(self):
        """Roundtrip datetimes to epoch16s and back"""
        dts = [datetime.datetime(2008, 12, 15, 3, 12, 5, 1000),
               datetime.datetime(1821, 1, 30, 2, 31, 5, 23000),
               datetime.datetime(2050, 6, 5, 15, 0, 5, 0),
               ]
        for dt in dts:
            self.assertEqual(dt, cdf.lib.epoch16_to_datetime(
                cdf.lib.datetime_to_epoch16(dt)))

    def testDatetimeEpochRT(self):
        """Roundtrip datetimes to epochs and back"""
        dts = [datetime.datetime(2008, 12, 15, 3, 12, 5, 1000),
               datetime.datetime(1821, 1, 30, 2, 31, 5, 23000),
               datetime.datetime(2050, 6, 5, 15, 0, 5, 0),
               ]
        for dt in dts:
            self.assertEqual(dt, cdf.lib.epoch_to_datetime(
                cdf.lib.datetime_to_epoch(dt)))

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
        self.assertTrue(cdf.lib.version[3] in (b'', b' ', b'a'))
        if cdf.lib.version == (3, 3, 0, ' '):
            self.assertTrue(cdf.lib._del_middle_rec_bug)
        elif cdf.lib.version == (3, 3, 1, ' '):
            self.assertFalse(cdf.lib._del_middle_rec_bug)

    def testTypeGuessing(self):
        """Guess CDF types based on input data"""
        samples = [[1, 2, 3, 4],
                   [[1.2, 1.3, 1.4], [2.2, 2.3, 2.4]],
                   ['hello', 'there', 'everybody'],
                   datetime.datetime(2009, 1, 1),
                   datetime.datetime(2009, 1, 1, 12, 15, 12, 1),
                   ]
        types = [([4], [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                        const.CDF_INT2, const.CDF_UINT2,
                        const.CDF_INT4, const.CDF_UINT4,
                        const.CDF_FLOAT, const.CDF_REAL4,
                        const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ([2, 3], [const.CDF_FLOAT, const.CDF_REAL4,
                           const.CDF_DOUBLE, const.CDF_REAL8], 1),
                 ([3], [const.CDF_CHAR, const.CDF_UCHAR], 9),
                 ([], [const.CDF_EPOCH, const.CDF_EPOCH16], 1),
                 ([], [const.CDF_EPOCH16, const.CDF_EPOCH], 1),
                 ]
        for (s, t) in zip(samples, types):
            t = (t[0], [i.value for i in t[1]], t[2])
            self.assertEqual(t, cdf._pycdf._Hyperslice.types(s))


class MakeCDF(unittest.TestCase):
    def setUp(self):
        self.testfspec='foo.cdf'
    
    def testOpenCDFNew(self):
        """Create a new CDF"""

        new = cdf.CDF(self.testfspec, '')
        self.assertTrue(os.path.isfile(self.testfspec))
        os.remove(self.testfspec)

    def testOpenCDFNonexistent(self):
        """Open a CDF which doesn't exist"""

        self.assertRaises(cdf.CDFError, cdf.CDF, self.testfspec)

    def testOpenCDFNoMaster(self):
        """Open a CDF from a master CDF which doesn't exist"""

        self.assertRaises(IOError, cdf.CDF, self.testfspec, 'nonexist.cdf')


class CreateVar(unittest.TestCase):
    def setUp(self):
        self.cdf = cdf.CDF('new.cdf', '')
        self.cdf.readonly(False)
        
    def tearDown(self):
        os.remove('new.cdf')

    def testCreateVarDefault(self):
        """Create a variable with default options in an empty CDF"""
        new = self.cdf._new_var('new_variable', cdf.const.CDF_UINT1)
        self.assertTrue(type(new).__name__ == 'Var')


class CDFTests(unittest.TestCase):
    """Tests that involve an existing CDF, read or write"""
    testmaster = 'po_l1_cam_test.cdf'
    testfile = 'test.cdf'

    def __init__(self, *args):
        self.expected_digest = 'b31bffefd7a63de6ce9856f08dbae43d'
        assert(self.calcDigest(self.testmaster) == self.expected_digest)
        super(CDFTests, self).__init__(*args)

    def calcDigest(self, file):
        """Calculate the MD5 digest from a file.

        @param file: path to the file
        @type file: string
        @return: hex digits of L{file}'s md5
        @rtype: string
        """
        m = hashlib.md5()
        with open(file, 'rb') as f:
            m.update(f.read())
        return m.hexdigest()


class ColCDFTests(unittest.TestCase):
    """Tests that involve an existing column-major CDF, read or write"""
    testmaster = 'po_l1_cam_testc.cdf'
    testfile = 'testc.cdf'

    def __init__(self, *args):
        self.expected_digest = '7728439e20bece4c0962a125373345bf'
        assert(self.calcDigest(self.testmaster) == self.expected_digest)
        super(ColCDFTests, self).__init__(*args)

    def calcDigest(self, file):
        """Calculate the MD5 digest from a file.

        @param file: path to the file
        @type file: string
        @return: hex digits of L{file}'s md5
        @rtype: string
        """
        m = hashlib.md5()
        with open(file, 'rb') as f:
            m.update(f.read())
        return m.hexdigest()


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
                    'MajorNumbers', 'MeanCharge', 'Epoch']
        with cdf.CDF(self.testfile) as f:
            names = list(f.keys())
        self.assertEqual(expected, names)
        self.assertRaises(cdf.CDFError, f.close)


class ReadCDF(CDFTests):
    """Tests that read an existing CDF, but do not modify it."""
    testfile = 'test_ro.cdf'

    def __init__(self, *args, **kwargs):
        super(ReadCDF, self).__init__(*args, **kwargs)
        #Unittest docs say 'the order in which the various test cases will be
        #run is determined by sorting the test function names with the built-in
        #cmp() function'
        testnames = [name for name in dir(self)
                     if name[0:4] == 'test' and callable(getattr(self,name))]
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
        expectedNames = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
                         'RateScalerNames', 'SectorRateScalerNames',
                         'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
                         'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
                         'MajorNumbers', 'MeanCharge', 'Epoch']
        names = [zVar.name() for zVar in self.cdf.values()]
        self.assertEqual(names, expectedNames)

    def testGetAllVarNames(self):
        """Getting a list of zVar names"""
        expectedNames = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
                         'RateScalerNames', 'SectorRateScalerNames',
                         'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
                         'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
                         'MajorNumbers', 'MeanCharge', 'Epoch']
        names = list(self.cdf.keys())
        self.assertEqual(expectedNames, names)

    def testGetVarNum(self):
        self.assertEqual(0, self.cdf['ATC']._num())

    def testCDFIterator(self):
        expected = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
                    'RateScalerNames', 'SectorRateScalerNames',
                    'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
                    'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
                    'MajorNumbers', 'MeanCharge', 'Epoch']
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
        self.assertEqual(len(self.cdf['ATC']), 747)
        self.assertEqual(len(self.cdf['MeanCharge']), 100)
        self.assertEqual(len(self.cdf['SpinNumbers']), 1)

    def testMajority(self):
        """Get majority of the CDF"""
        self.assertEqual(self.cdf._majority().value, cdf.const.ROW_MAJOR.value)

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
            self.assertEqual(self.cdf[i]._rec_vary(), expected[i])

    def testHyperslices(self):
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
            sliced = cdf._pycdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts, sliced.intervals, sliced.degen,
                      sliced.rev)
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)

    def testHyperslices2(self):
        """Additional checks: converting python slices to CDF counts, etc."""
        slices = {'ATC': Ellipsis,
                  } #Slice objects indexed by variable
        #Expected results [dims, dimsizes, starts, counts, intervals, degen, rev]
        #indexed by variable
        expected = {'ATC': [1, [747], [0], [747], [1], [False], [False]],
                    }
        for i in expected:
            zvar = self.cdf[i]
            sliced = cdf._pycdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts, sliced.intervals, sliced.degen,
                      sliced.rev)
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)

    def testHypersliceExpand(self):
        """Expand a slice to store the data passed in"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, None, 1))
        self.assertEqual(100, sliced.counts[0])
        sliced.expand(list(range(110)))
        self.assertEqual(110, sliced.counts[0])
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, 100, 2))
        sliced.expand(list(range(110)))
        self.assertEqual(50, sliced.counts[0])

    def testHypersliceExpectedDims(self):
        """Find dimensions expected by a slice"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, None, 1))
        self.assertEqual([100], sliced.expected_dims())
        sliced.expand(list(range(110)))
        self.assertEqual([110], sliced.expected_dims())
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, 100, 2))
        sliced.expand(list(range(110)))
        self.assertEqual([50], sliced.expected_dims())

        zvar = self.cdf['SpinRateScalersCounts']
        sliced = cdf._pycdf._Hyperslice(zvar, (slice(None, None, None),
                                               slice(None, None, 2),
                                               slice(0, None, 3)))
        self.assertEqual([100, 9, 6], sliced.expected_dims())

    def testPackBuffer(self):
        """Pack a buffer with data"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, None, 1))
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, list(range(100)))
        result = [buff[i] for i in range(100)]
        self.assertEqual(list(range(100)), result)
        
        sliced = cdf._pycdf._Hyperslice(zvar, slice(None, None, -1))
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, list(range(100)))
        result = [buff[i] for i in range(100)]
        self.assertEqual(list(reversed(range(100))), result)

        zvar = self.cdf['SectorRateScalersCounts']
        sliced = cdf._pycdf._Hyperslice(zvar, (0, slice(0, 3, 2),
                                              slice(3, None, -1), 1))
        buff = sliced.create_buffer()
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8]]
        expected = [[4, 3, 2, 1],
                    [8, 7, 6, 5]]
        sliced.pack_buffer(buff, data)
        self.assertEqual(expected[0], list(buff[0]))
        self.assertEqual(expected[1], list(buff[1]))

        sliced = cdf._pycdf._Hyperslice(zvar, (0, 1, 1, 1))
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, 10)
        self.assertEqual(10, buff.value)

        zvar = self.cdf['ATC']
        sliced = cdf._pycdf._Hyperslice(zvar, 10)
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, datetime.datetime(2009, 1, 1))
        self.assertEqual(list(buff), [63397987200.0, 0.0])

        zvar = self.cdf['SpinNumbers']
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, 3, None))
        buff = sliced.create_buffer()
        expected = [b'99', b'98', b'97']
        sliced.pack_buffer(buff, expected)
        for i in range(3):
            self.assertEqual(expected[i], buff[i].value)

    def testArraytoList(self):
        """Converts a ctypes array to Python nested lists"""
        sliced = cdf._pycdf._Hyperslice(self.cdf['PhysRecNo'], slice(0, 4, None))
        array_in = (ctypes.c_int * 4) (1, 2, 3, 4)
        expected = [1, 2, 3, 4]
        self.assertEqual(expected, sliced.c_array_to_list(array_in))

        sliced = cdf._pycdf._Hyperslice(self.cdf['SectorRateScalersCounts'],
                                (slice(5,10,2), slice(0,3,1),
                                 2, Ellipsis))
        array_in = (ctypes.c_float * 9 * 3 * 3)(
            (ctypes.c_float * 9 * 3)(
            (ctypes.c_float * 9)(1, 2, 3, 4, 5, 6, 7, 8, 9),
            (ctypes.c_float * 9)(1, 2, 3, 4, 5, 6, 7, 8, 9),
            (ctypes.c_float * 9)(1, 2, 3, 4, 5, 6, 7, 8, 9)
            ),
            (ctypes.c_float * 9 * 3)(
            (ctypes.c_float * 9)(11, 12, 13, 14, 15, 16, 17, 18, 19),
            (ctypes.c_float * 9)(11, 12, 13, 14, 15, 16, 17, 18, 19),
            (ctypes.c_float * 9)(11, 12, 13, 14, 15, 16, 17, 18, 19)
            ),
            (ctypes.c_float * 9 * 3)(
            (ctypes.c_float * 9)(21, 22, 23, 24, 25, 26, 27, 28, 29),
            (ctypes.c_float * 9)(21, 22, 23, 24, 25, 26, 27, 28, 29),
            (ctypes.c_float * 9)(21, 22, 23, 24, 25, 26, 27, 28, 29)
            ),
            )
        expected = [
            [
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
            ],
            [
            [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0],
            [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0],
            [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0],
            ],
            [
            [21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0],
            [21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0],
            [21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0],
            ]
            ]
        self.assertEqual(expected, sliced.c_array_to_list(array_in))

        array_in = (ctypes.c_char * 27)\
                   (*tuple([i for i in '123456789abcdefghiABCDEFGHI']))
        array_in = ctypes.cast(array_in,
                               ctypes.POINTER(ctypes.c_char * 9 * 3)).contents
        sliced = cdf._pycdf._Hyperslice(self.cdf['SectorRateScalerNames'],
                                slice(0,3,1))
        array_out = sliced.c_array_to_list(array_in)
        if isinstance(array_out[0], bytes):
            self.assertEqual([b'123456789', b'abcdefghi', b'ABCDEFGHI'],
                             array_out)
        else:
            self.assertEqual(['123456789', 'abcdefghi', 'ABCDEFGHI'],
                             array_out)

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
                             self.cdf[i]._cdf_type())

    def testCTypes(self):
        """Look up ctype to match variable"""
        expected = {'ATC': (ctypes.c_double * 2),
                    'PhysRecNo': ctypes.c_int,
                    'SpinNumbers': ctypes.c_char * 2,
                    'MeanCharge': ctypes.c_float,
                    'Epoch': ctypes.c_double,
                    }
        for i in expected:
            self.assertEqual(expected[i],
                             self.cdf[i]._c_type())

    def testCreateReadArrays(self):
        """Create an array for reading from a slice"""
        slices = {'ATC': 1,
                  'PhysRecNo': slice(10, 2, -2),
                  'SpinNumbers': slice(2, None, 2),
                  'SectorRateScalersCounts': (slice(3, 6, None),
                                              slice(None, None, None),
                                              slice(None, None, None)),
                  'SpinRateScalersCounts': (Ellipsis, slice(-1, None, -1)),
                  } #Slice objects indexed by variable
        types = {'ATC': (ctypes.c_double * 2),
                 'PhysRecNo': (ctypes.c_int * 4),
                 'SpinNumbers': (ctypes.c_char * 2 * 8),
                 'SectorRateScalersCounts': (ctypes.c_float *
                                             9 * 32 * 3 * 100),
                 'SpinRateScalersCounts': (ctypes.c_float * 16 * 18 * 100),
                 } #constructors, build inside out
        for i in types:
            zvar = self.cdf[i]
            sliced = cdf._pycdf._Hyperslice(zvar, slices[i])
            actual = sliced.create_buffer().__class__
            self.assertEqual(types[i], actual)

    def testSubscriptVariable(self):
        """Refer to an array by subscript"""
        self.assertEqual([3, 25, 47],
                         self.cdf['PhysRecNo'][0:5:2])
        self.assertEqual([1094, 1083, 1072, 1061],
                         self.cdf['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         self.cdf['SpinRateScalersCounts'][41, 2, 15])

    def testIncompleteSubscript(self):
        """Get data from a variable with a less-than-complete specification"""
        chargedata = self.cdf['MeanCharge'][0] #Should be the first record
        self.assertEqual(len(chargedata), 16)
        SpinRateScalersCounts = self.cdf['SpinRateScalersCounts'][...]
        self.assertEqual(100, len(SpinRateScalersCounts))

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
        self.assertEqual(expected,
                         self.cdf['ATC'][4:8])

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
        self.assertEqual(expected,
                         self.cdf['Epoch'][4:8])

    def testnElems(self):
        """Read number of elements in a string variable"""
        self.assertEqual(2, self.cdf['SpinNumbers']._nelems())
        self.assertEqual(2, self.cdf['SectorNumbers']._nelems())

    def testSubscriptString(self):
        """Refer to a string array by subscript"""
        self.assertEqual(['0 ', '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ',
                          '8 ', '9 ', '10', '11', '12', '13', '14', '15',
                          '16', '17'],
                         self.cdf['SpinNumbers'][:])

    def testSubscriptIrregString(self):
        """Refer to a variable-length string array by subscript"""
        self.assertEqual(['H+', 'He+', 'He++', 'O<=+2', 'O>=+3', 'CN<=+2',
                          'H0', 'He0', 'CNO0', 'CN>=+3', 'Ne-Si', 'S-Ni',
                          '3He', 'D', 'Molecules', 'Others'],
                         self.cdf['RateScalerNames'][:])

    def testGetAllNRV(self):
        """Get an entire non record varying variable"""
        self.assertEqual(['0 ', '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ',
                          '8 ', '9 ', '10', '11', '12', '13', '14', '15',
                          '16', '17'],
                         self.cdf['SpinNumbers'][...])

    def testGetsingleNRV(self):
        """Get single element of non record varying variable"""
        self.assertEqual('0 ',
                         self.cdf['SpinNumbers'][0])

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
        expected = ['ATC', 'Epoch', 'MajorNumbers', 'MeanCharge',
                    'PhysRecNo',
                    'RateScalerNames', 'SectorNumbers',
                    'SectorRateScalerNames',
                    'SectorRateScalersCounts',
                    'SectorRateScalersCountsSigma',
                    'SpinNumbers',
                    'SpinRateScalersCounts',
                    'SpinRateScalersCountsSigma',
                    ]
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
        self.assertEqual(13, result)

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

    def testzEntryType(self):
        """Get the type of a zEntry"""
        names = ['DEPEND_0', 'VALIDMAX', ]
        numbers = [1, 0, ]
        types = [cdf.const.CDF_CHAR, cdf.const.CDF_EPOCH16, ]
        for (name, number, cdf_type) in zip(names, numbers, types):
            attribute = cdf.zAttr(self.cdf, name)
            actual_type = attribute.entry_type(number)
            self.assertEqual(actual_type, cdf_type.value,
                             'zAttr ' + name + ' zEntry ' + str(number) +
                             ' ' + str(cdf_type.value) + ' != ' + 
                             str(actual_type))

    def testgEntryType(self):
        """Get the type of a gEntry"""
        names = ['PI_name', 'Project', ]
        numbers = [0, 0, ]
        types = [cdf.const.CDF_CHAR, cdf.const.CDF_CHAR, ]
        for (name, number, cdf_type) in zip(names, numbers, types):
            attribute = cdf.gAttr(self.cdf, name)
            actual_type = attribute.entry_type(number)
            self.assertEqual(actual_type, cdf_type.value,
                             'gAttr ' + name + ' gEntry ' + str(number) +
                             ' ' + str(cdf_type.value) + ' != ' +
                             str(actual_type))

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
        self.assertEqual([b'CATDESC', b'DEPEND_0', b'FIELDNAM', b'FILLVAL',
                          b'FORMAT', b'VALIDMIN', b'VALIDMAX', b'VAR_TYPE'],
                         list(attrlist))
        
    def testgAttribListIt(self):
        """Iterate over keys in a gAttrList"""
        attrlist = cdf.gAttrList(self.cdf)
        self.assertEqual([b'Project', b'Source_name', b'Discipline',
                          b'Data_type', b'Descriptor',
                          b'File_naming_convention', b'Data_version',
                          b'PI_name', b'PI_affiliation', b'TEXT',
                          b'Instrument_type', b'Mission_group',
                          b'Logical_source',
                          b'Logical_file_id', b'Logical_source_description',
                          b'Time_resolution', b'Rules_of_use', b'Generated_by',
                          b'Generation_date', b'Acknowledgement', b'MODS',
                          b'ADID_ref', b'LINK_TEXT', b'LINK_TITLE',
                          b'HTTP_LINK',
                          ],
                         list(attrlist))

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

    def testzVarCopy(self):
        """Make a copy of an entire zVar"""
        zvar = self.cdf['PhysRecNo']
        zvarcopy = zvar.copy()
        self.assertNotEqual(zvar, zvarcopy)
        for i in range(len(zvar)):
            self.assertEqual(zvar[i], zvarcopy[i])
        for i in zvarcopy.attrs:
            self.assertEqual(zvar.attrs[i], zvarcopy.attrs[i])
        self.assertEqual(zvar[...], zvarcopy[...])

    def testCDFCopy(self):
        """Make a copy of an entire CDF"""
        cdfcopy = self.cdf.copy()
        self.assertNotEqual(cdfcopy, self.cdf)
        for key in self.cdf:
            self.assertEqual(self.cdf[key][...], cdfcopy[key])
            self.assertNotEqual(self.cdf[key], cdfcopy[key])
        for key in self.cdf.attrs:
            self.assertEqual(self.cdf.attrs[key][:], cdfcopy.attrs[key])
            self.assertNotEqual(self.cdf.attrs[key], cdfcopy.attrs[key])

    def testSliceCDFCopy(self):
        """Slice a copy of a CDF"""
        cdfcopy = self.cdf.copy()
        self.assertEqual([3, 25, 47],
                         cdfcopy['PhysRecNo'][0:5:2])
        self.assertEqual([1094, 1083, 1072, 1061],
                         cdfcopy['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         cdfcopy['SpinRateScalersCounts'][41, 2, 15])

    def testVarString(self):
        """Convert a variable to a string representation"""
        for varname in self.cdf:
            var = self.cdf[varname]
            varcopy = var.copy()
            self.assertEqual(str(varcopy), str(var))

    def testCDFString(self):
        """Convert a CDF to a string representation"""
        #a bit funky because the order of keys may be different
        #self.assertEqual(eval(str(self.cdf)), self.cdf.copy())
        #a much cheaper way to do it; lengths different in py3k
        self.assertTrue(len(str(self.cdf)) in (7523076, 7614589))

    def testgAttrListString(self):
        """Convert a list of gattributes to a string"""
        #self.assertEqual(self.cdf.copy().attrs, eval(str(self.cdf.attrs)))
        self.assertTrue(len(str(self.cdf.attrs)) in (983, 1008))

    def testzAttrListString(self):
        """Convert a list of zAttributes to a string"""
        #for varname in self.cdf:
        #    var = self.cdf[varname]
        #    varcopy = var.copy()
        #    self.assertEqual(varcopy.attrs, eval(str(var.attrs)))
        expected = [[373, 210, 114, 131, 174, 192, 695, 561, 597, 493,
                     128, 396, 317],
                    [383, 218, 118, 135, 179, 197, 713, 573, 613, 504,
                     132, 407, 325],
                    [373, 210, 114, 131, 174, 192, 694, 560, 596,
                     492, 128, 395, 317]
                    ]
        self.assertTrue([len(str(self.cdf[varname].attrs))
                          for varname in self.cdf] in expected)


class ReadColCDF(ColCDFTests):
    """Tests that read a column-major CDF, but do not modify it."""
    testfile = 'testc_ro.cdf'

    def __init__(self, *args, **kwargs):
        super(ReadColCDF, self).__init__(*args, **kwargs)
        #Unittest docs say 'the order in which the various test cases will be
        #run is determined by sorting the test function names with the built-in
        #cmp() function'
        testnames = [name for name in dir(self)
                     if name[0:4] == 'test' and callable(getattr(self,name))]
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
        self.assertEqual(self.cdf._majority().value, cdf.const.COLUMN_MAJOR.value)

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
            sliced = cdf._pycdf._Hyperslice(zvar, slices[i])
            actual = (sliced.dims, sliced.dimsizes, sliced.starts,
                      sliced.counts, sliced.intervals, sliced.degen,
                      sliced.rev)
            self.assertEqual(tuple(expected[i]), actual,
                             '\n' + str(tuple(expected[i])) + '!=\n' +
                             str(actual) + ' variable ' + i)

    def testColCreateReadArrays(self):
        """Create an array for reading from a slice"""
        slices = {'ATC': 1,
                  'PhysRecNo': slice(10, 2, -2),
                  'SpinNumbers': slice(2, None, 2),
                  'SectorRateScalersCounts': (slice(3, 6, None),
                                              slice(None, None, None),
                                              slice(None, None, None)),
                  'SpinRateScalersCounts': (Ellipsis, slice(-1, None, -1)),
                  } #Slice objects indexed by variable
        types = {'ATC': (ctypes.c_double * 2),
                 'PhysRecNo': (ctypes.c_int * 4),
                 'SpinNumbers': (ctypes.c_char * 2 * 8),
                 'SectorRateScalersCounts': (ctypes.c_float *
                                             3 * 32 * 9 * 100),
                 'SpinRateScalersCounts': (ctypes.c_float * 18 * 16 * 100),
                 } #constructors, build inside out
        for i in types:
            zvar = self.cdf[i]
            sliced = cdf._pycdf._Hyperslice(zvar, slices[i])
            actual = sliced.create_buffer().__class__
            self.assertEqual(types[i], actual)

    def testColSubscriptVariable(self):
        """Refer to an column-major array by subscript"""
        #NB: Should be in SAME order as row-major,
        #since converted in convert_array
        self.assertEqual([3, 25, 47],
                         self.cdf['PhysRecNo'][0:5:2])
        self.assertEqual([1094, 1083, 1072, 1061],
                         self.cdf['PhysRecNo'][-1:-5:-1])
        self.assertEqual(1.0,
                         self.cdf['SpinRateScalersCounts'][41, 2, 15])

    def testColSubscriptString(self):
        """Refer to a string array by subscript"""
        self.assertEqual(['0 ', '1 ', '2 ', '3 ', '4 ', '5 ', '6 ', '7 ',
                          '8 ', '9 ', '10', '11', '12', '13', '14', '15',
                          '16', '17'],
                         self.cdf['SpinNumbers'][:])

    def testColSubscriptIrregString(self):
        """Refer to a variable-length string array by subscript"""
        self.assertEqual(['H+', 'He+', 'He++', 'O<=+2', 'O>=+3', 'CN<=+2',
                          'H0', 'He0', 'CNO0', 'CN>=+3', 'Ne-Si', 'S-Ni',
                          '3He', 'D', 'Molecules', 'Others'],
                         self.cdf['RateScalerNames'][:])

    def testContains(self):
        """See if variable exists in CDF"""
        self.assertTrue('ATC' in self.cdf)
        self.assertFalse('notthere' in self.cdf)

    def testColPackBuffer(self):
        """Pack a buffer with data"""
        zvar = self.cdf['PhysRecNo']
        sliced = cdf._pycdf._Hyperslice(zvar, slice(0, None, 1))
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, list(range(100)))
        result = [buff[i] for i in range(100)]
        self.assertEqual(list(range(100)), result)
        
        sliced = cdf._pycdf._Hyperslice(zvar, slice(None, None, -1))
        buff = sliced.create_buffer()
        sliced.pack_buffer(buff, list(range(100)))
        result = [buff[i] for i in range(100)]
        self.assertEqual(list(reversed(range(100))), result)

        zvar = self.cdf['SectorRateScalersCounts']
        sliced = cdf._pycdf._Hyperslice(zvar, (0, slice(0, 3, 2),
                                              slice(3, None, -1), 1))
        buff = sliced.create_buffer()
        data = [[1, 2, 3, 4],
                [5, 6, 7, 8]]
        expected = [[4, 8], [3, 7], [2, 6], [1, 5]]
        sliced.pack_buffer(buff, data)
        for i in range(4):
            self.assertEqual(expected[i], list(buff[i]))


class ChangeCDF(CDFTests):
    """Tests that modify an existing CDF"""
    def __init__(self, *args):
        super(ChangeCDF, self).__init__(*args)
        
    def setUp(self):
        super(ChangeCDF, self).setUp()
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)
        self.cdf.readonly(False)

    def tearDown(self):
        self.cdf.close()
        del self.cdf
        os.remove(self.testfile)
        super(ChangeCDF, self).tearDown()

    def testDeletezVar(self):
        """Delete a zVar"""
        self.cdf['PhysRecNo']._delete()
        self.assertRaises(KeyError, self.cdf.__getitem__, 'PhysRecNo')

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

    def testWriteSubscripted(self):
        """Write data to a slice of a zVar"""
        expected = ['0 ', '1 ', '99', '3 ', '98', '5 ', '97', '7 ',
                    '8 ', '9 ']
        self.cdf['SpinNumbers'][2:7:2] = ['99', '98', '97']
        self.assertEqual(expected, self.cdf['SpinNumbers'][0:10])

        expected = self.cdf['SectorRateScalersCounts'][...]
        expected[4][5][5][8:3:-1] = [101.0, 102.0, 103.0, 104.0, 105.0]
        self.cdf['SectorRateScalersCounts'][4, 5, 5, 8:3:-1] = \
            [101.0, 102.0, 103.0, 104.0, 105.0]
        self.assertEqual(expected, self.cdf['SectorRateScalersCounts'][...])

    def testWriteExtend(self):
        """Write off the end of the variable"""
        additional = [2000 + i for i in range(20)]
        expected = self.cdf['PhysRecNo'][0:95] + additional
        self.cdf['PhysRecNo'][95:] = additional
        self.assertEqual(115, len(self.cdf['PhysRecNo']))
        self.assertEqual(expected, self.cdf['PhysRecNo'][:])

    def testInsertRecord(self):
        """Insert a record into the middle of a variable"""
        PhysRecNoData = self.cdf['PhysRecNo'][...]
        PhysRecNoData[5:6] = [-1, -2, -3, -4]
        self.cdf['PhysRecNo'][5:6] = [-1, -2, -3, -4]
        self.assertEqual(103, len(self.cdf['PhysRecNo']))
        self.assertEqual(PhysRecNoData, self.cdf['PhysRecNo'][...])

    def testWriteAndTruncate(self):
        """Write with insufficient data to fill all existing records"""
        expected = [-1 * i for i in range(20)] 
        self.cdf['PhysRecNo'][:] = expected
        self.assertEqual(expected, self.cdf['PhysRecNo'][:])

    def testWriteWrongSizeData(self):
        """Write with data sized or shaped differently from expected"""
        message = 'attempt to assign data of dimensions [3] ' + \
                  'to slice of dimensions [5]'
        try:
            self.cdf['SpinNumbers'][0:5] = [b'99', b'98', b'97']
        except ValueError:
            (t, v, tb) = sys.exc_info()
            self.assertEqual(message, str(v))
        else:
            self.fail('Should have raised ValueError: ' + message)

        message = 'attempt to assign data of dimensions [2, 6, 2] ' + \
                  'to slice of dimensions [3, 6, 2]'
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
        del PhysRecCopy[5]
        self.assertEqual(oldlen - 1, len(self.cdf['PhysRecNo']))
        self.assertEqual(PhysRecCopy[0:15], self.cdf['PhysRecNo'][0:15])

        oldlen = len(self.cdf['ATC'])
        ATCCopy = self.cdf['ATC'].copy()
        del self.cdf['ATC'][0::2]
        self.assertEqual(int(oldlen / 2), len(self.cdf['ATC']))
        self.assertEqual(ATCCopy[1::2], self.cdf['ATC'][...])

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
        del SectorRateScalersCountsCopy[-1:-5:-1]
        del self.cdf['SectorRateScalersCounts'][-1:-5:-1]
        self.assertEqual(oldlen - 4, len(self.cdf['SectorRateScalersCounts']))
        self.assertEqual(SectorRateScalersCountsCopy[...],
                         self.cdf['SectorRateScalersCounts'][...])

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

    def testRenameVar(self):
        """Rename a variable"""
        zvar = self.cdf['PhysRecNo']
        zvardata = zvar[...]
        zvar.rename('foobar')
        self.assertEqual(zvardata, self.cdf['foobar'][...])
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

    def testChangezEntry(self):
        """Write new or changed zEntry"""
        zvar = self.cdf['PhysRecNo']
        zvar.attrs['DEPEND_0'] = 'foobar'
        self.assertEqual('foobar', zvar.attrs['DEPEND_0'])
        self.assertEqual(const.CDF_CHAR.value,
                         cdf.zAttr(self.cdf,
                                   'DEPEND_0').entry_type(zvar._num()))

        zvar.attrs['FILLVAL'] = [0, 1]
        self.assertEqual([0,1], zvar.attrs['FILLVAL'])
        self.assertEqual(const.CDF_INT4.value,
                         cdf.zAttr(self.cdf,
                                   'FILLVAL').entry_type(zvar._num()))

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
                                   'NEW_ATTRIBUTE').entry_type(zvar._num()))
        
        zvar.attrs['NEW_ATTRIBUTE2'] = [1, 2]
        self.assertEqual([1, 2], zvar.attrs['NEW_ATTRIBUTE2'])
        self.assertEqual(const.CDF_INT4.value,
                         cdf.zAttr(self.cdf,
                                   'NEW_ATTRIBUTE2').entry_type(zvar._num()))

        zvar = self.cdf['SpinNumbers']
        zvar.attrs['NEW_ATTRIBUTE3'] = 1
        self.assertEqual(1, zvar.attrs['NEW_ATTRIBUTE3'])
        self.assertEqual(const.CDF_BYTE.value,
                         cdf.zAttr(self.cdf,
                                   'NEW_ATTRIBUTE3').entry_type(zvar._num()))

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

        self.cdf.attrs['Source_name'][0] = datetime.datetime(2009, 1, 1)
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
        expected = ['0 ', '1 ', '99', '3 ', '98', '5 ', '97', '7 ',
                    '8 ', '9 ']
        self.cdf['SpinNumbers'][2:7:2] = ['99', '98', '97']
        self.assertEqual(expected, self.cdf['SpinNumbers'][0:10])

        expected = self.cdf['SectorRateScalersCounts'][...]
        expected[4][5][5][8:3:-1] = [101.0, 102.0, 103.0, 104.0, 105.0]
        self.cdf['SectorRateScalersCounts'][4, 5, 5, 8:3:-1] = \
            [101.0, 102.0, 103.0, 104.0, 105.0]
        self.assertEqual(expected, self.cdf['SectorRateScalersCounts'][...])


if __name__ == '__main__':
    unittest.main()
