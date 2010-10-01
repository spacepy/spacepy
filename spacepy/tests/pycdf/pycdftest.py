#!/usr/bin/env python

"""Unit test suite for pycdf"""

import ctypes
import datetime
import hashlib
import os, os.path
import shutil
import sys
import unittest

import pycdf as cdf


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
        self.cdf._readonly(False)
        
    def tearDown(self):
        os.remove('new.cdf')

    def testCreateVarDefault(self):
        """Create a variable with default options in an empty CDF"""
        new = self.cdf._new_var('new_variable', cdf.const.CDF_UINT1)
        self.assertTrue(type(new).__name__ == 'Var')


class CDFTests(unittest.TestCase):
    """Tests that involve an existing CDF, read or write"""
    def __init__(self, *args):
        self.testmaster = 'po_l1_cam_test.cdf'
        self.testfile = 'test.cdf'
        self.expected_digest = '4e78627e257ad2743387d3a7653b5e30'
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
    def __init__(self, *args):
        self.testmaster = 'po_l1_cam_testc.cdf'
        self.testfile = 'testc.cdf'
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
                    'MajorNumbers', 'MeanCharge']
        with cdf.CDF(self.testfile) as f:
            names = list(f.keys())
        self.assertEqual(expected, names)
        self.assertRaises(cdf.CDFError, f.close)


class ReadCDF(CDFTests):
    """Tests that read an existing CDF, but do not modify it."""
    def __init__(self, *args):
        super(ReadCDF, self).__init__(*args)
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)

    def __del__(self):
        del self.cdf
        os.remove(self.testfile)
    
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
        expectedNames = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers', 'RateScalerNames', 'SectorRateScalerNames', 'SectorRateScalersCounts', 'SectorRateScalersCountsSigma', 'SpinRateScalersCounts', 'SpinRateScalersCountsSigma', 'MajorNumbers', 'MeanCharge']
        names = [zVar.name() for zVar in self.cdf.values()]
        self.assertEqual(names, expectedNames)

    def testGetAllVarNames(self):
        """Getting a list of zVar names"""
        expectedNames = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers', 'RateScalerNames', 'SectorRateScalerNames', 'SectorRateScalersCounts', 'SectorRateScalersCountsSigma', 'SpinRateScalersCounts', 'SpinRateScalersCountsSigma', 'MajorNumbers', 'MeanCharge']
        names = list(self.cdf.keys())
        self.assertEqual(names, expectedNames)

    def testGetVarNum(self):
        self.assertEqual(0, self.cdf['ATC']._num())

    def testCDFIterator(self):
        expected = ['ATC', 'PhysRecNo', 'SpinNumbers', 'SectorNumbers',
                    'RateScalerNames', 'SectorRateScalerNames',
                    'SectorRateScalersCounts', 'SectorRateScalersCountsSigma',
                    'SpinRateScalersCounts', 'SpinRateScalersCountsSigma',
                    'MajorNumbers', 'MeanCharge']
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
                 'SpinNumbers': (ctypes.c_char * 16), #degenerate
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
        data = self.cdf.all_data()
        expected = ['ATC', 'MajorNumbers', 'MeanCharge',
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
        self.assertEqual(12, result)

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


class ReadColCDF(ColCDFTests):
    """Tests that read a column-major CDF, but do not modify it."""
    def setUp(self):
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)

    def tearDown(self):
        del self.cdf
        os.remove(self.testfile)
    
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
                 'SpinNumbers': (ctypes.c_char * 16), #degenerate
                 'SectorRateScalersCounts': (ctypes.c_float *
                                             9 * 32 * 3 * 100),
                 'SpinRateScalersCounts': (ctypes.c_float * 16 * 18 * 100),
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


class ChangeCDF(CDFTests):
    """Tests that modify an existing CDF"""
    def __init__(self, *args):
        super(ChangeCDF, self).__init__(*args)
        
    def setUp(self):
        shutil.copy(self.testmaster, self.testfile)
        self.cdf = cdf.CDF(self.testfile)
        self.cdf._readonly(False)

    def tearDown(self):
        del self.cdf
        os.remove(self.testfile)

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

    def testDeleteCount(self):
        """Verify correct number of records removed after deletion from zVar"""
        start_count = len(self.cdf['MeanCharge'])
        self.cdf['MeanCharge']._del_recs(50, 20)
        self.assertEqual(len(self.cdf['MeanCharge']), start_count - 20)

    def testReadonlySettable(self):
        """Readonly mode should prevent changes"""
        self.cdf._readonly(True)
        self.assertTrue(self.cdf._readonly())
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
        self.cdf._readonly(True)
        self.assertTrue(self.cdf._readonly())
        self.cdf._readonly(False)
        try:
            self.cdf['PhysRecNo']._delete()
        except:
            (type, val, traceback) = sys.exc_info()
            self.fail('Raised exception ' + str(val))


if __name__ == '__main__':
    unittest.main()
