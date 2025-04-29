#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for SeaPy

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import datetime as dt
import matplotlib.axes
import unittest
import warnings
import numpy as np
import numpy.testing as ntest

import spacepy_testing
from spacepy import seapy
import spacepy.datamodel as dm

__all__ = ['SEATestsUniform', 'SEATests2dUniform',
           'SEATestsUniWithBad', 'SeaClassExceptions']


class SEATestsUniform(unittest.TestCase):
    """Tests of the sea method using uniform input"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATestsUniform, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        time = list(range(200))
        self.epochs = [20, 40, 60, 80, 100, 120, 140, 160, 180]
        with warnings.catch_warnings():
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea(self.unidata, time,
                                 self.epochs, verbose=False)
            self.obj.sea()

    def testMeanUniform(self):
        """Check superposed means on uniform input"""
        ntest.assert_array_equal(self.obj.semean,
                                 [self.testval]*(int(self.obj.window)*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input"""
        ntest.assert_array_equal(self.obj.semedian,
                                 [self.testval]*(int(self.obj.window)*2+1))

    def testMeanMedEquality(self):
        """For uniform input mean and median are same"""
        ntest.assert_array_equal(self.obj.semedian, self.obj.semean)

    def testDatetimeEquiv(self):
        """Test of equivalence of serial and datetime handling"""
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs,
                                window=window, delta=delta, verbose=False)
            compobj.sea()
        ntest.assert_array_equal(self.obj.semedian, compobj.semedian)
        ntest.assert_array_equal(self.obj.semean, compobj.semean)

    def testSeaLen(self):
        """len should return number of epochs"""
        self.assertEqual(len(self.obj), len(self.epochs))

    def testRandomEpochsNoArgs(self):
        """Random epochs should have requested number"""
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea = self.obj.random()
        self.assertEqual(len(newsea), len(self.epochs))

    def testRandomEpochsArgs(self):
        """Random epochs should have requested number"""
        n_req = 27
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea = self.obj.random(n_req)
        self.assertEqual(len(newsea), n_req)

    def testRandomType(self):
        """Random epochs should have requested number"""
        newsea = self.obj.random()
        self.assertEqual(type(newsea), type(self.obj))

    def testRandomBeforeSea(self):
        n_req = 27
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs,
                                window=window, delta=delta, verbose=False)
            newsea = compobj.random(n_req)
        self.assertEqual(len(newsea), n_req)

    def testRandomBoundType(self):
        """Test bound_type attriubte is set correctly """
        newsea = self.obj.random()
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea.sea(mad=True, quartiles=False)
        self.assertEqual(newsea.bound_type, 'mad')
        newsea = self.obj.random()
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            newsea.sea(ci=True, ci_quan='mean')
        self.assertEqual(newsea.bound_type, 'ci')

    def testSeaCIFunc(self):
        '''Use Sea object to test mean and function passing'''
        time = list(range(200))
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj1 = seapy.Sea(self.unidata, time, self.epochs,
                                 window=window, delta=delta, verbose=False)
            compobj1.sea(ci=95, ci_quan=np.mean)
            compobj2 = seapy.Sea(self.unidata, time, self.epochs,
                                 window=window, delta=delta, verbose=False)
            compobj2.sea(ci=95, ci_quan='mean')
        ntest.assert_allclose(compobj1.bound_low, compobj2.bound_low)

    def testSeaDict(self):
        '''Test seadict grouping'''
        namelist = ['O1', 'O2']
        sd = seapy.seadict([self.obj, self.obj], namelist)
        for key in sd.keys():
            self.assertTrue(key in namelist)

    def testSeaDictFail(self):
        '''Test seadict fails with bad inputs'''
        namelist = ['O1']
        namelist2 = ['O2']
        with self.assertRaises(ValueError):
            sd = seapy.seadict([self.obj, self.obj], namelist)
        with self.assertRaises(ValueError):
            sd = seapy.seadict({'O1': self.obj, 'O2': self.obj}, namelist2)

    def testSeaPlotShowFalse(self):
        '''Test that plot method (show=False) returns axes'''
        ax = self.obj.plot(show=False)
        self.assertTrue(isinstance(ax, matplotlib.axes.SubplotBase))


class SEATestsUniWithBad(unittest.TestCase):
    """Tests of sea method's handling of badvals"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATestsUniWithBad, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        # insert badvals
        for ind in range(30, 180, 16):
            self.unidata[ind] = -99
        time = list(range(200))
        self.epochs = [20, 40, 60, 80, 100, 120, 140, 160, 180]
        self.obj = seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        self.obj.sea(badval=-99)

    def testMeanUniform(self):
        """Check superposed means on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semean,
                                 [self.testval]*(int(self.obj.window)*2+1))

    def testMedianUniform(self):
        """Check superposed medians on uniform input with bad data"""
        ntest.assert_array_equal(self.obj.semedian,
                                 [self.testval]*(int(self.obj.window)*2+1))

    def testMeanMedEquality(self):
        """For uniform input mean and median are same with bad data"""
        ntest.assert_array_equal(self.obj.semedian, self.obj.semean)

    def testDatetimeEquiv(self):
        """Test of equivalence of serial and datetime handling"""
        sttime = dt.datetime(2010, 1, 1)
        time = [sttime + dt.timedelta(minutes=x) for x in range(200)]
        epochs = [sttime + dt.timedelta(minutes=x) for x in self.epochs]
        window = dt.timedelta(minutes=3)
        delta = dt.timedelta(minutes=1)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            compobj = seapy.Sea(self.unidata, time, epochs,
                                window=window, delta=delta, verbose=False)
        compobj.sea(badval=-99)

        ntest.assert_array_equal(self.obj.semedian, compobj.semedian)
        ntest.assert_array_equal(self.obj.semean, compobj.semean)


class SeaClassExceptions(unittest.TestCase):
    """Tests of the exception handling in Sea class"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SeaClassExceptions, self).setUp()

        self.testval = 5
        self.unidata = [self.testval]*200
        self.time = list(range(200))
        self.epochs = [20, 40, 60, 80, 100, 120, 140, 160, 180]

    def testRestoreEpochs(self):
        """Check that restoreepochs fails with no bad epochs"""
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea(self.unidata, self.time,
                                 self.epochs, verbose=False)
        re_fun = self.obj.restoreepochs
        self.assertRaises(AttributeError, re_fun)

    def testWarnNonContiguous(self):
        """Warn if time inputs appear non-contiguous/non-monotonic"""
        with spacepy_testing.assertWarns(self, message='Input time not monotonic;'
                                         ' results are unlikely to be valid.'):
            seapy.Sea(self.unidata,
                      self.time[::-1], self.epochs, verbose=False)
        time = [dt.datetime(2020, 1, 1) + dt.timedelta(seconds=60 * i)
                for i in range(202)]
        del time[2:4]  # Introduce a gap
        with spacepy_testing.assertWarns(self, message='Input time not contiguous;'
                                         ' results are unlikely to be valid.'):
            seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        time = [dt.datetime(2020, 1, 1) + dt.timedelta(seconds=60 * i)
                for i in range(200)]
        time[2] = dt.datetime(2019, 12, 31, 23, 59, 59)  # Go backwards
        with spacepy_testing.assertWarns(self, message='Input time not monotonic or contiguous;'
                                         ' results are unlikely to be valid.'):
            seapy.Sea(self.unidata, time, self.epochs, verbose=False)
        time = [dt.datetime(2020, 1, 1) - dt.timedelta(seconds=60 * i)
                for i in range(200)]
        with spacepy_testing.assertWarns(self, message='Input time not monotonic;'
                                         ' results are unlikely to be valid.'):
            seapy.Sea(self.unidata, time, self.epochs, verbose=False)

    def testWarnTooManyEpochs(self):
        """Warn if there are too many time epochs"""
        with spacepy_testing.assertWarns(self, message='Too many epochs;'
                                         ' results are unlikely to be valid.'):
            seapy.Sea(self.unidata, self.time, list(
                range(0, 202, 2)), verbose=False)


class SEATests2dUniform(unittest.TestCase):
    """Tests of the sea method using uniform input"""

    def setUp(self):
        """Setup block executed for each test in this class"""
        super(SEATests2dUniform, self).setUp()

        self.testval = 5
        self.unidata = np.ones([200, 200])
        self.unidata.fill(self.testval)
        time = list(range(200))
        self.epochs = [20, 40, 60, 80, 100, 120, 140, 160, 180]
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj = seapy.Sea2d(
                self.unidata, time, self.epochs, verbose=False)
            self.obj.sea()

    def testSea2dLen(self):
        '''Test that object length is as expected'''
        self.assertEqual(len(self.obj), len(self.epochs))

    def testSea2dMean(self):
        '''Mean of semean attr should be same as uniform test val'''
        self.assertEqual(np.mean(self.obj.semean), self.testval)

    def test2dRandomType(self):
        """Random object should have same type as parent"""
        newsea = self.obj.random()
        self.assertEqual(type(newsea), type(self.obj))


class SEAGetWindowsTest(unittest.TestCase):
    """Tests for timecube and get_windows helper"""

    def setUp(self):
        super(SEAGetWindowsTest, self).setUp()
        self.n_epochs = 50
        self.cadence = 1            # minutes
        self.window = 90           # minutes (±90)
        start = dt.datetime(2025, 1, 1)

        # build synthetic series
        self.times = [start + dt.timedelta(minutes=i*self.cadence)
                      for i in range(10_000)]
        self.data = np.sin(np.linspace(0, 20*np.pi, len(self.times)))
        self.epochs = self.times[100::20][:self.n_epochs]

        self.obj = seapy.Sea(self.data, self.times, self.epochs,
                             window=dt.timedelta(minutes=self.window),
                             delta=dt.timedelta(minutes=self.cadence),
                             verbose=False)
        self.obj.sea(storedata=True)

    def test_cubes_exist_and_match(self):
        """timecube exists and shares mask with datacube"""
        self.assertTrue(hasattr(self.obj, 'timecube'))
        self.assertTrue(hasattr(self.obj, 'datacube'))
        self.assertEqual(self.obj.timecube.shape, self.obj.datacube.shape)
        self.assertTrue(np.array_equal(self.obj.timecube.mask,
                                       self.obj.datacube.mask))

    def test_get_windows_compress_true(self):
        """get_windows default (compress=True) returns trimmed lists"""
        xs, ys = self.obj.get_windows()      # default compress=True
        self.assertEqual(len(xs), len(self.obj))
        self.assertTrue(all(isinstance(x, np.ndarray) for x in xs))
        # each trimmed window should be exactly 2*window+1 long
        exp_len = 2*self.window//self.cadence + 1
        self.assertTrue(all(len(x) == exp_len for x in xs))

    def test_get_windows_compress_false(self):
        xs, ys = self.obj.get_windows(compress=False)
        exp_len = 2*self.window // self.cadence + 1
        self.assertEqual(xs[0].shape, (exp_len,))


class SEARerunTests(unittest.TestCase):
    """Verify that .sea() can be called multiple times safely"""

    def setUp(self):
        super(SEARerunTests, self).setUp()
        # simple deterministic data
        self.data = np.arange(200, dtype=float)
        self.times = list(range(200))
        self.epochs = [50, 100, 150]

        # first analysis: half-window = 3 points
        self.obj = seapy.Sea(self.data, self.times, self.epochs,
                             window=3, delta=1, verbose=False)

        self.obj.sea(storedata=True)
        self.obj.sea(storedata=True)        # run it twice in a row

    def test_rerun_with_new_window(self):
        """second call to sea() rebuilds all arrays with new size"""
        old_len = len(self.obj.semedian)        # should be 7
        old_shape = self.obj.datacube.shape

        # change parameter and re-run
        self.obj.window = 5                       # new half-window
        self.obj.sea(storedata=True)              # re-run with new window

        new_len = len(self.obj.semedian)          # expect 11
        new_shape = self.obj.datacube.shape

        # the length must match 2*window+1
        self.assertEqual(new_len, 2 * 5 + 1)
        # and must differ from the first run
        self.assertNotEqual(new_len, old_len)

        # datacube columns must track the new length
        self.assertEqual(new_shape[1], new_len)
        # ensure datacube really changed
        self.assertNotEqual(new_shape, old_shape)


class SEASignifTest(unittest.TestCase):
    """Tests for the sea_signif function with results object return"""

    def setUp(self):
        """Setup test objects for significance testing"""
        super(SEASignifTest, self).setUp()

        # Create two different datasets to enable statistical testing
        self.testval1 = 5
        self.testval2 = 5
        self.unidata1 = [self.testval1]*200
        self.unidata2 = [self.testval2]*200

        # Add significant variation to second dataset for U test to work
        for i in range(40, 160, 10):
            self.unidata2[i] = self.testval2 + 3

        self.time = list(range(200))
        self.epochs = [20, 40, 60, 80, 100, 120, 140, 160, 180]

        # Create two Sea objects for comparison
        with warnings.catch_warnings():
            warnings.filterwarnings('always', 'Window size changed .*',
                                    UserWarning, '^spacepy\\.seapy$')
            self.obj1 = seapy.Sea(self.unidata1, self.time,
                                  self.epochs, verbose=False)
            self.obj2 = seapy.Sea(self.unidata2, self.time,
                                  self.epochs, verbose=False)

        # Run the sea analysis with storedata=True for significance testing
        self.obj1.sea(storedata=True)
        self.obj2.sea(storedata=True)

    def testResultsReturned(self):
        """Test that sea_signif always returns a results object"""
        # Test with show=False first (doesn't call plt.show)
        results = seapy.sea_signif(self.obj1, self.obj2, show=False)
        self.assertIsInstance(results, dm.SpaceData)
        # results = seapy.sea_signif(self.obj1, self.obj2, show=True)
        # self.assertIsInstance(results, dm.SpaceData)

    def testResultsAttributes(self):
        """Test that results object has the expected attributes"""
        results = seapy.sea_signif(self.obj1, self.obj2, show=False)

        # Test for expected keys
        self.assertIn('x_values', results)
        self.assertIn('stat_values', results)
        self.assertIn('prob_values', results)

        # Test for expected attributes
        self.assertIn('test_name', results.attrs)
        self.assertIn('obj1_name', results.attrs)
        self.assertIn('obj2_name', results.attrs)

        # Test values
        self.assertEqual(results.attrs['test_name'], 'KS')
        self.assertEqual(results.attrs['obj1_name'], 'Sea')
        self.assertEqual(results.attrs['obj2_name'], 'Sea')

    def testResultsLength(self):
        """Test that arrays in results have the expected length"""
        results = seapy.sea_signif(self.obj1, self.obj2, show=False)

        # All arrays should match the length of the window
        expected_length = int(self.obj1.window)*2 + 1
        self.assertEqual(len(results['x_values']), expected_length)
        self.assertEqual(len(results['stat_values']), expected_length)
        self.assertEqual(len(results['prob_values']), expected_length)

    def testWindowAlignment(self):
        """Test that x_values match the window range"""
        results = seapy.sea_signif(self.obj1, self.obj2, show=False)

        # Check that x_values match obj1.x
        ntest.assert_array_equal(results['x_values'], self.obj1.x)

    def testAlternateTest(self):
        """Test results with Mann-Whitney U test"""
        try:
            results = seapy.sea_signif(self.obj1, self.obj2,
                                       test='U', show=False)

            # Test should be 'U'
            self.assertEqual(results.attrs['test_name'], 'U')
            self.assertEqual(results['stat_values'].attrs['test'], 'U')
        except ValueError as e:
            # Skip test if datasets are too similar for mannwhitneyu
            if "All numbers are identical in mannwhitneyu" in str(e):
                self.skipTest("Datasets too similar for Mann-Whitney U test")
            else:
                raise


if __name__ == '__main__':
    unittest.main()
