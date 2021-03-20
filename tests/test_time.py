# -*- coding: utf-8 -*-

"""
Test suite for time module

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import datetime
import itertools

try:
    import StringIO
except:
    import io as StringIO
try:
    from itertools import izip as zip
except ImportError:
    pass  # just use system zip. In python3 itertools.izip is just python zip
import unittest
import pickle
import time
import sys
import warnings

try:
    import astropy.time
    HAVE_ASTROPY = True
except: # Don't bring down whole test suite
    HAVE_ASTROPY = False
import numpy

import spacepy_testing
import spacepy
import spacepy.pycdf
import spacepy.time as t

__all__ = ['TimeFunctionTests', 'TimeClassTests']


class TimeFunctionTests(unittest.TestCase):

    def test_doy2dateconvert(self):
        """doy2date should return a known value for known input"""
        inval = [(2000, 1),
                 (2001, 34),
                 (2006, 34),
                 (2008, 60),
                 (2008, 366),
                 ([2008], [366])]
        real_ans = [(1, 1),
                    (2, 3),
                    (2, 3),
                    (2, 29),
                    (12, 31),
                    ([12], [31])]
        for i, val in enumerate(inval):
            ans = t.doy2date(*val)
            ans2 = t.doy2date(*val, dtobj=True)
            self.assertEqual(real_ans[i], ans)
            try:
                self.assertEqual(real_ans[i],
                                 ([ans2[0].month], [ans2[0].day]))
            except TypeError:
                self.assertEqual(real_ans[i], (ans2.month, ans2.day))
        self.assertRaises(ValueError, t.doy2date, (2000, 2000, 2000), (5, 4))
        numpy.testing.assert_array_equal(t.doy2date((2000, 2000, 2000), (5, 4, 3.3), flAns=True),
                                         [[1, 1, 1],[5, 4, 3]])

    def test_doy2datefail(self):
        '''doy2date should fail for bad input'''
        inval = ([[2007], [0.5]],
                 [2007, 0.5])
        for val in inval:
            func = lambda: t.doy2date(*val)
            self.assertRaises(ValueError, func)

    def test_doy2datefloat(self):
        '''doy2date should work with floats'''
        ans = (datetime.datetime(2000, 1, 2, 2, 58, 33, 600000),
               datetime.datetime(2000, 1, 2, 0, 0))
        inval = [(2000, 2.124, True, True),
                 (2000, 2, True, True)]
        for i, val in enumerate(ans):
            self.assertEqual(val, t.doy2date(*inval[i]))

    def test_tickrange(self):
        """tickrange should return a known value for known input"""
        inval = (('2002-02-01T00:00:00', '2002-02-04T00:00:00', 1),
                 ('2002-02-01T00:00:00', '2002-02-04T00:00:00', 0.5))
        strarray_dtype = numpy.array('x' * 19).dtype
        real_ans = (numpy.array(['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00',
                                 '2002-02-04T00:00:00'], dtype=strarray_dtype),
                    numpy.array(['2002-02-01T00:00:00', '2002-02-01T12:00:00', '2002-02-02T00:00:00',
                                 '2002-02-02T12:00:00', '2002-02-03T00:00:00', '2002-02-03T12:00:00',
                                 '2002-02-04T00:00:00'], dtype=strarray_dtype))
        for i, val in enumerate(inval):
            ans = t.tickrange(*val)
            numpy.testing.assert_equal(real_ans[i], ans.ISO)

    def test_tickrange2(self):
        """tickrange should return a known value for known input (timedelta)"""
        inval = (('2002-02-01T00:00:00', '2002-02-04T00:00:00', datetime.timedelta(days=1)),
                 ('2002-02-01T00:00:00', '2002-02-04T00:00:00', datetime.timedelta(hours=12)))
        strarray_dtype = numpy.array('x' * 19).dtype
        real_ans = (numpy.array(['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00',
                                 '2002-02-04T00:00:00'], dtype=strarray_dtype),
                    numpy.array(['2002-02-01T00:00:00', '2002-02-01T12:00:00', '2002-02-02T00:00:00',
                                 '2002-02-02T12:00:00', '2002-02-03T00:00:00', '2002-02-03T12:00:00',
                                 '2002-02-04T00:00:00'], dtype=strarray_dtype))
        for i, val in enumerate(inval):
            ans = t.tickrange(*val)
            numpy.testing.assert_equal(real_ans[i], ans.ISO)

    def test_tickrange3(self):
        """tickrange should return a known value for known input (test for bug #64 on SF tracker)"""
        inval = ('2009-01-01', '2010-12-31 23:00', datetime.timedelta(hours=1))
        real_ans = datetime.datetime(2010, 12, 31, 23)
        ans = t.tickrange(*inval)
        numpy.testing.assert_equal(real_ans, ans.UTC[-1])

    def test_sec2hms(self):
        """sec2hms should return a known value for known input"""
        inval = ((30, False, False),
                 (86401, False, False),
                 (86401, False, True),
                 (30.3, True, False),
                 (3599, False, False))
        real_ans = ([0, 0, 30],
                    [24, 0, 1],
                    [0, 0, 1],
                    [0, 0, 30],
                    [0, 59, 59])
        with spacepy_testing.assertWarns(
                self, 'always',
                r'Number of seconds > seconds in day\. Try days keyword\.$',
                UserWarning, r'spacepy\.time$'):
            for i, val in enumerate(inval):
                ans = t.sec2hms(*val)
                self.assertEqual(real_ans[i], ans)
        self.assertEqual(t.sec2hms(12, False, False, True), datetime.timedelta(seconds=12))

    def test_no_tzinfo(self):
        """no_tzinfo should have known output"""
        dt = datetime.datetime(2000, 1, 1, tzinfo=datetime.tzinfo('MST'))
        self.assertEqual(dt.replace(tzinfo=None), t.no_tzinfo(dt))
        ans = [datetime.datetime(2000, 1, 1)] * 10
        self.assertEqual(ans, t.no_tzinfo([dt] * 10))

    def test_no_tzinfo_attrs(self):
        """no_tzinfo should preserve attributes"""
        dt = spacepy.dmarray([
            datetime.datetime(2001, 1, 1),
            datetime.datetime(2001, 1, 2)
        ], attrs={'dtype': 'UTC'})
        dt2 = t.no_tzinfo(dt)
        self.assertEqual(
            {'dtype': 'UTC'}, dt2.attrs)

    def test_leapyear(self):
        """leapyear should give known output for known input"""
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
                self.assertTrue(t.leapyear(i))
            else:
                self.assertFalse(t.leapyear(i))

        data = (1993 + numpy.array(range(10)), 1900, [1993 + val for val in range(10)])
        real_ans = (
            numpy.array([365, 365, 365, 366, 365, 365, 365, 366, 365, 365]),
            365,
            numpy.array([365, 365, 365, 366, 365, 365, 365, 366, 365, 365]))
        for i, val in enumerate(real_ans):
            if i == 0:
                numpy.testing.assert_array_equal(val.tolist(), t.leapyear(data[i], True))
            else:
                numpy.testing.assert_array_equal(val, t.leapyear(data[i], True))
        real_ans = (numpy.array([False, False, False, True, False, False, False, True, False, False]),
                    False,
                    [False, False, False, True, False, False, False, True, False, False])
        for i, val in enumerate(real_ans):
            if i == 0:
                numpy.testing.assert_array_equal(val.tolist(), t.leapyear(data[i], False))
            else:
                numpy.testing.assert_array_equal(val, t.leapyear(data[i], False))

    def test_randomDate(self):
        """randomDate should give known result"""
        dt1 = datetime.datetime(2000, 1, 1)
        dt2 = datetime.datetime(2000, 2, 1)
        numpy.random.seed(8675309)
        ans = numpy.array([datetime.datetime(2000, 1, 26, 4, 28, 10, 500070),
                           datetime.datetime(2000, 1, 24, 6, 46, 39, 156905),
                           datetime.datetime(2000, 1, 12, 1, 52, 50, 481431),
                           datetime.datetime(2000, 1, 7, 6, 30, 26, 331312),
                           datetime.datetime(2000, 1, 13, 16, 17, 48, 619577)])
        ntests = len(ans)
        res = t.randomDate(dt1, dt2, ntests, sorted=False)
        #results are different at microsecond level between Python2 and Python3
        #one likely cause is the difference in behavior of round() between versions
        #so, we'll round off all the microseconds fields here
        for ii in range(ntests):
            ans[ii] = ans[ii].replace(microsecond=100*(ans[ii].microsecond//100))
            res[ii] = res[ii].replace(microsecond=100*(res[ii].microsecond//100))
        #TODO: improve testing for randomDate
        numpy.testing.assert_array_equal(ans, res)
        # check the exception
        class junktzinfo(datetime.tzinfo):
            def utcoffset(self, dt):
                return 0
            def dst(self, dt):
                return 0
            def tzname(self, dt):
                return 'junk'
        dt11 = dt1.replace(tzinfo=junktzinfo())
        self.assertRaises(ValueError, t.randomDate, dt11, dt2)
        ans.sort()

        numpy.random.seed(8675309)
        res = t.randomDate(dt1, dt2, ntests, sorted=True)
        for ii in range(ntests):
            res[ii] = res[ii].replace(microsecond=100*(res[ii].microsecond//100))
        numpy.testing.assert_array_equal(ans, res)

        numpy.random.seed(8675309)
        res = t.randomDate(dt1, dt2, ntests, sorted=True, tzinfo='MDT')
        for ii in range(ntests):
            res[ii] = res[ii].replace(microsecond=100*(res[ii].microsecond//100))
        numpy.testing.assert_array_equal(ans, res)

    def test_extract_YYYYMMDD(self):
        """extract_YYYYMMDD() should give known results"""
        filenames = ['rbspa_rel02_ect-hope-PA-L3_20130906_v4.0.0.cdf',
                     'rbspa_def_MagEphem_OP77Q_20150202_v1.0.0.h5',
                     '20150204_firebird-2-fu3_T89D_MagEphem.h5',
                     '20150202_firebird-2-fu3_T89D_MagEphem.h5',
                     'I_am_a_file_with_no_date.h5']
        ans = [datetime.datetime(2013, 9, 6),
               datetime.datetime(2015, 2, 2),
               datetime.datetime(2015, 2, 4),
               datetime.datetime(2015, 2, 2),
               None]
        for tst, ans in zip(filenames, ans):
            self.assertEqual(ans, t.extract_YYYYMMDD(tst))

    def test_valid_YYYYMMDD(self):
        """valid_YYYYMMDD() should give known results"""
        filenames = ['rbspa_rel02_ect-hope-PA-L3_20130906_v4.0.0.cdf',
                     'rbspa_def_MagEphem_OP77Q_20150202_v1.0.0.h5',
                     '20150204_firebird-2-fu3_T89D_MagEphem.h5',
                     '20150202_firebird-2-fu3_T89D_MagEphem.h5',
                     'I_am_a_file_with_no_date.h5']
        ans = [True, True, True, True, False]
        for tst, ans in zip(filenames, ans):
            self.assertEqual(ans, t.valid_YYYYMMDD(tst))

    def test_dtstr2iso(self):
        """convert datetime string to ISO + UTC"""
        inputs = [
            '2001-01-01T23:59:59',
            '2001-01-02',
            '2005-12-31T23:59:60',
            '2005-12-31T23:59:60.123',
            ]
        expectediso = [
            '2001-01-01T23:59:59',
            '2001-01-02T00:00:00',
            '2005-12-31T23:59:60',
            '2005-12-31T23:59:60',
            ]
        expectedUTC = [
            (2001, 1, 1, 23, 59, 59),
            (2001, 1, 2),
            (2005, 12, 31, 23, 59, 59, 999999),
            (2005, 12, 31, 23, 59, 59, 999999),
            ]
        expectedoffset = [
            0.,
            0,
            1,
            123001,
        ]
        expectedUTC = [datetime.datetime(*e) for e in expectedUTC]
        actualiso, actualUTC, actualoffset = t.dtstr2iso(inputs)
        # Make sure not overwritten.
        self.assertEqual('2005-12-31T23:59:60', inputs[-2])
        numpy.testing.assert_equal(expectedUTC, actualUTC)
        numpy.testing.assert_equal(expectediso, actualiso)
        numpy.testing.assert_equal(expectedoffset, actualoffset)
        self.assertEqual(
            ('2005-12-31T23:59:60.123000',
             datetime.datetime(2005, 12, 31, 23, 59, 59, 999999),
             123001),
            t.dtstr2iso('2005-12-31T23:59:60.123',
                        fmt='%Y-%m-%dT%H:%M:%S.%f'))
        # Make inputs numpy (more likely to be overwritten).
        inputs = numpy.array(inputs)
        actualiso, actualUTC, actualoffset = t.dtstr2iso(inputs)
        # Make sure not overwritten.
        self.assertEqual('2005-12-31T23:59:60', inputs[-2])

    def test_dtstr2isomicrosecond(self):
        """Test dtstr2iso with microseconds input"""
        inputs = [
            '2001-01-02T00:01:02.000600',
            ]
        expectediso = [
            '2001-01-02T00:01:02',
            ]
        expectedUTC = [
            (2001, 1, 2, 0, 1, 2, 600),
            ]
        expectedoffset = [
            0.,
        ]
        expectedUTC = [datetime.datetime(*e) for e in expectedUTC]
        actualiso, actualUTC, actualoffset = t.dtstr2iso(inputs)
        numpy.testing.assert_equal(expectedUTC, actualUTC)
        numpy.testing.assert_equal(expectediso, actualiso)
        numpy.testing.assert_equal(expectedoffset, actualoffset)
        actualiso, actualUTC, actualoffset = t.dtstr2iso(inputs[0])
        numpy.testing.assert_equal(expectedUTC[0], actualUTC)
        numpy.testing.assert_equal(expectediso[0], actualiso)
        numpy.testing.assert_equal(expectedoffset[0], actualoffset)

    def test_dtstr2isoearly(self):
        """convert datetime string to ISO + UTC before 1900"""
        inputs = ['1890-01-1', '1858-11-18',
                  '1858-11-17', '1858-11-16',
                  '1066-10-14',]
        expectedUTC = [(1890, 1, 1), (1858, 11, 18),
                       (1858, 11, 17), (1858, 11, 16),
                       (1066, 10, 14)]
        expectedUTC = [datetime.datetime(*e) for e in expectedUTC]
        expectediso = ['1890-01-01T00:00:00', '1858-11-18T00:00:00',
                       '1858-11-17T00:00:00', '1858-11-16T00:00:00',
                       '1066-10-14T00:00:00',]
        expectedoffset = [0, 0,
                          0, 0,
                          0]
        actualiso, actualUTC, actualoffset = t.dtstr2iso(inputs)
        numpy.testing.assert_equal(expectedUTC, actualUTC)
        numpy.testing.assert_equal(expectediso, actualiso)
        numpy.testing.assert_equal(expectedoffset, actualoffset)

    def test_dtstr2isobadleap(self):
        """Convert a string with bad leap second to UTC"""
        inputs = ['2008-12-31T23:59:60.123', '2009-12-31T23:59:60.100']
        with self.assertRaises(ValueError) as cm:
            t.dtstr2iso(inputs)
        self.assertEqual('2009-12-31T23:59:60.100 is not a valid leapsecond.',
                         str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            t.dtstr2iso(inputs[-1])
        self.assertEqual('2009-12-31T23:59:60.100 is not a valid leapsecond.',
                         str(cm.exception))
        more_bad = [
            # Potentially a valid leapsecond, but after the last known one,
            # so not valid for now (and I'll be dead before it is valid.)
            '2999-12-31T23:59:60',
            '2999-12-31T23:59:60.5']
        for b in more_bad:
            with self.assertRaises(ValueError) as cm:
                t.dtstr2iso(b)
            self.assertEqual('{} is not a valid leapsecond.'.format(b),
                             str(cm.exception))
            # Do same thing in an array
            with self.assertRaises(ValueError) as cm:
                t.dtstr2iso([b])
            self.assertEqual('{} is not a valid leapsecond.'.format(b),
                             str(cm.exception))

    def test_days1958(self):
        """Test fractional days since 1958"""
        self.assertRaises(ValueError, t._days1958, [0.], leaps='foo')
        inputs = [
            -31514400.,# 1957-01-01 06:00:00
            0.,        # 1958-01-01 00:00:00
            43200,     # 1958-01-01 12:00:00
            94694400,  # 1960-12-31 23:59:59
            94694401,  # 1960-12-31 23:59:60
            94694402,  # 1961-01-01 00:00:00
            126230400, # 1961-12-31 23:59:58
            126230401, # 1961-12-31 23:59:59
            126230402, # 1962-01-01 00:00:00
            126316802, # 1962-01-02 00:00:00
            1609459232,# 2008-12-31 23:59:59
            1609459233,# 2008-12-31 23:59:60
            1609459234,# 2009-01-01 00:00:00
            1609459235,# 2009-01-01 00:00:01
            ]

        expected = [
            -365.25,
            -0.5,
            0,
            1095 + 43199. / 86401,
            1095 + 43200. / 86401,
            1095 + 43201. / 86401,
            1460 + 43198. / 86400,
            1460 + 43199. / 86400,
            1460.5,
            1461.5,
            18627. + 43199. / 86401,
            18627. + 43200. / 86401,
            18627. + 43201. / 86401,
            18627. + 43202. / 86401,
            ]
        actual = t._days1958(inputs, leaps='rubber')
        numpy.testing.assert_equal(expected, actual)

        expected = [
            -365.25,
            -0.5,
            0,
            1095 + 43199. / 86400,
            1095 + 43199.999999 / 86400,
            1095.5,
            1460 + 43198. / 86400,
            1460 + 43199. / 86400,
            1460.5,
            1461.5,
            18627. + 43199. / 86400,
            18627. + 43199.999999 / 86400,
            18627.5,
            18627. + 43201. / 86400,
            ]
        actual = t._days1958(inputs, leaps='drop')
        numpy.testing.assert_equal(expected, actual)

        expected = [
            -365.25,
            -0.5,
            0,
            1095.5,
            1095.5 + 1. / 86400,
            1095.5 + 2. / 86400,
            1460.5,
            1460.5 + 1. / 86400,
            1460.5 + 2. / 86400,
            1461.5 + 2. / 86400,
            18627.5 + 32. / 86400,
            18627.5 + 33. / 86400,
            18627.5 + 34. / 86400,
            18627.5 + 35. / 86400,
            ]
        actual = t._days1958(inputs, leaps='continuous')
        numpy.testing.assert_equal(expected, actual)

    def test_days1958_scalar(self):
        """Fractional days since 1958, scalar input"""
        self.assertEqual(
            18627. + 43200. / 86401,
            t._days1958(1609459233, leaps='rubber'))
        self.assertEqual(
            18627. + 43199.999999 / 86400,
            t._days1958(1609459233, leaps='drop'))
        self.assertEqual(
            18627.5 + 33. / 86400,
            t._days1958(1609459233, leaps='continuous'))
        for handler in ('rubber', 'drop', 'continuous'):
            self.assertEqual(
                (), t._days1958(1609459233, leaps=handler).shape)

    def test_days1958_midnight(self):
        """Test fractional days since 1958, using midnight start time"""
        inputs = [
            -31514400.,# 1957-01-01 06:00:00
            0.,        # 1958-01-01 00:00:00
            43200,     # 1958-01-01 12:00:00
            94694400,  # 1960-12-31 23:59:59
            94694401,  # 1960-12-31 23:59:60
            94694402,  # 1961-01-01 00:00:00
            126230400, # 1961-12-31 23:59:58
            126230401, # 1961-12-31 23:59:59
            126230402, # 1962-01-01 00:00:00
            126316802, # 1962-01-02 00:00:00
            1609459232,# 2008-12-31 23:59:59
            1609459233,# 2008-12-31 23:59:60
            1609459234,# 2009-01-01 00:00:00
            1609459235,# 2009-01-01 00:00:01
            ]

        expected = [
            -364.75,
            0.,
            0.5,
            1095 + 86399. / 86401,
            1095 + 86400. / 86401,
            1095 + 86401. / 86401,
            1460 + 86398. / 86400,
            1460 + 86399. / 86400,
            1461,
            1462,
            18627. + 86399. / 86401,
            18627. + 86400. / 86401,
            18628.,
            18628. + 1. / 86400,
            ]
        actual = t._days1958(inputs, leaps='rubber', midnight=True)
        numpy.testing.assert_equal(expected, actual)
        self.assertEqual(
            expected[7],
            t._days1958(inputs[7], leaps='rubber', midnight=True))

        expected = [
            -364.75,
            0.,
            0.5,
            1095 + 86399. / 86400,
            1095 + 86399.999999 / 86400,
            1096,
            1460 + 86398. / 86400,
            1460 + 86399. / 86400,
            1461,
            1462,
            18627. + 86399. / 86400,
            18627. + 86399.999999 / 86400,
            18628.,
            18628. + 1. / 86400,
            ]
        actual = t._days1958(inputs, leaps='drop', midnight=True)
        numpy.testing.assert_equal(expected, actual)
        self.assertEqual(
            expected[7],
            t._days1958(inputs[7], leaps='drop', midnight=True))

        expected = [
            -364.75,
            0.,
            0.5,
            1096,
            1096 + 1. / 86400,
            1096 + 2. / 86400,
            1461,
            1461 + 1. / 86400,
            1461 + 2. / 86400,
            1462 + 2. / 86400,
            18628. + 32. / 86400,
            18628. + 33. / 86400,
            18628. + 34. / 86400,
            18628. + 35. / 86400,
            ]
        actual = t._days1958(inputs, leaps='continuous', midnight=True)
        numpy.testing.assert_equal(expected, actual)
        self.assertEqual(
            expected[6],
            t._days1958(inputs[6], leaps='continuous', midnight=True))

    def test_days1958toTAI(self):
        """Test fractional days since 1958 to TAI"""
        self.assertRaises(ValueError, t._days1958totai, [0.], leaps='foo')
        expected = [
            -31514400.,# 1957-01-01 06:00:00
            0.,        # 1958-01-01 00:00:00
            43200,     # 1958-01-01 12:00:00
            94694400,  # 1960-12-31 23:59:59
            94694401,  # 1960-12-31 23:59:60
            94694402,  # 1961-01-01 00:00:00
            126230400, # 1961-12-31 23:59:58
            126230401, # 1961-12-31 23:59:59
            126230402, # 1962-01-01 00:00:00
            126316802, # 1962-01-02 00:00:00
            1609459232,# 2008-12-31 23:59:59
            1609459233,# 2008-12-31 23:59:60
            1609459234,# 2009-01-01 00:00:00
            1609459235,# 2009-01-01 00:00:01
            ]

        inputs = [
            -365.25,
            -0.5,
            0,
            1095 + 43199. / 86401,
            1095 + 43200. / 86401,
            1095 + 43201. / 86401,
            1460 + 43198. / 86400,
            1460 + 43199. / 86400,
            1460.5,
            1461.5,
            18627. + 43199. / 86401,
            18627. + 43200. / 86401,
            18627. + 43201. / 86401,
            18627. + 43202. / 86401,
            ]
        actual = t._days1958totai(inputs, leaps='rubber')
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)

        inputs = [
            -365.25,
            -0.5,
            0,
            1095 + 43199. / 86400,
            1095 + 43199.999999 / 86400,
            1095.5,
            1460 + 43198. / 86400,
            1460 + 43199. / 86400,
            1460.5,
            1461.5,
            18627. + 43199. / 86400,
            18627. + 43199.999999 / 86400,
            18627.5,
            18627. + 43201. / 86400,
            ]
        actual = t._days1958totai(inputs, leaps='drop')
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)

        inputs = [
            -365.25,
            -0.5,
            0,
            1095.5,
            1095.5 + 1. / 86400,
            1095.5 + 2. / 86400,
            1460.5,
            1460.5 + 1. / 86400,
            1460.5 + 2. / 86400,
            1461.5 + 2. / 86400,
            18627.5 + 32. / 86400,
            18627.5 + 33. / 86400,
            18627.5 + 34. / 86400,
            18627.5 + 35. / 86400,
            ]
        actual = t._days1958totai(inputs, leaps='continuous')
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)

    def test_days1958toTAI_scalar(self):
        """Fractional days since 1958 to TAI, scalar input"""
        self.assertAlmostEqual(
            1609459233,
            t._days1958totai(18627. + 43200. / 86401, leaps='rubber'),
            places=6)
        self.assertAlmostEqual(
            1609459232.999999,
            t._days1958totai(18627. + 43199.999999 / 86400, leaps='drop'),
            places=6)
        self.assertAlmostEqual(
            1609459233,
            t._days1958totai(18627.5 + 33. / 86400, leaps='continuous'),
            places=6)
        for handler in ('rubber', 'drop', 'continuous'):
            self.assertEqual(
                (), t._days1958totai(18627.5, leaps=handler).shape)

    def test_days1958toTAI_midnight(self):
        """Test fractional days since 1958, using midnight start time"""
        expected = [
            -31514400.,# 1957-01-01 06:00:00
            0.,        # 1958-01-01 00:00:00
            43200,     # 1958-01-01 12:00:00
            94694400,  # 1960-12-31 23:59:59
            94694401,  # 1960-12-31 23:59:60
            94694402,  # 1961-01-01 00:00:00
            126230400, # 1961-12-31 23:59:58
            126230401, # 1961-12-31 23:59:59
            126230402, # 1962-01-01 00:00:00
            126316802, # 1962-01-02 00:00:00
            1609459232,# 2008-12-31 23:59:59
            1609459233,# 2008-12-31 23:59:60
            1609459234,# 2009-01-01 00:00:00
            1609459235,# 2009-01-01 00:00:01
            ]

        inputs = [
           -364.75,
            0.,
            0.5,
            1095 + 86399. / 86401,
            1095 + 86400. / 86401,
            1095 + 86401. / 86401,
            1460 + 86398. / 86400,
            1460 + 86399. / 86400,
            1461,
            1462,
            18627. + 86399. / 86401,
            18627. + 86400. / 86401,
            18628.,
            18628. + 1. / 86400,
             ]
        actual = t._days1958totai(inputs, leaps='rubber', midnight=True)
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)
        self.assertEqual(
            expected[7],
            t._days1958totai(inputs[7], leaps='rubber', midnight=True))

        inputs = [
            -364.75,
            0.,
            0.5,
            1095 + 86399. / 86400,
            1095 + 86399.999999 / 86400,
            1096,
            1460 + 86398. / 86400,
            1460 + 86399. / 86400,
            1461,
            1462,
            18627. + 86399. / 86400,
            18627. + 86399.999999 / 86400,
            18628.,
            18628. + 1. / 86400,
            ]
        actual = t._days1958totai(inputs, leaps='drop', midnight=True)
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)
        self.assertEqual(
            expected[7],
            t._days1958totai(inputs[7], leaps='drop', midnight=True))

        inputs = [
            -364.75,
            0.,
            0.5,
            1096,
            1096 + 1. / 86400,
            1096 + 2. / 86400,
            1461,
            1461 + 1. / 86400,
            1461 + 2. / 86400,
            1462 + 2. / 86400,
            18628. + 32. / 86400,
            18628. + 33. / 86400,
            18628. + 34. / 86400,
            18628. + 35. / 86400,
            ]
        actual = t._days1958totai(inputs, leaps='continuous', midnight=True)
        numpy.testing.assert_almost_equal(expected, actual, decimal=6)
        self.assertEqual(
            expected[7],
            t._days1958totai(inputs[7], leaps='continuous', midnight=True))


class TimeClassTests(unittest.TestCase):

    longMessage = True

    def test_TAIinit(self):
        """test that Ticktock can be made from TAI input"""
        t0 = 1663236947
        range_ex = list(range(t0, t0 + 5000, 500))
        tt = t.Ticktock(range_ex, 'TAI')
        ans = ['2010-09-15T10:15:13', '2010-09-15T10:23:33', '2010-09-15T10:31:53']
        numpy.testing.assert_equal(tt.ISO[0:3], ans)

    def test_TAIinit_leapsecond(self):
        """Make from TAI input across a leap second"""
        tt = t.Ticktock([1609459234, 1609459233, 1609459232, 1609459233.1],
                        dtype='TAI')
        numpy.testing.assert_equal(
            [1609459234, 1609459233, 1609459232, 1609459233.1],
            tt.TAI)
        numpy.testing.assert_equal(
            [34, 33, 33, 33], tt.leaps)
        numpy.testing.assert_equal(
            [datetime.datetime(2009, 1, 1),
             datetime.datetime(2008, 12, 31, 23, 59, 59, 999999),
             datetime.datetime(2008, 12, 31, 23, 59, 59),
             datetime.datetime(2008, 12, 31, 23, 59, 59, 999999)],
            tt.UTC)
        numpy.testing.assert_equal(
            ['2009-01-01T00:00:00',
             '2008-12-31T23:59:60',
             '2008-12-31T23:59:59',
             '2008-12-31T23:59:60',
            ],
            tt.ISO
            )

    def test_initRaises(self):
        """Ticktock init has a raise or two"""
        self.assertRaises(ValueError, t.Ticktock, 12345, 'BOGUS')

    def test_notypeguess(self):
        """Raise a reasonable error if can't guess type"""
        self.assertRaises(ValueError, t.Ticktock,
                          [6.36485472e+10, 3.89393000e+11])

    def test_sliceTicktock(self):
        """a ticktock sliced returns a ticktock"""
        n1 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11',
                         '2002-03-01T13:23:11'], 'ISO')
        self.assertTrue(isinstance(n1[:2], t.Ticktock))

    def test_subTicktock(self):
        """a ticktock minus a ticktock is a timedelta"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        n2 = t.Ticktock('2002-02-01T00:00:00', 'ISO')
        self.assertTrue(isinstance(n2 - n1, list))
        self.assertTrue(isinstance((n2 - n1)[0], datetime.timedelta))
        self.assertEqual(28, (n1 - n2)[0].days)
        self.assertEqual(40991, (n1 - n2)[0].seconds)

    def test_subTicktock2(self):
        """a ticktock minus a ticktock is a timedelta iterables too"""
        n1 = t.Ticktock(['2002-03-01T11:23:11'] * 3, 'ISO')
        n2 = t.Ticktock(['2002-02-01T00:00:00'] * 3, 'ISO')
        self.assertEqual(len(n1), 3)
        self.assertEqual(len(n2), 3)
        self.assertTrue(isinstance(n2 - n1, list))
        self.assertTrue(isinstance((n2 - n1)[0], datetime.timedelta))
        self.assertEqual(28, (n1 - n2)[0].days)
        self.assertEqual(40991, (n1 - n2)[0].seconds)
        numpy.testing.assert_equal((n1 - n2), [datetime.timedelta(days=28, seconds=40991)] * 3)

    def test_subTicktockRaises(self):
        """a ticktock minus a bad type errors"""
        n1 = t.Ticktock(['2002-03-01T11:23:11'] * 3, 'ISO')
        self.assertRaises(TypeError, n1.__sub__, ['bob'])

    def test_subTicktockRaises2(self):
        """a ticktock minus a bad shape errors"""
        n1 = t.Ticktock(['2002-03-01T11:23:11'] * 2, 'ISO')
        self.assertRaises(TypeError, n1.__sub__, [datetime.timedelta(seconds=1)] * 3)

    def test_subRaises(self):
        """subtracting ticktocks can raise"""
        n1 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11',
                         '2002-03-01T13:23:11'], 'ISO')
        n2 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11'], 'ISO')
        self.assertRaises(ValueError, n1.__sub__, n2)
        self.assertRaises(TypeError, n1.__sub__, 4)

    def test_subtimedeltalist(self):
        """a ticktock minus a list of timedeltas is a ticktock"""
        n1 = t.Ticktock(['2002-03-01T11:23:11', '2002-03-01T11:23:12'])
        diff = datetime.timedelta(hours=11, minutes=23)
        de = [diff, diff]
        res = t.Ticktock(['2002-03-01T00:00:11', '2002-03-01T00:00:12'])
        numpy.testing.assert_equal(res.UTC, (n1 - de).UTC)
        n1 = t.Ticktock(['2002-03-01T11:23:11', '2002-03-01T11:23:12'])
        diff = datetime.timedelta(hours=11, minutes=23)
        de = diff
        res = t.Ticktock(['2002-03-01T00:00:11', '2002-03-01T00:00:12'])
        numpy.testing.assert_equal(res.UTC, (n1 - de).UTC)

    def test_subtimedelta(self):
        """a ticktock minus a timedelta is a ticktock"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        de = datetime.timedelta(hours=12, seconds=2)
        self.assertEqual(t.Ticktock('2002-02-28T23:23:09', 'ISO'), n1 - de)

    def test_TickTock_with_xrange(self):
        try:
            xrange
        except NameError:
            return  # No xrange in Python 3, so this test is pointless
        t0 = 1663236947
        iter_ex = xrange(t0, t0 + 5000, 500)
        range_ex = list(range(t0, t0 + 5000, 500))
        numpy.testing.assert_equal(t.Ticktock(iter_ex, 'TAI').TAI, t.Ticktock(range_ex, 'TAI').TAI)

    def test_append(self):
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-03', '2002-01-04'])
        expected = t.Ticktock(['2002-01-01', '2002-01-02', '2002-01-03', '2002-01-04'])
        actual_1 = t1.append(t2)
        actual_2 = t1.append(t2.convert('UTC'))
        numpy.testing.assert_equal(expected.RDT, actual_1.RDT)
        numpy.testing.assert_equal(expected.RDT, actual_2.RDT)

    def test_ticktock_ticktock(self):
        """ticktocks are allowed inputs"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(t1)
        expected = t.Ticktock(['2002-01-01', '2002-01-02'])
        numpy.testing.assert_equal(expected.JD, t2.JD)

    def test_ticktock_cdf(self):
        """tests of the CDF time"""
        t1 = t.Ticktock(numpy.array([6.31770624e+13, 6.31771488e+13]))
        expected = ['2002-01-01T00:00:00', '2002-01-02T00:00:00']
        numpy.testing.assert_equal(expected, t1.ISO)
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        expected = [6.31770624e+13, 6.31771488e+13]
        numpy.testing.assert_equal(expected, t1.CDF)

    def test_ticktock_unix(self):
        """Create ticktock with Unix input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([946684800], dtype='UNX')
        self.assertEqual(
            [datetime.datetime(2000, 1, 1)], tt.UTC)

    def testCDFAgainstpycdf(self):
        """Compare CDF time to pycdf calculated time"""
        intimes = [
            datetime.datetime(1, 1, 1),
            datetime.datetime(1066, 1, 2),
            datetime.datetime(1850, 2, 3),
            datetime.datetime(1958, 3, 4),
            datetime.datetime(2002, 1, 1),
            datetime.datetime(2002, 1, 2),
            datetime.datetime(2005, 4, 5),
            datetime.datetime(2030, 5, 6),
        ]
        tt = t.Ticktock(intimes)
        expected = spacepy.pycdf.lib.v_datetime_to_epoch(intimes)
        actual = tt.CDF
        numpy.testing.assert_equal(expected, actual)
        # Go the other way
        tt = t.Ticktock(expected, dtype='CDF')
        actual = tt.UTC
        numpy.testing.assert_equal(intimes, actual)

    def test_setitem(self):
        """setitem should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        expected = ['1999-01-01T00:00:00', '2002-01-02T00:00:00']
        t1[0] = '1999-01-01T00:00:00'
        numpy.testing.assert_equal(expected, t1.ISO)
        t1[:] = ['1999-01-01T00:00:00', '2111-01-02T00:00:00']
        expected = ['1999-01-01T00:00:00', '2111-01-02T00:00:00']
        numpy.testing.assert_equal(expected, t1.ISO)

    def test_delitem(self):
        """delitem should remove item"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        expected = ['2002-01-02T00:00:00']
        del t1[0]
        numpy.testing.assert_equal(expected, t1.ISO)

    def test_removeitem(self):
        """remove should remove item"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        expected = ['2002-01-02T00:00:00']
        t1.remove(0)
        numpy.testing.assert_equal(expected, t1.ISO)

    def test_eq_ne(self):
        """the boolean operations should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t3 = t.Ticktock(['1999-01-01', '1999-01-02'])
        numpy.testing.assert_equal(t1 == t2, [True, True])
        numpy.testing.assert_equal(t1 == t3, [False, False])
        self.assertTrue(t1[0] == datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0] == datetime.datetime(1999, 1, 1))
        self.assertTrue(t1[0] != datetime.datetime(1999, 1, 1))

    def test_le_lt(self):
        """the boolean operations should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-01', '2002-01-02'])
        numpy.testing.assert_equal(t1 <= t2, [True, True])
        numpy.testing.assert_equal(t1 < t2, [False, False])
        self.assertTrue(t1[0] <= datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0] < datetime.datetime(2002, 1, 1))

    def test_ge_gt(self):
        """the boolean operations should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-01', '2002-01-02'])
        numpy.testing.assert_equal(t1 >= t2, [True, True])
        numpy.testing.assert_equal(t1 > t2, [False, False])
        self.assertTrue(t1[0] >= datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0] > datetime.datetime(2002, 1, 1))

    def test_sort(self):
        """sort should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02', '2001-12-12'])
        numpy.testing.assert_equal(t1.argsort(), [2, 0, 1])
        t1.sort()
        expected = ['2001-12-12T00:00:00', '2002-01-01T00:00:00', '2002-01-02T00:00:00']
        numpy.testing.assert_equal(t1.ISO, expected)
        t1 = t.Ticktock(['2002-01-01', '2002-01-01', '2001-12-12'])  # argsort is stable by defualt
        numpy.testing.assert_equal(t1.argsort(), [2, 0, 1])

    def test_sort_leapsecond(self):
        """sort should get the leapsecond right"""
        # Last second of 2008, first second of 2009
        tt = t.Ticktock([1609459234, 1609459233], dtype='TAI')
        # Use a stable sort so doesn't change if don't have to, i.e. this
        # will fail if the leapsecond is ignored.
        tt.sort(kind='mergesort')
        numpy.testing.assert_equal(
            [1609459233, 1609459234],
            tt.data)

    def test_isoformat1(self):
        """can change the iso format '%Y-%m-%dT%H:%M:%S'"""
        t1 = t.Ticktock(['2002-01-01T01:02:12', '2002-01-02T02:04:12', '2001-12-12T23:56:23'])
        t1.isoformat('microseconds')
        expected = ['2002-01-01T01:02:12.000000', '2002-01-02T02:04:12.000000', '2001-12-12T23:56:23.000000']
        numpy.testing.assert_equal(t1.ISO, expected)
        self.assertRaises(ValueError, t1.isoformat, 'badval')

    def test_isoformat2(self):
        """can change the iso format '%Y-%m-%dT%H:%M:%SZ'"""
        t1 = t.Ticktock(['2002-01-01T01:02:12Z', '2002-01-02T02:04:12Z', '2001-12-12T23:56:23Z'])
        t1.isoformat('microseconds')
        expected = ['2002-01-01T01:02:12.000000', '2002-01-02T02:04:12.000000', '2001-12-12T23:56:23.000000']
        numpy.testing.assert_equal(t1.ISO, expected)
        self.assertRaises(ValueError, t1.isoformat, 'badval')

    def test_isoformat3(self):
        """can change the iso format '%Y-%m-%d'"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02', '2001-12-12'])
        t1.isoformat('microseconds')
        expected = ['2002-01-01T00:00:00.000000', '2002-01-02T00:00:00.000000', '2001-12-12T00:00:00.000000']
        numpy.testing.assert_equal(t1.ISO, expected)
        self.assertRaises(ValueError, t1.isoformat, 'badval')

    def test_isoformat3(self):
        """can change the iso format other"""
        t1 = t.Ticktock(['2002-01-01T12', '2002-01-02T23', '2001-12-12T13'])
        t1.isoformat('microseconds')
        expected = ['2002-01-01T12:00:00.000000', '2002-01-02T23:00:00.000000', '2001-12-12T13:00:00.000000']
        numpy.testing.assert_equal(t1.ISO, expected)
        self.assertRaises(ValueError, t1.isoformat, 'badval')

    def test_isoformat4(self):
        """gives a message with nothing specified in isoformat"""
        t1 = t.Ticktock(['2002-01-01T12', '2002-01-02T23', '2001-12-12T13'])
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        self.assertTrue(t1.isoformat() is None)
        result = output.getvalue()
        output.close()
        self.assertTrue(result.startswith("Current ISO"))
        sys.stdout = realstdout

    def test_isoformat_input(self):
        """Supports ISO input format"""
        iso = ['2008-12-31T23:59:60', '2008-12-31T23:59:00']
        tai = [1609459233, 1609459173]
        utc = [datetime.datetime(2008, 12, 31, 23, 59, 59, 999999),
               datetime.datetime(2008, 12, 31, 23, 59)]
        t1 = t.Ticktock(iso)
        numpy.testing.assert_equal(utc, t1.UTC)
        actual = t1.TAI
        numpy.testing.assert_equal(tai, actual)
        t1 = t.Ticktock(tai, dtype='TAI')
        numpy.testing.assert_equal(iso, t1.ISO)

    def test_isoformat_doy(self):
        """Support ISO input in day-of-year, and special format"""
        iso = ['2008-010T12:00:00', '2008-123T10:00:12']
        utc = [(2008, 1, 10, 12),
               (2008, 5, 2, 10, 0, 12)]
        utc = [datetime.datetime(*u) for u in utc]
        t1 = t.Ticktock(iso, isoformat='%Y-%jT%H:%M:%S')
        numpy.testing.assert_equal(utc, t1.UTC)
        # Check the fallback constructor
        isoymd = ['2008-01-10T12:00', '2008-5-2T10:00:12']
        t1 = t.Ticktock(isoymd, isoformat='%Y-%jT%H:%M:%S')
        numpy.testing.assert_equal(utc, t1.UTC)
        # And that it's rendered as desired
        numpy.testing.assert_equal(iso, t1.ISO)

    def test_ISOtoTAIfractional(self):
        """Get fractional TAI seconds from fractional ISO input"""
        iso = ['2001-12-01T01:23:12.123',
               '2016-12-31T23:59:60.100']
        tai = [1385861024.123,
               1861920036.100]
        t1 = t.Ticktock(iso)
        numpy.testing.assert_equal(tai, t1.TAI)

    def test_isoinput_early(self):
        """Input ISO format before 1900."""
        t1 = t.Ticktock(['1890-01-1', '1858-11-18', '1858-11-17',
                         '1858-11-16', '1066-10-14',])
        expected = [(1890, 1, 1), (1858, 11, 18), (1858, 11, 17),
                    (1858, 11, 16), (1066, 10, 14)]
        expected = [datetime.datetime(*e) for e in expected]
        numpy.testing.assert_equal(expected, t1.UTC)

    def test_ISO(self):
        """converting to ISO format should work"""
        t0 = 1663236947
        range_ex = list(numpy.linspace(t0, t0 + 4000, 4))
        # make a TAI that is a leapsecond time
        with spacepy_testing.assertWarns(
                self, 'always',
                r'now\(\) returns UTC time as of 0\.2\.2\.$',
                DeprecationWarning, r'spacepy\.time$'):
            tt2 = t.Ticktock.now()
        tt2tai = tt2.TAI
        taileaps = tt2.TAIleaps
        # Get the TAI at the 2015 and 2012 leap seconds
        range_ex.append(taileaps[35] - 1)
        range_ex.append(taileaps[35])
        range_ex.append(taileaps[35] + 1)
        range_ex.append(taileaps[34])
        tt = t.Ticktock(range_ex, 'TAI')
        ans = ['2010-09-15T10:15:13', '2010-09-15T10:37:26', '2010-09-15T10:59:39',
               '2010-09-15T11:21:53',
               '2015-06-30T23:59:59', '2015-06-30T23:59:60',
               '2015-07-01T00:00:00', '2012-06-30T23:59:60']
        numpy.testing.assert_equal(tt.ISO, ans)

    def test_ISOearly(self):
        """Converting to ISO before 1900 should work"""
        tt = t.Ticktock([datetime.datetime(1800, 1, 1)])
        numpy.testing.assert_equal(
            ['1800-01-01T00:00:00'], tt.ISO)

    def test_DOY(self):
        """DOY conversion should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = [1., 2.]
        numpy.testing.assert_equal(expected, t1.DOY)

    def test_eDOY(self):
        """eDOY conversio should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = [0.04166667, 1.]
        numpy.testing.assert_almost_equal(expected, t1.eDOY)

    def test_astropy_input(self):
        """AstroPy time inputs"""
        apt = astropy.time.Time([
            '2010-01-01',
            '2010-01-01T00:01:00'])
        t1 = t.Ticktock(apt, dtype='APT')
        self.assertEqual('APT', t1.data.attrs['dtype'])
        numpy.testing.assert_array_equal(
            [datetime.datetime(2010, 1, 1),
             datetime.datetime(2010, 1, 1, 0, 1)],
            t1.UTC)
        numpy.testing.assert_array_equal(
            [1640995234.0, 1640995294.0], t1.TAI)
        apt = astropy.time.Time('2010-01-01')
        t1 = t.Ticktock(apt, dtype='APT')
        numpy.testing.assert_array_equal(
            [datetime.datetime(2010, 1, 1)],
            t1.UTC)

    def test_astropy_pickle(self):
        """Pickle with AstroPy time inputs"""
        apt = astropy.time.Time([
            '2010-01-01',
            '2010-01-01T00:01:00'])
        t1 = t.Ticktock(apt, dtype='APT')
        pkl = pickle.dumps(t1)
        t2 = pickle.loads(pkl)
        self.assertEqual('APT', t2.data.attrs['dtype'])
        numpy.testing.assert_array_equal(
            [datetime.datetime(2010, 1, 1),
             datetime.datetime(2010, 1, 1, 0, 1)],
            t2.UTC)
        numpy.testing.assert_array_equal(
            [1640995234.0, 1640995294.0], t2.TAI)

    def test_astropy_output(self):
        """AstroPy time outputs"""
        t1 = t.Ticktock(['2010-01-01',
                        '2010-01-01T06:00:00'])
        # Convert to UTC scale (from TAI)
        apt = t1.APT.utc
        numpy.testing.assert_array_equal(
            [2455197.5, 2455197.75],
            apt.jd)
        t1 = t.Ticktock('2010-01-01')
        apt = t1.APT.utc
        numpy.testing.assert_array_equal(
            [2455197.5],
            apt.jd)

    def test_str(self):
        """TickTock __str__ should give known results"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertEqual(str(t1), "Ticktock( ['2002-01-01T01:00:00' '2002-01-02'], dtype=ISO)")

    def test_TAIGregorian(self):
        """Test TAI across the Gregorian-Julian change"""
        t2 = t.Ticktock([datetime.datetime(1582, 10, 15)])
        t1 = t.Ticktock([datetime.datetime(1582, 10, 4)])
        #1582-10-15 was the day after 1582-10-4
        self.assertEqual(86400, t2.TAI - t1.TAI)

    def test_UTCfromTAIGregorian(self):
        """Test figuring UTC from TAI across calendar change"""
        tai = [
            -11840774400.0, -11840688000.0,
            -11840601600.0, -11840515200.0,
            ]
        utc = [datetime.datetime(*dt) for dt in [
            (1582, 10, 3), (1582, 10, 4),
            (1582, 10, 15), (1582, 10, 16)]]
        t1 = t.Ticktock(tai, dtype='TAI')
        numpy.testing.assert_equal(utc, t1.UTC)

    def test_pickle(self):
        """TickTock objects should pickle"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        pkl = pickle.dumps(t1)
        t2 = pickle.loads(pkl)
        self.assertTrue((t1 == t2).all())

    def test_add_list(self):
        """TickTocks should add properly"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = t.Ticktock(["2002-01-01T01:45:00", "2002-01-02T00:45:00"], dtype='UTC')
        addme = [datetime.timedelta(minutes=45), datetime.timedelta(minutes=45)]
        self.assertTrue((t1 + addme == expected).all())
        self.assertTrue((addme + t1 == expected).all())

    def test_add(self):
        """TickTocks should add properly"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = t.Ticktock(["2002-01-01T01:45:00", "2002-01-02T00:45:00"], dtype='UTC')
        self.assertTrue((t1 + datetime.timedelta(minutes=45) == expected).all())
        self.assertTrue((datetime.timedelta(minutes=45) + t1 == expected).all())

    def test_addRaises(self):
        """adding ticktocks can raise"""
        n1 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11',
                         '2002-03-01T13:23:11'], 'ISO')
        n2 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11'], 'ISO')
        self.assertRaises(TypeError, n1.__add__, n2)  # can't add Ticktocks
        self.assertRaises(TypeError, n1.__add__, [datetime.timedelta(seconds=5)] * 8)
        self.assertRaises(TypeError, n1.__add__, 345)

    def test_insert(self):
        """you can insert to a TickTock"""
        t1 = t.Ticktock(['2002-01-01T01:00:00'])
        expected = t.Ticktock(["2002-01-01T01:00:00", "2002-01-02T00:00:00"], dtype='UTC')
        t1.insert(1, '2002-01-02')
        numpy.testing.assert_array_equal(t1, expected)
        t1 = t.Ticktock(['2002-01-01T01:00:00'])
        t1.insert(1, '2002-01-02', dtype='ISO')
        numpy.testing.assert_array_equal(t1, expected)

    def test_MJD(self):
        """conversions to MJD should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02',
                         '1858-01-01', '1582-10-4', '1582-10-15',
                         '1961-1-1', '1971-12-31', '1972-1-1'])
        expected = numpy.asarray([52275.04166667, 52276.,
                                  -320., -100841., -100840.,
                                  37299.5 + 43201. / 86401,
                                  41316.,
                                  41317])
        numpy.testing.assert_almost_equal(t1.MJD, expected)

    def test_MJDLeapsecond(self):
        """Fractional modified Julian Day on day with leapsecond"""
        t1 = t.Ticktock([
            '1979-12-31T12:00:00',
            '1980-01-01T00:00:00',
            '1980-01-01T12:00:00',
        ])
        expected = numpy.array([
            44238.5, # JD starts at noon, MJD always 0.5 offset.
            44238.5 + 43201. / 86401,
            44239.5
        ])
        numpy.testing.assert_almost_equal(t1.MJD, expected)

    def test_JD(self):
        """Conversion to Julian Day should work"""
        t = t.Ticktock([
            '2000-01-01T12:00:00',
            ])
        expected = numpy.array([
            2451545,
            ])
        numpy.testing.assert_almost_equal(t.JD, expected)

    def test_JDLeapsecond(self):
        """Fractional Julian Day on day with leapsecond"""
        t1 = t.Ticktock([
            '1979-12-31T12:00:00',
            '1980-01-01T00:00:00',
            '1980-01-01T12:00:00',
            '1980-01-02T00:00:00',
        ])
        expected = numpy.array([
            2444239.,
            2444239. + 43201. / 86401,
            2444240.,
            2444240.5,
        ])
        numpy.testing.assert_almost_equal(t1.JD, expected)

    def test_fromMJD(self):
        """conversions from MJD should work"""
        t1 = t.Ticktock([52275 + 1. / 24, 52276.], dtype='MJD')
        expected = [datetime.datetime(2002, 1, 1, 1),
                    datetime.datetime(2002, 1, 2)]
        numpy.testing.assert_almost_equal(
            [(t1.UTC[i] - expected[i]).total_seconds() for i in range(len(t1))],
            0., decimal=6)
        t2 = t.Ticktock([52275 + 1. / 24, 52276.,
                         -320., -100841., -100840.,
                         37300., 41316., 41317.], dtype='MJD')
        expected = [(2002, 1, 1, 1), (2002, 1, 2),
                    (1858, 1, 1), (1582, 10, 4), (1582, 10, 15),
                    (1961, 1, 1), (1971, 12, 31), (1972, 1, 1)]
        expected = [datetime.datetime(*e) for e in expected]
        numpy.testing.assert_almost_equal(
            [(t2.UTC[i] - expected[i]).total_seconds() for i in range(len(t1))],
            0., decimal=6)
        # If get MJD on a non-JD basis, should be able to directly compare.
#        numpy.testing.assert_equal(expected, t2.UTC)

    def testMJDTAIJulian(self):
        """Convert between MJD and TAI around calendar reform"""
        mjd = [-100840, -100841]
        tai = [-11840601600.0, -11840688000.0]
        t1 = t.Ticktock(tai, dtype='TAI')
        numpy.testing.assert_equal(mjd, t1.MJD)
        # The following test fails because TAI to MJD goes through UTC
        # and the UTC-TAI conversion doesn't account for the calendar
        t2 = t.Ticktock(mjd, dtype='MJD')
        numpy.testing.assert_equal(tai, t2.TAI)

    def testMJDTAIOther(self):
        """Convert back and forth between MJD and TAI"""
        mjd = [
            36204.,  #1958-01-01T00
            36205.,  #1958-01-02T00
            36205.5, #1958-01-02T12
            44244.0, #1980-01-06T00
            54831.5 + 43199. / 86401, # 2008-12-31T23:59:59
            54831.5 + 43200. / 86401, # 2008-12-31T23:59:60
            54831.5 + 43201. / 86401, # 2009-01-01T00:00:00
        ]
        tai = [
            0.,         # 1958-01-01T00
            86400.,     # 1958-01-02T00
            129600.,    # 1958-01-02T12
            694656019., # 1980-01-06T00
            1609459232.0, # 2008-12-31T23:59:59
            1609459233.0, # 2008-12-31T23:59:60
            1609459234.0, # 2009-01-01T00:00:00
        ]
        t1 = t.Ticktock(tai, dtype='TAI')
        numpy.testing.assert_almost_equal(mjd, t1.MJD)
        t2 = t.Ticktock(mjd, dtype='MJD')
        numpy.testing.assert_almost_equal(tai, t2.TAI, decimal=6)
        # Round-trip a leapsecond
        self.assertEqual(
            1609459233.,
            t.Ticktock(
                t.Ticktock(1609459233., 'TAI').MJD,
                'MJD').TAI)

    def testRDTtofroTAI(self):
        """Convert RDT to and from TAI"""
        TAI = [
            -61756041600.0, # 0001-01-01T00:00:00
            -61755177600.0, # 0001-11-01T00:00:00
            1609459232.0, # 2008-12-31T23:59:59
            1609459233.0, # 2008-12-31T23:59:60
            1609459234.0, # 2009-01-01T00:00:00
            ]
        RDT = [
            1., # 0001-01-01T00:00:00
            11., # 0001-11-01T00:00:00
            733407. + numpy.float64(86399.) / 86400, # 2008-12-31T23:59:59
            # 2008-12-31T23:59.999999
            733407. + numpy.float64(86399.999999) / 86400,
            733408.0, # 2009-01-01T00:00:00
            ]
        tt1 = t.Ticktock(TAI, 'TAI')
        numpy.testing.assert_equal(RDT, tt1.RDT)
        # RDT can't represent the actual leapsecond, this is the closest
        TAI[-2] = 1609459232.999999
        # RDT is rounded into next day (even on input), so this is
        # what is returned up to the precision available. Can
        # remove this if precision improves.
        TAI[-2] = 1609459234
        tt1 = t.Ticktock(RDT, 'RDT')
        numpy.testing.assert_almost_equal(TAI, tt1.TAI, decimal=5)

    def testRDTtofroUTC(self):
        """Convert RDT to and from UTC"""
        UTC = [
            (1, 1, 1),
            (1, 1, 10),
            (1, 1, 11),
            (100, 1, 1),
            (1582, 10, 4),
            (1582, 10, 15),
            (1858, 11, 17),
            (1958, 1, 1),
            (1999, 7, 1, 14),
            (1999, 7, 1, 14, 0, 3),
            (1999, 7, 1, 14, 3),
            (2009, 1, 1),
            ]
        UTC = [datetime.datetime(*u) for u in UTC]
        RDT = [
            1.,
            10.,
            11.,
            36160.0,
            577725.,
            577736.,
            678576.,
            714780.,
            729936.5833333334,
            729936.5833680555,
            729936.5854166667,
            733408.,
            ]
        tt1 = t.Ticktock(UTC, 'UTC')
        numpy.testing.assert_almost_equal(RDT, tt1.RDT, decimal=6)
        tt1 = t.Ticktock(RDT, 'RDT')
        numpy.testing.assert_almost_equal(
            [(tt1.UTC[i] - UTC[i]).total_seconds() for i in range(len(UTC))],
            0., decimal=5)

    def test_GPS(self):
        """conversions to GPS should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([6.93882013e+08, 6.93964813e+08])
        numpy.testing.assert_almost_equal(t1.GPS, expected)

    def test_now(self):
        """now() is at least deterministic"""
        with spacepy_testing.assertWarns(
                self, 'always',
                r'now\(\) returns UTC time as of 0\.2\.2\.$',
                DeprecationWarning, r'spacepy\.time$'):
            v1 = t.Ticktock.now()
        time.sleep(0.1)
        with spacepy_testing.assertWarns(
                self, 'always',
                r'now\(\) returns UTC time as of 0\.2\.2\.$',
                DeprecationWarning, r'spacepy\.time$'):
            v2 = t.Ticktock.now()
        self.assertTrue(v1 < v2)

    def test_today(self):
        """today() has 0 time"""
        with spacepy_testing.assertWarns(
                self, 'always',
                r'today\(\) returns UTC day as of 0\.2\.2\.$',
                DeprecationWarning, r'spacepy\.time$'):
            v1 = t.Ticktock.today()
        self.assertEqual(v1.UTC[0].hour, 0)
        self.assertEqual(v1.UTC[0].minute, 0)
        self.assertEqual(v1.UTC[0].second, 0)
        self.assertEqual(v1.UTC[0].microsecond, 0)

    def test_UTCGPS(self):
        """testing get UTC from GPS"""
        t1 = t.Ticktock([6.93882013e+08, 6.93964813e+08], 'GPS')
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        numpy.testing.assert_array_equal(t1, expected)

    def test_GPSinput(self):
        """Regressions on GPS input, correct TAI/UTC conversions"""
        gps = [1167264016.5, 1167264017.5, 1167264018., 1167264018.5]
        tai = [1861920035.5, 1861920036.5, 1861920037., 1861920037.5]
        utc = [(2016, 12, 31, 23, 59, 59, 500000),
               (2016, 12, 31, 23, 59, 59, 999999),
               (2017, 1, 1),
               (2017, 1, 1, 0, 0, 0, 500000)]
        utc = [datetime.datetime(*u) for u in utc]
        t1 = t.Ticktock(gps, dtype='GPS')
        numpy.testing.assert_array_equal(tai, t1.TAI)
        # See if straight from TAI works...
        numpy.testing.assert_array_equal(utc, t.Ticktock(tai, dtype='TAI').UTC)
        numpy.testing.assert_array_equal(utc, t1.UTC)

    def test_UTCUNX(self):
        """testing get UTC from UNX"""
        t1 = t.Ticktock([1.00984680e+09, 1.00992960e+09], 'UNX')
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertTrue((t1 == expected).all())

    def test_UTCJD(self):
        """testing get UTC from JD/MJD"""
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        t1 = t.Ticktock([52275.04166667, 52276.], 'MJD')
        #        numpy.testing.assert_allclose(t1.UNX, expected.UNX)
        numpy.testing.assert_almost_equal(t1.UNX, expected.UNX, decimal=-2)
        t1 = t.Ticktock([2452275.54166667, 2452276.5], 'JD')
        self.assertTrue((t1.ISO == expected.ISO).all())

    def test_JD(self):
        """test converting to JD"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([2452275.54166667, 2452276.5])
        numpy.testing.assert_almost_equal(t1.JD, expected)
        t2 = t.Ticktock([datetime.datetime(1582, 10, 4),
                         datetime.datetime(1582, 10, 15)])
        ans = t2.JD
        numpy.testing.assert_almost_equal(ans, [2299159.5, 2299160.5,])

    def test_UTCHasDtype(self):
        """Conversion to UTC has dtype"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertEqual('UTC', t1.UTC.attrs['dtype'])

    def test_getleapsecs(self):
        """preform tests on just getleapsecs"""
        t1 = t.Ticktock([datetime.datetime(1995, 3, 22, 4, 18, 14, 350699),
                         datetime.datetime(1997, 9, 24, 4, 46, 42, 764556),
                         datetime.datetime(1999, 12, 20, 23, 38, 18, 111738),
                         datetime.datetime(2003, 12, 12, 16, 40, 9, 348465),
                         datetime.datetime(2008, 2, 6, 23, 2, 55, 773692),
                         datetime.datetime(2009, 7, 30, 0, 11, 4, 235111),
                         datetime.datetime(2009, 12, 1, 0, 49, 43, 125943),
                         datetime.datetime(2010, 10, 30, 20, 2, 33, 58859),
                         datetime.datetime(2011, 3, 20, 6, 32, 30, 529110),
                         datetime.datetime(2014, 1, 8, 11, 58, 55, 40726)])
        ans = [29, 31, 32, 32, 33, 34, 34, 34, 34, 35]
        numpy.testing.assert_equal(ans, t1.getleapsecs())
        self.assertEqual(29, t.Ticktock(datetime.datetime(1995, 3, 22, 4, 18, 14, 350699)).getleapsecs())

    def test_getleapsecs_early(self):
        """Test of leapseconds in the fractional era"""
        t1 = t.Ticktock([
            datetime.datetime(1958, 1, 1),
            datetime.datetime(1960, 12, 31),
            datetime.datetime(1961, 1, 1),
            datetime.datetime(1964, 1, 1),
            datetime.datetime(1965, 3, 1)], dtype='UTC')
        numpy.testing.assert_equal(
            [0, 1, 2, 3, 4], t1.getleapsecs())

    def test_readleapsecs(self):
        """Test that the leap second file was properly read"""
        numpy.testing.assert_equal(
            numpy.arange(12) + 1, spacepy.time.secs[:12])
        # The date in the file (the moment after the leapsecond, i.e.
        # the first time where the TAI-UTC changes).
        expected = [(1959, 1, 1),
                    (1961, 1, 1),
                    (1963, 7, 1),
                    (1965, 1, 1),
                    (1966, 7, 1),
                    (1967, 7, 1),
                    (1968, 7, 1),
                    (1969, 7, 1),
                    (1970, 7, 1),
                    (1971, 7, 1),
                    (1972, 7, 1),
                    (1973, 1, 1)]
        actual = [(int(y), int(m), int(d))
                  for y, m, d in zip(t.year, t.mon, t.day)][:12]
        numpy.testing.assert_equal(expected, actual)

    def test_leapsgood(self):
        """Test the check for out-of-date leapseconds"""
        # Each case is current time, file mtime, lastleap, leapsgood (or not)
        # Last leap must be july or january
        cases = [[(2021, 1, 5), (2020, 12, 26), (2017, 1, 1), True],
                 [(2021, 7, 5), (2020, 12, 26), (2017, 1, 1), False],
                 [(2021, 7, 5), (2020, 12, 26), (2020, 7, 1), False],
                 [(2021, 7, 5), (2020, 12, 26), (2021, 7, 1), True],
                 [(2020, 12, 26), (2020, 12, 24), (2017, 1, 1), True],
                 [(2018, 6, 2), (2018, 5, 12), (2017, 1, 1), True],
                 [(2018, 6, 2), (2017, 5, 12), (2017, 1, 1), False],
                 [(2017, 6, 2), (2017, 5, 12), (2017, 1, 1), True],
                 ]
        for caseno, caseinfo in enumerate(cases):
            currtime, filetime, lastleap, isgood = caseinfo
            self.assertEqual(isgood, t._leapsgood(
                datetime.datetime(*currtime), datetime.datetime(*filetime),
                datetime.datetime(*lastleap)), 'Case {}'.format(caseno))

    def test_diffAcrossLeaps(self):
        """Do TAI differences across leapseconds"""
        t1 = t.Ticktock([
            datetime.datetime(1960, 12, 31, 23, 59, 58),
            # 1 normal second in between
            datetime.datetime(1960, 12, 31, 23, 59, 59),
            # 2 seconds, one normal and one leap
            datetime.datetime(1961, 1, 1),
            # 1 normal second
            datetime.datetime(1961, 1, 1, 0, 0, 1)], dtype='UTC')
        numpy.testing.assert_equal(
            [1, 2, 1], numpy.diff(t1.TAI))

    def testLeapCount(self):
        """Check leap-second counts (TAI-UTC) at various times"""
        t1 = t.Ticktock([datetime.datetime(*dt) for dt in (
            (1958, 1, 1), (1959, 1, 1), (1961, 1, 1),
            (1971, 12, 31), (1972, 1, 1), (1982, 7, 1),
        )], dtype='UTC')
        numpy.testing.assert_equal(
            [0, 1, 2, 10, 10, 21], t1.leaps)

    def testTAIBase(self):
        """Test the baseline of TAI"""
        t1 = t.Ticktock([
            datetime.datetime(1958, 1, 1),
            datetime.datetime(1961, 1, 1)], dtype='UTC')
        numpy.testing.assert_equal(
            [0, # Start epoch.
             # 1958, 1959, 1960 (leap year) + leap seconds 1959 and 1961
             (3 * 365 + 1) * 86400 + 2,
             ],
            t1.TAI)

    def test_callable_input(self):
        """can pass in a callable to convert to datetime"""
        times = ['2002-01-01T01:00:00', '2002-01-02T02:03:04']
        tt = t.Ticktock(times, dtype=lambda x: datetime.datetime.strptime(x, '%Y-%m-%dT%H:%M:%S'))
        ans = [datetime.datetime(2002, 1, 1, 1, 0, 0), datetime.datetime(2002, 1, 2, 2, 3, 4)]
        numpy.testing.assert_equal(ans, tt.UTC)

    def test_iso_nonstr(self):
        """ISO string works with types other than str"""
        for isostr in (b'2020-01-01T00:00:00',
                       b'2020-01-01T00:00:00Z',
                       b'2020-01-01',
                       b'20200101',
                       b'20200101 00:00:00',
                       b'2020 Jan 1'):
            if str is bytes: #Py2k
                isostr = isostr.decode('ascii')
            tt = t.Ticktock([isostr])
            self.assertEqual(
                datetime.datetime(2020, 1, 1),
                tt.UTC[0], isostr)
            # Do same thing with explicit type
            tt = t.Ticktock([isostr], dtype='ISO')
            self.assertEqual(
                datetime.datetime(2020, 1, 1),
                tt.UTC[0])

    def testUpdateItems(self):
        """Change data and call update"""
        tt = t.Ticktock(['2001-01-01'])
        # All the possible dtypes
        attrs = ['ISO', 'TAI', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'UTC',
                 'eDOY']
        preattrs = [a for a in attrs if hasattr(tt, a)]
        self.assertEqual(
            datetime.datetime(2001, 1, 1),
            tt.UTC[0])
        tt.data[0] = '2002-01-01'
        tt.update_items('data')
        postattrs = [a for a in attrs if hasattr(tt, a)]
        self.assertEqual(
            datetime.datetime(2002, 1, 1),
            tt.UTC[0])
        # Nothing new calculated
        self.assertEqual(preattrs, postattrs)

    def testUpdateItemseDOY(self):
        tt = t.Ticktock(['2001-01-01'])
        # All the possible dtypes
        attrs = ['ISO', 'TAI', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'UTC',
                 'eDOY']
        preattrs = [a for a in attrs if hasattr(tt, a)]
        self.assertEqual(
            0,
            tt.eDOY[0])
        tt.data[0] = '2002-02-02'
        tt.update_items('data')
        postattrs = [a for a in attrs if hasattr(tt, a)]
        self.assertEqual(
            32,
            tt.eDOY[0])
        # Nothing new calculated
        self.assertEqual(preattrs, postattrs)

    def testUpdateItemsNonDefault(self):
        """Change non-data attribute and call update"""
        tt = t.Ticktock([datetime.datetime(2001, 1, 1)])
        self.assertEqual(
            '2001-01-01T00:00:00',
            tt.ISO[0])
        # All the possible dtypes
        attrs = ['ISO', 'TAI', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'UTC',
                 'eDOY']
        preattrs = [a for a in attrs if hasattr(tt, a)]
        tt.ISO[0] = '2002-01-01'
        tt.update_items('ISO')
        postattrs = [a for a in attrs if hasattr(tt, a)]
        # Not smashed by update_items
        self.assertEqual(
            '2002-01-01',
            tt.ISO[0])
        # The update also changed the UTC in the main data
        self.assertEqual(
            datetime.datetime(2002, 1, 1),
            tt.UTC[0])
        self.assertEqual(
            datetime.datetime(2002, 1, 1),
            tt.data[0])
        # But dtype is the same
        self.assertEqual('UTC', tt.data.attrs['dtype'])
        # Nothing new calculated
        self.assertEqual(preattrs, postattrs)

    def testUpdateItemsGiveCls(self):
        """Change data and call update with a class"""
        tt = t.Ticktock(['2001-01-01'])
        self.assertEqual(
            datetime.datetime(2001, 1, 1),
            tt.UTC[0])
        tt.data[0] = '2002-01-01'
        with spacepy_testing.assertWarns(
                self, 'always',
                r'cls argument of update_items was deprecated in 0\.2\.2'
                r' and will be ignored\.$',
                DeprecationWarning, r'spacepy\.time$'):
            tt.update_items(type(tt), 'data')
        self.assertEqual(
            datetime.datetime(2002, 1, 1),
            tt.UTC[0])

    def testDataPersistsTAI(self):
        """Verify that the input data is returned for TAI input"""
        # Leap second at the end of 2008
        tt = t.Ticktock([1609459233], dtype='TAI')
        self.assertEqual([1609459233], tt.data)
        self.assertEqual([1609459233], tt.TAI)
        self.assertEqual([1609459233], tt.getTAI())
        self.assertEqual([1609459233], tt.TAI)

    def testDataPersistsGPS(self):
        """Verify that the input data is returned for GPS input"""
        tt = t.Ticktock([0, 86400], dtype='GPS')
        numpy.testing.assert_equal([0, 86400], tt.data)
        self.assertTrue(tt.GPS is tt.data)
        tt.getGPS()
        self.assertTrue(tt.GPS is tt.data)
        tt.getTAI()
        tt.TAI[0] = -1
        self.assertTrue(tt.GPS is tt.data)
        tt.getGPS()
        self.assertTrue(tt.GPS is tt.data)

    def testDataPersistsCDF(self):
        """Verify input data is returned for CDF input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([63113904000000.], dtype='CDF')
        # Calculate RDT, then munge it
        oldrdt = tt.RDT[0]
        tt.RDT[0] += 1
        # But data is untouched
        self.assertEqual([63113904000000.], tt.getCDF())
        self.assertTrue(tt.data is tt.CDF)
        # And if recalc RDT, it's corrected
        tt.update_items('data')
        self.assertEqual(oldrdt, tt.RDT[0])

    def testDataPersistsJD(self):
        """Verify input data is returned for JD input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([2451544.5], dtype='JD')
        # Calculate UTC, then munge it
        oldutc = tt.UTC[0]
        tt.UTC[0] = datetime.datetime(2001, 1, 1)
        # But data is untouched
        self.assertEqual([2451544.5], tt.getJD())
        self.assertTrue(tt.data is tt.JD)
        # And if recalc UTC, it's corrected
        tt.update_items('data')
        self.assertEqual(datetime.datetime(2000, 1, 1), tt.UTC[0])

    def testDataPersistsMJD(self):
        """Verify input data is returned for MJD input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([51544], dtype='MJD')
        # Calculate JD, then munge it
        oldjd = tt.JD[0]
        tt.JD[0] = 2451545.5
        # But data is untouched
        self.assertEqual([51544], tt.getMJD())
        self.assertTrue(tt.data is tt.MJD)
        # And if recalc JD, it's corrected
        tt.update_items('data')
        self.assertEqual([2451544.5], tt.JD)

    def testDataPersistsUNX(self):
        """Verify input data is returned for Unix input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([946684800], dtype='UNX')
        # Calculate UTC, then munge it
        oldutc = tt.UTC[0]
        tt.UTC[0] = datetime.datetime(2001, 1, 1)
        # But data is untouched
        self.assertEqual([946684800], tt.getUNX())
        self.assertTrue(tt.data is tt.UNX)
        # And if recalc UTC, it's corrected
        tt.update_items('data')
        self.assertEqual(datetime.datetime(2000, 1, 1), tt.UTC[0])

    def testDataPersistsRDT(self):
        """Verify input data is returned for RDT input"""
        # 2000-01-01 00:00:00
        tt = t.Ticktock([730120.], dtype='RDT')
        # Calculate UTC, then munge it
        oldutc = tt.UTC[0]
        tt.UTC[0] = datetime.datetime(2001, 1, 1)
        # But data is untouched
        self.assertEqual([730120.], tt.getRDT())
        self.assertTrue(tt.data is tt.RDT)
        # And if recalc UTC, it's corrected
        tt.update_items('data')
        self.assertEqual(datetime.datetime(2000, 1, 1), tt.UTC[0])

    def testDataPersistsISO(self):
        """Input is returned (transformed) for ISO input; data is untouched."""
        iniso = ['2010-1-1',
                 '2012-06-30T23:59:60',
                 '2012-2-3T23:59:42.123']
        expected = ([
            '2010-01-01T00:00:00',
            '2012-06-30T23:59:60',
            '2012-02-03T23:59:42'])
        tt = t.Ticktock(iniso)
        numpy.testing.assert_equal(expected, tt.ISO)
        numpy.testing.assert_equal(iniso, tt.data)
        numpy.testing.assert_equal(expected, tt.getISO())
        numpy.testing.assert_equal(iniso, tt.data)


if __name__ == "__main__":
    unittest.main()
