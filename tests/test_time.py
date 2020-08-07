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

import numpy

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
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings(
                'always', 'Number of seconds > seconds in day.*',
                UserWarning, '^spacepy\\.time')
            for i, val in enumerate(inval):
                ans = t.sec2hms(*val)
                self.assertEqual(real_ans[i], ans)
        self.assertEqual(1, len(w))
        self.assertEqual(
            'Number of seconds > seconds in day. Try days keyword.',
            str(w[0].message))
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
        try:
            from matplotlib.dates import date2num, num2date
        except ImportError:
            return  # don't even do the test
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
        dt11 = num2date(date2num(dt1))
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
            (2006, 1, 1),
            (2006, 1, 1),
            ]
        expectedUTC = [datetime.datetime(*e) for e in expectedUTC]
        actualiso, actualUTC = t.dtstr2iso(inputs)
        # Make sure not overwritten.
        self.assertEqual('2005-12-31T23:59:60', inputs[-2])
        numpy.testing.assert_equal(expectedUTC, actualUTC)
        numpy.testing.assert_equal(expectediso, actualiso)
        self.assertEqual(
            ('2005-12-31T23:59:60.123000',
             datetime.datetime(2006, 1, 1)),
            t.dtstr2iso('2005-12-31T23:59:60.123',
                        fmt='%Y-%m-%dT%H:%M:%S.%f'))
        # Make inputs numpy (more likely to be overwritten).
        inputs = numpy.array(inputs)
        actualiso, actualUTC = t.dtstr2iso(inputs)
        # Make sure not overwritten.
        self.assertEqual('2005-12-31T23:59:60', inputs[-2])


class TimeClassTests(unittest.TestCase):

    def test_TAIinit(self):
        """test that Ticktock can be made from TAI input"""
        t0 = 1663236947
        range_ex = list(range(t0, t0 + 5000, 500))
        tt = t.Ticktock(range_ex, 'TAI')
        ans = ['2010-09-15T10:15:13', '2010-09-15T10:23:33', '2010-09-15T10:31:53']
        numpy.testing.assert_equal(tt.ISO[0:3], ans)

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
            datetime.datetime(2005, 4, 5),
            datetime.datetime(2030, 5, 6),
        ]
        tt = t.Ticktock(intimes)
        expected = spacepy.pycdf.lib.v_datetime_to_epoch(intimes)
        actual = tt.CDF
        numpy.testing.assert_equal(expected, actual)

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

    @unittest.expectedFailure
    def test_isoformat_input(self):
        """Supports ISO input format"""
        t1 = t.Ticktock(['2008-12-31T23:59:60', '2008-12-31T23:59:00'])
        expected = [datetime.datetime(2009, 1, 1),
                    datetime.datetime(2008, 12, 31, 23, 59)]
        numpy.testing.assert_equal(expected, t1.UTC)

    def test_ISO(self):
        """converting to ISO format should work"""
        t0 = 1663236947
        range_ex = list(numpy.linspace(t0, t0 + 4000, 4))
        # make a TAI that is a leapsecond time
        tt2 = t.Ticktock.now()
        tt2tai = tt2.TAI
        taileaps = tt2.TAIleaps
        range_ex.append(taileaps[39] - 1)
        range_ex.append(taileaps[39])
        range_ex.append(taileaps[39] + 1)
        range_ex.append(taileaps[38])
        tt = t.Ticktock(range_ex, 'TAI')
        ans = ['2010-09-15T10:15:13', '2010-09-15T10:37:26', '2010-09-15T10:59:39',
               '2010-09-15T11:21:53',
               '2015-06-30T23:59:59', '2015-06-30T23:59:60',
               '2015-07-01T00:00:00', '2012-06-30T23:59:60']
        numpy.testing.assert_equal(tt.ISO, ans)

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

    def test_str(self):
        """TickTock __str__ should give known results"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertEqual(str(t1), "Ticktock( ['2002-01-01T01:00:00' '2002-01-02'], dtype=ISO)")

    @unittest.expectedFailure
    def test_TAIGregorian(self):
        """Test TAI across the Gregorian-Julian change"""
        t2 = t.Ticktock([datetime.datetime(1582, 10, 15)])
        t1 = t.Ticktock([datetime.datetime(1582, 10, 4)])
        #1582-10-15 was the day after 1582-10-4
        self.assertEqual(86400, t2.TAI - t1.TAI)

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
        self.assertTrue((t1 == expected).all())
        t1 = t.Ticktock(['2002-01-01T01:00:00'])
        t1.insert(1, '2002-01-02', dtype='ISO')
        self.assertTrue((t1 == expected).all())

    def test_MJD(self):
        """conversions to MJD should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([52275.04166667, 52276.])
        numpy.testing.assert_almost_equal(t1.MJD, expected)

    @unittest.expectedFailure
    def test_MJDLeapsecond(self):
        """Fractional modified Julian Day on day with leapsecond"""
        t1 = t.Ticktock([
            '1979-12-31T12:00:00',
            '1980-01-01T12:00:00',
        ])
        expected = numpy.array([
            44238 + 43200. / 86491,
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

    @unittest.expectedFailure
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
            2444239. + 43201. / 86491,
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
            0., decimal=5)

    def test_GPS(self):
        """conversions to GPS should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([6.93882013e+08, 6.93964813e+08])
        numpy.testing.assert_almost_equal(t1.GPS, expected)

    def test_now(self):
        """now() is at least deterministic"""
        v1 = t.Ticktock.now()
        time.sleep(0.1)
        v2 = t.Ticktock.now()
        self.assertTrue(v1 < v2)

    def test_today(self):
        """today() has 0 time"""
        v1 = t.Ticktock.today()
        self.assertEqual(v1.UTC[0].hour, 0)
        self.assertEqual(v1.UTC[0].minute, 0)
        self.assertEqual(v1.UTC[0].second, 0)
        self.assertEqual(v1.UTC[0].microsecond, 0)

    def test_UTCGPS(self):
        """testing get UTC from GPS"""
        t1 = t.Ticktock([6.93882013e+08, 6.93964813e+08], 'GPS')
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        #        numpy.testing.assert_array_equal(t1, expected)
        self.assertTrue((t1 == expected).all())

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
        t2 = t.Ticktock(datetime.datetime(1582, 10, 14))
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings(
                'always', 'Calendar date before the switch from Julian.*',
                UserWarning, '^spacepy\\.time')
            ans = t2.JD
        self.assertEqual(1, len(w))
        self.assertEqual(UserWarning, w[0].category)
        numpy.testing.assert_almost_equal(ans, [2299169.5])

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
            datetime.datetime(1965, 3, 1)], dtype='UTC')
        numpy.testing.assert_equal(
            [0, 0, 1, 3], t1.getleapsecs())

    def test_readleapsecs(self):
        """Test that the leap second file was properly read"""
        # Tickle the system to force the globals to be read
        t.Ticktock([datetime.datetime(1958, 1, 1)]).TAI
        numpy.testing.assert_equal(
            [1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 4, 4, 10],
            spacepy.time.secs[:14])
        # The date in the file (the moment after the leapsecond, i.e.
        # the first time where the TAI-UTC changes).
        expected = [(1961, 1, 1),
                    (1961, 8, 1),
                    (1962, 1, 1),
                    (1963, 11, 1),
                    (1964, 1, 1),
                    (1964, 4, 1),
                    (1964, 9, 1),
                    (1965, 1, 1),
                    (1965, 3, 1),
                    (1965, 7, 1),
                    (1965, 9, 1),
                    (1966, 1, 1),
                    (1968, 2, 1),
                    (1972, 1, 1)]
        actual = [(int(y), int(m), int(d))
                  for y, m, d in zip(t.year, t.mon, t.day)][:14]
        numpy.testing.assert_equal(expected, actual)

    def test_diffAcrossLeaps(self):
        """Do TAI differences across the first leapsecond"""
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

    def testTAIBase(self):
        """Test the baseline of TAI"""
        t1 = t.Ticktock([
            datetime.datetime(1958, 1, 1),
            datetime.datetime(1961, 1, 1)], dtype='UTC')
        numpy.testing.assert_equal(
            [0, # Start epoch.
             (3 * 365 + 1) * 86400 + 1, # 1958, 1959, 1960 (leap) + 1 second
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
        with warnings.catch_warnings(record=True) as w:
            tt.update_items(type(tt), 'data')
        self.assertEqual(1, len(w))
        self.assertEqual(w[0].category, DeprecationWarning)
        self.assertEqual(
            'cls argument of update_items was deprecated in 0.2.2'
            ' and will be ignored.',
            str(w[0].message))
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
                 '2011-12-31T23:59:60',
                 '2012-2-3T23:59:42.123']
        expected = ([
            '2010-01-01T00:00:00',
            '2011-12-31T23:59:60',
            '2012-02-03T23:59:42'])
        tt = t.Ticktock(iniso)
        numpy.testing.assert_equal(expected, tt.ISO)
        numpy.testing.assert_equal(iniso, tt.data)
        numpy.testing.assert_equal(expected, tt.getISO())
        numpy.testing.assert_equal(iniso, tt.data)


if __name__ == "__main__":
    unittest.main()
