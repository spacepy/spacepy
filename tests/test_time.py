# -*- coding: utf-8 -*-

"""
Test suite for time module

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import datetime
import unittest
import pickle
import warnings

import numpy

import spacepy.time as t

__all__ = ['TimeFunctionTests', 'TimeClassTests']


class TimeFunctionTests(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_doy2dateconvert(self):
        """doy2date should return a known value for known input"""
        inval = [ (2000, 1),
                  (2001, 34),
                  (2006, 34),
                  (2008, 60),
                  (2008, 366),
                  ([2008],[366])]
        real_ans = [ (1, 1),
                     (2, 3),
                     (2, 3),
                     (2, 29),
                     (12, 31),
                     ([12], [31])]
        for i, val in enumerate(inval):
            ans = t.doy2date(*val)
            ans2 = t.doy2date(*val, dtobj = True)
            self.assertEqual(real_ans[i] , ans)
            try:
                self.assertEqual(real_ans[i],
                                 ([ans2[0].month], [ans2[0].day]))
            except TypeError:
                self.assertEqual(real_ans[i], (ans2.month , ans2.day))

    def test_doy2datefail(self):
        '''doy2date should fail for bad input'''
        inval = ([[2007],[0.5]],
                [2007, 0.5])
        for val in inval:
            func = lambda:t.doy2date(*val)
            self.assertRaises(ValueError, func)

    def test_doy2datefloat(self):
        '''doy2date should work with floats'''
        ans = ( datetime.datetime(2000, 1, 2, 2, 58, 33, 600000),
                datetime.datetime(2000, 1, 2, 0, 0) )
        inval = [ (2000, 2.124, True, True),
                  (2000, 2, True, True)]
        for i, val in enumerate(ans):
            self.assertEqual(val, t.doy2date(*inval[i]))

    def test_tickrange(self):
        """tickrange should return a known value for known input"""
        inval = ( ('2002-02-01T00:00:00', '2002-02-04T00:00:00', 1),
                  ('2002-02-01T00:00:00', '2002-02-04T00:00:00', 0.5) )
        real_ans = ( numpy.array(['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00',
                      '2002-02-04T00:00:00'], dtype='|S19'),
                     numpy.array(['2002-02-01T00:00:00', '2002-02-01T12:00:00', '2002-02-02T00:00:00',
                      '2002-02-02T12:00:00', '2002-02-03T00:00:00', '2002-02-03T12:00:00',
                      '2002-02-04T00:00:00'], dtype='|S19'))
        for i, val in enumerate(inval):
            ans = t.tickrange(*val)
            numpy.testing.assert_equal(real_ans[i], ans.ISO)

    def test_tickrange2(self):
        """tickrange should return a known value for known input (timedelta)"""
        inval = ( ('2002-02-01T00:00:00', '2002-02-04T00:00:00', datetime.timedelta(days=1)),
                  ('2002-02-01T00:00:00', '2002-02-04T00:00:00', datetime.timedelta(hours=12)) )
        real_ans = ( numpy.array(['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00',
                      '2002-02-04T00:00:00'], dtype='|S19'),
                     numpy.array(['2002-02-01T00:00:00', '2002-02-01T12:00:00', '2002-02-02T00:00:00',
                      '2002-02-02T12:00:00', '2002-02-03T00:00:00', '2002-02-03T12:00:00',
                      '2002-02-04T00:00:00'], dtype='|S19'))
        for i, val in enumerate(inval):
            ans = t.tickrange(*val)
            numpy.testing.assert_equal(real_ans[i], ans.ISO)

    def test_sec2hms(self):
        """sec2hms should return a known value for known input"""
        inval = ( (30, False, False),
                  (86401, False, False),
                  (86401, False, True),
                  (30.3, True, False),
                  (3599, False, False) )
        real_ans = ( [0, 0, 30],
                     [24, 0, 1],
                     [0, 0, 1],
                     [0, 0, 30],
                     [0, 59, 59] )
        with warnings.catch_warnings(record=True) as w:
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
        ans = [datetime.datetime(2000, 1, 1)]*10
        self.assertEqual(ans, t.no_tzinfo([dt]*10))

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

        data = ( 1993 + numpy.array(range(10)), 1900, [1993 + val for val in range(10)] )
        real_ans = (
            numpy.array([365, 365, 365,  366, 365, 365, 365,  366, 365, 365]),
            365,
            [365, 365, 365,  366, 365, 365, 365,  366, 365, 365] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), t.leapyear(data[i], True))
            else:
                self.assertEqual(val, t.leapyear(data[i], True))
        real_ans = ( numpy.array([False, False, False,  True, False, False, False, True, False, False]),
                     False,
                     [False, False, False, True, False, False, False,  True, False, False] )
        for i, val in enumerate(real_ans):
            if i == 0:
                self.assertEqual(val.tolist(), t.leapyear(data[i], False))
            else:
                self.assertEqual(val, t.leapyear(data[i], False))

    def test_randomDate(self):
        """randomDate should give known result"""
        try:
            from matplotlib.dates import date2num, num2date
        except ImportError:
            return # don't even do the test
        dt1 = datetime.datetime(2000, 1, 1)
        dt2 = datetime.datetime(2000, 2, 1)
        numpy.random.seed(8675309)
        ans = numpy.array([datetime.datetime(2000,01,26,04,28,10,500070),
                           datetime.datetime(2000,01,24,06,46,39,156905),
                           datetime.datetime(2000,01,12,01,52,50,481431),
                           datetime.datetime(2000,01,07,06,30,26,331312),
                           datetime.datetime(2000,01,13,16,17,48,619577)])
        numpy.testing.assert_array_equal(ans, t.randomDate(dt1, dt2, 5, sorted=False))
        # check the exception
        dt11 = num2date(date2num(dt1))
        self.assertRaises(ValueError, t.randomDate, dt11, dt2)
        ans.sort()
        numpy.random.seed(8675309)
        numpy.testing.assert_array_equal(ans, t.randomDate(dt1, dt2, 5, sorted=True))


class TimeClassTests(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

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
        self.assertEqual(28, (n1-n2)[0].days)
        self.assertEqual(40991, (n1-n2)[0].seconds)

    def test_subtimedelta(self):
        """a ticktock minus a timedelta is a ticktock"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        de = datetime.timedelta(hours=12, seconds=2)
        self.assertEqual(t.Ticktock( '2002-02-28T23:23:09', 'ISO'), n1-de)

    def test_TickTock_with_xrange(self):
        t0 = 1663236947
        iter_ex = range(t0, t0+5000, 500)
        range_ex = list(range(t0, t0+5000, 500))
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
        t1 = t.Ticktock(numpy.array([  6.31770624e+13,   6.31771488e+13]))
        expected = ['2002-01-01T00:00:00', '2002-01-02T00:00:00']
        numpy.testing.assert_equal(expected, t1.ISO)
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        expected = [  6.31770624e+13,   6.31771488e+13]
        numpy.testing.assert_equal(expected, t1.CDF)

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
        numpy.testing.assert_equal(t1==t2, [ True,  True])
        numpy.testing.assert_equal(t1==t3, [ False,  False])
        self.assertTrue(t1[0]==datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0]==datetime.datetime(1999, 1, 1))
        self.assertTrue(t1[0]!=datetime.datetime(1999, 1, 1))

    def test_le_lt(self):
        """the boolean operations should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-01', '2002-01-02'])
        numpy.testing.assert_equal(t1 <= t2, [ True,  True])
        numpy.testing.assert_equal(t1 < t2, [ False,  False])
        self.assertTrue(t1[0] <= datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0] < datetime.datetime(2002, 1, 1))

    def test_ge_gt(self):
        """the boolean operations should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02'])
        t2 = t.Ticktock(['2002-01-01', '2002-01-02'])
        numpy.testing.assert_equal(t1 >= t2, [ True,  True])
        numpy.testing.assert_equal(t1 > t2, [ False,  False])
        self.assertTrue(t1[0] >= datetime.datetime(2002, 1, 1))
        self.assertFalse(t1[0] > datetime.datetime(2002, 1, 1))

    def test_sort(self):
        """sort should work"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02', '2001-12-12'])
        numpy.testing.assert_equal(t1.argsort(), [2, 0, 1])
        t1.sort()
        expected = ['2001-12-12T00:00:00', '2002-01-01T00:00:00', '2002-01-02T00:00:00']
        numpy.testing.assert_equal(t1.ISO, expected)

    def test_isoformat(self):
        """can change the iso format"""
        t1 = t.Ticktock(['2002-01-01', '2002-01-02', '2001-12-12'])
        t1.isoformat('microseconds')
        expected = ['2002-01-01T00:00:00.000000', '2002-01-02T00:00:00.000000', '2001-12-12T00:00:00.000000']
        numpy.testing.assert_equal(t1.ISO, expected)
        self.assertRaises(ValueError, t1.isoformat, 'badval')

    def test_DOY(self):
        """DOY conversion should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = [ 1.,  2.]
        numpy.testing.assert_equal(expected, t1.DOY)

    def test_eDOY(self):
        """eDOY conversio should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = [ 0.04166667,  1.        ]
        numpy.testing.assert_allclose(expected, t1.eDOY)

    def test_str(self):
        """TickTock __str__ should give known results"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertEqual(str(t1), "Ticktock( ['2002-01-01T01:00:00' '2002-01-02'], dtype=ISO)")

    def test_pickle(self):
        """TickTock objects should pickle"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        pkl = pickle.dumps(t1)
        t2 = pickle.loads(pkl)
        self.assertTrue((t1 == t2).all())

    def test_add(self):
        """TickTocks should add properly"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = t.Ticktock( ["2002-01-01T01:45:00",  "2002-01-02T00:45:00"], dtype='UTC')
        self.assertTrue((t1 + datetime.timedelta(minutes=45) == expected).all())
        self.assertTrue((datetime.timedelta(minutes=45) + t1 == expected).all())

    def test_insert(self):
        """you can insert to a TickTock"""
        t1 = t.Ticktock(['2002-01-01T01:00:00'])
        expected = t.Ticktock( ["2002-01-01T01:00:00",  "2002-01-02T00:00:00"], dtype='UTC')
        t1.insert(1, '2002-01-02')
        self.assertTrue((t1 == expected).all())
        t1 = t.Ticktock(['2002-01-01T01:00:00'])
        t1.insert(1, '2002-01-02', dtype='ISO')
        self.assertTrue((t1 == expected).all())

    def test_MJD(self):
        """conversions to MJD should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([ 52275.04166667,  52276.        ])
        numpy.testing.assert_allclose(t1.MJD, expected)

    def test_GPS(self):
        """conversions to GPS should work"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([  6.93882013e+08,   6.93964813e+08])
        numpy.testing.assert_allclose(t1.GPS, expected)

    def test_now(self):
        """now() is at least deterministic"""
        v1 = t.Ticktock.now()
        v2 = t.Ticktock.now()
        self.assertTrue(v1 < v2)

    def test_UTCGPS(self):
        """testing get UTC from GPS"""
        t1 = t.Ticktock([  6.93882013e+08,   6.93964813e+08], 'GPS')
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertTrue((t1 == expected).all())

    def test_UTCUNX(self):
        """testing get UTC from UNX"""
        t1 = t.Ticktock([  1.00984680e+09,   1.00992960e+09], 'UNX')
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        self.assertTrue((t1 == expected).all())

    def test_UTCJD(self):
        """testing get UTC from JD/MJD"""
        expected = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        t1 = t.Ticktock([ 52275.04166667,  52276.        ], 'MJD')
        numpy.testing.assert_allclose(t1.UNX, expected.UNX)
        t1 = t.Ticktock([ 2452275.54166667,  2452276.5       ], 'JD')
        self.assertTrue((t1.ISO == expected.ISO).all())

    def test_JD(self):
        """test converting to JD"""
        t1 = t.Ticktock(['2002-01-01T01:00:00', '2002-01-02'])
        expected = numpy.asarray([ 2452275.54166667,  2452276.5       ])
        numpy.testing.assert_allclose(t1.JD, expected)
        with warnings.catch_warnings(record=True) as w:
            t2 = t.Ticktock(datetime.datetime(1582,10,14))
            ans = t2.JD
            self.assertEqual(1, len(w))
            self.assertEqual(UserWarning, w[0].category)
            numpy.testing.assert_allclose(ans, [2299169.5])




if __name__ == "__main__":
    unittest.main()
