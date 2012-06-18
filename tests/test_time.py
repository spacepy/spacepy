# -*- coding: utf-8 -*-

"""
Test suite for time module

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import datetime
import unittest
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

    def test_no_tzinfo(self):
        """no_tzinfo should have known output"""
        dt = datetime.datetime(2000, 1, 1, tzinfo=datetime.tzinfo('MST'))
        self.assertEqual(dt.replace(tzinfo=None), t.no_tzinfo(dt))
        ans = [datetime.datetime(2000, 1, 1)]*10
        self.assertEqual(ans, t.no_tzinfo([dt]*10))

    def test_num2date_rt(self):
        """round-trip num2date should have same output as input"""
        indate = datetime.datetime(2000, 11, 12, 1, 0)
        self.assertEqual(indate,
                         t.num2date(t.date2num(indate)))


class TimeClassTests(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_Tickdelta(self):
        """Tickdelta should function"""
        tst = t.Tickdelta(hours=12)
        self.assertEqual(43200.0, tst.seconds)

    def test_sliceTicktock(self):
        """a ticktock sliced returns a ticktock"""
        n1 = t.Ticktock(['2002-03-01T11:23:11',
                         '2002-03-01T12:23:11',
                         '2002-03-01T13:23:11'], 'ISO')
        self.assertTrue(isinstance(n1[:2], t.Ticktock))

    def test_subTicktock(self):
        """a ticktock minus a ticktock is a tickdelta"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        n2 = t.Ticktock('2002-02-01T00:00:00', 'ISO')
        self.assertTrue(isinstance(n2 - n1, list))
        self.assertTrue(isinstance((n2 - n1)[0], t.Tickdelta))
        self.assertAlmostEqual(28.47443287, (n1-n2)[0].days, places=8)

    def test_subTickdelta(self):
        """a ticktock minus a Tickdelta is a ticktock"""
        n1 = t.Ticktock('2002-03-01T11:23:11', 'ISO')
        de = datetime.timedelta(hours=12, seconds=2) #t.Tickdelta(hours=12, seconds=2)
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



if __name__ == "__main__":
    unittest.main()
