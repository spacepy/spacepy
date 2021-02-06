#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for plot.utils

Copyright 2010-2015 Los Alamos National Security, LLC.
"""

try:
    import cStringIO as StringIO
    sio = StringIO.StringIO
except ImportError:
    import io
    sio = io.BytesIO
import datetime
import unittest
import warnings

import matplotlib.dates
import matplotlib.pyplot as plt
import numpy
import numpy.testing
import spacepy_testing
import spacepy.plot.utils
import spacepy.time as st
import spacepy.toolbox as tb

__all__ = ['PlotUtilFunctionTests']

class PlotUtilFunctionTests(unittest.TestCase):
    """Tests for plot.utils functions"""

    def test_applySmartTimeTicks(self):
        """applySmartTimeTicks should have known behaviour"""
        ticks = st.tickrange('2002-02-01T00:00:00', '2002-02-07T00:00:00', deltadays=1)
        y = list(range(len(ticks)))
        fig = plt.figure()
        ax = fig.add_subplot(111)
        line = ax.plot(ticks.UTC, y)
        spacepy.plot.utils.applySmartTimeTicks(ax, ticks.UTC)
        fp = sio()
        fig.savefig(fp, format='png')
        fp.close()
        # should not have moved the ticks
        real_ans = matplotlib.dates.date2num([
            datetime.datetime(2002, 2, i) for i in range (1, 8)])
        numpy.testing.assert_almost_equal(real_ans, ax.get_xticks())
        # should have named them 01 Feb, 02 Feb etc
        try:
            real_ans = ['{0:02d} Feb'.format(i+1).decode() for i in range(7)]
        except AttributeError: #Py3k
            real_ans = ['{0:02d} Feb'.format(i+1) for i in range(7)]
        ans = [t.get_text()
               for t in ax.xaxis.get_majorticklabels()]
        numpy.testing.assert_array_equal(real_ans, ans)

    def test_smartTimeTicksSubDay(self):
        """smartTimeTicks should give known output (regression)"""
        # hits all the different cases
        # else
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 10), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%d %b', fmt.fmt)
        # elif nHours < 4:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 1), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nHours < 24:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 13), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nHours < 12:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 11), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # if nHours < 1:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 1, 0, 30), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)

    def test_smartTimeTicksLonger(self):
        """smartTimeTicks should give known output (regression)"""
        # elif nHours < 48:
        t1 = tb.linspace(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 2, 0, 30), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%H:%M UT', fmt.fmt)
        # elif nDays < 32:
        t1 = tb.linspace(datetime.datetime(2000, 1, 2), datetime.datetime(2000, 2, 1, 0, 30), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%d %b', fmt.fmt)
        # elif nDays < 731:
        t1 = tb.linspace(datetime.datetime(2000, 1, 2), datetime.datetime(2001, 2, 1, 0, 30), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%b %Y', fmt.fmt)
        # elif nDays >= 731:
        t1 = tb.linspace(datetime.datetime(2000, 1, 2), datetime.datetime(2010, 2, 1, 0, 30), 20)
        Mtick, mtick, fmt = spacepy.plot.utils.smartTimeTicks(t1)
        self.assertEqual('%Y', fmt.fmt)

    def test_set_target_figureIn(self):
        '''Test that set_target returns expected objects and types'''
        testfig = plt.figure()
        retfig, retax = spacepy.plot.utils.set_target(testfig)
        self.assertTrue(testfig is retfig)
        self.assertTrue(retax is retfig.axes[0])
        plt.close()

    def test_set_target_axesIn(self):
        '''Test that set_target returns expected objects and types'''
        testfig = plt.figure()
        testax = testfig.add_subplot(111)
        retfig, retax = spacepy.plot.utils.set_target(testax)
        self.assertTrue(testfig is retfig)
        self.assertTrue(retax is retfig.axes[0])
        plt.close()

    def test_timestamp(self):
        '''Test that timestamp behaves as expected'''
        ann = spacepy.plot.utils.timestamp(strnow='I am a test string')
        self.assertEqual(ann.get_text(), 'I am a test string')
        self.assertEqual(ann.get_position(), (1.003, 0.01))
        plt.close()

    def test_timestamp2(self):
        '''Test that timestamp behaves as expected2'''
        ann = spacepy.plot.utils.timestamp()
        now = datetime.datetime.now().strftime("%d%b%Y %H:%M")
        self.assertTrue(ann.get_text().startswith(now[:7]))
        plt.close()

    def test_add_arrows(self):
        '''Test different cases for adding arrows to lines.'''
        from matplotlib.collections import LineCollection
        
        fig = plt.figure()
        ax  = fig.add_subplot(111)
        # Basic line2d and tuple of line2Ds:
        x, y = numpy.arange(10), numpy.arange(10)
        line = ax.plot(x,y)[0]
        lines= ax.plot(x,y+5, x,y+10, x,y+15)
        # Gin up a LineCollection:
        xs = numpy.arange(100)
        ys = xs[:10, numpy.newaxis] + xs[numpy.newaxis, :]
        segs = numpy.zeros((10, 100, 2), float)
        segs[:, :, 1] = ys/5
        segs[:, :, 0] = xs/10.
        collect = LineCollection(segs)
        ax.add_collection(collect)

        #### Should finish successfully ####
        # Add arrows to all:
        spacepy.plot.utils.add_arrows(line)
        spacepy.plot.utils.add_arrows(lines)
        spacepy.plot.utils.add_arrows(collect)

        #### Should fail ####
        with self.assertRaises(ValueError):
            spacepy.plot.utils.add_arrows( (line, 1) )
        with self.assertRaises(ValueError):
            spacepy.plot.utils.add_arrows( 'not a line' )
        
        plt.close('all')
        
if __name__ == "__main__":
    unittest.main()
