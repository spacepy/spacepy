#!/usr/bin/env python

"""
Test suite for plot
"""

import datetime
import unittest
import os

import matplotlib as mpl
import matplotlib.collections
import matplotlib.pyplot as plt
import numpy
import numpy.testing

import spacepy_testing
import spacepy.plot


class PlotFunctionTests(spacepy_testing.TestPlot):
    """Tests for plot functions"""

    def test_levelPlot(self):
        """basic test of levelPlot execution and output"""
        time = numpy.array([datetime.datetime(2010, 1, i) for i in range(1, 7)])
        values = numpy.array([1., 5.1, 2, 3.2, 1, 7])
        ax = spacepy.plot.levelPlot(values, time=time)
        pc = [c for c in ax.get_children()
              if isinstance(c, matplotlib.collections.PolyCollection)]
        self.assertEqual(3, len(pc[0].get_paths()))  # green
        self.assertEqual(1, len(pc[1].get_paths()))  # yellow
        self.assertEqual(2, len(pc[2].get_paths()))  # red, one bar two patches
        numpy.testing.assert_array_equal(
            pc[0].get_fc(),
            [[0., 1., 0., 0.75]])
        numpy.testing.assert_array_equal(
            pc[1].get_fc(),
            [[1., 1., 0., 0.75]])
        numpy.testing.assert_allclose(
            pc[2].get_fc(),
            [[0.8627, 0.07843, .235290, 0.75]], rtol=1e-3)

    def test_plot_STT(self):
        """test STT correctly replaces x tick labels"""
        time = numpy.array([datetime.datetime(2010, 1, i) for i in range(1, 7)])
        values = numpy.array([1., 5.1, 2, 3.2, 1, 7])
        line = spacepy.plot.plot(time, values, smartTimeTicks=True)
        ax = plt.gca()
        xticks = ax.get_xticklabels()
        self.assertEqual('01 Jan', xticks[0].get_text())

    def test_available(self):
        """test available returns the expected plot styles"""
        expected = ['default', 'spacepy', 'spacepy_altgrid', 'altgrid',
                    'spacepy_polar', 'polar']
        actual = spacepy.plot.available()
        self.assertEqual(expected, actual)

    def test_available_retvals(self):
        """test available returns dict of expected plot styles"""
        from spacepy import __path__ as basepath
        expected = {
            'default': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy.mplstyle'),
            'spacepy': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy.mplstyle'),
            'spacepy_altgrid': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_altgrid.mplstyle'),
            'altgrid': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_altgrid.mplstyle'),
            'spacepy_polar': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_polar.mplstyle'),
            'polar': os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_polar.mplstyle')
        }
        actual = spacepy.plot.available(returnvals=True)
        self.assertEqual(expected, actual)

    def test_style(self):
        """test style updates rcParams"""
        spacepy.plot.style(look='spacepy')
        # It seems impossible to test the style:
        # https://stackoverflow.com/questions/44968099/
        self.assertEqual('plasma', mpl.rcParams['image.cmap'])

    def test_revert_style(self):
        """test revert_style removes new keys in rcParams"""
        original = mpl.rcParams.get('image.cmap')
        spacepy.plot.style(look='spacepy')
        spacepy.plot.revert_style()
        self.assertEqual(original, mpl.rcParams.get('image.cmap'))

    def test_dual_half_circle(self):
        """basic test of dual_half_circle execution and output"""
        wedges = spacepy.plot.dual_half_circle(colors=('y', 'k'))
        self.assertEqual((0.75, 0.75, 0.0, 1), wedges[0].get_facecolor())
        self.assertEqual((0.0, 0.0, 0.0, 1), wedges[1].get_facecolor())
        self.assertEqual(71, len(wedges[0].properties()['verts']))
        self.assertEqual(71, len(wedges[1].properties()['verts']))


if __name__ == "__main__":
    unittest.main()
