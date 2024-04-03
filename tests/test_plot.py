#!/usr/bin/env python

"""
Test suite for plot
"""

import datetime
import unittest

import matplotlib.collections
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


if __name__ == "__main__":
    unittest.main()
