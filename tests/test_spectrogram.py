# -*- coding: utf-8 -*-

"""
Test suite for spectrogram module

Copyright ©2010 Los Alamos National Security, LLC.
"""

import unittest

import numpy as np
import datetime
import matplotlib.axes
import matplotlib.collections
import matplotlib.dates as mdates
import matplotlib.pyplot

import spacepy_testing
import spacepy.datamodel as dm
import spacepy.toolbox as tb

from spacepy.plot.spectrogram import simpleSpectrogram, Spectrogram
import spacepy.plot


__all__ = ['spectrogramTests', 'spectrogramDateTests', 'SimpleSpectrogramTests']

class spectrogramTests(unittest.TestCase):
    def setUp(self):
        super(spectrogramTests, self).setUp()
        self.kwargs = {}
        self.kwargs['variables'] = ['xval', 'yval', 'zval']
        np.random.seed(8675309)
        self.data = dm.SpaceData(xval=dm.dmarray(np.random.random_sample(200)),
                                 yval=dm.dmarray(np.random.random_sample(200)),
                                 zval=dm.dmarray(np.random.random_sample(200)))


    def tearDown(self):
        super(spectrogramTests, self).tearDown()

    def test_keywords(self):
        """there is some input checking"""
        self.assertRaises(KeyError, Spectrogram, self.data, variables=['bad'])
        self.assertRaises(KeyError, Spectrogram, self.data, bad_keyword=['bad'])

    def test_init_raise_len(self):
        """__init__ does some checking on data length"""
        self.data['zval'] = []
        self.assertRaises(ValueError, Spectrogram, self.data, variables=self.kwargs['variables'])
        self.data['yval'] = []
        self.assertRaises(ValueError, Spectrogram, self.data, variables=self.kwargs['variables'])
        self.data['xval'] = []
        self.assertRaises(ValueError, Spectrogram, self.data, variables=self.kwargs['variables'])

    def test_defaults(self):
        """run it and check that defaults were set correctly"""
        a = Spectrogram(self.data, variables=self.kwargs['variables'], extended_out=False)
        ans = {'bins': [dm.dmarray([0.00120857, 0.07751865, 0.15382872, 0.2301388, 0.30644887,
                                    0.38275895, 0.45906902, 0.5353791, 0.61168917, 0.68799925,
                                    0.76430932, 0.8406194, 0.91692947, 0.99323955]),
                        dm.dmarray([0.00169679, 0.07848775, 0.1552787, 0.23206965, 0.30886061,
                                    0.38565156, 0.46244251, 0.53923347, 0.61602442, 0.69281538,
                                    0.76960633, 0.84639728, 0.92318824, 0.99997919])],
               'variables': ['xval', 'yval', 'zval'],
               'xlim': (0.0012085702179961411, 0.99323954710300699),
               'ylim': (0.001696792515639145, 0.99997919064162388),
               'zlim': (0.012544022260691956, 0.99059103521121727)}
        for key in ans:
            if key == 'variables':
                self.assertEqual(a.specSettings[key], ans[key])
            else:
                np.testing.assert_allclose(a.specSettings[key], ans[key], rtol=1e-5)
        self.assertRaises(NotImplementedError, a.add_data, self.data)


    def test_defaults_extended(self):
        """run it and check that defaults were set correctly (extended_out)"""
        a = Spectrogram(self.data, variables=self.kwargs['variables'])
        ans = {'bins': [dm.dmarray([0.00120857, 0.07751865, 0.15382872, 0.2301388, 0.30644887,
                                    0.38275895, 0.45906902, 0.5353791, 0.61168917, 0.68799925,
                                    0.76430932, 0.8406194, 0.91692947, 0.99323955]),
                        dm.dmarray([0.00169679, 0.07848775, 0.1552787, 0.23206965, 0.30886061,
                                    0.38565156, 0.46244251, 0.53923347, 0.61602442, 0.69281538,
                                    0.76960633, 0.84639728, 0.92318824, 0.99997919])],
               'variables': ['xval', 'yval', 'zval'],
               'xlim': (0.0012085702179961411, 0.99323954710300699),
               'ylim': (0.001696792515639145, 0.99997919064162388),
               'zlim': (0.012544022260691956, 0.99059103521121727)}
        for key in ans:
            if key == 'variables':
                self.assertEqual(a.specSettings[key], ans[key])
            else:
                np.testing.assert_allclose(a.specSettings[key], ans[key], rtol=1e-5)

    def test_add_data(self):
        """run it and check that add_data correctly"""
        data = dm.SpaceData(xval = dm.dmarray(np.arange(3)), 
                            yval = dm.dmarray(np.arange(3)), 
                            zval = dm.dmarray(np.arange(3)))
        xbins = np.arange(-0.5, 3.5, 2.0)
        ybins = np.arange(-0.5, 3.5, 2.0)
        a = Spectrogram(self.data, variables=self.kwargs['variables'], extended_out=True)
        count = a['spectrogram']['count'][:].copy()
        sm = a['spectrogram']['sum'][:].copy()
        spect = a['spectrogram']['spectrogram'][:].copy()
        a.add_data(self.data) # add te same data back, sum, count will double, spectrogram stays the same
        np.testing.assert_almost_equal(a['spectrogram']['count'].filled(), (count*2).filled())
        np.testing.assert_almost_equal(a['spectrogram']['sum'], sm*2)
        np.testing.assert_almost_equal(a['spectrogram']['spectrogram'], spect)

class spectrogramDateTests(spacepy_testing.TestPlot):
    def setUp(self):
        super(spectrogramDateTests, self).setUp()
        self.kwargs = {'variables': ['xval', 'yval', 'zval']}
        np.random.seed(8675309)
        self.data = dm.SpaceData(
            xval=dm.dmarray([datetime.datetime(2000, 1, 1) + datetime.timedelta(days=nn) for nn in range(200)]),
            yval=dm.dmarray(np.random.random_sample(200)),
            zval=dm.dmarray(np.random.random_sample(200)))

    def test_defaults(self):
        """run it and check that defaults were set correctly"""
        a = Spectrogram(self.data, variables=self.kwargs['variables'])
        ans = {'bins': [dm.dmarray(np.linspace(mdates.date2num(self.data['xval'][0]),
                                               mdates.date2num(self.data['xval'][-1]), 14)),
                        dm.dmarray([0.00169679, 0.07848775, 0.1552787, 0.23206965, 0.30886061,
                                    0.38565156, 0.46244251, 0.53923347, 0.61602442, 0.69281538,
                                    0.76960633, 0.84639728, 0.92318824, 0.99997919])],
               'variables': ['xval', 'yval', 'zval'],
               'ylim': (0.0012085702179961411, 0.99323954710300699),
               'zlim': (0.001696792515639145, 0.99997919064162388)}
        for key in ans:
            if key == 'variables':
                self.assertEqual(a.specSettings[key], ans[key])
            else:
                if key == 'bins':
                    np.testing.assert_allclose(a.specSettings[key], ans[key], atol=1e-2, rtol=1e-3)
                else:
                    np.testing.assert_allclose(a.specSettings[key], ans[key], rtol=1e-5)
        ax = a.plot()
        self.assertTrue(isinstance(ax, matplotlib.axes.SubplotBase))


class SimpleSpectrogramTests(spacepy_testing.TestPlot):
    """Test simpleSpectrogram function"""

    def testBadInputs(self):
        """Pass invalid number of arrays"""
        z = np.full((12, 6), 1.)
        x = np.arange(12)
        with self.assertRaises(TypeError) as cm:
            simpleSpectrogram(x, z)
        self.assertEqual('simpleSpectrogram, takes Z or X, Y, Z', str(cm.exception))

    def testSimpleZ(self):
        """Simple, single input"""
        x = np.linspace(0, np.pi, 12)
        y = np.logspace(0, 2, 6)
        # Power-law in energy, sin in time
        z = 1e4 * np.sin(x)[:, None] * (y ** -2)[None, :] + 1
        ax = simpleSpectrogram(z, ylog=False, cbtitle='COLORBAR')
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        np.testing.assert_array_almost_equal(
            z, mesh.get_array().reshape(z.shape[::-1]).transpose())  # mesh swaps row/column
        self.assertEqual((0, 12), ax.get_xlim())
        self.assertEqual((0, 6), ax.get_ylim())
        fig = ax.get_figure()
        axes = fig.get_axes()
        self.assertEqual(2, len(axes))
        self.assertIs(ax, axes[0])
        cb = axes[1]
        zlim = cb.get_ylim()
        self.assertEqual(1., zlim[0])
        self.assertAlmostEqual(1e4, zlim[1], delta=1e3)
        self.assertEqual('COLORBAR', cb.get_ylabel())

    def testSimpleZGiveAxes(self):
        """Simple, single input, provide axes"""
        z = np.full((12, 6), 1.)
        fig = matplotlib.pyplot.figure()
        ax0 = fig.add_subplot(111)
        ax = simpleSpectrogram(z, ax=ax0)
        self.assertIs(ax0, ax)

    def testSimpleXYZ(self):
        """Simple, three inputs"""
        x = np.linspace(0, np.pi, 12)
        y = np.logspace(0, 2, 6)
        # Power-law in energy, sin in time
        z = 1e4 * np.sin(x)[:, None] * (y ** -2)[None, :] + 1
        ax = simpleSpectrogram(x, y, z)
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        np.testing.assert_array_almost_equal(
            z, mesh.get_array().reshape(z.shape[::-1]).transpose())  # mesh swaps row/column
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        self.assertGreater(xlim[0], -0.5)
        self.assertLess(xlim[0], x[1])
        self.assertGreater(xlim[1], x[-1])
        self.assertLess(xlim[1], 4)
        self.assertGreater(ylim[0], .1)
        self.assertLess(ylim[0], y[1])
        self.assertGreater(ylim[1], y[-1])
        self.assertLess(ylim[1], 250)

    def testSimpleXYZDates(self):
        """Simple, three inputs, x is time"""
        x = np.linspace(0, np.pi, 12)
        y = np.logspace(0, 2, 6)
        # Power-law in energy, sin in time
        z = 1e4 * np.sin(x)[:, None] * (y ** -2)[None, :] + 1
        dt = np.vectorize(lambda d: datetime.datetime(2010, 1, 1)
                          + datetime.timedelta(days=d))(x)
        ax = simpleSpectrogram(dt, y, z)
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        np.testing.assert_array_almost_equal(
            z, mesh.get_array().reshape(z.shape[::-1]).transpose())  # mesh swaps row/column
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        self.assertGreater(xlim[0], mdates.date2num(dt[0]) - .2)
        self.assertLess(xlim[0], mdates.date2num(dt[1]))
        self.assertGreater(xlim[1], mdates.date2num(dt[-1]))
        self.assertLess(xlim[1], mdates.date2num(dt[-1]) + .2)
        self.assertGreater(ylim[0], .1)
        self.assertLess(ylim[0], y[1])
        self.assertGreater(ylim[1], y[-1])
        self.assertLess(ylim[1], 250)

    def testLinearZ(self):
        """Linear Z axis"""
        z = np.arange(72).reshape((12, 6))
        x = np.arange(12)
        y = np.arange(6)
        ax = simpleSpectrogram(x, y, z, zlog=False, ylog=False)
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        np.testing.assert_array_almost_equal(
            z, mesh.get_array().reshape(z.shape[::-1]).transpose())  # mesh swaps row/column

    def testTimeDepY(self):
        """Time-dependent Y axis, linear Z"""
        z = np.full((12, 6), 1.)
        x = np.arange(12)
        # Values are all the same, but "time-dependent"
        y = np.tile(np.logspace(0, 2, 6), (12, 1))
        ax = simpleSpectrogram(x, y, z)
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        data = mesh.get_array()
        np.testing.assert_array_almost_equal(1., data)
        self.assertEqual(6 * 12, data.size)

    def testFillAndLow(self):
        """Distinguish between fill and below range cutoff"""
        z = np.tile(np.arange(1., 6), (10, 1))
        z[0, :2] = .5
        z[1, :2] = .1
        z[2, :2] = 0
        z[3, :2] = np.nan
        ax = simpleSpectrogram(z, ylog=False, vmin=0.5, vmax=6, zero_valid=True)
        mesh =  [c for c in ax.get_children() if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(mesh))
        mesh = mesh[0]
        matplotlib.pyplot.draw()  # Force rendering of mesh colors
        colors = np.rollaxis(mesh.get_facecolor().reshape(5, 10, 4), 1, 0)
        np.testing.assert_array_equal(  # the two .5 are same
            colors[0, 0, :], colors[0, 1, :])
        np.testing.assert_array_equal(  # top strips is same
            colors[0, 4, :], colors[8, 4, :])
        np.testing.assert_array_equal(  # .5 and .1 are the same
            colors[0, 0, :], colors[1, 0, :])
        np.testing.assert_array_equal(  # .1 and 0 are the same (under-value)
            colors[1, 0, :], colors[2, 0, :])
        # 0 and nan are different
        self.assertFalse(np.allclose(colors[0, 2, :], colors[0, 3, :]))


if __name__ == "__main__":
    unittest.main()
