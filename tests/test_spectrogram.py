# -*- coding: utf-8 -*-

"""
Test suite for spectrogram module

Copyright ©2010 Los Alamos National Security, LLC.
"""

import unittest

import numpy as np
import datetime
import matplotlib.dates as mdates

import spacepy_testing
import spacepy.datamodel as dm
import spacepy.toolbox as tb

from spacepy.plot.spectrogram import Spectrogram
import spacepy.plot


__all__ = ['spectrogramTests', 'spectrogramDateTests']

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

class spectrogramDateTests(unittest.TestCase):
    def setUp(self):
        super(spectrogramDateTests, self).setUp()
        self.kwargs = {'variables': ['xval', 'yval', 'zval']}
        np.random.seed(8675309)
        self.data = dm.SpaceData(
            xval=dm.dmarray([datetime.datetime(2000, 1, 1) + datetime.timedelta(days=nn) for nn in range(200)]),
            yval=dm.dmarray(np.random.random_sample(200)),
            zval=dm.dmarray(np.random.random_sample(200)))

    def tearDown(self):
        super(spectrogramDateTests, self).tearDown()

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

if __name__ == "__main__":
    unittest.main()
