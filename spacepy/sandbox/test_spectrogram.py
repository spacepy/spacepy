# -*- coding: utf-8 -*-

"""
Test suite for spectrogram module

Copyright ©2010 Los Alamos National Security, LLC.
"""

import unittest

import numpy as np

import spacepy.datamodel as dm
#from spcepy import spectrogram
import spectrogram


kwargs = {}
kwargs['variables'] = ['Epoch', 'lValue', 'mepOmni']
np.random.seed(8675309)
data = dm.SpaceData(xval = np.random.random_sample(200), 
                    yval = np.random.random_sample(200), 
                    zval = np.random.random_sample(200))


class spectrogramTests(unittest.TestCase):
    def setUp(self):
        super(spectrogramTests, self).setUp()

    def tearDown(self):
        super(spectrogramTests, self).tearDown()

    def test_keywords(self):
        """there is some input checking"""
        self.assertRaises(KeyError, spectrogram.spectrogram, data, variables=['bad'] )
        self.assertRaises(KeyError, spectrogram.spectrogram, data, bad_keyword=['bad'] )

    def test_defaults(self):
        """run it and check that defualts were set correctly"""
        a = spectrogram.spectrogram(data, variables=['xval', 'yval', 'zval'])
        ans = {'bins': [np.array([ 0.00120857,  0.07751865,  0.15382872,  0.2301388 ,  0.30644887,
                               0.38275895,  0.45906902,  0.5353791 ,  0.61168917,  0.68799925,
                               0.76430932,  0.8406194 ,  0.91692947,  0.99323955]),
                        np.array([ 0.00169679,  0.07848775,  0.1552787 ,  0.23206965,  0.30886061,
                               0.38565156,  0.46244251,  0.53923347,  0.61602442,  0.69281538,
                               0.76960633,  0.84639728,  0.92318824,  0.99997919])],
                'variables': ['xval', 'yval', 'zval'],
                'xlim': (0.0012085702179961411, 0.99323954710300699),
                'ylim': (0.001696792515639145, 0.99997919064162388),
                'zlim': (0.012544022260691956, 0.99059103521121727)}
        for key in ans:
            if key == 'variables':
                self.assertEqual(a.specSettings[key], ans[key])
            else:
                np.testing.assert_allclose(a.specSettings[key], ans[key], rtol=1e-5)

     
if __name__ == "__main__":
    unittest.main()
