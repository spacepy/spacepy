#!/usr/bin/env	python
#	-*-	coding:	utf-8	-*-
"""
testing	the	omni	module

Copyright (C)2010 Los Alamos National Security, LLC.
"""

import datetime
import unittest

from numpy import array
import numpy.testing

import spacepy
import spacepy.omni as om


__all__ = ['OmniBigTests', 'OmniOtherTests']

class OmniBigTests(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2001-02-02T12:00:00', '2001-02-02T12:10:00'], 'ISO')
        self.actual = om.get_omni(self.ticks, dbase='Test')
        self.expected = {'ByIMF': array([ 0.2       , -0.04999999]),
                    'Bz1': array([-0.1       ,  0.13333333]),
                    'Bz2': array([-0.1       ,  0.13333333]),
                    'Bz3': array([-0.1       ,  0.13333333]),
                    'Bz4': array([-0.1       ,  0.13333333]),
                    'Bz5': array([-0.1       ,  0.13333333]),
                    'Bz6': array([-0.1       ,  0.13333333]),
                    'BzIMF': array([-0.1       ,  0.13333333]),
                    'DOY': array([ 33.,  33.]),
                    'Dst': array([-9., -9.]),
                    'G': array([[ 0.01      ,  0.03      ,  0.01      ],
                           [ 0.01      ,  0.025     ,  0.00833333]]),
                    'Hr': array([ 12.        ,  12.16666667]),
                    'Kp': array([ 0.30000001,  0.30000001]),
                    'Pdyn': array([ 1.07000005,  1.05500005]),
                    'Qbits': {'ByIMF': array([ 1.,  1.,  1.,  2.]),
                     'BzIMF': array([ 0.,  2.,  2.,  2.]),
                     'G': array([[ 2.,  2.,  2.],
                            [ 2.,  2.,  2.],
                            [ 2.,  2.,  2.],
                            [ 2.,  2.,  2.]]),
                     'Pdyn': array([ 0.,  1.,  0.,  0.]),
                     'W': array([[ 1.,  0.,  0.,  2.,  2.,  2.],
                            [ 2.,  0.,  0.,  2.,  2.,  2.],
                            [ 1.,  0.,  0.,  2.,  2.,  2.],
                            [ 1.,  0.,  0.,  2.,  2.,  2.]]),
                     'dens': array([ 1.,  1.,  1.,  1.]),
                     'velo': array([ 2.,  2.,  2.,  2.])},
                    'RDT': array([ 730518.5       ,  730518.50694444]),
                    'UTC': array([datetime.datetime(2001, 2, 2, 12, 0),
                           datetime.datetime(2001, 2, 2, 12, 10)]),
                    'W': array([[ 0.026     ,  0.017     ,  0.31600001,  0.006     ,  0.017     ,
                             0.022     ],
                           [ 0.02466667,  0.01566667,  0.31433334,  0.0055    ,  0.015     ,
                             0.01983333]]),
                    'Year': array([2001, 2001]),
                    'akp3': array([ 1.02999997,  0.99499997]),
                    'dens': array([ 3.20000005,  3.15000006]),
                    'velo': array([ 396.,  396.])}

    def tearDown(self):
        pass

    def test_get_omni(self):
        keylist = list(self.expected.keys())
        keylist.remove('UTC')
        keylist.remove('Qbits')
        for key in keylist:
            numpy.testing.assert_almost_equal(self.expected[key], self.actual[key], decimal=6)
        for key in self.actual['Qbits']:
            numpy.testing.assert_almost_equal(self.expected['Qbits'][key], self.actual['Qbits'][key])

    def test_get_omni_timerange(self):
        self.assertTrue( (self.ticks.data == self.actual['ticks'].data).all() )

    def test_get_omni_outside_range(self):
        ticks = spacepy.time.Ticktock(['2525-01-01T12:00:00', '2525-01-03T12:10:00'], 'ISO')
        self.assertRaises(ValueError, om.get_omni, ticks, dbase='Test')

class OmniOtherTests(unittest.TestCase):
    def test_omnirange(self):
        expected = (datetime.datetime(1999, 7, 1, 14, 0), datetime.datetime(2001, 10, 11, 23, 0))
        actual = om.omnirange('Test')
        numpy.testing.assert_equal(expected, actual)
        


# -----------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()





