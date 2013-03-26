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


__all__ = ['OmniBigTests']

class OmniBigTests(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
        #self.loci = spacepy.coordinates.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')

    def tearDown(self):
        pass

    def test_get_omni(self):
        expected = {'ByIMF': array([-3.        , -3.41666667]),
           'Bz1': array([-5.        , -5.18333332]),
           'Bz2': array([-5.        , -5.18333332]),
           'Bz3': array([-5.        , -5.18333332]),
           'Bz4': array([-5.        , -5.18333332]),
           'Bz5': array([-5.        , -5.18333332]),
           'Bz6': array([-5.        , -5.18333332]),
           'BzIMF': array([-5.        , -5.18333332]),
           'DOY': array([ 33.,  33.]),
           'Dst': array([-75.        , -74.83333333]),
           'G1': array([ 7.71000004,  7.84500011]),
           'G2': array([ 10.61999989,  10.59499995]),
           'G3': array([ 11.43000031,  11.16166687]),
           'Hr': array([ 12.        ,  12.16666667]),
           'Kp': array([ 3.70000005,  3.70000005]),
           'Pdyn': array([ 2.4000001 ,  2.44166676]),
           'Qbits': {'ByIMF': array([ 2.,  2.]),
               'BzIMF': array([ 2.,  2.]),
               'G1': array([ 2.,  2.]),
               'G2': array([ 2.,  2.]),
               'G3': array([ 2.,  2.]),
               'Pdyn': array([ 2.,  2.]),
               'W1': array([ 0.,  0.]),
               'W2': array([ 2.,  2.]),
               'W3': array([ 1.,  1.]),
               'W4': array([ 2.,  2.]),
               'W5': array([ 2.,  2.]),
               'W6': array([ 2.,  2.]),
               'dens': array([ 2.,  2.]),
               'velo': array([ 2.,  2.])},
           'RDT': array([ 730883.5       ,  730883.50694444]),
           'UTC': [datetime.datetime(2002, 2, 2, 12, 0),
                   datetime.datetime(2002, 2, 2, 12, 10)],
           'W1': array([ 1.86000001,  1.84266669]),
           'W2': array([ 1.86699998,  1.84616665]),
           'W3': array([ 1.17499995,  1.17749997]),
           'W4': array([ 2.09100008,  2.06300006]),
           'W5': array([ 1.48500001,  1.48333335]),
           'W6': array([ 3.14599991,  3.09949994]),
           'Year': array([2002, 2002]),
           'akp3': array([ 3.95000005,  3.94333339]),
           'dens': array([ 9.5      ,  9.4666667]),
           'ticks': spacepy.time.Ticktock( ['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO'),
           'velo': array([ 359.        ,  363.50000004])}

        actual = om.get_omni(self.ticks)
        keylist = expected.keys()
        keylist.remove('ticks')
        keylist.remove('UTC')
        keylist.remove('Qbits')
        for key in keylist:
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
        for key in actual['Qbits']:
            numpy.testing.assert_almost_equal(expected['Qbits'][key], actual['Qbits'][key])
        self.assertTrue( (expected['ticks'].data == actual['ticks'].data).all() )

    def test_get_omni_outside_range(self):
        ticks = spacepy.time.Ticktock(['2525-01-01T12:00:00', '2525-01-03T12:10:00'], 'ISO')
        self.assertRaises(ValueError, om.get_omni, ticks)



# -----------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()





