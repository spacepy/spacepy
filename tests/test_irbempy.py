#!/usr/bin/env	python
#	-*-	coding:	utf-8	-*-
"""
testing	the	irbempy	module

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import unittest
import spacepy
import spacepy.time
import spacepy.coordinates
import spacepy.irbempy as ib
import glob
import os
import numpy as np
import numpy.testing
from numpy import array

__all__ = ['IRBEMBigTests']

class IRBEMBigTests(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
        self.loci = spacepy.coordinates.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')

    def test_prep_irbem(self):
        expected = {
            'badval': -1e31,
            'degalpha': [0.0] * 25,
            'idoysat': [33.0] * 2 + [0.0] * 99998,
            'ntime_max': 100000,
            'nalp_max': 25,
            'magin': np.zeros((25, 100000)),
            'sysaxes': 1,
            'kext': 10,
            'iyearsat': [2002.] * 2 + [0.0] * 99998,
            'xin3': 0.0 * 100000,
            'xin2': 0.0 * 100000,
            'xin1': [3., 2.] + [0.0] * 99998,
            'utsat': [43200., 43800.] + [0.0] * 99998,
            'options': [1, 0, 0, 0, 0],
            }
        expected['magin'][:, 0:2] = numpy.array(
            [[  37.00000048,   37.00000048],
             [ -75.        ,  -74.83333333],
             [   9.        ,    8.94999997],
             [ 358.        ,  362.66666671],
             [   2.20000005,    2.23666672],
             [  -3.        ,   -3.41666667],
             [  -5.        ,   -5.18333332],
             [   7.69000006,    7.82666675],
             [  10.59000015,   10.56666676],
             [  10.86999989,   10.59666649],
             [   1.83500004,    1.81683336],
             [   1.82799995,    1.80666663],
             [   1.16999996,    1.17216663],
             [   2.05200005,    2.02350004],
             [   1.45000005,    1.44783338],
             [   2.94400001,    2.89649999],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ],
             [   0.        ,    0.        ]]
            )
                            
        actual = ib.prep_irbem(self.ticks, self.loci)
        for key in expected:
            numpy.testing.assert_allclose(expected[key],
                                          actual[key],
                                          rtol=1e-6)        
        
    def test_get_dtype(self):
        sysaxes = 3
        expected = ('GSE', 'car')
        self.assertEqual(expected, ib.get_dtype(sysaxes))

    def test_get_sysaxes(self):
        dtype = 'GSE'
        carsph = 'car'
        expected = 3
        self.assertEqual(expected, ib.get_sysaxes(dtype, carsph))
        
    def test_sph2car(self):
        loc = [1,45,45]
        expected = array([ 0.5,  0.5,  0.70710678])	
        numpy.testing.assert_almost_equal(expected, ib.sph2car(loc))

    def test_car2sph(self):
        loc = [ 0.5,  0.5,  0.70710678]
        expected = [1,45,45]
        numpy.testing.assert_almost_equal(expected, ib.car2sph(loc))
        
    def test_coord_trans(self):
        self.loci.ticks = self.ticks
        expected = array([[ 2.86714166, -0.02178308,  0.88262348],
            [ 1.91462214,  0.06992421,  0.57387514]])
        numpy.testing.assert_almost_equal(expected, ib.coord_trans(self.loci, 'GSM', 'car'))
        
    def test_find_Bmirror(self):
        expected = {'Blocal': array([977.70182354,  3398.91281775]),
            'Bmirr': array([2366.26469914,  8226.02095512])}
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key], ib.find_Bmirror(self.ticks, self.loci, [40])[key], decimal=6)

    def test_find_magequator(self):
        expected = {'Bmin': array([976.95969476,  3390.67989733])}
        Bmin_loci =  [[ 2.99929568,  0.0061249 , -0.03561588],
                      [ 2.00296095, -0.00739215,  0.04580058]]
        actual = ib.find_magequator(self.ticks, self.loci)
        numpy.testing.assert_almost_equal(expected['Bmin'], actual['Bmin'], decimal=6)
        numpy.testing.assert_almost_equal(Bmin_loci, actual['loci'].data, decimal=6)
            
    def test_get_Bfield(self):
        """test get_Bfield"""	
        expected = {'Blocal': array([ 977.70182354,  3398.91281775]),
        'Bvec':	array([[6.56029545e-01, -1.65416734e+02, 9.63606626e+02],
                      [3.35038317e+02, -5.42656917e+02, 3.33854476e+03]])}
        actual = ib.get_Bfield(self.ticks, self.loci)
        for key in expected.keys():
            numpy.testing.assert_allclose(actual[key], expected[key], atol=1e-6)

    def test_get_AEP8(self):
        """test get_AEP8"""
        c=self.loci
        c.ticks = self.ticks
        E = 2.0 # energy in MeV
        expected = 99492.059080021136
        actual = ib.get_AEP8(E, c)
        numpy.testing.assert_almost_equal(expected, actual)
        
    def test_get_Lstar_T01(self):
        # test T01STORM
        expected = {'Xj': array([[ 0.00068834], [ 0.00279593]]), 
            'Lstar': array([[ 2.92660617], [ 2.0246814 ]]), 
            'Bmirr': array([[  977.70182354], [ 3398.91281775]]), 
            'Lm': array([[3.13365664], [2.06963969]]), 
            'Bmin': array([  976.95969476,  3390.67989733]), 
            'MLT': array([ 11.97222034,  12.13378624])}    
        actual = ib.get_Lstar(self.ticks, self.loci, [90])
        for key in expected.keys():
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
    
    def test_get_Lstari_T05(self):
        # test T05
        expected = {'Xj': array([[ 0.25023298], [ 0.1833702 ]]), 
                    'Lstar': array([[ 2.92777304], [ 2.02445757]]), 
                    'Bmirr': array([[ 1096.4439606 ], [ 3842.38188155]]), 
                    'Lm': array([[ 3.1275911 ], [ 2.06731527]]), 
                    'Bmin': array([  967.49112103,  3384.64545328]), 
                    'MLT': array([ 11.97222034,  12.13378624])}
        actual = ib.get_Lstar(self.ticks, self.loci, [70], extMag='T05')
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
        
        
# -----------------------------------------------------------------------


if	__name__	==	"__main__":
    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()





