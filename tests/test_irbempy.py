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
        expected['magin'][:, :2] = array(
             [[ 37.00000048,   37.00000048],
              [-75.,          -74.83333333],
              [ 9.5,             9.4666667],
              [ 359.,         363.50000004],
              [ 2.4000001,      2.44166676],
              [-3.,            -3.41666667],
              [-5.,            -5.18333332],
              [ 7.71000004,     7.84500011],
              [ 10.61999989,   10.59499995],
              [ 11.43000031,   11.16166687],
              [ 1.86000001,     1.84266669],
              [ 1.86699998,     1.84616665],
              [ 1.17499995,     1.17749997],
              [ 2.09100008,     2.06300006],
              [ 1.48500001,     1.48333335],
              [ 3.14599991,     3.09949994],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.],
              [ 0.,                     0.]]
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
        expected = {'Blocal': array([  978.7877728 ,  3399.56718463]),
            'Bmirr': array([ 2368.85605375,  8227.73648272])}
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key], ib.find_Bmirror(self.ticks, self.loci, [40])[key], decimal=6)

    def test_find_magequator(self):
        expected = {'Bmin': array([  978.04944134,  3391.33985425])}
        Bmin_loci = [[ 2.99929697,  0.00611795, -0.03557354],
                     [ 2.00296012, -0.00739183,  0.045796  ]] 
        actual = ib.find_magequator(self.ticks, self.loci)
        numpy.testing.assert_almost_equal(expected['Bmin'], actual['Bmin'], decimal=6)
        numpy.testing.assert_almost_equal(Bmin_loci, actual['loci'].data, decimal=6)
            
    def test_get_Bfield(self):
        """test get_Bfield"""	
        expected = {'Blocal': array([  978.7877728 ,  3399.56718463]),
        'Bvec':	array([[  6.88144933e-01,  -1.65610809e+02,   9.64675122e+02],
                       [  3.35026863e+02,  -5.42786970e+02,   3.33919097e+03]])}
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
        expected = {'Xj': array([[ 0.00068403], [ 0.00279439]]), 
            'Lstar': array([[ 2.92661098], [ 2.02468639]]), 
            'Bmirr': array([[  978.7877728 ], [ 3399.56718463]]), 
            'Lm': array([[ 3.13249547], [ 2.06950628]]), 
            'Bmin': array([  978.04944134,  3391.33985425]), 
            'MLT': array([ 11.97222034,  12.13378624])}    
        actual = ib.get_Lstar(self.ticks, self.loci, [90])
        for key in expected.keys():
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
    
    def test_get_Lstar_T05(self):
        # test T05
        expected = {'Xj': array([[0.25054612], [0.18343388]]), 
                    'Lstar': array([[2.92611998], [ 2.02463598]]), 
                    'Bmirr': array([[ 1097.45134429], [ 3843.28257415]]), 
                    'Lm': array([[ 3.12680641], [ 2.06718824]]), 
                    'Bmin': array([  968.39952101,  3385.42654697]), 
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





