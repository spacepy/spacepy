#!/usr/bin/env	python
#	-*-	coding:	utf-8	-*-
"""
testing	the	irbempy	module

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import unittest
import spacepy
import spacepy.omni
import spacepy.time
import spacepy.coordinates
try: #if IRBEM install fails, test suite should not break entirely...
    import spacepy.irbempy as ib
    ibflag = False
except ImportError:
    ibflag = True
import glob
import os
import numpy as np
import numpy.testing
from numpy import array

__all__ = ['IRBEMBigTests', 'IRBEMTestsWithoutOMNI']

@unittest.skipIf(ibflag, "Warning: import spacepy.irbempy failed, skipping tests")
class IRBEMBigTests(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2001-02-02T12:00:00', '2001-02-02T12:10:00'], 'ISO')
        self.loci = spacepy.coordinates.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
        self.omnivals = spacepy.omni.get_omni(self.ticks, dbase='Test')

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
            'iyearsat': [2001.] * 2 + [0.0] * 99998,
            'xin3': 0.0 * 100000,
            'xin2': 0.0 * 100000,
            'xin1': [3., 2.] + [0.0] * 99998,
            'utsat': [43200., 43800.] + [0.0] * 99998,
            'options': [1, 0, 0, 0, 0],
            }
        expected['magin'][:, :2] = array(
            [[  3.00000012e+00,   3.00000012e+00],
             [ -9.00000000e+00,  -9.00000000e+00],
             [  3.20000005e+00,   3.15000006e+00],
             [  3.96000000e+02,   3.96000000e+02],
             [  1.07000005e+00,   1.05500005e+00],
             [  2.00000003e-01,  -4.99999917e-02],
             [ -1.00000001e-01,   1.33333326e-01],
             [  9.99999978e-03,   9.99999978e-03],
             [  2.99999993e-02,   2.49999994e-02],
             [  9.99999978e-03,   8.33333313e-03],
             [  2.60000005e-02,   2.46666670e-02],
             [  1.70000009e-02,   1.56666674e-02],
             [  3.16000015e-01,   3.14333344e-01],
             [  6.00000005e-03,   5.50000004e-03],
             [  1.70000009e-02,   1.50000007e-02],
             [  2.19999999e-02,   1.98333332e-02],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00],
             [  0.00000000e+00,   0.00000000e+00]])
                            
        actual = ib.prep_irbem(self.ticks, self.loci, omnivals=self.omnivals)
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key],
                                              actual[key],
                                              decimal=5)        

    def test_find_Bmirror(self):
        expected = {'Blocal': array([1031.00899202,  3451.98936965]),
            'Bmirr': array([ 2495.24300379,  8354.35553522])}
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key], ib.find_Bmirror(self.ticks, self.loci, [40], omnivals=self.omnivals)[key], decimal=6)

    def test_find_magequator(self):
        expected = {'Bmin': array([1030.4563374 ,  3444.07701499])}
        Bmin_loci = [[ 2.999354  ,  0.00551125, -0.03235315],
                     [ 2.00289878, -0.0073488 ,  0.04538181]]
        actual = ib.find_magequator(self.ticks, self.loci, omnivals=self.omnivals)
        numpy.testing.assert_almost_equal(expected['Bmin'], actual['Bmin'], decimal=6)
        numpy.testing.assert_almost_equal(Bmin_loci, actual['loci'].data, decimal=6)
            
    def test_get_Bfield(self):
        """test get_Bfield"""	
        expected = {'Blocal': array([1031.00899202,  3451.98936965]),
        'Bvec': array([[    3.49178255,  -172.79037333,  1016.42059993],
                       [  335.09279969,  -553.03590378,  3390.88406067]])}
        actual = ib.get_Bfield(self.ticks, self.loci, omnivals=self.omnivals)
        for key in expected.keys():
            numpy.testing.assert_almost_equal(actual[key], expected[key], decimal=5)

    def test_get_Lstar_T01(self):
        # test T01STORM
        expected = {'Xj': array([[ 0.00040309], [ 0.00269014]]),
            'Lstar': array([[ 3.02588643], [ 2.05239571]]), 
            'Bmirr': array([[ 1031.00899202], [ 3451.98936965]]), 
            'Lm': array([[ 3.07915135], [ 2.05932648]]), 
            'Bmin': array([ 1030.4563374 ,  3444.07701499]), 
            'MLT': array([ 11.97159175,  12.13313906])}    
        actual = ib.get_Lstar(self.ticks, self.loci, [90], omnivals=self.omnivals)
        for key in expected.keys():
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
    
    def test_get_Lstar_T05(self):
        # test T05
        expected = {'Xj': array([[ 0.2661136 ], [ 0.18600774]]),
                    'Lstar': array([[ 3.01546075], [ 2.04304329]]),
                    'Bmirr': array([[ 1150.67044093], [ 3895.81080489]]), 
                    'Lm': array([[ 3.08702644], [ 2.05973379]]), 
                    'Bmin': array([ 1015.46803142,  3432.14690674]), 
                    'MLT': array([ 11.97159175,  12.13313906])}
        actual = ib.get_Lstar(self.ticks, self.loci, [70], extMag='T05', omnivals=self.omnivals)
        for key in expected:
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)

    def test_AlphaOfK(self):
        '''test calculation of eq. pitch angle from K (regression)'''
        t = spacepy.time.Ticktock(['2001-09-01T04:00:00'], 'ISO') 
        loci = spacepy.coordinates.Coords([-4,0,0], 'GSM', 'car') 
        ans = spacepy.irbempy.AlphaOfK(t, loci, 0.11, extMag='T89', omnivals=self.omnivals) 
        numpy.testing.assert_almost_equal(ans, 50.625, decimal=5)

    def test_find_footpoint(self):
        '''test computation of field line footpoint location/magnitude (regression)'''
        expected = {'Bfoot': numpy.array([47626.95330601,  47625.98974627]),
                    'loci': spacepy.coordinates.Coords([[ 99.2876423 ,  56.14644481, -10.29427227],
                                     [ 99.33380087,  56.14603124, -10.29737286]],
                                     dtype='GDZ', carsph='sph', units=['km', 'deg', 'deg'])}
        y = spacepy.coordinates.Coords([[3,0,0],[3,0,0]], 'GEO', 'car')
        ans = spacepy.irbempy.find_footpoint(self.ticks, y, omnivals=self.omnivals)
        numpy.testing.assert_almost_equal(expected['Bfoot'], ans['Bfoot'], decimal=5)
        numpy.testing.assert_almost_equal(expected['loci'].data, ans['loci'].data, decimal=5)


@unittest.skipIf(ibflag, "Warning: import spacepy.irbempy failed, skipping tests")
class IRBEMTestsWithoutOMNI(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
        self.loci = spacepy.coordinates.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')

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

    def test_get_AEP8(self):
        """test get_AEP8"""
        c=self.loci
        c.ticks = self.ticks
        E = 2.0 # energy in MeV
        expected = 99492.059080021136
        actual = ib.get_AEP8(E, c)
        numpy.testing.assert_almost_equal(expected, actual)
# -----------------------------------------------------------------------


if	__name__	==	"__main__":
    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()





