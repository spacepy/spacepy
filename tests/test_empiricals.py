# -*- coding: utf-8 -*-
#!/usr/bin/env python2.6

import unittest
import warnings

import dateutil.parser as dup
import numpy as np
import spacepy_testing
import spacepy.time as spt
import spacepy.toolbox as tb
import spacepy.empiricals as em
import scipy.integrate as integ
import spacepy.omni as om

__all__ = ['empFunctionTests', 'PAmodelTests']


class empFunctionTests(unittest.TestCase):
    def setUp(self):
        super(empFunctionTests, self).setUp()
        self.ticks = spt.tickrange('2001-01-01T12:00:00','2001-01-04T00:00:00',.25)
        self.omnivals = om.get_omni(self.ticks, dbase='Test')

    def tearDown(self):
        super(empFunctionTests, self).tearDown()

    def test_getPlasmaPause_regress1(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 5.07249999,  4.90250001,  4.64750002,  4.64750002,  4.64750002,
                              4.775     ,  4.22250002,  4.22250002,  4.22250002,  4.22250002,
                              4.22250002])
        ans = em.getPlasmaPause(self.ticks, LT=12, omnivals=self.omnivals)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getPlasmaPause_regress2(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 5.46199999,  5.27800001,  5.00200002,  5.00200002,  5.00200002,
                              5.00200002,  4.54200002,  4.54200002,  4.54200002,  4.54200002,
                              4.54200002])
        ans = em.getPlasmaPause(self.ticks, 'CA1992', LT='all', omnivals=self.omnivals)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getPlasmaPause_regress3(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 5.2754    ,  5.1226    ,  4.89340002,  4.89340002,  4.89340002,
                              5.008     ,  4.51140002,  4.51140002,  4.51140002,  4.51140002,
                              4.51140002])
        ans = em.getPlasmaPause(self.ticks, omnivals=self.omnivals)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getPlasmaPause_regress4(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 5.2127764 ,  4.98740518,  4.75066318,  5.64      ,  4.98740518,
                              4.86      ,  4.45707144,  4.45707144,  4.45707144,  4.98740518,
                              4.45707144])
        ans = em.getPlasmaPause(self.ticks, 'RT1970', omnivals=self.omnivals)
        np.testing.assert_almost_equal(real_ans, ans)
        

    def test_getPlasmaPauseErrors1(self):
        '''tests for exceptions in getPlasmaPause'''
        #check for fail on bad LT
        foo = lambda: em.getPlasmaPause(self.ticks, LT='bad')
        self.assertRaises(ValueError, foo)
        #check for fail on bad model
        bar = lambda: em.getPlasmaPause(self.ticks, model=3)
        self.assertRaises(ValueError, bar)
        #check for fail on LT out of legal range
        spam = lambda: em.getPlasmaPause(self.ticks, LT=25)
        self.assertRaises(IndexError, spam)

    def test_getPlasmaPauseErrors2(self):
        '''test more exceptions in getPlasmaPause'''
        #check for fail on omnivals of wrong type
        spam = lambda: em.getPlasmaPause(self.ticks, omnivals=['Kp', 2.7])
        self.assertRaises(TypeError, spam)
        #check for fail on omnivals without correct inputs
        spam1 = lambda: em.getPlasmaPause(self.ticks, omnivals={'Kp': [2.7]*10})
        self.assertRaises(KeyError, spam1)

    def test_getPlasmaPauseCA1992warn(self):
        with spacepy_testing.assertWarns(
                self, 'always',
                r'No LT dependence currently supported for CA1992 model',
                RuntimeWarning):
            em.getPlasmaPause(self.ticks, model='CA1992', LT=12, omnivals=self.omnivals)

    def test_getLmax(self):
        """getLmax should give known results (regression)"""
        real_ans = np.array([  9.5895175,   9.5895175,   9.2397463,  10.0376575,   9.7226848,
                              10.1287423,  10.1287423,   7.9973023,   8.8555408,   9.4136607,
                               8.9825167])
        ans = em.getLmax(self.ticks, dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)
        self.assertRaises(ValueError, em.getLmax, self.ticks, model='bad', dbase='Test')

    def test_getMagnetopause(self):
        """regression test for getMagnetopause"""
        real_ans = np.array([[[9.59273467e-16,  1.56661246e+01], [1.29836523e+00,  1.48403825e+01],
                              [2.47027293e+00,  1.40096140e+01], [3.52974970e+00,  1.31732052e+01],
                              [4.48803597e+00,  1.23307775e+01], [5.35421146e+00,  1.14821435e+01],
                              [6.13565949e+00,  1.06272740e+01], [6.83841599e+00,  9.76627016e+00],
                              [7.46743497e+00,  8.89934246e+00], [8.02679239e+00,  8.02679239e+00],
                              [8.51984398e+00,  7.14899794e+00], [8.94934819e+00,  6.26640106e+00],
                              [9.31756243e+00,  5.37949717e+00], [9.62631832e+00,  4.48882595e+00],
                              [9.87708053e+00,  3.59496331e+00], [1.00709922e+01,  2.69851422e+00],
                              [1.02089094e+01,  1.80010618e+00], [1.02914267e+01,  9.00383166e-01],
                              [1.03188939e+01,  1.43963780e-14], [1.02914267e+01, -9.00383166e-01],
                              [1.02089094e+01, -1.80010618e+00], [1.00709922e+01, -2.69851422e+00],
                              [9.87708053e+00, -3.59496331e+00], [9.62631832e+00, -4.48882595e+00],
                              [9.31756243e+00, -5.37949717e+00], [8.94934819e+00, -6.26640106e+00],
                              [8.51984398e+00, -7.14899794e+00], [8.02679239e+00, -8.02679239e+00],
                              [7.46743497e+00, -8.89934246e+00], [6.83841599e+00, -9.76627016e+00],
                              [6.13565949e+00, -1.06272740e+01], [5.35421146e+00, -1.14821435e+01],
                              [4.48803597e+00, -1.23307775e+01], [3.52974970e+00, -1.31732052e+01],
                              [2.47027293e+00, -1.40096140e+01], [1.29836523e+00, -1.48403825e+01],
                              [4.27022148e-14, -1.56661246e+01]]])
        ans = em.getMagnetopause(self.ticks[0], dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getMPstandoffTicks(self):
        """getMPstandoff should give known results (regression)"""
        real_ans = np.array([ 10.31889392,  10.95526872,  10.22524928,  10.29643815,
                              11.10576765,  10.07913129,   9.32063995,   9.72841668,
                               9.80611695,  10.10876463,  10.00481683])
        ans = em.ShueMP(self.ticks, dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getMPstandoffScalar(self):
        """getMPstandoff should cast scalars to singletons (regression)"""
        real_ans = [9.96096838]
        ans = em.ShueMP({'P': 2, 'Bz': -2.4})
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getMPstandoffError(self):
        """getMPstandoff should give known exception on bad input"""
        self.assertRaises(TypeError, em.ShueMP, 'bad')
        self.assertRaises(TypeError, em.ShueMP, {'P': 'bad', 'Bz': 'bad'})

    def test_getMPstandoffDict(self):
        """getMPstandoff should give known results (regression)"""
        data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
        real_ans = np.array([ 9.96096838,  8.96790412])
        ans = em.ShueMP(data)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststar(self):
        """getDststar should give known results (regression)"""
        real_ans = np.array([  6.39592274,   7.5959423 ,  -2.44532089,  16.2867821 ,
                              11.01399987,  17.01362653,  13.6810416 , -34.3206288 ,
                             -12.22368788,   1.03764047,  -9.03424456])
        ans = em.getDststar(self.ticks, dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststarTuple(self):
        """getDststar should give known results from tuple input (regression)"""
        real_ans = np.array([  6.39592274,   7.5959423 ,  -2.44532089,  16.2867821 ,
                              11.01399987,  17.01362653,  13.6810416 , -34.3206288 ,
                             -12.22368788,   1.03764047,  -9.03424456])
        ans = em.getDststar(self.ticks, model=(7.26, 11), dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststarError(self):
        """getDststar should give known exceptions on bad input"""
        real_ans = np.array([-30.80228229, -26.85289053, -11.2457748 , -17.98012397,
                             -16.1640001 , -13.64888467, -14.75155876, -10.43928609,
                             -21.22360883,  -8.49354146,  -3.29620967])
        ans = em.getDststar(self.ticks, dbase='Test')
        self.assertRaises(ValueError, em.getDststar, self.ticks, model='bad')
        self.assertRaises(ValueError, em.getDststar, self.ticks, model={'a':7.26, 'b':11})

    def test_getDststarOmnivals(self):
        """getDststar should give known result using omnival dict input"""
        dst, pdyn = -10.0, 5.0
        real_ans = dst - np.sqrt(pdyn)
        ans = em.getDststar({'Pdyn': pdyn, 'Dst': dst}, model=(1,0))
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststarOmniLists(self):
        """getDststar should give known result using omnival dict with list values"""
        real_ans = [-9.26719046, -9.57468886]
        ans = em.getDststar({'Pdyn': [2, 3], 'Dst': [-10, -8]}, dbase='Test')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getExpectedSWTemp_BS06(self):
        """getExpectedSWTemp should give known results with BS06 formula (regression)"""
        real_ans = [70671.539589, 137049.406001, 203427.272413]
        ans = em.getExpectedSWTemp([400, 500, 600], model='BS06')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getExpectedSWTemp_XB15(self):
        """getExpectedSWTemp should give known results with XB15 formula (regression)"""
        real_ans = [45443.09144407, 91022.49564712, 160560.96058384]
        ans = em.getExpectedSWTemp([400, 500, 600], model='XB15')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getExpectedSWTemp_L87(self):
        """getExpectedSWTemp should give known results with L87 formula (regression)"""
        real_ans = [53290, 110250, 156250]
        ans = em.getExpectedSWTemp([400, 500, 600], model='L87')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getExpectedSWTemp_bad_model(self):
        """getExpectedSWTemp should raise a KeyError when passed an invalid model"""
        with self.assertRaises(KeyError):
            em.getExpectedSWTemp([400, 500, 600], model='asdf')

    def test_getExpectedSWTemp_ev(self):
        """getExpectedSWTemp should give known results with ev units"""
        real_ans = np.array([70671.539589, 137049.406001, 203427.272413]) / 1.16045221e4
        ans = em.getExpectedSWTemp([400, 500, 600], model='BS06', units='ev')
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getExpectedSWTemp_bad_units(self):
        """getExpectedSWTemp should raise a ValueError when passed invalid units"""
        with self.assertRaises(ValueError):
            em.getExpectedSWTemp([400, 500, 600], units='asdf')

    def test_getSolarRotation_Carrington(self):
        """make sure getSolarRotation returns known values"""
        dates = spt.Ticktock([dup.parse(t) for t in ['1853-11-10T00:00:00','1973-08-22T18:00:00','2003-03-19T12:02:45']])
        real_ans = np.array([1,1605,2001]).astype(int)
        ans = em.getSolarRotation(dates, rtype='carrington')
        np.testing.assert_almost_equal(real_ans, ans)
        float_ans = np.array([1.003596660715006, 1605.009785410243, 2001.000000356448])
        ans = em.getSolarRotation(dates, rtype='carrington', fp=True)
        np.testing.assert_almost_equal(float_ans, ans)

    def test_getSolarRotation_Bartels(self):
        """make sure getSolarRotation returns known values"""
        dates = spt.Ticktock([dup.parse(t) for t in ['1832-02-08T00:00:00','2004-05-06T12:00:00','2012-12-12T00:00:00']])
        real_ans = np.array([1,2331,2447]).astype(int)
        ans = em.getSolarRotation(dates, rtype='bartels')
        np.testing.assert_almost_equal(real_ans, ans)
        float_ans = np.array([1.0, 2331.0185185185187, 2447.3703703703704])
        ans = em.getSolarRotation(dates, rtype='bartels', fp=True)
        np.testing.assert_almost_equal(float_ans, ans)

    def test_getSolarRotation_BartelsDateTime(self):
        """make sure getSolarRotation returns known values"""
        dates = [dup.parse(t) for t in ['1832-02-08T00:00:00','2004-05-06T12:00:00','2012-12-12T00:00:00']]
        real_ans = np.array([1,2331,2447]).astype(int)
        for dd, aa in zip(dates, real_ans):
            ans = em.getSolarRotation(dd, rtype='bartels')
            np.testing.assert_almost_equal(aa, ans)

class PAmodelTests(unittest.TestCase):
    def setUp(self):
        super(PAmodelTests, self).setUp()
        self.PA = tb.linspace(0.01,179.99,20000)

    def test_vampola_singleval(self):
        """sin^n model should have d_flux that integrates to input omniflux"""
        omniflux = 3000
        dnflux, alphas = em.vampolaPA(omniflux, order=2, alpha=self.PA)
        d_sum = em.omniFromDirectionalFlux(dnflux, alphas, norm=False)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)

        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = em.omniFromDirectionalFlux(dnflux, alphas, norm=False)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)

    def test_vampola_len1list(self):
        """sin^n model should have d_flux that integrates to input omniflux"""
        omniflux = [3000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = em.omniFromDirectionalFlux(dnflux, alphas, norm=False)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)

    def test_vampola_multival(self):
        """sin^n model should have d_flux that integrates to input omniflux"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = em.omniFromDirectionalFlux(dnflux[:,i], alphas, norm=False)
            np.testing.assert_almost_equal(d_sum, omniflux[i])

    def test_vampola_multi_n(self):
        """sin^n model should have d_flux that integrates to input omniflux"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=[2,4], alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = em.omniFromDirectionalFlux(dnflux[:,i], alphas, norm=False)
            np.testing.assert_almost_equal(d_sum, omniflux[i])

    def test_vampola_mismatched_order_len(self):
        """mismatch between lengths of input arrays should raise error"""
        omniflux = [3000, 4500, 6000]
        self.assertRaises(ValueError, em.vampolaPA, omniflux, order=[2,4], alpha=self.PA)

    def test_vampola_defaults(self):
        """vampolaPA should use defaults if no order or alpha is passed"""
        ans = 119.36249766276326
        omniflux = 3000
        dnflux, alphas = em.vampolaPA(omniflux)
        d_sum = em.omniFromDirectionalFlux(dnflux, alphas)
        # somewhat high tolerance since older numpy results vary slightly
        np.testing.assert_allclose(d_sum, ans, rtol=1e-4)

    def test_getSolarProtonSpectra_regress(self):
        """getSolarProtonSpectra() should return constant values"""
        dat = em.getSolarProtonSpectra()
        Eans = [0.1      ,  0.109185 ,  0.1192137,  0.1301636]
        Eflu = np.asarray([2.8990431649e08,   2.6628832715e08,   2.4458237105e08,
                2.2463193728e08], dtype=float)
        np.testing.assert_almost_equal(dat['Energy'][0:4], Eans)
        np.testing.assert_almost_equal(np.asarray(dat['Fluence'][0:4]), Eflu, decimal=2)

    def test_getVampolaOrder(self):
        """getVampolaOrder should give known results (regression)"""
        real_ans = np.array([ 3.095 ,  1.6402])
        ans = em.getVampolaOrder([4, 6.6])
        np.testing.assert_almost_equal(real_ans, ans)

    def test_omniFromDirectionalFlux_scalar_input(self):
        """omniFromDirectionalFlux should raise TypeError on non-iterable flux"""
        d_sum = em.omniFromDirectionalFlux(100, self.PA)
        np.testing.assert_almost_equal(d_sum, 100)

if __name__ == "__main__":
    unittest.main()
