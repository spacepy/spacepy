# -*- coding: utf-8 -*-
#!/usr/bin/env python2.6

import unittest
import datetime as dt
import dateutil.parser as dup
import numpy as np
import spacepy.time as spt
import spacepy.toolbox as tb
import spacepy.empiricals as em
import scipy.integrate as integ

__all__ = ['empFunctionTests', 'PAmodelTests']


class empFunctionTests(unittest.TestCase):
    def setUp(self):
        super(empFunctionTests, self).setUp()
        self.ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)

    def tearDown(self):
        super(empFunctionTests, self).tearDown()

    def test_getPlasmaPause_regress(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 4.05249998,  4.05249998,  4.05249998,  4.05249998,  4.05249998,
                            4.05249998,  4.05249998,  4.22250002,  4.775     ,  4.90250001,
                            4.90250001])
        ans = em.getPlasmaPause(self.ticks, LT=12)
        np.testing.assert_almost_equal(real_ans, ans)

        real_ans = np.array([ 3.76 ,  3.76 ,  4.358,  4.358,  4.358,  4.358,  4.358,  4.358,
                            4.358,  4.542,  5.14])
        ans = em.getPlasmaPause(self.ticks, 'CA1992')
        np.testing.assert_almost_equal(real_ans, ans)

        real_ans = np.array([ 4.35859998, 4.35859998, 4.35859998, 4.35859998, 4.35859998, 
                            4.35859998, 4.35859998, 4.51140002, 5.008, 5.1226, 5.1226])
        ans = em.getPlasmaPause(self.ticks)
        np.testing.assert_almost_equal(real_ans, ans)

        real_ans = np.array([4.7506632,  4.3583292,  4.5369134,  4.86     ,  4.3583292,
                             4.4570714,  4.86     ,  4.9874052,  5.2127764,  5.64     ,
                             4.9874052])
        ans = em.getPlasmaPause(self.ticks, 'RT1970')
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

    def test_getLmax(self):
        """getLmax should give known results (regression)"""
        real_ans = np.array([ 7.9973023,  8.11663  ,  8.7714972,  8.52228  ,  8.6463423,
        8.6463423,  8.6048668,  8.7714972,  8.3179375,  8.8134583,
        9.0677743])
        ans = em.getLmax(self.ticks)
        np.testing.assert_almost_equal(real_ans, ans)
        self.assertRaises(ValueError, em.getLmax, self.ticks, model='bad')

    def test_getMPstandoff(self):
        """getMPstandoff should give known results (regression)"""
        real_ans = np.array([ 10.52909163,  10.91327764,  10.71260773,  10.69958165,
                               9.75129057,  10.76640718,  11.18228247,  11.05199603,
                              11.42648675,  11.8202582 ,  11.18712131])
        ans = em.ShueMP(self.ticks)
        np.testing.assert_almost_equal(real_ans, ans)
        self.assertRaises(TypeError, em.ShueMP, 'bad')
        data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
        real_ans = np.array([ 9.96096838,  8.96790412])
        ans = em.ShueMP(data)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststar(self):
        """getDststar should give known results (regression)"""
        real_ans = np.array([-30.80228229, -26.85289053, -11.2457748 , -17.98012397,
                             -16.1640001 , -13.64888467, -14.75155876, -10.43928609,
                             -21.22360883,  -8.49354146,  -3.29620967])
        ans = em.getDststar(self.ticks)
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststarTuple(self):
        """getDststar should give known results from tuple input (regression)"""
        real_ans = np.array([-30.80228229, -26.85289053, -11.2457748 , -17.98012397,
                             -16.1640001 , -13.64888467, -14.75155876, -10.43928609,
                             -21.22360883,  -8.49354146,  -3.29620967])
        ans = em.getDststar(self.ticks, model=(7.26, 11))
        np.testing.assert_almost_equal(real_ans, ans)

    def test_getDststarError(self):
        """getDststar should give known exceptions on bad input"""
        real_ans = np.array([-30.80228229, -26.85289053, -11.2457748 , -17.98012397,
                             -16.1640001 , -13.64888467, -14.75155876, -10.43928609,
                             -21.22360883,  -8.49354146,  -3.29620967])
        ans = em.getDststar(self.ticks)
        self.assertRaises(ValueError, em.getDststar, self.ticks, model='bad')
        self.assertRaises(ValueError, em.getDststar, self.ticks, model={'a':7.26, 'b':11})

    def test_getDststarOmnivals(self):
        """getDststar should give known result using omnival dict input"""
        dst, pdyn = -10.0, 5.0
        real_ans = dst - np.sqrt(pdyn)
        ans = em.getDststar({'Pdyn': pdyn, 'Dst': dst}, model=(1,0))
        np.testing.assert_almost_equal(real_ans, ans)

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
        """sin^n model should have d_flux that integrates to omniflux divided by 4pi"""
        omniflux = 3000
        dnflux, alphas = em.vampolaPA(omniflux, order=2, alpha=self.PA)
        d_sum = 4*np.pi*integ.simps(dnflux, np.deg2rad(alphas))
#        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = 4*np.pi*integ.simps(dnflux, np.deg2rad(alphas))
#        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)

    def test_vampola_len1list(self):
        """sin^n model should have d_flux that integrates to omniflux divided by 4pi"""
        omniflux = [3000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = 4*np.pi*integ.simps(dnflux, np.deg2rad(alphas))
#        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)
        np.testing.assert_almost_equal(d_sum, omniflux, decimal=3)

    def test_vampola_multival(self):
        """sin^n model should have d_flux that integrates to omniflux divided by 4pi"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = 4*np.pi*integ.simps(dnflux[:,i], np.deg2rad(alphas))
            np.testing.assert_almost_equal(d_sum, omniflux[i])

    def test_vampola_multi_n(self):
        """sin^n model should have d_flux that integrates to omniflux divided by 4pi"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=[2,4], alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = 4*np.pi*integ.simps(dnflux[:,i], np.deg2rad(alphas))
            np.testing.assert_almost_equal(d_sum, omniflux[i])

    def test_vampola_mismatched_order_len(self):
        """mismatch between lengths of input arrays should raise error"""
        omniflux = [3000, 4500, 6000]
        self.assertRaises(ValueError, em.vampolaPA, omniflux, order=[2,4], alpha=self.PA)


if __name__ == "__main__":
    unittest.main()
