# -*- coding: utf-8 -*-
#!/usr/bin/env python2.6

import unittest
import numpy as np
import spacepy.time as spt
import spacepy.toolbox as tb
import spacepy.empiricals as em
import scipy.integrate as integ

class empFunctionTests(unittest.TestCase):
    def setUp(self):
        super(empFunctionTests, self).setUp()
        self.ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)

    def tearDown(self):
        super(empFunctionTests, self).tearDown()

    def test_getPlasmaPause_regress(self):
        """regression test for plasmapause location"""
        real_ans = np.array([ 6.3475,  6.3475,  6.3475,  6.3475,  6.3475,  6.3475,  6.3475,
                           6.1775,  5.625 ,  5.4975,  5.4975])
        ans = em.getPlasmaPause(self.ticks, LT=12)

        real_ans = np.array([ 3.76 ,  3.76 ,  4.358,  4.358,  4.358,  4.358,  4.358,  4.358,
                            4.358,  4.542,  5.14])
        ans = em.getPlasmaPause(self.ticks, 'CA1992')
        np.testing.assert_allclose(real_ans, ans)

        real_ans = np.array([ 6.4214,  6.4214,  6.4214,  6.4214,  6.4214,  6.4214,  6.4214,
                            6.2686,  5.772 ,  5.6574,  5.6574])
        ans = em.getPlasmaPause(self.ticks)
        np.testing.assert_allclose(real_ans, ans)

    def test_getPlasmaPause(self):
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

    def test_getLmax(self):
        """getLmax should give known results (regression)"""
        real_ans = np.array([ 7.9973023,  8.11663  ,  8.7714972,  8.52228  ,  8.6463423,
        8.6463423,  8.6048668,  8.7714972,  8.3179375,  8.8134583,
        9.0677743])
        ans = em.getLmax(self.ticks)
        np.testing.assert_allclose(real_ans, ans)
        self.assertRaises(ValueError, em.getLmax, self.ticks, model='bad')

    def test_getMPstandoff(self):
        """getMPstandoff should give known results (regression)"""
        real_ans = np.array([ 10.57319537,  10.91327764,  10.75086873,  10.77577207,
            9.78180261,  11.0374474 ,  11.4065    ,  11.27555451,
            11.47988573,  11.8202582 ,  11.23834814])
        ans = em.ShueMP(self.ticks)
        np.testing.assert_allclose(real_ans, ans)
        self.assertRaises(TypeError, em.ShueMP, 'bad')
        data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
        real_ans = np.array([ 9.96096838,  8.96790412])
        ans = em.ShueMP(data)
        np.testing.assert_allclose(real_ans, ans)

    def test_getDststar(self):
        """getDststar should give known results (regression)"""
        real_ans = np.array([-30.68169714, -26.85289053, -11.14932976, -17.77229149,
             -16.05975098, -13.07617265, -14.26      ,  -9.96354744,
            -21.11331813,  -8.49354146,  -3.18703339])
        ans = em.getDststar(self.ticks)
        np.testing.assert_allclose(real_ans, ans)
        self.assertRaises(ValueError, em.getDststar, self.ticks, model='bad')

class PAmodelTests(unittest.TestCase):
    def setUp(self):
        super(PAmodelTests, self).setUp()
        self.PA = tb.linspace(0.01,179.99,20000)

    def test_vampola_singleval(self):
        """sin^n model should have d_flux that integrates to omniflux"""
        omniflux = 3000
        dnflux, alphas = em.vampolaPA(omniflux, order=2, alpha=self.PA)
        d_sum = integ.simps(dnflux, np.deg2rad(alphas))
        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = integ.simps(dnflux, np.deg2rad(alphas))
        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)

    def test_vampola_len1list(self):
        """sin^n model should have d_flux that integrates to omniflux"""
        omniflux = [3000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        d_sum = integ.simps(dnflux, np.deg2rad(alphas))
        np.testing.assert_allclose(d_sum, omniflux, atol=0.001)

    def test_vampola_multival(self):
        """sin^n model should have d_flux that integrates to omniflux"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=4, alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = integ.simps(dnflux[:,i], np.deg2rad(alphas))
            np.testing.assert_allclose(d_sum, omniflux[i])

    def test_vampola_multi_n(self):
        """sin^n model should have d_flux that integrates to omniflux"""
        omniflux = [3000, 6000]
        dnflux, alphas = em.vampolaPA(omniflux, order=[2,4], alpha=self.PA)
        for i in range(len(omniflux)):
            d_sum = integ.simps(dnflux[:,i], np.deg2rad(alphas))
            np.testing.assert_allclose(d_sum, omniflux[i])

    def test_vampola_mismatched_order_len(self):
        """sin^n model should have d_flux that integrates to omniflux"""
        omniflux = [3000, 4500, 6000]
        self.assertRaises(ValueError, em.vampolaPA, omniflux, order=[2,4], alpha=self.PA)


if __name__ == "__main__":
    unittest.main()
