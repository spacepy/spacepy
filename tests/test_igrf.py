# -*- coding: utf-8 -*-

"""
Test suite for IGRF module

"""

import datetime as dt
import unittest
import warnings

import numpy

import spacepy_testing
from spacepy import igrf

__all__ = ['IgrfClassTests']


class IgrfClassTests(unittest.TestCase):
    def setUp(self):
        self.IGRF = igrf.IGRF()
        self.date1899 = dt.datetime(1899, 3, 17)
        self.date2010 = dt.datetime(2010, 1, 1)
        self.date2049 = dt.datetime(2049, 6, 21)

    def tearDown(self):
        del self.IGRF

    def test_initRaises(self):
        """IGRF init raises error on bad input"""
        self.assertRaises(TypeError, self.IGRF.initialize, 'Incorrect input')

    def test_initRaises_outofrange(self):
        """IGRF init raises error on time out of range"""
        self.assertRaises(ValueError, self.IGRF.initialize, self.date1899,
                          limits='error')
        self.assertRaises(ValueError, self.IGRF.initialize, self.date2049,
                          limits='error')

    def test_initWarns_outofrange(self):
        """IGRF init warns on time out of range"""
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            self.assertRaises(UserWarning, self.IGRF.initialize, self.date2049,
                              limits='warn')

    def test_loadCoefficients_default(self):
        """IGRF coefficients are loaded without error"""
        self.IGRF.initialize(self.date2010)
        self.assertTrue(hasattr(self.IGRF, '_IGRF__coeffs'))

    def test_CDlocation(self):
        """Test CD Northern pole location against SPENVIS docs

        See www.spenvis.oma.be/help/background/magfield/cd.html
        """
        epochs = [dt.datetime(1950, 1, 1), dt.datetime(1975, 1, 1),
                  dt.datetime(1985, 1, 1)]
        expected_lat = [78.47, 78.69, 78.97]
        expected_long = [291.15, 289.53, 289.10]
        for ep, la, lo in zip(epochs, expected_lat, expected_long):
            self.IGRF.initialize(ep)
            got_lon = (self.IGRF.dipole['cd_glon'] + 360) % 360
            got_lat = 90 - self.IGRF.dipole['cd_gcolat']
            numpy.testing.assert_almost_equal(got_lon, lo, decimal=2)
            numpy.testing.assert_almost_equal(got_lat, la, decimal=2)

    def test_IGRF13_SV_LGM(self):
        """Test projected dipole from IGRF14 calc using SV

        Values taken LANLGeoMag (revision 06116e71, PR#58)
        Test will need updating when IGRF15 or later is introduced
        """
        self.IGRF.initialize(dt.datetime(2027, 1, 19))
        numpy.testing.assert_allclose(self.IGRF.moment['cd'], 29700.2016, rtol=8e-7)
        self.IGRF.initialize(dt.datetime(2026, 9, 1, 12))
        numpy.testing.assert_allclose(self.IGRF.moment['cd'],  29706.3816, rtol=8e-7)
        numpy.testing.assert_almost_equal(self.IGRF.dipole['cd_gcolat'], 9.14254876, decimal=4)
        numpy.testing.assert_almost_equal(self.IGRF.dipole['cd_glon'], -72.8271716, decimal=4)

    def test_dipMoment_CCMC(self):
        """Another external validation, this time against CCMC Vitmo ModelWeb"""
        self.IGRF.initialize(dt.datetime(1999, 1, 1))
        # CCMC calculator only gives value to 1 decimal
        numpy.testing.assert_almost_equal(self.IGRF.moment['cd'], 30138.7, decimal=1)
        self.IGRF.initialize(dt.datetime(1988, 1, 1))
        numpy.testing.assert_almost_equal(self.IGRF.moment['cd'], 30364.8, decimal=1)

    def test_SVmoment_noDisco(self):
        """Test for no discontinuities across interpolate/extrapolate bound

        In increments of 10 minutes the dipole axis shouldn't change appreciably
        """
        epochs = [dt.datetime(2024, 12, 31, 23, 50), dt.datetime(2025, 1, 1),
                  dt.datetime(2025, 1, 1, 0, 10)]
        ans = []
        for ep in epochs:
            self.IGRF.initialize(ep)
            ans.append(self.IGRF.moment['cd'])
        center = numpy.repeat(ans[1], 3)
        numpy.testing.assert_allclose(ans, center, atol=1e-7)


if __name__ == "__main__":
    unittest.main()
