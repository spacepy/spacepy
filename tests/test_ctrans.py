# -*- coding: utf-8 -*-

"""
Test suite for ctrans module

"""

import copy
import datetime as dt
import unittest

import numpy
from scipy.constants import arcsec

import spacepy_testing
from spacepy import ctrans
import spacepy.time

__all__ = ['CTransClassTests', 'CTransRegressionTests', 'ModuleTests']


class CTransClassTests(unittest.TestCase):
    def setUp(self):
        self.t2k = dt.datetime(2000, 1, 1)
        self.CTrans2000 = ctrans.CTrans(self.t2k)

    def tearDown(self):
        del self.CTrans2000

    def test_initRaises(self):
        """CTrans init raises error on bad input"""
        self.assertRaises(TypeError, ctrans.CTrans, 'Incorrect input')

    def test_initRaises_ephmodel(self):
        """CTrans init raises not implemented error for setting ephmodel"""
        self.assertRaises(NotImplementedError, ctrans.CTrans, self.t2k, ephmodel='BAD')

    def test_initRaises_pnmodel(self):
        """CTrans init raises not implemented error for setting ephmodel"""
        self.assertRaises(NotImplementedError, ctrans.CTrans, self.t2k, pnmodel='BAD')

    def test_getEOP_raises(self):
        """getEOP method raises an error for requesting non-zero EOP"""
        # If EOP support is added, this test can be removed or
        # changed to an expected failure
        self.assertRaises(NotImplementedError, self.CTrans2000.getEOP, useEOP=True)

    def test_getEOP_raises2(self):
        """getEOP method raises an error for requesting non-zero EOP"""
        # If EOP support is added, this test can be removed or
        # changed to an expected failure
        self.assertRaises(NotImplementedError, ctrans.CTrans, self.t2k, eop=True)

    def test_calcCore_times(self):
        """calcCoreTransforms should set up required times"""
        self.CTrans2000.calcCoreTransforms()
        self.assertTrue('TT_JD' in self.CTrans2000)
        self.assertTrue('GMST' in self.CTrans2000)

    def test_convert_badsys_in(self):
        """Calling convert with undefined systems should raise an error"""
        self.CTrans2000.calcCoreTransforms()
        self.assertRaises(ValueError, self.CTrans2000.convert,
                          [1, 1, 1], 'ECIMOD', 'UNDEFINED')

    def test_convert_badsys_out(self):
        """Calling convert with undefined systems should raise an error"""
        self.CTrans2000.calcCoreTransforms()
        self.assertRaises(ValueError, self.CTrans2000.convert,
                          [1, 1, 1], 'UNDEFINED', 'ECIMOD')

    def test_convert_badsys_force(self):
        """Calling convert with undefined systems should raise an error"""
        self.CTrans2000.calcCoreTransforms()
        self.CTrans2000['Transform']['UNK_ECITOD'] = numpy.identity(3)
        self.assertRaises(ValueError, self.CTrans2000.convert,
                          [1, 1, 1], 'UNK', 'ECIMOD')

    def test_time_as_ticktock(self):
        """Construction with datetime and ticktock should give same answer"""
        self.CTrans2000.calcTimes()
        tick2k = spacepy.time.Ticktock(self.t2k)
        tt_test = ctrans.CTrans(tick2k)
        tt_test.calcTimes()
        numpy.testing.assert_equal(tt_test['UTC_JD'], self.CTrans2000['UTC_JD'])

    def test_time_length(self):
        """CTrans should raise error with non-singular time (list of datetimes)"""
        self.assertRaises(TypeError, ctrans.CTrans, [self.t2k]*2)

    def test_time_length2(self):
        """CTrans should raise error with non-singular time (spacepy.Ticktock)"""
        self.assertRaises(ValueError, ctrans.CTrans, spacepy.time.Ticktock([self.t2k]*2))

    def test_time_list_ISO(self):
        """Test construction with ISO string"""
        self.CTrans2000.calcTimes()
        tt_test = ctrans.CTrans(['2000-01-01T00:00:00'])
        tt_test.calcTimes()
        numpy.testing.assert_equal(tt_test['UTC_JD'], self.CTrans2000['UTC_JD'])

    def test_time_recalc(self):
        """Multiple calls to calcTimes should preserve state"""
        self.CTrans2000.calcTimes()
        saved = copy.copy(self.CTrans2000['UTC'])
        self.CTrans2000.calcTimes()
        self.assertEqual(saved, self.CTrans2000['UTC'])

    def test_orbit_recalc(self):
        """Multiple calls to calcOrbitParams should preserve state"""
        self.CTrans2000.calcOrbitParams()
        saved = copy.copy(self.CTrans2000['constants'])
        self.CTrans2000.calcOrbitParams()
        self.assertEqual(saved, self.CTrans2000['constants'])

    def test_MagTransforms_calls_Core(self):
        """calcMagTransforms should call calcCoreTransforms if they aren't set"""
        self.CTrans2000.calcMagTransforms()
        self.assertIn('DipoleTilt_rad', self.CTrans2000)

    def test_GMST_versions_different(self):
        """GMST versions should not be identical"""
        ct2k_82 = ctrans.CTrans(self.t2k, pnmodel='IAU82')
        ct2k_00 = ctrans.CTrans(self.t2k, pnmodel='IAU00')
        # IAU82
        ct2k_82.gmst()
        exp1 = ct2k_82['GMST']
        # IAU00
        ct2k_00.gmst()
        exp2 = ct2k_00['GMST']
        self.assertNotAlmostEqual(exp1, exp2, places=10)

    def test_GMST_versions_similar(self):
        """GMST versions should all be roughly similar (not a regression test)"""
        ct2k_82 = ctrans.CTrans(self.t2k, pnmodel='IAU82')
        ct2k_00 = ctrans.CTrans(self.t2k, pnmodel='IAU00')
        # IAU82
        ct2k_82.gmst()
        exp1 = ct2k_82['GMST']
        # IAU00
        ct2k_00.gmst()
        exp2 = ct2k_00['GMST']
        numpy.testing.assert_approx_equal(exp1, exp2, significant=3)

    def test_convert_2D_same_as_1D(self):
        """2D input should give same answer as 1D input repeated"""
        self.CTrans2000.calcCoreTransforms()
        pos = [1.1, 2.2, 3.3]
        pos_1d = numpy.array(pos)
        pos_2d = numpy.array([pos, pos])
        got_1d = self.CTrans2000.convert(pos_1d, 'GSE', 'GEO')
        got_2d = self.CTrans2000.convert(pos_2d, 'GSE', 'GEO')
        numpy.testing.assert_allclose(got_1d, got_2d[0])
        numpy.testing.assert_allclose(got_1d, got_2d[1])

    def test_convert_GEO_RLL(self):
        """ """
        self.CTrans2000.calcCoreTransforms()
        got_eq = self.CTrans2000.convert([ctrans.WGS84['A'], 0, 0], 'GEO', 'RLL')
        exp_eq = [ctrans.WGS84['A'], 0, 0]
        got_po = self.CTrans2000.convert([0, 0, -ctrans.WGS84['A'] + 101], 'GEO', 'RLL')
        exp_po = [ctrans.WGS84['A'] - 101, -90, 0]
        numpy.testing.assert_allclose(exp_eq, got_eq)
        numpy.testing.assert_allclose(exp_po, got_po)

    def test_convert_RLL_GEO(self):
        """ """
        self.CTrans2000.calcCoreTransforms()
        got_eq = self.CTrans2000.convert([ctrans.WGS84['A'], 0, 0], 'RLL', 'GEO')
        exp_eq = [ctrans.WGS84['A'], 0, 0]
        got_po = self.CTrans2000.convert([ctrans.WGS84['A'] - 101, -90, 0], 'RLL', 'GEO')
        exp_po = [0, 0, -ctrans.WGS84['A'] + 101]
        numpy.testing.assert_allclose(exp_eq, got_eq, rtol=0, atol=1e-9)
        numpy.testing.assert_allclose(exp_po, got_po, rtol=0, atol=1e-9)


class CTransRegressionTests(unittest.TestCase):
    def setUp(self):
        self.CTrans2000 = ctrans.CTrans(dt.datetime(2000, 1, 1))
        self.CTrans2014 = ctrans.CTrans(dt.datetime(2014, 8, 29, 15, 32, 13, 811285),
                                        pnmodel='IAU82')

    def tearDown(self):
        del self.CTrans2000
        del self.CTrans2014

    def test_times2014_LGM(self):
        """Test that time conversions agree with LGM"""
        testdict = dict()
        testdict['UTC_JD'] = 2456899.147382074967
        testdict['UT1_JD'] = 2456899.147382074967
        testdict['TT_JD'] = 2456899.148159667850
        testdict['UT1_JC'] = 0.146588566244
        testdict['TT_JC'] = 0.146588587534
        self.CTrans2014.calcTimes()
        for key in testdict:
            numpy.testing.assert_almost_equal(testdict[key], self.CTrans2014[key], decimal=7)

    def test_gmst2014_LGM(self):
        """Test that GMST agrees with LGM"""
        self.CTrans2014.calcTimes()
        numpy.testing.assert_almost_equal(14.0546293, self.CTrans2014['GMST'], decimal=7)

    def test_gmst2014_astropy(self):
        """Test that GMST agrees with astropy (includes DUT1)"""
        ct14 = self.CTrans2014
        # getEOP method currently doesn't support EOP
        # to set, we need to call the factory method as the generated values are set
        # to be immutable
        ct14['EarthOrientationParameters'] = ct14._factory['eop'](DUT1=-0.32591566, xp=0,
                                                                  yp=0, ddPsi=0, ddEps=0)
        ct14._CTrans__status['useEOP'] = True
        self.CTrans2014.calcTimes()
        numpy.testing.assert_almost_equal(self.CTrans2014['GMST'], 14.05453852, decimal=7)

    def test_prec_orbit2014_LGM(self):
        """Test precession/orbital quantities agree with LGM"""
        testdict = dict()
        testdict['PrecessionZeta'] = 0.0939088
        testdict['PrecessionZee'] = 0.0939136
        testdict['PrecessionTheta'] = 0.0816111
        testdict['Eccentricity'] = 0.01670295
        self.CTrans2014.calcOrbitParams()
        for key in testdict:
            numpy.testing.assert_almost_equal(testdict[key], self.CTrans2014[key], decimal=7)

    def test_true_obliquity2014_LGM(self):
        """ """
        obl_true_rad = numpy.deg2rad(23.43715838)
        self.CTrans2014.calcOrbitParams()
        numpy.testing.assert_allclose(obl_true_rad,
                                      self.CTrans2014['ObliquityTrue_rad'],
                                      atol=5e-5)

    def test_mean_obliquity2014_LGM(self):
        """ """
        obl_mean_rad = numpy.deg2rad(23.43738485)
        self.CTrans2014.calcOrbitParams()
        numpy.testing.assert_allclose(obl_mean_rad,
                                      self.CTrans2014['ObliquityEcliptic_rad'],
                                      atol=5e-5)

    def test_dPsi2014_LGM(self):
        """ """
        dPsi = numpy.deg2rad(0.00197768)
        self.CTrans2014.calcOrbitParams()
        numpy.testing.assert_allclose(self.CTrans2014['dPsi']*arcsec, dPsi, atol=1e-7)

    def test_dEps2014_LGM(self):
        """ """
        dEps = numpy.deg2rad(-0.00227433)
        self.CTrans2014.calcOrbitParams()
        numpy.testing.assert_allclose(self.CTrans2014['dEps']*arcsec, dEps, atol=1e-7)

    def test_coreTransforms_MOD_J2000_2014_LGM(self):
        """ """
        exp = dict()
        exp['ECIMOD_GSE'] = [[-0.91509195, 0.36997560, 0.16038943],
                             [-0.40324523, -0.83959256, -0.36397474],
                             [0.00000000, -0.39774663, 0.91749530]]
        exp['ECI2000_ECIMOD'] = [[0.99999361, -0.00327811, -0.00142438],
                                 [0.00327811, 0.99999463, -0.00000233],
                                 [0.00142438, -0.00000233, 0.99999899]]
        self.CTrans2014.calcCoreTransforms()
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(got[tra], exp[tra], atol=1e-7)

    def test_coreTransforms_TOD_2014_LGM(self):
        """ """
        exp = dict()
        exp['ECIMOD_ECITOD'] = [[1.00000000, -0.00003167, -0.00001373],
                                [0.00003167, 1.00000000, 0.00003969],
                                [0.00001373, -0.00000233, 1.00000000]]
        self.CTrans2014.calcCoreTransforms()
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(got[tra], exp[tra], atol=5e-5)

    def test_coreTransforms_TEME_2019_astropy(self):
        """Compare TEME to GEO conversion to expected astropy result

        See example at https://docs.astropy.org/en/stable/coordinates/satellites.htm
        """
        tt19 = spacepy.time.Ticktock('2019-12-09T20:42:09.072')
        ct19 = ctrans.CTrans(tt19)
        # UTC-UT1 is required to get these to agree to tolerance
        ct19['EarthOrientationParameters'] = ct19._factory['eop'](DUT1=-0.17249231, xp=0,
                                                                  yp=0, ddPsi=0, ddEps=0)
        ct19._CTrans__status['useEOP'] = True
        ct19.calcCoreTransforms()
        uteme = numpy.array([-6102.443276428913, -986.3320160861297, -2820.3130707199225])  # kilometers
        expect_gdz = numpy.array([420.17927591, -24.6609379, 160.34199789])  # altitude, lat, long
        gotgdz = ct19.convert(uteme, 'TEME', 'GDZ')
        numpy.testing.assert_almost_equal(gotgdz[0], expect_gdz[0], decimal=4)  # alt same to cm
        numpy.testing.assert_allclose(gotgdz[1:], expect_gdz[1:], rtol=2.78e-5)  # lat/lon same to 0.1 arcsec

    def test_dipoleValues_LGM(self):
        """ """
        cd_moment = 29872.9290856547  # nT
        cd_gcolat = 9.707223  # (deg.)  (  09° 42′ 26″.002 )
        cd_glon = -72.584807  # (deg.)  ( -72° 35′ 05″.305 )
        self.CTrans2014.calcMagTransforms()
        numpy.testing.assert_allclose(cd_moment, self.CTrans2014['IGRF'].moment['cd'], rtol=1e-6)
        numpy.testing.assert_allclose(cd_gcolat, self.CTrans2014['IGRF'].dipole['cd_gcolat'],
                                      rtol=1e-5)
        numpy.testing.assert_allclose(cd_glon, self.CTrans2014['IGRF'].dipole['cd_glon'],
                                      rtol=1e-5)

    def test_dipoleTilt2014_LGM(self):
        """ """
        psi = 18.344635  # degrees
        self.CTrans2014.calcMagTransforms()
        numpy.testing.assert_allclose(self.CTrans2014['DipoleTilt'], psi, atol=1e-5)

    def test_magTransforms_MOD_GSM_2014_LGM(self):
        """ """
        exp = dict()
        exp['GSM_ECIMOD'] = [[-0.91509195, -0.36522551, 0.17092994],
                             [0.36997560, -0.92903127, -0.00435396],
                             [0.16038943, 0.05925564, 0.98527357]]
        self.CTrans2014.calcMagTransforms()
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(got[tra], exp[tra], atol=1e-7)

    def test_convert_viaMOD2014_LGM(self):
        """Test that convert constructs transformation via MOD correctly"""
        ugsm = numpy.array([-6.6, 3.4, -2.3])
        ugse = numpy.array([-6.6, 4.05437055, -0.64193415])
        self.CTrans2014.calcMagTransforms()
        gotgse = self.CTrans2014.convert(ugsm, 'GSM', 'GSE')
        numpy.testing.assert_almost_equal(gotgse, ugse, decimal=6)

    def test_convert_ECI2000_GEO_2014_LGM(self):
        """Test conversion to GEO (compare to LGM)"""
        exp = dict()
        exp['ECI2000_GEO'] = [[-0.86046014, -0.50951631, 0.00121589],
                              [0.50951573, -0.86046100, -0.00076888],
                              [0.00143798, -0.00004207, 0.99999897]]
        self.CTrans2014.calcCoreTransforms()
        self.CTrans2014.convert([2, 2, 2], 'ECI2000', 'GEO')
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(got[tra], exp[tra], atol=5e-5)

    def test_geodetic_func_1D(self):
        """Test GEO to geodetic with 1D input (regression test against IRBEM)"""
        pos_geo = [-6374.2283916, -869.65264893, -4263.36088969]
        expected = [1346.134508626895, -33.6790290525164, -172.2309539745096]
        got = ctrans.geo_to_gdz(pos_geo)
        numpy.testing.assert_allclose(got, expected, rtol=1e-5)

    def test_convert_GEO_GDZ_IRBEM(self):
        """Test GEO to geodetic using CTrans convert (regression test against IRBEM)"""
        tt = spacepy.time.Ticktock(2459218.5, 'JD')
        pos_geo = [-6374.2283916, -869.65264893, -4263.36088969]
        expected = [1346.134508626895, -33.6790290525164, -172.2309539745096]
        ct = ctrans.CTrans(tt)
        got = ct.convert(pos_geo, 'GEO', 'GDZ')
        numpy.testing.assert_allclose(got, expected, rtol=1e-5)

    def test_convert_ECITOD_GDZ_IRBEM(self):
        """Test inertial to geodetic using CTrans convert (regression test against IRBEM)"""
        tt = spacepy.time.Ticktock(2459218.5, 'JD')
        pos_eci = [2367.83158, -5981.75882, -4263.24591]
        expected = [1346.134508626895, -33.6790290525164, -172.2309539745096]
        ct = ctrans.CTrans(tt)
        ct.calcCoreTransforms()
        got = ct.convert(pos_eci, 'ECITOD', 'GDZ')
        numpy.testing.assert_allclose(got, expected, rtol=5e-5)

    def test_geodetic_func_2D(self):
        """Test GEO to geodetic with 2D input (regression test against IRBEM)"""
        pos_geo1 = [-6374.2283916, -869.65264893, -4263.36088969]
        pos_geo = [pos_geo1, pos_geo1]
        expected1 = [1346.134508626895, -33.6790290525164, -172.2309539745096]
        expected = [expected1, expected1]
        got = ctrans.geo_to_gdz(pos_geo)
        numpy.testing.assert_allclose(got, expected, rtol=1e-5)


class ModuleTests(unittest.TestCase):
    def test_geodetic_roundtrip_1D(self):
        """Test GEO->geodetic round trip with 1D input"""
        pos_geo = [-6374.2283916, -869.65264893, -4263.36088969]
        gdz = ctrans.geo_to_gdz(pos_geo)
        got = ctrans.gdz_to_geo(gdz)
        numpy.testing.assert_allclose(got, pos_geo, rtol=1e-9)

    def test_geodetic_1D_Re_and_km(self):
        """Test GEO->geodetic round trip with 1D Re input"""
        pos_geo = numpy.array([-6374.2283916, -869.65264893, -4263.36088969])
        gdz_re = ctrans.geo_to_gdz(pos_geo/ctrans.WGS84['A'], units='Re')
        gdz_km = ctrans.geo_to_gdz(pos_geo, units='km')
        numpy.testing.assert_allclose(gdz_re[1:], gdz_km[1:], rtol=1e-9)
        numpy.testing.assert_almost_equal(gdz_re[0]*ctrans.WGS84['A'], gdz_km[0], decimal=7)

    def test_geo_gdz_longitude(self):
        """Geodetic longitude should be same as geographic"""
        geo = numpy.array([7000, 4500, 5500])
        gdz = ctrans.geo_to_gdz(geo)
        gdzlon = gdz[2]
        geolon = numpy.rad2deg(numpy.arctan2(geo[1], geo[0]))
        numpy.testing.assert_allclose(gdzlon, geolon, rtol=1e-9)

    def test_gdz_geo_longitude(self):
        """Geographic longitude should be same as geodetic"""
        gdz = numpy.array([100, 45, 64])
        geo = ctrans.gdz_to_geo(gdz)
        gdzlon = gdz[2]
        geolon = numpy.rad2deg(numpy.arctan2(geo[1], geo[0]))
        numpy.testing.assert_allclose(geolon, gdzlon, rtol=1e-7)

    def test_geo_gdz_equator(self):
        """Geodetic conversion at equator gives expected result"""
        geo = numpy.array([ctrans.WGS84['A'], 0, 0])
        gdz = ctrans.geo_to_gdz(geo)
        expected = [0, 0, 0]
        numpy.testing.assert_allclose(gdz, expected, atol=1e-9, rtol=0)

    def test_geo_gdz_pole(self):
        """Geodetic conversion at pole gives expected result"""
        geo = numpy.array([0, 0, ctrans.WGS84['B']])
        gdz = ctrans.geo_to_gdz(geo)
        expected = [0, 90, 0]
        numpy.testing.assert_allclose(gdz, expected, atol=1e-9, rtol=0)

    def test_geo_gdz_pole_S(self):
        """Geodetic conversion at south pole gives expected result"""
        alt = 100
        # x, y, z - so altitude is all in z
        sp = [0, 0, -(ctrans.WGS84['B'] + alt)]
        np = [0, 0, ctrans.WGS84['B']]
        geo = numpy.array([sp, np, sp])
        gdz = ctrans.geo_to_gdz(geo)
        exp_s = [alt, -90, 0]
        exp_n = [0, 90, 0]
        expected = [exp_s, exp_n, exp_s]
        numpy.testing.assert_allclose(gdz, expected, atol=1e-9, rtol=0)

    def test_gdz_geo_equator(self):
        """Geodetic/geo conversion at equator gives expected result"""
        gdz = [0, 0, 0]
        expected = numpy.array([ctrans.WGS84['A'], 0, 0])
        got = ctrans.gdz_to_geo(gdz)
        numpy.testing.assert_allclose(got, expected, atol=1e-9, rtol=0)

    def test_gdz_geo_pole(self):
        """Geodetic/geo conversion at pole gives expected result"""
        gdz = [0, 90, 0]
        expected = numpy.array([0, 0, ctrans.WGS84['B']])
        got = ctrans.gdz_to_geo(gdz)
        numpy.testing.assert_allclose(got, expected, atol=1e-9, rtol=0)

    def test_geo_rll_equator(self):
        """geo/RLL conversion at equator gives expected result"""
        geo = [7500, 0, 0]
        expected = numpy.array([7500, 0, 0])
        got = ctrans.geo_to_rll(geo)
        numpy.testing.assert_allclose(got, expected, atol=1e-9, rtol=0)

    def test_rll_geo_equator(self):
        """RLL/geo conversion at equator gives expected result"""
        alt = 100
        rll = [alt + ctrans.WGS84['A'], 0, 0]
        expected = numpy.array([alt + ctrans.WGS84['A'], 0, 0])
        got = ctrans.rll_to_geo(rll)
        numpy.testing.assert_allclose(got, expected, atol=1e-9, rtol=0)

    def test_rll_geo_roundtrip_pole(self):
        """RLL/geo roundtrip at poles gives inital result"""
        alt = 100
        # test cases have negative in sqrt term on geo_to_gdz conversion
        rll = [[alt + ctrans.WGS84['A'], 90, 0], [alt + ctrans.WGS84['B'], -90, 0]]
        stage1 = ctrans.rll_to_geo(rll)
        got = ctrans.geo_to_rll(stage1)
        numpy.testing.assert_allclose(got, rll, atol=1e-9, rtol=0)


if __name__ == "__main__":
    unittest.main()
