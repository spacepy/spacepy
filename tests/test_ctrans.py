# -*- coding: utf-8 -*-

"""
Test suite for ctrans module

"""

import datetime as dt
import unittest
import warnings

import numpy
from scipy.constants import arcsec
from spacepy import ctrans

__all__ = ['CTransClassTests']


class CTransClassTests(unittest.TestCase):
    def setUp(self):
        self.CTrans2000 = ctrans.CTrans(dt.datetime(2000, 1, 1))
        self.CTrans2014 = ctrans.CTrans(dt.datetime(2014, 8, 29, 15, 32, 13, 811285),
                                        pnmodel='IAU82')

    def tearDown(self):
        del self.CTrans2000
        del self.CTrans2014

    def test_initRaises(self):
        """CTrans init raises error on bad input"""
        self.assertRaises(TypeError, ctrans.CTrans, 'Incorrect input')

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

    def test_gmst2014_LMG_IAU(self):
        """Test that LGM GMST formulation agrees with alternate IAU form"""
        ctlgm = ctrans.CTrans(dt.datetime(2014, 8, 29, 15, 32, 13, 811285),
                              pnmodel='LGMDEFAULT')
        ctlgm.calcTimes()
        self.CTrans2014.calcTimes()
        numpy.testing.assert_almost_equal(ctlgm['GMST'], self.CTrans2014['GMST'], decimal=7)

    def test_gmst2014_astropy(self):
        """Test that GMST agrees with astropy (includes DUT1)"""
        ct14 = self.CTrans2014
        ct14['EarthOrientationParameters'] = ct14._factory['eop'](DUT1=-0.32591566, xp=0,
                                                                  yp=0, ddPsi=0, ddEps=0)
        self.CTrans2014.calcTimes()
        numpy.testing.assert_almost_equal(14.05453852, self.CTrans2014['GMST'], decimal=7)

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
        numpy.testing.assert_allclose(dPsi, self.CTrans2014['dPsi']*arcsec, atol=1e-7)

    def test_dEps2014_LGM(self):
        """ """
        dEps = numpy.deg2rad(-0.00227433)
        self.CTrans2014.calcOrbitParams()
        numpy.testing.assert_allclose(dEps, self.CTrans2014['dEps']*arcsec, atol=1e-7)

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
            numpy.testing.assert_allclose(exp[tra], got[tra], atol=1e-7)

    def test_coreTransforms_TOD_2014_LGM(self):
        """ """
        exp = dict()
        exp['ECIMOD_ECITOD'] = [[1.00000000, -0.00003167, -0.00001373],
                                [0.00003167, 1.00000000, 0.00003969],
                                [0.00001373, -0.00000233, 1.00000000]]
        self.CTrans2014.calcCoreTransforms()
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(exp[tra], got[tra], atol=5e-5)

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
        numpy.testing.assert_allclose(psi, self.CTrans2014['DipoleTilt'], atol=1e-5)

    def test_magTransforms_MOD_GSM_2014_LGM(self):
        """ """
        exp = dict()
        exp['GSM_ECIMOD'] = [[-0.91509195, -0.36522551, 0.17092994],
                             [0.36997560, -0.92903127, -0.00435396],
                             [0.16038943, 0.05925564, 0.98527357]]
        self.CTrans2014.calcMagTransforms()
        got = self.CTrans2014['Transform']
        for tra in exp:
            numpy.testing.assert_allclose(exp[tra], got[tra], atol=1e-7)

    def test_convert_viaMOD2014_LGM(self):
        """Test that convert constructs transformation via MOD correctly"""
        ugsm = numpy.array([-6.6, 3.4, -2.3])
        ugse = numpy.array([-6.6, 4.05437055, -0.64193415])
        self.CTrans2014.calcMagTransforms()
        gotgse = self.CTrans2014.convert(ugsm, 'GSM', 'GSE')
        numpy.testing.assert_almost_equal(ugse, gotgse, decimal=6)

'''
Time Quantitites:
    fYear (UTC)       = 2014.659308
    UTC               = 20140829  15.53716980138889 (  15ʰ 32ᵐ 13ˢ.81128500 )
    TT = TAI+32.184s  = 20140829  15.55583202361111 (  15ʰ 33ᵐ 20ˢ.99528500 )
    TDB               = 20140829  15.55583165196660 (  15ʰ 33ᵐ 20ˢ.99394708 )
    DUT1 = UT1-UTC    = 0.0000000 seconds
    DAT  = TAI-UTC    = 35.0000000 seconds
    gmst (hours)      = 14.0546293 (  14ʰ 03ᵐ 16ˢ.665 )
    gmst (degrees)    = 210.8194396 (  210° 49′ 09″.982 )
    gast (hours)      = 14.0547503 (  14ʰ 03ᵐ 17ˢ.101 )
    gast (degrees)    = 210.8212538 (  210° 49′ 16″.514 )
Eccentricity and Obliquity:
    eccentricity                      = 0.01670295
    epsilon mean (obliq. of ecliptic) = 23.43738485 (  23° 26′ 14″.585 )
    epsilon true (obliq. of ecliptic) = 23.43511052 (  23° 26′ 06″.398 )
Precession Quantities:
    Zeta              = 0.0939088 (  00° 05′ 38″.072 )
    Zee               = 0.0939136 (  00° 05′ 38″.089 )
    Theta             = 0.0816111 (  00° 04′ 53″.800 )
Nutation Quantities:
    dPsi (w.o. corrections)           = 0.00197768 (  00° 00′ 07″.120 )
    dEps (w.o. corrections)           = -0.00227433 ( -00° 00′ 08″.188 )
    ddPsi (EOP correction)            = 0.00000000 (  00° 00′ 00″.000 )
    ddEps (EOP correction)            = 0.00000000 (  00° 00′ 00″.000 )
    dPsi (w. corrections)             = 0.00197768 (  00° 00′ 07″.120 )
    dEps (w. corrections)             = -0.00227433 ( -00° 00′ 08″.188 )
    epsilon true (obliq. of ecliptic) = 23.43511052 (  23° 26′ 06″.398 )
    Equation of the Equinox           = 0.00181425 (  00° 00′ 06″.531 )
Low Accuracy Position of Sun:
    lambda_sun      =      156.218789  (  156° 13′ 07″.641 )
    earth_sun_dist  =    23683.928241 Re
Sun vector and Ecliptic Pole in GEI2000:
    Sun               = (-0.915092, 0.369976, 0.160389)
    EcPole            = (0.000000, -0.397747, 0.917495)
Geo-dipole tilt angle:
    psi                      = 18.344635  (  18° 20′ 40″.684 )
    sin_psi                  = 0.314732
    cos_psi                  = 0.949181
    tan_psi                  = 0.331583
IGRF-derived quantities:
    M_cd              = 29872.9290856547 nT
    M_cd_McIlwain    = 31165.3000000000 nT
    M_cd_2010         = 29950.1686985232 nT
    CD_gcolat         = 9.707223 (deg.)  (  09° 42′ 26″.002 )
    CD_glon           = -72.584807 (deg.)  ( -72° 35′ 05″.305 )
    ED_x0             = -0.062767  Re  (-400.334621 km)
    ED_y0             = 0.055013  Re  (350.879290 km)
    ED_z0             = 0.034686  Re  (221.232022 km)
Transformation Matrices:

                        [     -0.86046014      -0.50951631       0.00121589 ]
    Agei_to_wgs84     = [      0.50951573      -0.86046100      -0.00076888 ]
                        [      0.00143798      -0.00004207       0.99999897 ]

                        [     -0.91509195      -0.40324523       0.00000000 ]
    Agse_to_mod       = [      0.36997560      -0.83959256      -0.39774663 ]
                        [      0.16038943      -0.36397474       0.91749530 ]

                        [      1.00000000       0.00000000      -0.00000000 ]
    Agse_to_gsm       = [     -0.00000000       0.90571563       0.42388582 ]
                        [      0.00000000      -0.42388582       0.90571563 ]

                        [     -0.86046014       0.50951573       0.00143798 ]
    Awgs84_to_gei     = [     -0.50951631      -0.86046100      -0.00004207 ]
                        [      0.00121589      -0.00076888       0.99999897 ]

                        [     -0.91509195      -0.36522551       0.17092994 ]
    Agsm_to_mod       = [      0.36997560      -0.92903127      -0.00435396 ]
                        [      0.16038943       0.05925564       0.98527357 ]

                        [      0.94918058       0.00000000      -0.31473198 ]
    Agsm_to_sm        = [      0.00000000       1.00000000       0.00000000 ]
                        [      0.31473198       0.00000000       0.94918058 ]

                        [      1.00000000      -0.00000000       0.00000000 ]
    Agsm_to_gse       = [      0.00000000       0.90571563      -0.42388582 ]
                        [     -0.00000000       0.42388582       0.90571563 ]

                        [      0.94918058       0.00000000       0.31473198 ]
    Asm_to_gsm        = [      0.00000000       1.00000000       0.00000000 ]
                        [     -0.31473198       0.00000000       0.94918058 ]

                        [      0.99999361      -0.00327811      -0.00142438 ]
    Agei_to_mod       = [      0.00327811       0.99999463      -0.00000233 ]
                        [      0.00142438      -0.00000233       0.99999899 ]

                        [      0.99999361       0.00327811       0.00142438 ]
    Amod_to_gei       = [     -0.00327811       0.99999463      -0.00000233 ]
                        [     -0.00142438      -0.00000233       0.99999899 ]

                        [      1.00000000      -0.00003167      -0.00001373 ]
    Amod_to_tod       = [      0.00003167       1.00000000       0.00003969 ]
                        [      0.00001373      -0.00000233       1.00000000 ]

                        [      1.00000000       0.00003167       0.00001373 ]
    Atod_to_mod       = [     -0.00003167       1.00000000      -0.00003969 ]
                        [     -0.00001373       0.00003969       1.00000000 ]

                        [     -0.85876990      -0.51236146       0.00000000 ]
    Atod_to_pef       = [      0.51236146      -0.85876990       0.00000000 ]
                        [      0.00000000       0.00000000       1.00000000 ]

                        [     -0.85876990       0.51236146       0.00000000 ]
    Apef_to_tod       = [     -0.51236146      -0.85876990       0.00000000 ]
                        [      0.00000000       0.00000000       1.00000000 ]

                        [ -8.58786119e-01  -5.12334267e-01   0.00000000e+00 ]
    Ateme_to_pef      = [  5.12334267e-01  -8.58786119e-01   0.00000000e+00 ]
                        [  0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [ -8.58786119e-01   5.12334267e-01   0.00000000e+00 ]
    Apef_to_teme      = [ -5.12334267e-01  -8.58786119e-01   0.00000000e+00 ]
                        [  0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [  1.00000000e+00   0.00000000e+00  -0.00000000e+00 ]
    Awgs84_to_pef     = [  0.00000000e+00   1.00000000e+00   0.00000000e+00 ]
                        [  0.00000000e+00  -0.00000000e+00   1.00000000e+00 ]

                        [  1.00000000e+00   0.00000000e+00   0.00000000e+00 ]
    Apef_to_wgs84     = [  0.00000000e+00   1.00000000e+00  -0.00000000e+00 ]
                        [ -0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [     -0.91364483      -0.40651337      -0.00000985 ]
    Agse2000_to_gei   = [      0.37297301      -0.83825275      -0.39777314 ]
                        [      0.16169184      -0.36342704       0.91748381 ]

                        [     -0.91364483       0.37297301       0.16169184 ]
    Agei_to_gse2000   = [     -0.40651337      -0.83825275      -0.36342704 ]
                        [     -0.00000985      -0.39777314       0.91748381 ]

Date = 20140829
UTC  = 15.537170
Ugsm = -6.60000000 3.40000000 -2.30000000 Re
Usm  = -5.54070829 3.40000000 -4.26034642 Re
Umod = 4.40470133 -5.59053118 -3.12323028 Re

Going to GSE from SM and GSM
Ugse  = -6.60000000 4.05437055 -0.64193415 Re
Ugse  = -6.60000000 4.05437055 -0.64193415 Re
They are different by
x:0.000000 y:0.000000 z:0.000000 mag:0.000000
'''

if __name__ == "__main__":
    unittest.main()
