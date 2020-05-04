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
        exp = dict()
        psi = 18.344635  # degrees
        cd_moment = 29872.9290856547  # nT
        cd_gcolat = 9.707223  # (deg.)  (  09° 42′ 26″.002 )
        cd_glon = -72.584807  # (deg.)  ( -72° 35′ 05″.305 )
        self.CTrans2014.calcMagTransforms()
        numpy.testing.assert_allclose(cd_moment, self.CTrans2014['IGRF'].moment['cd'], rtol=1e-6)
        numpy.testing.assert_allclose(cd_gcolat, self.CTrans2014['IGRF'].dipole['cd_gcolat'],
                                      rtol=1e-5)
        numpy.testing.assert_allclose(cd_glon, self.CTrans2014['IGRF'].dipole['cd_glon'],
                                      rtol=1e-5)
        # numpy.testing.assert_allclose(psi, self.CTrans2014['DipoleTilt'], rtol=1e-5)

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
                        [     -0.91509195       0.36997560       0.16038943 ]
    Amod_to_gsm       = [     -0.36522551      -0.92903127       0.05925564 ]
                        [      0.17092994      -0.00435396       0.98527357 ]

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


'''
Going from GSM to GSE2000
Ugse2000  = -6.60000000 4.05439084 -0.64180598 Re
Difference between GSE and GSE2000
x:0.000000 y:-0.000020 z:-0.000128 mag:0.000130
and in km:      7.27596e-12 -0.129282 -0.816614 0.826784
Ang. diff. between LGM GSE and GSE2000 = 0.00095662 (3.44383)

The ground track point (geocentric):
Lat:-23.690776 Lon:97.414295


Time Quantitites:
    fYear (UTC)       = 1999.000394
    UTC               = 19990101  3.45055555555556 (  03ʰ 27ᵐ 02ˢ.00000000 )
    UT1               = 19990101  3.45055555555556 (  03ʰ 27ᵐ 02ˢ.00000000 )
    TAI               = 19990101  3.45944444444444 (  03ʰ 27ᵐ 34ˢ.00000000 )
    TT = TAI+32.184s  = 19990101  3.46838444444444 (  03ʰ 28ᵐ 06ˢ.18400000 )
    TDB               = 19990101  3.46838442348996 (  03ʰ 28ᵐ 06ˢ.18392456 )
    DUT1 = UT1-UTC    = 0.0000000 seconds
    DAT  = TAI-UTC    = 32.0000000 seconds
    JD (UTC)          = 2451179.643773148302
    JD (UT1)          = 2451179.643773148302
    JD (TT)           = 2451179.644516018685
    T  (UTC)          = -0.010002908333 Julian Centuries of UTC
    T  (UT1)          = -0.010002908333 Julian Centuries of UT1
    T  (TT)           = -0.010002887994 Julian Centuries of TT
    Year   (UTC)      = 1999
    Month  (UTC)      = 1
    Day    (UTC)      = 1
    Doy    (UTC)      = 1
    Dow    (UTC)      = 5
    Dowstr (UTC)      = Fri
    gmst (hours)      = -13.8595634 ( -13ʰ 51ᵐ 34ˢ.428 )
    gmst (degrees)    = -207.8934510 ( -207° 53′ 36″.423 )
    gast (hours)      = -13.8597291 ( -13ʰ 51ᵐ 35ˢ.025 )
    gast (degrees)    = -207.8959371 ( -207° 53′ 45″.374 )

Eccentricity and Obliquity:
    eccentricity                      = 0.01670953
    epsilon mean (obliq. of ecliptic) = 23.43942119 (  23° 26′ 21″.916 )
    epsilon true (obliq. of ecliptic) = 23.43715838 (  23° 26′ 13″.770 )

Precession Quantities:
    Zeta              = -0.006408 ( -00° 00′ 23″.069 )
    Zee               = -0.00640798 ( -00° 00′ 23″.069 )
    Theta             = -0.00556915 ( -00° 00′ 20″.049 )

Nutation Quantities:
    dPsi (w.o. corrections)           = -0.00271024 ( -00° 00′ 09″.757 )
    dEps (w.o. corrections)           = -0.00226281 ( -00° 00′ 08″.146 )
    ddPsi (EOP correction)            = 0.00000000 (  00° 00′ 00″.000 )
    ddEps (EOP correction)            = 0.00000000 (  00° 00′ 00″.000 )
    dPsi (w. corrections)             = -0.00271024 ( -00° 00′ 09″.757 )
    dEps (w. corrections)             = -0.00226281 ( -00° 00′ 08″.146 )
    epsilon true (obliq. of ecliptic) = 23.43715838 (  23° 26′ 13″.770 )
    Equation of the Equinox           = -0.00248618 ( -00° 00′ 08″.950 )

DE421 Position of Sun:
    lambda_sun      =      280.260712  (  280° 15′ 38″.561 )
    earth_sun_dist  =    23063.091904 Re
    beta_sun        =   -0.000174591   ( -00° 00′ 00″.629 )
    RA_sun  (MOD)   =      281.161332  (  18ʰ 44ᵐ 38ˢ.720 )
    DEC_sun (MOD)   =      -23.042917  ( -23° 02′ 34″.501 )
    RA_sun  (TOD)   =      281.158209  (  18ʰ 44ᵐ 37ˢ.970 )
    DEC_sun (TOD)   =      -23.040906  ( -23° 02′ 27″.260 )
    RA_sun  (J2000) =      281.176472  (  18ʰ 44ᵐ 42ˢ.353 )
    DEC_sun (J2000) =      -23.041838  ( -23° 02′ 30″.618 )

Sun vector and Ecliptic Pole in GEI2000:
    Sun               = (0.178128, -0.902807, -0.391421)
    EcPole            = (0.000000, -0.397779, 0.917481)

Geo-dipole tilt angle:
    psi                      = -32.810092  ( -32° 48′ 36″.332 )
    sin_psi                  = -0.541856
    cos_psi                  = 0.840471
    tan_psi                  = -0.644705

Position of Moon:
   RA_moon                      = 86.916668  (  05ʰ 47ᵐ 40ˢ.000 )
   DEC_moon                     = 19.151396  (  19° 09′ 05″.027 )
   EarthMoonDistance            = 57.253380
   MoonPhase                    = -9999999999999999635896294965248.000000

IGRF-derived quantities:
    M_cd              = 30138.6640162994 nT
    M_cd_McIlwain    = 31165.3000000000 nT
    M_cd_2010         = 29950.1686985232 nT
    CD_gcolat         = 10.500783 (deg.)  (  10° 30′ 02″.817 )
    CD_glon           = -71.538674 (deg.)  ( -71° 32′ 19″.228 )
    ED_x0             = -0.062943  Re  (-401.456816 km)
    ED_y0             = 0.046670  Re  (297.668292 km)
    ED_z0             = 0.031255  Re  (199.350760 km)

Transformation Matrices:
                        [      0.17812751      -0.90280705      -0.39142052 ]
    Amod_to_gse       = [      0.98400741       0.16342863       0.07085543 ]
                        [      0.00000054      -0.39778199       0.91747997 ]

                        [      0.17812751      -0.90280705      -0.39142052 ]
    Amod_to_gsm       = [      0.97244938       0.22229767      -0.07018517 ]
                        [      0.15037553      -0.36813473       0.91753148 ]

                        [     -0.88392373       0.46763108      -0.00008406 ]
    Agei_to_wgs84     = [     -0.46763108      -0.88392373      -0.00008917 ]
                        [     -0.00011601      -0.00003951       0.99999999 ]

                        [      0.17812751       0.98400741       0.00000054 ]
    Agse_to_mod       = [     -0.90280705       0.16342863      -0.39778199 ]
                        [     -0.39142052       0.07085543       0.91747997 ]

                        [      1.00000000       0.00000000      -0.00000000 ]
    Agse_to_gsm       = [     -0.00000000       0.98825420      -0.15281897 ]
                        [     -0.00000000       0.15281897       0.98825420 ]

                        [     -0.88392373      -0.46763108      -0.00011601 ]
    Awgs84_to_gei     = [      0.46763108      -0.88392373      -0.00003951 ]
                        [     -0.00008406      -0.00008917       0.99999999 ]

                        [      0.17812751       0.97244938       0.15037553 ]
    Agsm_to_mod       = [     -0.90280705       0.22229767      -0.36813473 ]
                        [     -0.39142052      -0.07018517       0.91753148 ]

                        [      0.84047117       0.00000000       0.54185626 ]
    Agsm_to_sm        = [      0.00000000       1.00000000       0.00000000 ]
                        [     -0.54185626       0.00000000       0.84047117 ]

                        [      1.00000000      -0.00000000      -0.00000000 ]
    Agsm_to_gse       = [      0.00000000       0.98825420       0.15281897 ]
                        [     -0.00000000      -0.15281897       0.98825420 ]

                        [      0.84047117       0.00000000      -0.54185626 ]
    Asm_to_gsm        = [      0.00000000       1.00000000       0.00000000 ]
                        [      0.54185626       0.00000000       0.84047117 ]

                        [      0.99999997       0.00022368       0.00009720 ]
    Agei_to_mod       = [     -0.00022368       0.99999997      -0.00000001 ]
                        [     -0.00009720      -0.00000001       1.00000000 ]

                        [      0.99999997      -0.00022368      -0.00009720 ]
    Amod_to_gei       = [      0.00022368       0.99999997      -0.00000001 ]
                        [      0.00009720      -0.00000001       1.00000000 ]

                        [      1.00000000       0.00004340       0.00001882 ]
    Amod_to_tod       = [     -0.00004340       1.00000000       0.00003949 ]
                        [     -0.00001881      -0.00000001       1.00000000 ]

                        [      1.00000000      -0.00004340      -0.00001881 ]
    Atod_to_mod       = [      0.00004340       1.00000000      -0.00003949 ]
                        [      0.00001882       0.00003949       1.00000000 ]

                        [     -0.88379881       0.46786714       0.00000000 ]
    Atod_to_pef       = [     -0.46786714      -0.88379881       0.00000000 ]
                        [      0.00000000       0.00000000       1.00000000 ]

                        [     -0.88379881      -0.46786714       0.00000000 ]
    Apef_to_tod       = [      0.46786714      -0.88379881       0.00000000 ]
                        [      0.00000000       0.00000000       1.00000000 ]

                        [ -8.83819110e-01   4.67828795e-01   0.00000000e+00 ]
    Ateme_to_pef      = [ -4.67828795e-01  -8.83819110e-01   0.00000000e+00 ]
                        [  0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [ -8.83819110e-01  -4.67828795e-01   0.00000000e+00 ]
    Apef_to_teme      = [  4.67828795e-01  -8.83819110e-01   0.00000000e+00 ]
                        [  0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [  1.00000000e+00   0.00000000e+00  -0.00000000e+00 ]
    Awgs84_to_pef     = [  0.00000000e+00   1.00000000e+00   0.00000000e+00 ]
                        [  0.00000000e+00  -0.00000000e+00   1.00000000e+00 ]

                        [  1.00000000e+00   0.00000000e+00   0.00000000e+00 ]
    Apef_to_wgs84     = [  0.00000000e+00   1.00000000e+00  -0.00000000e+00 ]
                        [ -0.00000000e+00   0.00000000e+00   1.00000000e+00 ]

                        [      0.17836749       0.98396394       0.00000094 ]
    Agse2000_to_gei   = [     -0.90276718       0.16364897      -0.39778189 ]
                        [     -0.39140319       0.07095051       0.91748001 ]

                        [      0.17836749      -0.90276718      -0.39140319 ]
    Agei_to_gse2000   = [      0.98396394       0.16364897       0.07095051 ]
                        [      0.00000094      -0.39778189       0.91748001 ]

Date = 19990101
UTC  = 3.450556
Ugsm = -6.60000000 3.40000000 -2.30000000 Re
Usm  = -6.79337914 3.40000000 1.64316763 Re
Umod = 1.78482258 7.56104847 0.23442341 Re

Going to GSE from SM and GSM
Ugse  = -6.60000000 3.00858065 -2.79256915 Re
Ugse  = -6.60000000 3.00858065 -2.79256915 Re
They are different by
x:-0.000000 y:0.000000 z:0.000000 mag:0.000000


Going from GSM to GSE2000
Ugse2000  = -6.60000000 3.00858234 -2.79256733 Re
Difference between GSE and GSE2000
x:-0.000000 y:-0.000002 z:-0.000002 mag:0.000002
and in km:      -7.27596e-12 -0.0107676 -0.0116005 0.0158276
Ang. diff. between LGM GSE and GSE2000 = 1.83114e-05 (0.0659211)

The ground track point (geocentric):
Lat:1.725913 Lon:-75.388382
'''

if __name__ == "__main__":
    unittest.main()
