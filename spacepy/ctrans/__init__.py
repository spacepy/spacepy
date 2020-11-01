'''
CTrans: Module for backend coordinate transformations in SpacePy
'''
import datetime as dt
import collections
from math import fmod

import numpy as np
import scipy.constants

from spacepy import datamodel as dm
from spacepy import time as spt
from spacepy import igrf


class CTrans(dm.SpaceData):
    def __init__(self, ctime, ephmodel=None, pnmodel=None, eop=False):
        """
        """
        super(CTrans, self).__init__()
        if ephmodel is not None:
            self._raiseErr(NotImplementedError, 'ephmodel')
        else:
            self.attrs['ephmodel'] = 'LGMDEFAULT'
        if pnmodel is not None:
            if pnmodel not in ['LGMDEFAULT', 'IAU82', 'IAU00']:
                self._raiseErr(NotImplementedError, 'pnmodel')
            self.attrs['pnmodel'] = pnmodel
        else:
            self.attrs['pnmodel'] = 'IAU82'

        if isinstance(ctime, (spt.Ticktock)):
            try:
                if len(ctime) == 1:
                    ctime = ctime[0]
                else:  # Input time is Ticktock, but has a length > 1
                    self._raiseErr(ValueError, 'time_in')
            except:
                pass  # Nothing to do here, Ticktock is what we want
        elif isinstance(ctime, (dt.datetime)):
            # Input time is datetime
            ctime = spt.Ticktock(ctime, dtype='UTC')
        else:
            try:
                if len(ctime) == 1:
                    ctime = ctime[0]
                    ctime = spt.Ticktock(ctime)  # Guess dtype
                else:
                    self._raiseErr(TypeError, 3)
            except TypeError:
                self._raiseErr(TypeError, 3)
        self.attrs['time'] = ctime

        # Make key information immutable, but referenced by name
        self._factory = dict()
        self._factory['constants'] = collections.namedtuple('Constants',
                                                            'twopi, j2000_jd, j1900_jd, ' +
                                                            'j1990_jd, daysec, daycentury, ' +
                                                            'arcsec, AU, Re',
                                                            )
        self._factory['eop'] = collections.namedtuple('EarthOrientationParameters',
                                                      'DUT1, xp, yp, ddPsi, ddEps',
                                                      )
        self._setconstants()
        self.getEOP(useEOP=eop)
        self.__status = {'time': False,
                         'orbital': False,
                         'transformCore': False,
                         'transformMag': False,
                         'timeInProgress': False,
                         'coreInProgress': False,
                         }

    def _setconstants(self):
        self['constants'] = self._factory['constants'](twopi=scipy.constants.pi*2,
                                                       j2000_jd=2451545.0,
                                                       j1990_jd=2447891.5,
                                                       j1900_jd=2415020.0,
                                                       daysec=86400,
                                                       daycentury=36525.0,
                                                       arcsec=scipy.constants.arcsec,
                                                       AU=scipy.constants.au/1e3,  # 1AU in km
                                                       Re=6378.137,  # WGS84_A radius in km
                                                       )

    def getEOP(self, useEOP=False):
        """Get/set Earth Orientation Parameters"""
        if not useEOP:
            # DUT1 in seconds
            # xp, yp in arcseconds
            # ddPsi and ddEps in arcseconds
            self['EarthOrientationParameters'] = self._factory['eop'](DUT1=0,
                                                                      xp=0,
                                                                      yp=0,
                                                                      ddPsi=0,
                                                                      ddEps=0)
        else:
            self._raiseErr(NotImplementedError, 'eop')

    def calcTimes(self, recalc=False, **kwargs):
        """Calculate time in systems required to set up coordinate transforms"""
        if self.__status['time'] and not recalc:
            return
        # Set up various time systems required and variable shortcuts
        eop = self['EarthOrientationParameters']
        const = self['constants']
        ctime = self.attrs['time']
        nleaps = ctime.getleapsecs()[0]

        # UT1
        self['UTC'] = ctime
        self['UTC_JD'] = ctime.JD[0]
        self['UTC_JC'] = (self['UTC_JD'] - const.j2000_jd)/const.daycentury
        self['TAI'] = spt.Ticktock(ctime.UTC[0]+dt.timedelta(seconds=float(nleaps)), dtype='UTC')
        self['TAI_JD'] = self['TAI'].JD[0]
        self['UT1'] = spt.Ticktock(ctime.UTC[0] + dt.timedelta(seconds=eop.DUT1), dtype='UTC')
        self['UT1_JD'] = self['UT1'].JD[0]
        self['UT1_JC'] = (self['UT1_JD'] - const.j2000_jd)/const.daycentury

        # Terrestrial Time in Julian centuries since J2000
        self['TT_JD'] = self['TAI_JD'] + 32.184/const.daysec
        self['TT_JC'] = (self['TT_JD'] - const.j2000_jd)/const.daycentury

        # Greenwich Mean Sidereal Time
        if 'from_gmst' not in kwargs:
            self.__status['timeInProgress'] = True
            self.gmst()
            self.__status['timeInProgress'] = False
        self.__status['time'] = True

    def calcOrbitParams(self, recalc=False):
        """Calculate Earth orbit parameters needed for coordinate transforms"""
        if not (self.__status['time'] or recalc):
            self.calcTimes(recalc=recalc)
        if self.__status['orbital'] and not recalc:
            # Don't recalculate unless it's asked for
            return
        const = self['constants']

        # Julian centuries with 1900 reference epoch
        self['JulianCenturies'] = (self['TT_JD'] - const.j1900_jd)/const.daycentury
        TU = self['JulianCenturies']  # since 1900
        TU2 = TU*TU

        # Orbital parameters required for later use
        # self['MeanEclipticLongitude'] = np.deg2rad(279.6966778 + 36000.76892*TU
        #                                            + 0.0003025*TU2)
        # self['EclipticLongitudePerigee'] = np.deg2rad(281.2208444 + 1.719175*TU
        #                                               + 0.000452778*TU2)
        self['Eccentricity'] = 0.01675104 - 0.0000418*TU - 0.000000126*TU2
        TT = self['TT_JC']  # in Julian century (J2000) format
        TT2 = TT*TT
        self['ObliquityEcliptic'] = 84381.448 - 46.8150*TT\
            - 0.00059*TT2 + 0.001813*TT2*TT
        self['ObliquityEcliptic_rad'] = self['ObliquityEcliptic']*const.arcsec
        self['ObliquityEcliptic'] = self['ObliquityEcliptic']

        # Anomalies
        days_since_1990 = self['TT_JD'] - const.j1990_jd
        varep90 = np.deg2rad(279.403303)
        varpi90 = np.deg2rad(282.768422)
        M = fmod((const.twopi/365.242191*days_since_1990), const.twopi)
        self['MeanAnomaly'] = fmod((M + varep90 - varpi90), const.twopi)
        ecc_anom = self._kepler(self['MeanAnomaly'], self['Eccentricity'])
        self['TrueAnomaly'] = 2.0*np.arctan(np.sqrt((1.0+self['Eccentricity'])
                                            / (1.0-self['Eccentricity']))*np.tan(ecc_anom/2.0))
        self['LambdaSun'] = fmod((self['TrueAnomaly'] + varpi90), const.twopi)
        self['OmegaMoon'] = fmod((125.04455501 - (5*360 + 134.1361851)*TT
                                 + 0.0020756*TT2 + 2.139e-6*TT2*TT), 360.0)  # degrees

        # Precession angles in degrees
        self['PrecessionZeta'] = 0.6406161*TT + 0.0000839*TT2 + 5.0e-6*TT2*TT
        self['PrecessionZee'] = 0.6406161*TT + 0.0003041*TT2 + 5.1e-6*TT2*TT
        self['PrecessionTheta'] = 0.5567530*TT - 0.0001185*TT2 - 1.16e-5*TT2*TT

        # Nutation terms
        self._calcNutation()
        self.__status['orbital'] = True

    def _calcNutation(self, nTerms=106):
        # Get dPSi and dEps in arcsec
        dPsi, dEps = self._nutation(nTerms=nTerms)  # default 106-term nutation series
        # Now apply EOP corrections
        dPsi += self['EarthOrientationParameters'].ddPsi
        dEps += self['EarthOrientationParameters'].ddEps
        # Equation of the equinoxes
        omegamoon_rad = np.deg2rad(self['OmegaMoon'])
        EquatEquin = dPsi*np.cos(self['ObliquityEcliptic_rad']) \
            + 0.00264*np.sin(omegamoon_rad) \
            + 0.000063*np.sin(2.0*omegamoon_rad)  # arcsec
        self['EquatEquin'] = EquatEquin
        # true obliquity of ecliptic
        self['ObliquityTrue'] = dEps + self['ObliquityEcliptic']  # arcsec
        self['ObliquityTrue_rad'] = self['ObliquityTrue']*scipy.constants.arcsec
        self['dPsi'] = dPsi
        self['dEps'] = dEps

    @staticmethod
    def _normVec(vec):
        tmpmag = np.linalg.norm(vec)
        vec /= tmpmag

    def _nutation(self, nTerms=106):
        pnmodel = self.attrs['pnmodel'].upper()
        if not pnmodel == 'IAU82' or pnmodel == 'LGMDEFAULT':
            import warnings
            warnings.warn('Only the IAU1980 nutation model is currently implemented. '
                          + 'This will be used with the selected precession model.')
        else:
            from . import iau80n
            dPsi, dEps = iau80n.nutation(self['TT_JC'], self['constants'], nTerms)
        # TODO: look at implementing IAU2000 as well
        #       See p88 onwards of USNO circular 179, for the coefficients.
        #       That has 1365 terms, and we'll likely want to have the default number evaluated
        #       set much lower.
        #       The equations are p45-48 of the same reference
        return dPsi, dEps

    def calcCoreTransforms(self, recalc=False):
        """
        Calculate core coordinate transform matrices

        Excludes coordinate systems requiring magnetic field information
        """
        if not self.__status['orbital'] or recalc:
            self.calcOrbitParams(recalc=recalc)
        if not recalc and self.__status['transformCore']:
            return
        self['Transform'] = dm.SpaceData()
        # ECI J2000 to ECI Mean of Date, and inverse
        zeta = np.deg2rad(self['PrecessionZeta'])
        cosZeta, sinZeta = np.cos(zeta), np.sin(zeta)
        zee = np.deg2rad(self['PrecessionZee'])
        cosZee, sinZee = np.cos(zee), np.sin(zee)
        theta = np.deg2rad(self['PrecessionTheta'])
        cosTheta, sinTheta = np.cos(theta), np.sin(theta)
        eci2000_ecimod = np.empty((3, 3))
        eci2000_ecimod[0, 0] = cosZeta*cosTheta*cosZee-sinZeta*sinZee
        eci2000_ecimod[0, 1] = -sinZeta*cosTheta*cosZee-cosZeta*sinZee
        eci2000_ecimod[0, 2] = -sinTheta*cosZee
        eci2000_ecimod[1, 0] = cosZeta*cosTheta*sinZee+sinZeta*cosZee
        eci2000_ecimod[1, 1] = -sinZeta*cosTheta*sinZee+cosZeta*cosZee
        eci2000_ecimod[1, 2] = -sinTheta*sinZee
        eci2000_ecimod[2, 0] = cosZeta*sinTheta
        eci2000_ecimod[2, 1] = -sinZeta*sinTheta
        eci2000_ecimod[2, 2] = cosTheta
        self['Transform']['ECI2000_ECIMOD'] = eci2000_ecimod
        self['Transform']['ECIMOD_ECI2000'] = eci2000_ecimod.T

        # ECI Mean of Date to ECI True of Date, and inverse
        ecitod_ecimod = np.empty((3, 3))
        dpsi_rad = self['dPsi']*self['constants'].arcsec
        sdps, cdps = np.sin(dpsi_rad), np.cos(dpsi_rad)
        sobl, cobl = np.sin(self['ObliquityEcliptic_rad']), np.cos(self['ObliquityEcliptic_rad'])
        sobt, cobt = np.sin(self['ObliquityTrue_rad']), np.cos(self['ObliquityTrue_rad'])
        ecitod_ecimod[0, 0] = cdps
        ecitod_ecimod[0, 1] = sdps*cobt
        ecitod_ecimod[0, 2] = sobt*sdps
        ecitod_ecimod[1, 0] = -sdps*cobl
        ecitod_ecimod[1, 1] = cobt*cdps*cobl + sobt*sobl
        ecitod_ecimod[1, 2] = sobt*cdps*cobl - sobl*cobt
        ecitod_ecimod[2, 0] = -sdps*sobl
        ecitod_ecimod[2, 1] = cobt*cdps*sobl - sobt*cobl
        ecitod_ecimod[2, 2] = sobt*sobl*cdps + cobt*cobl
        # TODO: if scipy version constraints allow, consider scipy.spatial.transform rotations
        self['Transform']['ECITOD_ECIMOD'] = ecitod_ecimod
        self['Transform']['ECIMOD_ECITOD'] = ecitod_ecimod.T

        # Need to get Sun vector for subsequent transforms, but Sun vector requires
        # J2000 <-> MOD transforms to be set up already
        self.__status['coreInProgress'] = True
        self._calcSun()
        self.__status['coreInProgress'] = False

        # ECI Mean of Date to Geocentric Solar Ecliptic, and reverse
        # First define axes for GSE coordinate system
        self['EclipticPole'] = np.array([0,  # x-hat
                                         -np.sin(self['ObliquityEcliptic_rad']),  # y-hat
                                         np.cos(self['ObliquityEcliptic_rad'])])
        gse_y_vec = np.cross(self['EclipticPole'], self['SunVector_MOD'])
        self._normVec(gse_y_vec)
        gse_z_vec = np.cross(self['SunVector_MOD'], gse_y_vec)
        # Set up transformation matrices for MOD <-> GSE
        ecimod_gse = np.empty((3, 3))
        ecimod_gse[0, :] = self['SunVector_MOD'][...]
        ecimod_gse[1, :] = gse_y_vec[...]
        ecimod_gse[2, :] = gse_z_vec[...]

        self['Transform']['ECIMOD_GSE'] = ecimod_gse
        self['Transform']['GSE_ECIMOD'] = ecimod_gse.T

        # Pseudo Earth Fixed (PEF) to/from ECI TOD
        # Use PEF as a pass-through to get to GEO
        # See "Revisiting Spacetrack Report #3" [Vallado]

        # First, get apparent sidereal time
        # Uses GMST in degrees and Equation of Equinoxes in degrees (IAU82? Or general?)
        gast_deg = self['GMST']*15.0 + self['EquatEquin']/3600.0
        gast_rad = np.deg2rad(gast_deg)
        csga = np.cos(gast_rad)
        snga = np.sin(gast_rad)
        ecitod_pef = np.zeros((3, 3))
        ecitod_pef[0, 0] = csga
        ecitod_pef[0, 1] = snga
        ecitod_pef[1, 0] = -snga
        ecitod_pef[1, 1] = csga
        ecitod_pef[2, 2] = 1.0
        pef_ecitod = ecitod_pef.T

        # TEME (used by WGS84)
        sn = np.sin(self['GMST_rad'])
        cs = np.cos(self['GMST_rad'])
        teme_pef = np.zeros((3, 3))
        teme_pef[0, 0] = cs
        teme_pef[0, 1] = sn
        teme_pef[1, 0] = -sn
        teme_pef[1, 1] = cs
        teme_pef[2, 2] = 1.0
        pef_teme = teme_pef.T

        # GEO (Geocentric geographic) to/from PEF (Psuedo Earth Fixed)
        xprad = self['EarthOrientationParameters'].xp*self['constants'].arcsec
        yprad = self['EarthOrientationParameters'].yp*self['constants'].arcsec
        sxp = np.sin(xprad)
        cxp = np.cos(xprad)
        syp = np.sin(yprad)
        cyp = np.cos(yprad)
        pef_geo = np.zeros((3, 3))
        pef_geo[0, 0] = cxp
        pef_geo[0, 1] = sxp*syp
        pef_geo[0, 2] = sxp*cyp
        pef_geo[1, 1] = cyp
        pef_geo[1, 2] = -syp
        pef_geo[2, 0] = -sxp
        pef_geo[2, 1] = cxp*syp
        pef_geo[2, 2] = cxp*cyp
        geo_pef = pef_geo.T

        # make/store ECIMOD<->GEO transforms
        ecitod_geo = pef_geo.dot(ecitod_pef)
        ecimod_geo = ecitod_geo.dot(self['Transform']['ECIMOD_ECITOD'])
        self['Transform']['ECIMOD_GEO'] = ecimod_geo
        self['Transform']['GEO_ECIMOD'] = ecimod_geo.T

        self.__status['transformCore'] = True

    def calcMagTransforms(self, recalc=False):
        """Calculate geophysical coordinate systems"""
        if not self.__status['transformCore'] or recalc:
            self.calcCoreTransforms()

        # Call IGRF to calculate centered dipole
        self['IGRF'] = igrf.IGRF()
        self['IGRF'].initialize(self['UTC'].UTC[0])

        # Get dipole axis unit vector
        glon = self['IGRF'].dipole['cd_glon_rad']
        gcolat = self['IGRF'].dipole['cd_gcolat_rad']
        dip = np.empty(3)
        dip[0] = np.cos(glon)*np.sin(gcolat)
        dip[1] = np.sin(glon)*np.sin(gcolat)
        dip[2] = np.cos(gcolat)

        # Convert dipole axis to MOD
        dipmod = np.empty(3)
        dipmod[0] = np.cos(self['GMST_rad'])*dip[0] - np.sin(self['GMST_rad'])*dip[1]
        dipmod[1] = np.sin(self['GMST_rad'])*dip[0] + np.cos(self['GMST_rad'])*dip[1]
        dipmod[2] = dip[2]

        # Compute Ygsm axis in MOD
        gsm_y_vec = np.cross(dipmod, self['SunVector_MOD'])
        self._normVec(gsm_y_vec)

        # Compute Zgsm axis in MOD system
        gsm_z_vec = np.cross(self['SunVector_MOD'], gsm_y_vec)

        # Set up transformation matrices for MOD <-> GSM
        ecimod_gsm = np.empty((3, 3))
        ecimod_gsm[0, :] = self['SunVector_MOD'][...]
        ecimod_gsm[1, :] = gsm_y_vec[...]
        ecimod_gsm[2, :] = gsm_z_vec[...]
        self['Transform']['ECIMOD_GSM'] = ecimod_gsm
        self['Transform']['GSM_ECIMOD'] = ecimod_gsm.T

        # Compute the Dipole Tilt angle (psi). First convert dipmod into GSM.
        # Then psi is just the angle between dipgsm and Zgsm
        dipgsm = self.convert(dipmod, 'ECIMOD', 'GSM')
        psidot = np.array([0, 0, 1]).dot(dipgsm)
        self._normVec(psidot)
        psi = np.arccos(psidot)
        psi *= -1 if (dipgsm[0] < 0) else 1
        self['DipoleTilt_rad'] = psi
        self['DipoleTilt'] = np.rad2deg(psi)
        cos_psi = np.cos(psi)
        sin_psi = np.sin(psi)

        # Set up transformation matrices for GSM <-> SM
        gsm_sm = np.empty((3, 3))
        gsm_sm[:, 0] = np.array([cos_psi, 0, sin_psi])
        gsm_sm[:, 1] = np.array([0, 1, 0], dtype=np.float)
        gsm_sm[:, 2] = np.array([-1*sin_psi, 0, cos_psi])
        self['Transform']['GSM_SM'] = gsm_sm
        self['Transform']['SM_GSM'] = gsm_sm.T
        # Get conversion to MOD for completeness
        self['Transform']['ECIMOD_SM'] = gsm_sm.dot(self['Transform']['ECIMOD_GSM'])
        self['Transform']['SM_ECIMOD'] = self['Transform']['ECIMOD_SM'].T

    def convert(self, vec, sys_in, sys_out):
        """Convert an input vector between two coordinate systems
        """
        transform = '{0}_{1}'.format(sys_in, sys_out)
        if transform not in self['Transform']:
            try:
                trans1 = '{0}_{1}'.format(sys_in, 'ECIMOD')
                trans2 = '{0}_{1}'.format('ECIMOD', sys_out)
                assert trans1 in self['Transform']
                assert trans1 in self['Transform']
                # If we have A->MOD and MOD->B then we can just
                # construct the conversion
                tmatr = self['Transform'][trans2].dot(self['Transform'][trans1])
            except (KeyError, AssertionError):
                self._raiseError(ValueError, 'transform')
        else:
            tmatr = self['Transform'][transform]
        return tmatr.dot(vec)

    def gmst(self):
        """Calculate Greenwich Mean Sidereal Time
        """
        if not (self.__status['time'] or self.__status['timeInProgress']):
            self.calcTimes(from_gmst=True)
        pnmodel = self.attrs['pnmodel']
        const = self['constants']
        # Greenwich Mean Sidereal Time
        UT1_P1 = self['UT1_JC']
        UT1_P2 = UT1_P1*UT1_P1  # UT1 squared
        UT1_P3 = UT1_P2*UT1_P1  # UT1 cubed
        coeffB = 8640184.812866
        coeffC = 0.093104
        coeffD = 6.2e-6
        if pnmodel.upper() == 'LGMDEFAULT':
            self['GMST'] = fmod((67310.54841 + (876600*3600 + coeffB)*UT1_P1 +
                                coeffC*UT1_P2 - coeffD*UT1_P3)/3600, 24)
            self['GMST_rad'] = np.deg2rad(self['GMST']*15)
        elif pnmodel.upper() == 'IAU82':
            # total fractional part of UT1 Julian day
            fts = const.daysec * (self['UT1_JD'] % 1 + const.j2000_jd % 1)
            gmst_rad = ((const.twopi/const.daysec) * ((-19089.45159 +
                        (coeffB + (coeffC + coeffD*UT1_P1) * UT1_P1)
                        * UT1_P1) + fts)) % const.twopi
            self['GMST'] = np.rad2deg(gmst_rad)/15
            self['GMST_rad'] = gmst_rad
        elif pnmodel.upper() == 'IAU00':
            du = self['UT1_JD'] - const.j2000_jd
            theta = const.twopi*(0.7790572732640 + 0.00273781191135448*du + du % 1) % const.twopi
            t = self['TT_JC']
            angle = (0.014506 +
                     (4612.15739966 +
                      (1.39667721 +
                       (- 0.00009344 +
                        (0.00001882)
                        * t) * t) * t) * t)
            self['GMST_rad'] = (theta + angle*const.arcsec) % const.twopi
            self['GMST'] = np.rad2deg(self['GMST_rad'])/15
        elif pnmodel.upper() == 'P03':
            raise NotImplementedError
            du = self['UT1_JD'] - const.j2000_jd
            theta = const.twopi*(0.7790572732640 + 0.00273781191135448*du + du % 1)
            t = self['TT_JC']
            # angle = 0.014506 + 4612.156534* + 1.3915817*t*t - 0.00000044*t*t*t \
            #         - 0.000029956*t**4 - 0.0000000368*t**5
            # self['GMST'] = (86400*theta) + angle/15)/3600.0 %24
            self['GMST'] = (UT1_P1 + 24110.5493771
                            + 8640184.79447825*UT1_P1 + 307.4771013*(t-UT1_P1)
                            + 0.0927721*t*t - 0.0000002926*t*t*t - 0.00000199708*t**4
                            - 0.000000002454*t**5)/3600.0 % 24
            self['GMST_rad'] = np.deg2rad(self['GMST']*15)
        else:
            raise Exception

    def _calcSun(self):
        """
        Calculate the Sun vector.

        This should only be called from within calcCoreTransforms, after
        the ECI2000 <-> ECIMOD transforms have been defined
        """
        validEntry = self.__status['time'] and self.__status['orbital']
        validEntry = validEntry and (self.__status['coreInProgress']
                                     or self.__status['transformCore'])
        if not validEntry:
            self._raiseErr(RuntimeError, 'sun')
        # TODO: check units of everything in this section
        const = self['constants']
        eccen = self['Eccentricity']
        cos_epsilon = np.cos(self['ObliquityEcliptic_rad'])
        sin_epsilon = np.sin(self['ObliquityEcliptic_rad'])
        cos_lambda = np.cos(self['LambdaSun'])
        sin_lambda = np.sin(self['LambdaSun'])

        if self.attrs['ephmodel'] == 'LGMDEFAULT':
            earth_sun_distance = const.AU*((1.0 - eccen * eccen)
                                           / (1.0 + eccen * np.cos(self['TrueAnomaly']))
                                           / const.Re)
            self['EarthSunDistance_Re'] = earth_sun_distance
            # Direction of the Sun in MOD coords
            sunVec = np.empty(3)
            sunVec[0] = cos_lambda
            sunVec[1] = cos_epsilon*sin_lambda
            sunVec[2] = sin_epsilon*sin_lambda
            self['SunVector_MOD'] = sunVec
            # Convert to ECI2000
            self['SunVector_J2000'] = self.convert(sunVec, 'ECIMOD', 'ECI2000')
        elif self.attrs['ephmodel'] == 'DE421':
            raise NotImplementedError
            # TODO: get Sun vector from DE421 (read HDF5 of coeffs from LGM github?)
            #       and use numpy chebyshev module for polynomials from coeffs
            # TODO: get Earth-Sun distance as magnitude of Sun vector (put in R_E)
            # TODO: normalize Sun vector (this is in ECIJ2000)

    def _raiseErr(self, errtype, code):
        if code == 'ephmodel':
            err = 'No alternate models of ephemerides at this time'
        elif code == 'pnmodel':
            err = 'No alternate models of precession/nutation at this time'
        elif code == 'time_in':
            err = 'Input time must be single-valued. Allowed types are datetime.datetime, ' + \
                  'spacepy.time.Ticktock, or any valid input type for spacepy.time.Ticktock'
        elif code == 'eop':
            err = 'Earth Oriention Parameters not found'
        elif code == 'sun':
            err = 'Invalid entry into Sun vector routine. Not intended for direct use.'
        elif code == 'transform':
            err = 'Requested coordinate transform not defined.'
        else:
            err = 'Operation not defined'

        raise errtype('CTrans: {0}'.format(err))

    def _kepler(self, mean_anom, ecc, tol=1e-8, maxiter=100):
        """
        Solve Kepler's equation for eccentric anomaly

        mean_anom: float
            Mean anomaly in radians
        ecc: float
            Eccentricity of orbit

        Returns
        -------
        ecc_anom: float
            Eccentric anomaly in radians
        """
        ecc_anom = mean_anom + ecc * np.sin(mean_anom)

        count = 0
        error = 1
        while (error > tol) and (count < maxiter):
            ecc_anom_old = ecc_anom
            ecc_anom = ecc_anom_old + (mean_anom - ecc_anom_old + ecc * np.sin(ecc_anom_old)) \
                / (1 - ecc * np.cos(ecc_anom_old))
            count += 1
            error = np.abs(ecc_anom - ecc_anom_old)
        if count >= maxiter:
            import warnings
            warnings.warn('Estimation of eccentric anomaly failed to converge')
        return ecc_anom
