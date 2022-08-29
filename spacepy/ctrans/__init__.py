"""CTrans: Module for backend coordinate transformations in SpacePy

This module is primarily intended to provide a backend for the standard
:class:`~spacepy.coordinates.Coords` class rather than direct use, and the
interface is subject to change.

The :class:`CTrans` class calculates all of the necessary information
to convert between different coordinate systems *at a single time*. By
using :class:`~spacepy.coordinates.Coords` the handling of multiple
times is built in, and the calling syntax is backwards compatible with
the legacy IRBEM-backed coordinate transforms.

Coordinate systems supported by this module can broadly be described
by two categories. The first category is a broad set of Earth-centered
coordinate systems that are specified by astronomical parameters.
If we consider the International Celestial Reference Frame to be our
starting point, then taking the origin as the center of the Earth
instead of the solar barycenter gives us the Geocentric Celestial
Reference Frame (GCRF). All coordinate systems described here are
right-handed Cartesian systems, except geodetic.

Systems and their relationships:

- ECI2000: Earth-Centered Inertial, J2000 epoch
    This system can be considered equivalent to the GCRF, to within 10s
    of milliarcseconds. The z-axis is aligned with the mean celestial pole
    at the J2000 epoch. The x-axis is aligned with the mean equinox at the
    J2000 epoch. The y-axis completes and lies in the plane of the
    celestial equator.
- ECIMOD: Earth-Centered Inertial, Mean-of-Date
    This system accounts for precession between the J2000 epoch and the
    date of interest: The coordinate system is time-dependent. The system
    is defined similarly to ECI2000, but uses the mean equinox and mean
    equator of the date of interest.
- ECITOD: Earth-Centered Inertial, True-of-Date
    This system builds on ECIMOD and accounts for the nutation (the short-
    period perturbations on the precession). This system is therefore
    considered to use the true equator and true equinox of date.
- TEME: Earth-Centered Inertial, True Equator Mean Equinox
    This system is poorly defined in the literature, despite being used in
    the SGP4 orbital propagator (note that multiple versions of SGP4 exist,
    see e.g. Vallado et al. 2006; AIAA 2006-6753-Rev2). The mean equinox here
    is not the same as the mean equinox used in, e.g., ECIMOD, but lies along
    the true equator between the origin of the Pseudo Earth Fixed and ECITOD
    frames. It is highly recommended that TEME coordinates are converted to a
    standard system (e.g., ECI2000) before passing to other users or to
    different software.
- GSE: Geocentric Solar Ecliptic
    This system is not inertial. It is Earth-centered, with the x-axis
    pointing towards the Sun. The y-axis lies in the mean ecliptic plane
    of date, pointing in the anti-orbit direction. The z-axis is parallel
    to the mean ecliptic pole.
- GEO: Geocentric Geographic
    This system is not inertial. It is Earth-Centered and Earth-Fixed (also
    called ECEF), so that the coordinates of a point fixed on (or relative
    to) the surface of the Earth do not change. The x-axis lies in the
    Earth's equatorial plane (zero latitude) and intersects the Prime
    Meridian (zero longitude; Greenwich, UK). The z-axis points to True
    North (which is roughly aligned with the instantaneous rotation axis).
- GDZ: Geodetic
    This system is not inertial and is defined in terms of altitude above
    a reference ellipsoid, the geodetic latitude, and geodetic longitude.
    Geodetic longitude is identical to GEO longitude. Both the altitude
    and latitude depend on the ellipsoid used. While geodetic latitude is
    close to geographic latitude, they are not the same. The default here is
    to use the WGS84 reference ellipsoid.

The remaining coordinate systems are also reference to Earth's magnetic field.
Different versions of these systems exist, but the most common (and those given
here) use a centered dipole axis.

- GSM: Geocentric Solar Magnetospheric
    This system is similar to GSE, but is defined such that the centered dipole
    lies in the x-z plane. As in all of these systems, z is positive northward.
    The y-axis is thus perpendicular to both the Sun-Earth line and the
    centered dipole axis (of date, defined using the first 3 coefficients of the
    IGRF/DGRF). GSM is therefore a rotation about the x-axis from the GSE system.
- SM: Solar Magnetic
    The z-axis is aligned with the centered dipole axis of date (positive
    northward), and the y-axis is perpendicular to both the Sun-Earth line and
    the dipole axis. As with GSE and GSM, y is positive in the anti-orbit
    direction. The x-axis therefore is not aligned with the Sun-Earth line and
    SM is a rotation about the y-axis from the GSM system.
- CDMAG: Geomagnetic
    This is a geomagnetic analog of the GEO system. The z-axis is aligned with
    the centered dipole axis of date. The y-axis is perpendicular to
    to both the dipole axis and True North, i.e., y is the cross product of
    the z-axis of the GEO system with the dipole axis. The x-axis completes.

Classes
-------
.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    CTrans
    Ellipsoid

Functions
---------
.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    convert_multitime
    gdz_to_geo
    geo_to_gdz
    geo_to_rll
    rll_to_geo

Submodules
----------
.. autosummary::
    :template: clean_module.rst
    :toctree: autosummary

    iau80n
"""

__contact__ = 'Steve Morley, smorley@lanl.gov'


import datetime as dt
import collections
from math import fmod

import numpy as np
import scipy.constants

from spacepy import datamodel as dm
from spacepy import time as spt
from spacepy import igrf


class Ellipsoid(dm.SpaceData):
    """Ellipsoid definition class for geodetic coordinates

    Other Parameters
    ----------------
    name : str
        Name for ellipsoid, stored in attrs of returned Ellipsoid instance.
        Default is 'WGS84'
    A : float
        Semi-major axis (equatorial radius) of ellipsoid in km.
        Default is 6378.137km (WGS84_A)
    iFlat : float
        Inverse flattening of ellipsoid. Default is WGS84 value
        of 298.257223563.

    Returns
    -------
    out : Ellipsoid
        Ellipsoid instance storing all relevant paramters for geodetic conversion

    Notes
    -----

    .. versionadded:: 0.3.0
    """
    def __init__(self, name='WGS84', A=6378.137, iFlat=298.257223563):
        super(Ellipsoid, self).__init__(A=A, iFlat=iFlat, attrs={'name': name})
        self['A2'] = self['A']**2  # [km^2]
        self['Flat'] = 1/self['iFlat']
        self['E2'] = 2*self['Flat'] - self['Flat']**2
        self['E4'] = self['E2']*self['E2']
        self['1mE2'] = 1 - self['E2']
        self['B'] = self['A']*np.sqrt(self['1mE2'])
        self['B2'] = self['B']**2  # [km^2]
        self['A2mB2'] = self['A2'] - self['B2']  # [km^2]
        self['EP2'] = self['A2mB2']/self['B2']  # 2nd eccentricity squared


# World Geodetic System 1984 (WGS84) parameters
WGS84 = Ellipsoid()

magsys = ['GSM', 'SM', 'CDMAG']


class CTrans(dm.SpaceData):
    """Coordinate transformation class for a single instance in time

    A general coordinate conversion routine, which takes a numpy array (Nx3)
    of Cartesian vectors along with the names of the input and output
    coordinate systems and returns an array of the converted coordinates.

    Parameters
    ----------
    ctime : (spacepy.time.Ticktock, datetime, float, string)
        Input time stamp. Must have one time only. Accepted input formats

    Returns
    -------
    out : CTrans
        instance with self.convert, etc.

    Other Parameters
    ----------------
    ephmodel : str, optional
        Select ephemerides model (e.g., for determining Sun direction).
        Currently only 'LGMDEFAULT' is supported, for consistency with
        LANLGeoMag implementation.
    pnmodel : str, optional
        Select precession/nutation model set. Options are: 'IAU82' (default),
        and 'IAU00'.
    eop : bool, optional
        Use Earth Orientation Parameters

    See Also
    --------
    spacepy.coordinates.Coords

    Notes
    -----

    .. versionadded:: 0.3.0


    .. rubric:: Methods

    .. autosummary::

        ~CTrans.calcTimes
        ~CTrans.calcOrbitParams
        ~CTrans.calcCoreTransforms
        ~CTrans.calcMagTransforms
        ~CTrans.convert
        ~CTrans.getEOP
        ~CTrans.gmst

    .. automethod:: calcTimes
    .. automethod:: calcOrbitParams
    .. automethod:: calcCoreTransforms
    .. automethod:: calcMagTransforms
    .. automethod:: convert
    .. automethod:: getEOP
    .. automethod:: gmst
    """
    def __init__(self, ctime, ephmodel=None, pnmodel=None, eop=False):
        super(CTrans, self).__init__()
        if ephmodel is not None:
            self._raiseErr(NotImplementedError, 'ephmodel')
        else:
            self.attrs['ephmodel'] = 'LGMDEFAULT'
        if pnmodel is not None:
            if pnmodel not in ['IAU82', 'IAU00']:
                self._raiseErr(NotImplementedError, 'pnmodel')
            self.attrs['pnmodel'] = pnmodel
        else:
            self.attrs['pnmodel'] = 'IAU82'

        if isinstance(ctime, spt.Ticktock):
            # input time is ticktock
            if len(ctime.data) != 1:
                # Input time is Ticktock, but has a length > 1
                self._raiseErr(ValueError, 'time_in')
        elif isinstance(ctime, dt.datetime):
            # Input time is datetime
            ctime = spt.Ticktock(ctime, dtype='UTC')
        else:
            try:
                if len(ctime) == 1:
                    ctime = ctime[0]
                    ctime = spt.Ticktock(ctime)  # Guess dtype
                else:
                    self._raiseErr(TypeError, 'time_in')
            except TypeError:
                self._raiseErr(TypeError, 'time_in')
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
                         'useEOP': eop
                         }

    def _setconstants(self):
        """Set constants to be used in calculations

        Constants set through this method are immutable to ensure
        consistency in calculations.
        """
        self['constants'] = self._factory['constants'](twopi=scipy.constants.pi*2,
                                                       j2000_jd=2451545.0,
                                                       j1990_jd=2447891.5,
                                                       j1900_jd=2415020.0,
                                                       daysec=86400,
                                                       daycentury=36525.0,
                                                       arcsec=scipy.constants.arcsec,
                                                       AU=scipy.constants.au/1e3,  # 1AU in km
                                                       Re=WGS84['A'],  # WGS84_A radius in km
                                                       )

    def getEOP(self, useEOP=False):
        """Get/set Earth Orientation Parameters

        Parameters
        ----------
        useEOP : bool
            If True, use Earth Orientation Parameters. Default False.

        Notes
        -----
        Currently Earth Orientation Parameters are all set to zero.
        Use is not yet supported.
        """
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
            # Note: These case still be set manually by directly calling the
            # factory method.
            self._raiseErr(NotImplementedError, 'eop')

    def calcTimes(self, recalc=False, **kwargs):
        """Calculate time in systems required to set up coordinate transforms

        Sets Julian Date and Julian centuries in UTC, TAI, UT1, and TT systems.

        Parameters
        ----------
        recalc : bool, optional
            If True, recalculate the times for coordinate transformation. Default
            is False.
        """
        if self.__status['time'] and not recalc:
            return
        # Set up various time systems required and variable shortcuts
        eop = self['EarthOrientationParameters']
        const = self['constants']
        ctime = self.attrs['time']

        # UT1
        self['UTC'] = ctime
        self['UTC_JD'] = ctime.JD[0]
        self['UTC_JC'] = (self['UTC_JD'] - const.j2000_jd)/const.daycentury
        self['TAI_JD'] = spt._days1958(self['UTC'].TAI[0], leaps='continuous') + 2436205.0
        if self.__status['useEOP']:
            self['UT1'] = spt.Ticktock(ctime.UTC[0] + dt.timedelta(seconds=eop.DUT1), dtype='UTC')
            self['UT1_JD'] = self['UT1'].JD[0]
        else:
            # no EOP, so UT1 = UTC and we can save time recalculating times
            self['UT1'] = self['UTC']
            self['UT1_JD'] = self['UTC_JD']
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
        """Calculate Earth orbit parameters needed for coordinate transforms

        Calculates  Earth's orbital parameters required for defining coordinate
        system transformations, such as orbital eccentricity, the obliquity of
        the ecliptic, anomalies, and precession angles.

        Parameters
        ----------
        recalc : bool, optional
            If True, recalculate the orbital parameters for coordinate transformation.
            Default is False.
        """
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
        """Calculate nutation-related quantities

        Calculates nutation parameters, equation of equinoxes,
        true obliquity of the ecliptic, and updates object values
        """
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
        """Normalize input vector"""
        tmpmag = np.linalg.norm(vec)
        vec /= tmpmag

    def _nutation(self, nTerms=106):
        """Calculate nutation terms dPsi and dEps"""
        pnmodel = self.attrs['pnmodel'].upper()
        if not pnmodel == 'IAU82':
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

        These coordinate systems do not require information about
        Earth's magnetic field.
        The systems are:
        Earth-Centered Inertial, J2000 (ECI2000)
        Earth-Centered Inertial, Mean-of-date (ECIMOD)
        Earth-Centered Inertial, True-of-date (ECITOD)
        Geocentric Solar Ecliptic (GSE)
        Geocentric Geographic (GEO)

        Parameters
        ----------
        recalc : bool, optional
            If True, recalculate the core (non-magnetic) coordinate transformations.
            Default is False.
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

        # TEME (SGP4 coordinate system)
        sn = np.sin(self['GMST_rad'])
        cs = np.cos(self['GMST_rad'])
        teme_pef = np.zeros((3, 3))
        teme_pef[0, 0] = cs
        teme_pef[0, 1] = sn
        teme_pef[1, 0] = -sn
        teme_pef[1, 1] = cs
        teme_pef[2, 2] = 1.0

        # add TEME to transforms
        pef_teme = teme_pef.T
        ecitod_teme = pef_teme.dot(ecitod_pef)
        ecimod_teme = ecitod_teme.dot(self['Transform']['ECIMOD_ECITOD'])
        self['Transform']['ECIMOD_TEME'] = ecimod_teme
        self['Transform']['TEME_ECIMOD'] = ecimod_teme.T

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

        # make/store ECIMOD<->GEO transforms
        ecitod_geo = pef_geo.dot(ecitod_pef)
        ecimod_geo = ecitod_geo.dot(self['Transform']['ECIMOD_ECITOD'])
        self['Transform']['ECIMOD_GEO'] = ecimod_geo
        self['Transform']['GEO_ECIMOD'] = ecimod_geo.T

        self.__status['transformCore'] = True

    def calcMagTransforms(self, recalc=False):
        """Calculate geophysical coordinate systems

        Calculate transforms for coordinate systems requiring magnetic field
        information.

        These are:
        Solar Magnetic (SM)
        Geocentric Solar Magnetospheric (GSM)
        Geomagnetic, centered dipole (CDMAG)

        Parameters
        ----------
        recalc : bool, optional
            If True, recalculate the core (non-magnetic) coordinate transformations.
            Default is False.
        """
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
        gsm_sm[:, 1] = np.array([0, 1, 0], dtype=float)
        gsm_sm[:, 2] = np.array([-1*sin_psi, 0, cos_psi])
        self['Transform']['GSM_SM'] = gsm_sm
        self['Transform']['SM_GSM'] = gsm_sm.T
        # Get conversion to MOD for completeness
        self['Transform']['ECIMOD_SM'] = gsm_sm.dot(self['Transform']['ECIMOD_GSM'])
        self['Transform']['SM_ECIMOD'] = self['Transform']['ECIMOD_SM'].T

        # Set up transformation matrices for GEO <-> CDMAG
        # CDMAG is centered dipole geomagnetic
        #    Z: parallel to dipole axis.
        #    Y: perpendicular to both Z and rotation axis
        #    X: completes
        zhat = np.array([0, 0, 1])  # GEO Z-axis in GEO
        mag_y_vec = np.cross(zhat, dip)
        self._normVec(mag_y_vec)
        mag_x_vec = np.cross(mag_y_vec, dip)
        self._normVec(mag_x_vec)
        geo_cdmag = np.empty((3, 3))
        geo_cdmag[0, ...] = mag_x_vec
        geo_cdmag[1, ...] = mag_y_vec
        geo_cdmag[2, ...] = dip
        self['Transform']['GEO_CDMAG'] = geo_cdmag
        self['Transform']['CDMAG_GEO'] = geo_cdmag.T
        # Get conversion to MOD for completeness
        self['Transform']['ECIMOD_CDMAG'] = geo_cdmag.dot(self['Transform']['ECIMOD_GEO'])
        self['Transform']['CDMAG_ECIMOD'] = self['Transform']['ECIMOD_CDMAG'].T

    def convert(self, vec, sys_in, sys_out, defaults=None):
        """Convert an input vector between two coordinate systems

        Parameters
        ----------
        vec : array-like
            Input 3-vector (can be an array of input 3-vectors) to convert.
            E.g., for 2D input the array must be like [[x1, y1, z1], [x2, y2, z2]]
        sys_in : str
            String name for initial coordinate system. For supported systems,
            see module level documentation.
        sys_out : str
            String name for target coordinate system. For supported systems,
            see module level documentation.

        Other Parameters
        ----------------
        defaults : namedtuple or None
            Named tuple containing default settings passed from Coordinates module
        """
        # must be at least 2D for conversion methods
        gvec = np.atleast_2d(vec)
        # Needs 3xN vectors for a broadcast dot product
        trvec = gvec.T

        # Special case geodetic transforms
        to_geodetic = False
        to_rll = False
        if sys_out == 'GDZ':
            sys_out = 'GEO'
            to_geodetic = True
        elif sys_out == 'RLL':
            sys_out = 'GEO'
            to_rll = True
        if sys_in == 'GDZ':
            trvec = gdz_to_geo(vec) if defaults is None else gdz_to_geo(vec, geoid=defaults.ellipsoid)
            sys_in = 'GEO'
        elif sys_in == 'RLL':
            trvec = rll_to_geo(vec) if defaults is None else rll_to_geo(vec, geoid=defaults.ellipsoid)
            sys_in = 'GEO'

        if sys_in != sys_out:
            transform = '{0}_{1}'.format(sys_in, sys_out)
            if transform not in self['Transform']:
                try:
                    # Construct requested transform via ECIMOD
                    trans1 = '{0}_{1}'.format(sys_in, 'ECIMOD')
                    trans2 = '{0}_{1}'.format('ECIMOD', sys_out)
                    assert trans1 in self['Transform']
                    assert trans2 in self['Transform']
                    tmatr = self['Transform'][trans2].dot(self['Transform'][trans1])
                    self['Transform'][transform] = tmatr
                except (KeyError, AssertionError):
                    # Can't construct the transform requested
                    self._raiseErr(ValueError, 'transform')
            else:
                # Required transform stored already, just use it
                tmatr = self['Transform'][transform]

            converted = tmatr.dot(trvec).T
        else:
            converted = trvec
        # squeeze to fix return of 1D input vectors
        # squeeze returns either input array or a view, so this isn't wasteful
        converted_squeezed = converted.squeeze()

        if to_geodetic:
            if defaults is not None:
                converted_squeezed = geo_to_gdz(converted_squeezed.T, geoid=defaults.ellipsoid)
            else:
                converted_squeezed = geo_to_gdz(converted_squeezed.T)
        elif to_rll:
            if defaults is not None:
                converted_squeezed = geo_to_rll(converted_squeezed.T, geoid=defaults.ellipsoid)
            else:
                converted_squeezed = geo_to_rll(converted_squeezed.T)
        return converted_squeezed

    def gmst(self):
        """Calculate Greenwich Mean Sidereal Time

        Notes
        -----
        The formulation used to calculate GMST is selected using the
        status of the 'pnmodel' variable in the CTrans object attributes.
        """
        if not (self.__status['time'] or self.__status['timeInProgress']):
            self.calcTimes(from_gmst=True)
        pnmodel = self.attrs['pnmodel']
        const = self['constants']
        # Greenwich Mean Sidereal Time
        UT1_P1 = self['UT1_JC']
        coeffB = 8640184.812866
        coeffC = 0.093104
        coeffD = 6.2e-6
        if pnmodel.upper() == 'IAU82':
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
        """Common error raising method

        Parameters
        ----------
        errtype : Exception
            Exception to raise
        code : str
            Error category. Defines error message.
        """
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

        Parameters
        ----------
        mean_anom : float
            Mean anomaly in radians
        ecc : float
            Eccentricity of orbit

        Returns
        -------
        ecc_anom : float
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


def geo_to_gdz(geovec, units='km', geoid=WGS84):
    """
    Convert geocentric geographic (cartesian GEO) to geodetic (spherical GDZ)

    Uses Heikkinen's exact solution [#Heikkinen]_, see Zhu et al. [#Zhu] for
    details.

    Parameters
    ----------
    geovec : array-like
        Nx3 array (or array-like) of geocentric geographic [x, y, z] coordinates

    Returns
    -------
    out : numpy.ndarray
        Nx3 array of geodetic altitude, latitude, and longitude

    Notes
    -----
    .. versionadded:: 0.3.0

    References
    ----------
    .. [#Heikkinen] Heikkinen, M., "Geschlossene formeln zur berechnung raumlicher geodatischer
            koordinaten aus rechtwinkligen koordinaten", Z. Vermess., vol. 107, pp. 207-211,
            1982.
    .. [#Zhu] J. Zhu, "Conversion of Earth-centered Earth-fixed coordinates to geodetic
            coordinates," in IEEE Transactions on Aerospace and Electronic Systems, vol. 30,
            no. 3, pp. 957-961, July 1994, doi: 10.1109/7.303772.
    """
    posarr = np.atleast_2d(geovec).T
    x_geo = posarr[0, :]
    y_geo = posarr[1, :]
    z_geo = posarr[2, :]
    if units == 'Re':
        # Make sure positions are in km
        rx = x_geo*geoid['A']
        ry = y_geo*geoid['A']
        rz = z_geo*geoid['A']
    else:
        rx = x_geo
        ry = y_geo
        rz = z_geo

    rad2 = rx*rx + ry*ry
    rad = np.sqrt(rad2)
    z2 = rz*rz
    # Heikkinen's method, as written in Zhu et al.
    # Each equation not explicitly in Zhu is marked with a trailing comment
    F = 54.0*geoid['B2']*z2
    G = rad2 + geoid['1mE2']*z2 - geoid['E2']*geoid['A2mB2']
    G2 = G*G  # Square G for convenience
    c = (geoid['E4']*F*rad2)/(G*G*G)
    s = np.cbrt(1 + c + np.sqrt(c*c + 2*c))
    denom_p = s + 1/s + 1  # Shorthand for term in P's denominator
    P = F/(3*denom_p*denom_p*G2)
    Q = np.sqrt(1 + 2*geoid['E4']*P)
    # The second term in r0 can evaluate as negative...
    sqrt_for_r0 = np.abs((0.5*geoid['A2'])*(1 + 1/Q) -
                         (geoid['1mE2']*P*z2)/(Q*(1 + Q)) -
                         0.5*P*rad2)
    # Back to Heikinnen's method per Zhu
    r0 = -(geoid['E2']*P*rad)/(1 + Q) + np.sqrt(sqrt_for_r0)
    brac_u = (rad - geoid['E2']*r0)  # Shorthand for term in U and V denominator
    U = np.sqrt(brac_u*brac_u + z2)
    V = np.sqrt(brac_u*brac_u + geoid['1mE2']*z2)
    z0 = (geoid['B2']*rz)/(geoid['A']*V)
    # Heikkinen's solution has a singularity at the pole, so let's trap that here
    polemask = rad <= 1e-12
    lati_gdz = np.empty_like(ry)
    if polemask.any():
        # Make sure we assign the correct pole
        zmask = rz > 0
        nhem = np.logical_and(polemask, zmask)
        shem = np.logical_and(polemask, ~zmask)
        lati_gdz[nhem] = 90
        lati_gdz[shem] = -90
    lati_gdz[~polemask] = np.rad2deg(np.arctan((rz[~polemask] +
                                                geoid['EP2']*z0[~polemask]) /
                                                rad[~polemask]))
                                                # Geodetic latitude (phi in Zhu paper)
    alti_gdz = U*(1 - geoid['B2']/(geoid['A']*V))  # Geodetic altitude [km] (h in Zhu paper)
    long_gdz = np.rad2deg(np.arctan2(ry, rx))  # Geodetic longitude (same as GEO)

    out = np.c_[alti_gdz, lati_gdz, long_gdz]
    if units == 'Re':
        # Return in Re
        out[:, 0] /= geoid['A']
    return out.squeeze()


def gdz_to_geo(gdzvec, units='km', geoid=WGS84):
    """
    Convert geodetic (GDZ) coordinates to geocentric geographic

    Parameters
    ----------
    gdzvec : array-like
        Nx3 array of geodetic altitude, latitude, longitude (in specified units)

    Returns
    -------
    out : numpy.ndarray
        Nx3 array of geocentric geographic x, y, z coordinates

    Other Parameters
    ----------------
    units : str
        Units for input geodetic altitude. Options are 'km' or 'Re'. Default is 'km'.
        Output units will be the same as input units.
    geoid : spacepy.ctrans.Ellipsoid
        Instance of a reference ellipsoid to use for geodetic conversion.
        Default is WGS84.

    Notes
    -----
    .. versionadded:: 0.3.0
    """
    posarr = np.atleast_2d(gdzvec)
    h = posarr[:, 0] if units == 'km' else posarr[:, 0]*geoid['A']  # convert to km
    lam = np.deg2rad(posarr[:, 1])
    phi = np.deg2rad(posarr[:, 2])
    chi = np.sqrt(1 - geoid['E2']*np.sin(lam)*np.sin(lam))

    # Convert to GEO [km]
    x = (geoid['A']/chi + h)*np.cos(lam)*np.cos(phi)
    y = (geoid['A']/chi + h)*np.cos(lam)*np.sin(phi)
    z = (geoid['A']*(1 - geoid['E2'])/chi + h)*np.sin(lam)

    out = np.c_[x, y, z]
    if units == 'km':
        return out.squeeze()
    else:
        # Return in Re
        return out.squeeze()/geoid['A']


def geo_to_rll(geovec, units='km', geoid=WGS84):
    """Calculate RLL from geocentric geographic (GEO) coordinates

    Parameters
    ----------
    geovec : array-like
        Nx3 array of geographic radius, latitude, longitude (in specified units)

    Returns
    -------
    rllvec : numpy.ndarray
        Nx3 array of [distance from Earth's center, geodetic latitude, geodetic longitude]

    Other Parameters
    ----------------
    units : str
        Units for input geodetic altitude. Options are 'km' or 'Re'. Default is 'km'.
        Output units will be the same as input units.
    geoid : spacepy.ctrans.Ellipsoid
        Instance of a reference ellipsoid to use for geodetic conversion.
        Default is WGS84.

    Notes
    -----
    .. versionadded:: 0.3.0
    """
    rllvec = np.atleast_2d(geo_to_gdz(geovec, units=units, geoid=geoid))
    # Reolace altitude with norm of Cartesian position
    rllvec[:, 0] = np.linalg.norm(geovec, axis=-1)

    return rllvec.squeeze()


def rll_to_geo(rllvec, units='km', geoid=WGS84):
    """Calculate geocentric geographic (GEO) from RLL coordinates

    Parameters
    ----------
    rllvec : array-like
        Nx3 array of geocentric radius, geodetic latitude, geodetic longitude
        (in specified units)

    Returns
    -------
    geoarr : numpy.ndarray
        Nx3 array of [altitude, geodetic latitude, geodetic longitude]

    Other Parameters
    ----------------
    units : str
        Units for input geocentric radii. Options are 'km' or 'Re'. Default is 'km'.
        Output units will be the same as input units.
    geoid : spacepy.ctrans.Ellipsoid
        Instance of a reference ellipsoid to use for geodetic conversion.
        Default is WGS84.

    Notes
    -----
    .. versionadded:: 0.3.0
    """
    posarr = np.atleast_2d(rllvec).astype(float)
    surf = np.zeros_like(posarr)
    surf[:, 1:] = posarr[:, 1:]
    geoid_at_pos = np.atleast_2d(gdz_to_geo(surf, units=units, geoid=geoid))
    gdz = posarr.copy()
    # remove geoid height from geocentric radius at each location
    gdz[:, 0] -= np.linalg.norm(geoid_at_pos, axis=-1)
    geoarr = gdz_to_geo(gdz, units=units, geoid=geoid)

    return geoarr


def convert_multitime(coords, ticks, sys_in, sys_out, defaults=None, itol=None):
    '''
    Convert coordinates for N times, where N >= 1

    Parameters
    ----------
    coords : array-like
        Coordinates as Nx3 array. Cartesian assumed unless input system is geodetic.
    ticks : spacepy.time.Ticktock
        Times for each element of coords. Must contain either N times or 1 time.
    sys_in : str
        Name of input coordinate system.
    sys_out : str
        Name of output coordinate system.

    Other Parameters
    ----------------
    defaults : namedtuple or None
        Named tuple with parameters from coordinates module
    itol : float
        Time tolerance, in seconds, for using a unique set of conversions.
        Default is 1. Supplying a defaults namedtuple (i.e., if routine is called
        by spacepy.cooordinates.Coords.convert) will override this value.
    '''
    itol = itol if defaults is None else defaults.itol
    # Calculate a CTrans for each unique time
    ctdict = dict()
    ttused = set()
    tais = ticks.TAI
    for idx, tt in enumerate(tais):
        if tt not in ctdict:
            # for speed, let's not recalculate transforms
            # for times very close to each other (defined
            # by the itol kwarg)
            if ttused:
                usedlist = list(ttused)
                closest = np.argmin(np.abs(tt - usedlist))
                tdiff = np.abs(tt - usedlist[closest])
            else:
                tdiff = itol + 1
            if tdiff < itol and ctdict:
                ctdict[tt] = ctdict[usedlist[closest]]
            else:
                ctdict[tt] = CTrans(ticks[idx])
                ctdict[tt].calcCoreTransforms()
                if sys_in in magsys or sys_out in magsys:
                    ctdict[tt].calcMagTransforms()
                ttused.add(tt)
    # Unless speed becomes an issue, let's just loop over each value
    # except the spacial case of only 1 unique time
    loopover = np.atleast_2d(coords)
    if len(ctdict) == 1:
        newcoords = ctdict[tais[0]].convert(loopover, sys_in, sys_out, defaults=defaults)
    else:
        newcoords = np.empty_like(loopover)
        for idx, cc in enumerate(loopover):
            newcoords[idx] = ctdict[tais[idx]].convert(cc, sys_in, sys_out, defaults=defaults)
    return newcoords.squeeze()
