#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Implementation of Coords class functions for coordinate transformations

The coordinate systems supported by this module cover the most commonly
used geophysical and magnetospheric systems. The naming conventions can
follow the names used by the popular IRBEM library, but for inertial
systems we use a more consistent, fine-grained naming convention that
clarifies the different systems.

Earth-centered Inertial Systems
-------------------------------
    * **ECI2000** Earth-centered Inertial, J2000 epoch
    * **ECIMOD** Earth-centered Inertial, mean-of-date
    * **ECITOD** Earth-centered Inertial, true-of-date
    * **GEI** Geocentric Equatorial Inertial (IRBEM approximation of TOD)

Magnetospheric Systems
----------------------
    * **GSM** Geocentric Solar Magnetospheric
    * **GSE** Geocentric Solar Ecliptic
    * **SM** Solar Magnetic
    * **MAG** Geomagnetic Coordinate System (aka CDMAG)

Earth-fixed Systems
-------------------
    * **GEO** Geocentric geographic, aka Earth-centered Earth-fixed
    * **GDZ** Geodetic coordinates

By convention *all* systems are treated as natively Cartesian except
geodetic (GDZ), which is defined in [altitude, latitude, longitude]
where altitude is relative to a reference ellipsoid. Similarly, distance
units are assumed to be Earth radii (Re) in all systems except GDZ, where
altitude is given in km. Conversions to GDZ will output altitude in km
regardless of the input distance units and conversions from GDZ will
output in Re regardless of input units. In all other cases, the distance
units will be preserved.

    .. versionchanged:: 0.3.0

    The new CTrans backend was added, which includes support for the names
    ``ECI2000``, ``ECIMOD``, ``ECITOD``, and ``CDMAG``. With the exception
    of ``ECIMOD``, these can be used with the existing IRBEM backend, and
    will be converted to their closest equivalents.

    .. versionchanged:: 0.4.0

    The default backend for coordinate transformations was changed from IRBEM
    to the CTrans-based SpacePy backend.

Notes on differences between representations
--------------------------------------------
IRBEM's coordinate transformations are low-accuracy and were written for
a library with a driving philosophy of speed and robustness as priorities.
The coordinate transformations are therefore approximate. Further, most of
the geophysical systems (e.g., GSE, SM) are derived from an inertial
system. It is standard practice to use ECIMOD as this system. However,
IRBEM does not currently make ECIMOD available as one of its inertial
systems. IRBEM's default inertial system (called GEI) is consistent with
an approximation of ECITOD. Hence there will be small differences between
IRBEM's transformations and those using SpacePy's CTrans backend.
Further sources of difference include: IRBEM uses a low-order approximation
to the sidereal time and other parameters; the calculation of the Earth-Sun
vector differs between the representations; the definitions of an Earth
radius differ (SpacePy = 6378.137km; IRBEM = 6371.2 km). SpacePy's in-built
representation is higher accuracy and is comprehensively tested, including
tests for consistency with other high accuracy packages such as LANLGeoMag
and AstroPy. However, for use cases where the required precision is of order
1 percent the output can be considered equivalent.

Setting options for coordinate transformation
---------------------------------------------
The backend for coordinate transformations can be provided at
instantiation of a :class:`Coords` object using a keyword
argument. However, for convenience and flexibility the options can be
set at the module level. Configurable options include the backend used
(:mod:`~spacepy.irbempy` or SpacePy's :mod:`~spacepy.ctrans`) and the
reference ellipsoid (only configurable for the SpacePy backend). A
warning will be raised if the backend is not set (either through the
defaults or the keyword argument). The final configurable option
(``itol``) is the maximum separation, in seconds, for which the
coordinate transformations will not be recalculated. To force all
transformations to use an exact transform for the time, set ``itol``
to zero. Values between 10s and 60s are recommended for speed while
also preserving accuracy, though different applications will require
different accuracies.  For example, assuming this module has been
imported as ``spc``, to set the SpacePy backend as the default and set
``itol`` to 5 seconds:

    >>> spc.DEFAULTS.set_values(use_irbem=False, itol=5)

Authors: Steven Morley and Josef Koller
Institution: Los ALamos National Laboratory
Contact: smorley@lanl.gov

Copyright 2010-2016 Los Alamos National Security, LLC.

"""

from collections import namedtuple
import numbers
import numpy as np
from spacepy import help
import spacepy
import spacepy.datamodel as dm
import spacepy.time
from spacepy import ctrans

__contact__ = 'Steven Morley, smorley@lanl.gov'


# -----------------------------------------------
# space coordinate class
# -----------------------------------------------

SYSAXES_TYPES = {'GDZ': {'sph': 0, 'car': None},
                 'GEO': {'sph': None, 'car': 1}, 'GSM': {'sph': None, 'car': 2},
                 'GSE': {'sph': None, 'car': 3}, 'SM': {'sph': None, 'car': 4},
                 'GEI': {'sph': None, 'car': 5}, 'CDMAG': {'sph': None, 'car': 6},
                 'MAG': {'sph': None, 'car': 6}, 'SPH': {'sph': 7, 'car': None},
                 'RLL': {'sph': 8, 'car': None}, 'TOD': {'sph': None, 'car': 12},
                 'ECITOD': {'sph': None, 'car': 12}, 'J2000': {'sph': None, 'car': 13},
                 'ECI2000': {'sph': None, 'car': 13}, 'TEME': {'sph': None, 'car': 14},
                 'ECIMOD': {'sph': None, 'car': 15}}

SYS_EQUIV = {'SPH': 'GEO', 'GEI': 'ECITOD', 'TOD': 'ECITOD',
             'J2000': 'ECI2000', 'MAG': 'CDMAG'}

IRBEM_RE = 6371.2  # kilometers


class __Defaults(object):
    '''Configuration factory for coordinates module default settings
    '''
    def __init__(self):
        self._fac = namedtuple('values', 'use_irbem, ellipsoid, itol')
        self.set_values()

    def __repr__(self):
        full = '{}'.format(self.values)
        b1 = full.index('(')
        return 'DEFAULTS{}'.format(full[b1:])

    def set_values(self, use_irbem=False, ellipsoid=ctrans.WGS84, itol=30):
        '''Set options for coordinate transforms

        Parameters
        ----------
        use_irbem : bool
            If True, use IRBEM as coordinate transform backend. If False, use SpacePy's
            implementations. Default is to use SpacePy.
        ellipsoid : spacepy.ctrans.Ellipsoid
            A reference ellipsoid to use for the Earth radius definition and for geodetic
            coordinate conversion. SpacePy defaults to the WGS84 ellipsoid and the semi-major
            axis as the Earth radius (6738.137 km). The IRBEM backend uses WGS84, but expresses
            an Earth radius using 6371.2 km.
        '''
        self.values = self._fac(use_irbem=use_irbem,
                                ellipsoid=ellipsoid,
                                itol=itol)

DEFAULTS = __Defaults()


class Coords(object):
    '''
    a = Coords( data, dtype, carsph, [units, ticks, use_irbem])

    A class holding spatial coordinates and enabling transformation between
    coordinate systems. Coordinates can be stored as Cartesian or spherical
    and units are assumed to be Re (distance) and degrees (angle)

    .. note:: Although other units may be specified and will be carried
         through, most functions throughout SpacePy assume distances
         in Re and angles in degrees, regardless of specified units.

    By default, coordinate transforms are based on the SpacePy library's
    high-accuracy coordinates backend.
    The legacy transforms provided by the IRBEM library can also be used
    by setting the `use_irbem` flag to True; its manual
    <http://svn.code.sf.net/p/irbem/code/trunk/manual/user_guide.html>`_
    may prove useful.
    For a good reference on heliospheric and magnetospheric
    coordinate systems, see Franz & Harper, "Heliospheric Coordinate Systems",
    Planet. Space Sci., 50, pp 217-233, 2002
    (https://doi.org/10.1016/S0032-0633(01)00119-2).

    Parameters
    ----------
    data : list or ndarray, dim = (n,3)
        coordinate points [X,Y,Z] or [rad, lat, lon]
    dtype : string
        coordinate system; supported systems are defined in
        module-level documentation. Common systems include
        GEO, GSE, GSM, SM, MAG, ECIMOD
    carsph : string
        Cartesian or spherical, 'car' or 'sph'
    units : list of strings, optional
        standard are  ['Re', 'Re', 'Re'] or ['Re', 'deg', 'deg'] depending
        on the carsph content. See note.
    ticks : Ticktock instance, optional
        used for coordinate transformations (see a.convert)
    use_irbem : bool

        .. versionadded:: 0.3.0

        Set to True to use IRBEM for coordinate transforms.
        Otherwise use SpacePy's coordinate transform library.

    Returns
    -------
    out : Coords instance
        instance with a.data, a.carsph, etc.

    See Also
    --------
    spacepy.coordinates.DEFAULTS
    spacepy.time.Ticktock

    Examples
    --------
    >>> from spacepy import coordinates as coord
    >>> cvals = coord.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
    >>> cvals.x  # returns all x coordinates
    array([1, 1])
    >>> from spacepy.time import Ticktock
    >>> cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
    >>> newcoord = cvals.convert('GSM', 'sph')
    >>> newcoord

    .. currentmodule:: spacepy.coordinates
    .. autosummary::
        ~Coords.append
        ~Coords.convert
        ~Coords.from_skycoord
        ~Coords.to_skycoord
    .. automethod:: append
    .. automethod:: convert
    .. automethod:: from_skycoord
    .. automethod:: to_skycoord
    '''
    def __init__(self, data, dtype, carsph, units=None, ticks=None, use_irbem=None):
        if use_irbem is None:
            use_irbem = DEFAULTS.values.use_irbem
        if use_irbem:
            from . import irbempy as op
        self.use_irbem = use_irbem
        # setup units
        self.Re = IRBEM_RE if use_irbem else DEFAULTS.values.ellipsoid['A']  # kilometers

        # Make sure that inputs are all formed correctly
        self.data = np.atleast_2d(data).astype(np.float_)
        if len(self.data.shape) != 2 or self.data.shape[-1] != 3:
            raise ValueError('Input position vectors must be Nx3')

        dtype = dtype.upper()
        carsph = carsph.lower()

        if dtype not in list(SYSAXES_TYPES.keys()):
            raise ValueError('This dtype=' + dtype + ' is not supported. Only ' +
                             str(list(SYSAXES_TYPES.keys()))
                             )
        if carsph not in ['car', 'sph']:
            raise ValueError('This carsph=' + str(carsph) +
                             ' is not supported. Only "car" or "sph"')

        if dtype in ['GDZ', 'RLL'] and carsph == 'car':
            raise ValueError("Geodetic coordinates are referenced to an ellipsoid "
                             + "and can not be properly represented as Cartesian")

        # add ticks
        if ticks:
            assert len(ticks) == len(self.data), 'Ticktock dimensions seem off'
        self.ticks = ticks

        # GEO and SPH are cartesian/spherical versions of GEO, so use GEO for everything
        if dtype == 'SPH':
            dtype = 'GEO'
        # Make sure we're using the SpacePy name if we're using the SpacePy backend
        self.dtype = SYS_EQUIV[dtype] if (not use_irbem and dtype in SYS_EQUIV) else dtype
        self.sysaxes = SYSAXES_TYPES[dtype][carsph]

        self.carsph = carsph
        # setup units
        if units is None and carsph == 'car':
            # use standard units
            self.units = ['Re', 'Re', 'Re']
        elif units is None and carsph == 'sph':
            self.units = ['Re', 'deg', 'deg']
        else:
            self.units = units
        if dtype == 'GDZ' and carsph == 'sph':
            self.units = ['km', 'deg', 'deg']
        # setup x,y,z etc
        if carsph == 'car':
            self.x = spacepy.datamodel.dmcopy(self.data[:, 0])
            self.y = spacepy.datamodel.dmcopy(self.data[:, 1])
            self.z = spacepy.datamodel.dmcopy(self.data[:, 2])
        else:
            self.radi = spacepy.datamodel.dmcopy(self.data[:, 0])
            self.lati = spacepy.datamodel.dmcopy(self.data[:, 1])
            self.long = spacepy.datamodel.dmcopy(self.data[:, 2])
        self.shape = np.shape(self.data)
        return None

    # -----------------------------------------------
    def __str__(self):
        '''
        a.__str__() or a

        Will be called when printing Coords instance a

        Returns
        -------
        out : string
            string represenation of the instance

        Examples
        --------
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y
        Coords( [[1 2 4]
         [1 2 2]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']
        '''
        rstr = "Coords( {0} , '{1}', '{2}')".format(self.data.tolist(), self.dtype, self.carsph)
        return rstr
    __repr__ = __str__

    # -----------------------------------------------
    def __getitem__(self, idx):
        '''
        a.__getitem__(idx) or a[idx]

        Will be called when requesting items in this instance

        Parameters
        ----------
        idx : int
            integer numbers as index

        Returns
        -------
        out : numpy array
            new values

        Examples
        --------
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[0]
        Coords( [[1 2 4]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']
        '''
        arr = np.array(self.data)
        t_select = self.ticks[idx] if self.ticks else self.ticks

        return Coords(arr[idx].tolist(), self.dtype, self.carsph, self.units,
                      t_select, use_irbem=self.use_irbem)

    # -----------------------------------------------
    def __setitem__(self, idx, vals):
        '''
        a.__setitem__(idx, vals) or a[idx] = vals

        Will be called setting items in this instance

        Parameters
        ----------
        idx : int
            integer numbers as index
        vals : numpy array or list
            new values

        Examples
        --------
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[1] = [9,9,9]
        >>> y
        Coords( [[1 2 4]
            [9 9 9]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']
        '''
        self.data[idx] = vals
        return

    # -----------------------------------------------
    def __len__(self):
        '''
        a.__len__() or len(a)

        Will be called when requesting the length, i.e. number of items

        Returns
        -------
        out : int
            length

        Examples
        --------
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> len(y)
        2

        '''
        if isinstance(self.data, (list, np.ndarray)):
            return len(self.data)
        else:
            return 1
        return

    # -----------------------------------------------
    def convert(self, returntype, returncarsph):
        '''Create a new Coords instance with new coordinate types

        Parameters
        ----------
        returntype : string
            coordinate system, see module level documentation for
            supported systems
        returncarsph : string
            coordinate type, possible 'car' for Cartesian and 'sph' for spherical

        Returns
        -------
        out : Coords object
            Coords object in the new coordinate system

        Examples
        --------
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> from spacepy.time import Ticktock
        >>> y.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> x = y.convert('SM','car')
        >>> x
        Coords( [[ 0.81134097  2.6493305   3.6500375 ]
         [ 0.92060408  2.30678864  1.68262126]] ), dtype=SM,car, units=['Re', 'Re', 'Re']
        '''
        # Check whether the name is an equivalent IRBEM name
        returnname = SYS_EQUIV[returntype] if returntype in SYS_EQUIV else returntype
        is_dtype_returntype = (self.dtype == returntype) or (self.dtype == returnname)

        if returntype in ['GDZ', 'RLL'] and returncarsph == 'car':
            raise ValueError("Geodetic coordinates are referenced to an ellipsoid "
                             + "and can not be properly represented as Cartesian")

        # Then check whether a full coordinate transformation is required
        if is_dtype_returntype:
            if self.carsph == returncarsph:
                # no change necessary
                return self

            if self.carsph != returncarsph:
                # only car2sph or sph2car is needed
                if returncarsph == 'car':
                    carsph = 'car'
                    units = [self.units[0]]*3
                    data = sph2car(self.data)
                else:
                    carsph = 'sph'
                    units = [self.units[0], 'deg', 'deg']
                    data = car2sph(self.data)
                return Coords(data, self.dtype, carsph, units, self.ticks,
                              use_irbem=self.use_irbem)

        # check the length of ticks and do the more complex conversions
        if self.ticks and (len(self.ticks) != len(self)):
            raise ValueError('Ticktock dimension does not match Coords dimensions')

        # check if car2sph is needed first for oneralib compatibility
        if (self.sysaxes is None):  # a car2sph or sph2car is needed
            if self.carsph == 'sph':
                carsph = 'car'
                units = [self.units[0]]*3
                data = sph2car(self.data)
            else:
                carsph = 'sph'
                units = [self.units[0], 'deg', 'deg']
                data = car2sph(self.data)
        else:
            data = self.data.copy()
            units = dm.dmcopy(self.units)
            carsph = dm.dmcopy(self.carsph)

        if self.use_irbem:
            from . import irbempy as op
            # irbempy.coord_trans needs the passed Coords object to have the "from" dtype
            NewCoords = Coords(data, self.dtype, carsph, units, self.ticks,
                               use_irbem=self.use_irbem)
        else:
            # SpacePy conversion takes input/output type separately, so set to desired output
            NewCoords = Coords(data, returntype if self.use_irbem else returnname,
                               returncarsph, units, self.ticks, use_irbem=self.use_irbem)

        # now convert to other coordinate system
        if not is_dtype_returntype:
            assert self.ticks, "Time information required; add a .ticks attribute"
            if self.use_irbem:
                if returnname in ['GDZ', 'RLL'] and units[0] == 'km':
                    NewCoords.data /= self.Re  # IRBEM geodetic calculation requires Re input
                NewCoords.data = op.coord_trans(NewCoords, returntype, returncarsph)
                NewCoords.dtype = returntype
                if returnname == 'RLL' and units[0] == 'km':
                    NewCoords.data[:, 0] *= self.Re  # Preserve units
            else:
                if returnname in ['GDZ', 'RLL'] and self.units[0] == 'Re':
                    # converting from a Cartesian referenced system to an
                    # ellipsoid referenced system, input must be in km
                    data *= self.Re  # geodetic is defined in [km, deg, deg]
                NewCoords.data = np.atleast_2d(ctrans.convert_multitime(data, self.ticks,
                                                                        self.dtype, returnname,
                                                                        DEFAULTS.values))
                if returnname == 'RLL' and units[0] == 'Re':
                    # RLL is defined in Re, deg, deg, but all geodetic calcs are in km
                    NewCoords.data[:, 0] = NewCoords.data[:, 0]/self.Re
                if (self.dtype in ['RLL', 'GDZ'] and
                    returnname not in ['RLL', 'GDZ'] and
                    units[0] == 'km'):
                    # input was in km, needs to be Earth radii
                    NewCoords.data /= self.Re
                if returncarsph == 'sph' and returnname not in ['GDZ', 'RLL']:
                    NewCoords.data = car2sph(NewCoords.data)
            if self.dtype == 'GDZ':  # GDZ is defined in km and everything else is strictly Re
                units[0] = 'Re'  # Converting back from GDZ set units to Re
            NewCoords.carsph = returncarsph
            NewCoords.sysaxes = SYSAXES_TYPES[NewCoords.dtype][NewCoords.carsph]

        # Correct corresponding attributes and units
        if returncarsph == 'sph':
            if returnname == 'GDZ':
                NewCoords.units = ['km', 'deg', 'deg']
            else:
                NewCoords.units = [units[0], 'deg', 'deg']
            for k in ('x', 'y', 'z'):
                if hasattr(NewCoords, k):
                    delattr(NewCoords, k)
            NewCoords.radi = NewCoords.data[:, 0]
            NewCoords.lati = NewCoords.data[:, 1]
            NewCoords.long = NewCoords.data[:, 2]
        else:  # 'car'
            NewCoords.units = [units[0]]*3
            for k in ('radi', 'lati', 'long'):
                if hasattr(NewCoords, k):
                    delattr(NewCoords, k)
            NewCoords.x = NewCoords.data[:, 0]
            NewCoords.y = NewCoords.data[:, 1]
            NewCoords.z = NewCoords.data[:, 2]

        return NewCoords

    # -----------------------------------------------
    def __getstate__(self):
        '''
        Is called when pickling
        See Also http://docs.python.org/library/pickle.html
        '''
        odict = self.__dict__.copy()  # copy the dict since we change it
        return odict

    def __setstate__(self, dict):
        '''
        Is called when unpickling
        See Also http://docs.python.org/library/pickle.html
        '''
        self.__dict__.update(dict)
        return

    # -----------------------------------------------
    def append(self, other):
        '''Append another Coords instance to the current one

        Parameters
        ----------
        other : Coords instance
            Coords instance to append

        Examples
        --------
        '''
        data = list(self.data)
        otherdata = other.convert(self.dtype, self.carsph)
        data.extend(list(otherdata.data))
        newobj = Coords(data, dtype=self.dtype, carsph=self.carsph, use_irbem=self.use_irbem)
        return newobj

    # -----------------------------------------------
    def to_skycoord(self):
        '''
        Create an Astropy SkyCoord instance based on this instance

        Returns
        =======
        out : `astropy.coordinates.SkyCoord`
            This coordinate as an Astropy SkyCoord

        Notes
        =====
        This method requires Astropy to be installed.

        This method uses the GEO coordinate frame as the common frame between the two libraries.
        '''
        if self.ticks is None:
            raise ValueError("This method requires the attribute `ticks` to be specified.")

        try:
            import astropy.coordinates
            import astropy.time
        except ImportError as e:
            raise e.__class__("This method requires Astropy to be installed.")

        coord = self.convert('GEO', 'car')  # GEO is Astropy ITRS
        data = coord.data * self.Re * 1000 # convert to meters
        obstime = coord.ticks.APT

        return astropy.coordinates.SkyCoord(x=data[:, 0], y=data[:, 1], z=data[:, 2],
                                            unit='m', frame='itrs', obstime=obstime)

    # -----------------------------------------------
    @classmethod
    def from_skycoord(cls, skycoord, use_irbem=None):
        '''
        Create a Coords instance from an Astropy SkyCoord instance

        Parameters
        ==========
        skycoord : `astropy.coordinates.SkyCoord`
            The coordinate to be converted

        Returns
        =======
        out : Coords instance
            The converted coordinate

        Notes
        =====
        This method requires Astropy to be installed.

        This method uses the GEO coordinate frame as the common frame between the two libraries.
        '''
        try:
            import astropy.coordinates
        except ImportError as e:
            raise e.__class__("This method requires Astropy to be installed.")

        skycoord = astropy.coordinates.SkyCoord(skycoord).itrs  # Astropy ITRS is GEO
        # Ticktock constructor warns if use_irbem is None; just do right thing
        use_Re = IRBEM_RE if use_irbem else DEFAULTS.values.ellipsoid['A']  # kilometers
        data = (skycoord.cartesian.xyz.to('m').value / (1000 * use_Re)).T  # convert Cartesian to Re units
        ticks = spacepy.time.Ticktock(skycoord.obstime)

        return cls(data, 'GEO', 'car', ticks=ticks, use_irbem=use_irbem)


def car2sph(car_in):
    """
    Coordinate transformation from Cartesian to spherical

    Parameters
    ----------
        - car_in (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
            units of [Re, Re, Re] = [x,y,z]

    Returns
    -------
        - results (ndarray) : values after conversion to spherical coordinates in
            radius, latitude, longitude in units of [Re, deg, deg]

    Examples
    --------
    >>> sph2car([1,45,0])
    array([ 0.70710678,  0.        ,  0.70710678])

    See Also
    --------
    sph2car
    """
    if isinstance(car_in[0], numbers.Number):
        car = np.array([car_in])
    else:
        car = np.asanyarray(car_in)

    res = np.zeros(np.shape(car))
    for i in np.arange(len(car)):
        x, y, z = car[i, 0], car[i, 1], car[i, 2]
        r = np.sqrt(x*x+y*y+z*z)
        sq = np.sqrt(x*x+y*y)
        if (x == 0) & (y == 0):  # on the poles
            longi = 0.
            if z < 0:
                lati = -90.
            else:
                lati = 90.0
        else:
            longi = np.arctan2(y, x)*180./np.pi
            lati = 90. - np.arctan2(sq, z)*180./np.pi
        res[i, :] = [r, lati, longi]

    if isinstance(car_in[0], numbers.Number):
        return res[0]
    else:
        return res


def sph2car(sph_in):
    """
    Coordinate transformation from spherical to Cartesian

    Parameters
    ----------
        - sph_in (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
            units of [Re, deg, deg] = [r, latitude, longitude]

    Returns
    -------
        - results (ndarray) : values after conversion to cartesian coordinates x,y,z

    Examples
    --------
    >>> sph2car([1,45,45])
    array([ 0.5       ,  0.5       ,  0.70710678])

    See Also
    --------
    car2sph
    """
    if isinstance(sph_in[0], numbers.Number):
        sph = np.array([sph_in])
    else:
        sph = np.asanyarray(sph_in)

    res = np.zeros(np.shape(sph))
    for i in np.arange(len(sph)):
        r, lati, longi = sph[i, 0], sph[i, 1], sph[i, 2]
        colat = np.pi/2. - lati*np.pi/180.
        x = r*np.sin(colat)*np.cos(longi*np.pi/180.)
        y = r*np.sin(colat)*np.sin(longi*np.pi/180.)
        z = r*np.cos(colat)
        res[i, :] = [x, y, z]

    if isinstance(sph_in[0], numbers.Number):
        return res[0]
    else:
        return res


def quaternionNormalize(Qin, scalarPos='last'):
    '''
    Given an input quaternion (or array of quaternions), return the unit quaternion

    Parameters
    ----------
    vec : array_like
        input quaternion to normalize

    Returns
    -------
    out : array_like
        normalized quaternion

    Examples
    --------
    >>> import spacepy.coordinates
    >>> spacepy.coordinates.quaternionNormalize([0.707, 0, 0.707, 0.2])
    array([ 0.69337122,  0.        ,  0.69337122,  0.19614462])
    '''
    w = {'first': 0, 'last': 3}.get(scalarPos.lower())
    if w is None:
        errmsg = 'quaternionNormalize: scalarPos must be set to "First" or "Last"'
        raise NotImplementedError(errmsg)

    Quse = np.asanyarray(Qin).astype(float)
    squeeze = Quse.ndim < 2
    if squeeze:
        Quse = np.expand_dims(Quse, axis=0)
    magn = np.linalg.norm(Quse, axis=-1)
    # Find places where the magnitude is tiny, convert to unit real
    idx = np.where(magn <= 1e-12)
    magn[idx] = 1
    outarr = Quse / np.expand_dims(magn, axis=magn.ndim)
    outarr[idx + (None,)] = 0
    outarr[idx + (w,)] = 1
    if squeeze:
        outarr = outarr[0, ...]
    return outarr


def quaternionRotateVector(Qin, Vin, scalarPos='last', normalize=True):
    '''
    Given quaternions and vectors, return the vectors rotated by the quaternions

    Parameters
    ----------
    Qin : array_like
        input quaternion to rotate by
    Vin : array-like
        input vector to rotate

    Returns
    -------
    out : array_like
        rotated vector

    Examples
    --------
    >>> import spacepy.coordinates
    >>> import numpy as np
    >>> vec = [1, 0, 0]
    >>> quat_wijk = [np.sin(np.pi/4), 0, np.sin(np.pi/4), 0.0]
    >>> quat_ijkw = [0.0, np.sin(np.pi/4), 0, np.sin(np.pi/4)]
    >>> spacepy.coordinates.quaternionRotateVector(quat_ijkw, vec)
    array([ 0.,  0., -1.])
    >>> spacepy.coordinates.quaternionRotateVector(
    ...     quat_wijk, vec, scalarPos='first')
    array([ 0.,  0., -1.])

    See Also
    --------
    quaternionMultiply
    '''
    if scalarPos.lower() == 'last':
        i, j, k = 0, 1, 2
        w = 3
    elif scalarPos.lower() == 'first':
        i, j, k = 1, 2, 3
        w = 0
    else:
        errmsg = 'quaternionRotateVector: scalarPos must be set to "First" or "Last"'
        raise NotImplementedError(errmsg)

    Quse = np.asanyarray(Qin).astype(float)
    if normalize:
        Quse = quaternionNormalize(Quse, scalarPos=scalarPos)
    Vuse = np.asanyarray(Vin).astype(float)
    try:
        Quse.shape[1]
    except IndexError:
        Quse = np.asanyarray([Quse])
    try:
        Vuse.shape[1]
    except IndexError:
        Vuse = np.asanyarray([Vuse])
    try:
        assert Vuse.shape[0] == Quse.shape[0]
    except AssertionError:
        errmsg = ' '.join(['quaternionRotateVector:',
                           'Input vector array must have same length',
                           'as input quaternion array'])
        raise ValueError(errmsg)

    outarr = np.empty((Quse.shape[0], 3))
    for idx, row in enumerate(Quse):
        Vrow = Vuse[idx]
        ii, jj, kk, ww = row[i]**2, row[j]**2, row[k]**2, row[w]**2
        iw, jw, kw = row[i]*row[w], row[j]*row[w], row[k]*row[w]
        ij, ik, jk = row[i]*row[j], row[i]*row[k], row[j]*row[k]

        outarr[idx, 0] = ww*Vrow[0] + 2*jw*Vrow[2] - 2*kw*Vrow[1] + \
                         ii*Vrow[0] + 2*ij*Vrow[1] + 2*ik*Vrow[2] - \
                         kk*Vrow[0] - jj*Vrow[0]
        outarr[idx, 1] = 2*ij*Vrow[0] + jj*Vrow[1] + 2*jk*Vrow[2] + \
                         2*kw*Vrow[0] - kk*Vrow[1] + ww*Vrow[1] - \
                         2*iw*Vrow[2] - ii*Vrow[1]
        outarr[idx, 2] = 2*ik*Vrow[0] + 2*jk*Vrow[1] + kk*Vrow[2] - \
                         2.0*jw*Vrow[0] - jj*Vrow[2] + 2.0*iw*Vrow[1] - \
                         ii*Vrow[2] + ww*Vrow[2]
    return outarr.squeeze()


def quaternionMultiply(Qin1, Qin2, scalarPos='last'):
    '''
    Given quaternions, return the product, i.e. Qin1*Qin2

    Parameters
    ----------
    Qin1 : array_like
        input quaternion, first position
    Qin2 : array-like
        input quaternion, second position

    Returns
    -------
    out : array_like
        quaternion product

    Examples
    --------
    >>> import spacepy.coordinates
    >>> import numpy as np
    >>> vecX = [1, 0, 0] #shared X-axis
    >>> vecZ = [0, 0, 1] #unshared, but similar, Z-axis
    >>> quat_eci_to_gsm = [-0.05395384,  0.07589845, -0.15172533,  0.98402634]
    >>> quat_eci_to_gse = [ 0.20016056,  0.03445775, -0.16611386,  0.96496352]
    >>> quat_gsm_to_eci = spacepy.coordinates.quaternionConjugate(
    ...     quat_eci_to_gsm)
    >>> quat_gse_to_gsm = spacepy.coordinates.quaternionMultiply(
    ...     quat_gsm_to_eci, quat_eci_to_gse)
    >>> spacepy.coordinates.quaternionRotateVector(quat_gse_to_gsm, vecX)
    array([  1.00000000e+00,   1.06536725e-09,  -1.16892107e-08])
    >>> spacepy.coordinates.quaternionRotateVector(quat_gse_to_gsm, vecZ)
    array([  1.06802834e-08,  -4.95669027e-01,   8.68511494e-01])
    '''
    if scalarPos.lower() == 'last':
        i, j, k = 0, 1, 2
        w = 3
    elif scalarPos.lower() == 'first':
        i, j, k = 1, 2, 3
        w = 0
    else:
        raise NotImplementedError('quaternionMultiply: scalarPos must be set to "First" or "Last"')

    Quse1 = np.asanyarray(Qin1).astype(float)
    Quse2 = np.asanyarray(Qin2).astype(float)
    try:
        Quse1.shape[1]
    except IndexError:
        Quse1 = np.asanyarray([Quse1])
    try:
        Quse2.shape[1]
    except IndexError:
        Quse2 = np.asanyarray([Quse2])
    try:
        assert Quse2.shape == Quse1.shape
    except AssertionError:
        raise ValueError('quaternionMultiply: Input quaternion arrays must have same length')

    outarr = np.empty_like(Quse1)
    for idx, row in enumerate(Quse1):
        row1 = row
        row2 = Quse2[idx]
        # vector components
        outarr[idx, i] = row1[w]*row2[i] + row1[i]*row2[w] + row1[j]*row2[k] - row1[k]*row2[j]
        outarr[idx, j] = row1[w]*row2[j] - row1[i]*row2[k] + row1[j]*row2[w] + row1[k]*row2[i]
        outarr[idx, k] = row1[w]*row2[k] + row1[i]*row2[j] - row1[j]*row2[i] + row1[k]*row2[w]
        # real part
        outarr[idx, w] = row1[w]*row2[w] - row1[i]*row2[i] - row1[j]*row2[j] - row1[k]*row2[k]
    return outarr.squeeze()


def quaternionConjugate(Qin, scalarPos='last'):
    '''
    Given an input quaternion (or array of quaternions), return the conjugate

    Parameters
    ----------
    Qin : array_like
        input quaternion to conjugate

    Returns
    -------
    out : array_like
        conjugate quaternion

    Examples
    --------
    >>> import spacepy.coordinates
    >>> spacepy.coordinates.quaternionConjugate(
    ...     [0.707, 0, 0.707, 0.2], scalarPos='last')
    array([-0.707, -0.   , -0.707,  0.2  ])

    See Also
    --------
    quaternionMultiply
    '''
    if scalarPos.lower() == 'last':
        i, j, k = 0, 1, 2
        w = 3
    elif scalarPos.lower() == 'first':
        i, j, k = 1, 2, 3
        w = 0
    else:
        errmsg = 'quaternionConjugate: scalarPos must be set to "First" or "Last"'
        raise NotImplementedError(errmsg)

    Quse = np.asanyarray(Qin).astype(float)
    try:
        Quse.shape[1]
    except IndexError:
        Quse = np.asanyarray([Quse])
    outarr = np.empty_like(Quse)
    for idx, row in enumerate(Quse):
        outarr[idx, i] = -row[i]
        outarr[idx, j] = -row[j]
        outarr[idx, k] = -row[k]
        outarr[idx, w] = row[w]

    return outarr.squeeze()


def quaternionFromMatrix(matrix, scalarPos='last'):
    '''
    Given an input rotation matrix, return the equivalent quaternion

    The output has one fewer axis than the input (the last axis) and the
    shape is otherwise unchanged, allowing multi-dimensional matrix input.

    Parameters
    ----------
    matrix : array_like
        input rotation matrix or array of matrices

    Returns
    -------
    out : array_like
        Quaternions representing the same rotation as the input rotation
        matrices.

    Other Parameters
    ----------------
    scalarPos : str
        Location of the scalar component of the output quaternion, either
        'last' (default) or 'first'.

    Raises
    ------
    NotImplementedError
        for invalid values of ``scalarPos``

    ValueError
        for inputs which are obviously not valid 3D rotation matrices or
        arrays thereof: if the size doesn't end in (3, 3), if the matrix is
        not orthogonal, or not a proper rotation.

    See Also
    --------
    quaternionToMatrix

    Notes
    -----
    .. versionadded:: 0.2.2

    No attempt is made to resolve the sign ambiguity; in particular,
    conversions of very similar matrices may result in equivalent quaternions
    with the opposite sign. This may have implications for interpolating a
    sequence of quaternions.

    The conversion of a rotation matrix to a quaternion suffers from some
    of the same disadvantages inherent to rotation matrices, such as potential
    numerical instabilities. Working in quaternion space as much as possible
    is recommended.

    There are several algorithms; the most well-known algorithm for this
    conversion is Shepperd's [#Shepperd]_, although the many "rediscoveries"
    indicate
    it is not sufficiently well-known. This function uses the method of
    Bar-Itzhack [#BarItzhack]_ (version 3), which should be resistant to small
    errors in the rotation matrix. As a result, the input checking is quite
    coarse and will likely accept many matrices that do not represent valid
    rotations.

    Also potentially of interest, although not implemented here, is Sarabandi
    and Thomas [#Sarabandi]_.

    References
    ----------
    .. [#Shepperd] S.W. Shepperd, "Quaternion from rotation matrix," Journal of
            Guidance and Control, Vol. 1, No. 3, pp. 223-224, 1978,
            `doi:10.2514/3.55767b <https://doi.org/10.2514/3.55767b>`_

    .. [#BarItzhack] I. Y. Bar-Itzhack, "New method for extracting the
            quaternion from a rotation matrix", AIAA Journal of Guidance,
            Control and Dynamics, 23 (6): 1085–1087,
            `doi:10.2514/2.4654 <https://doi.org/10.2514/2.4654>`_

    .. [#Sarabandi] S. Sarabandi and F. Thomas, "Accurate Computation of
            Quaternions from Rotation Matrices", In: Lenarcic J.,
            Parenti-Castelli V. (eds) Advances in Robot Kinematics 2018,
            Springer.
            `doi:10.1007/978-3-319-93188-3_5
            <https://doi.org/10.1007/978-3-319-93188-3_5>`_

    Examples
    --------
    >>> import spacepy.coordinates
    >>> spacepy.coordinates.quaternionFromMatrix(
    ...     [[ 0.,  0.,  1.],
    ...      [ 1.,  0.,  0.],
    ...      [ 0.,  1.,  0.]])
    array([0.5, 0.5, 0.5, 0.5])
    '''
    if scalarPos.lower() not in ('last', 'first'):
        raise NotImplementedError(
            'quaternionFromMatrix: scalarPos must be set to "First" or "Last"')
    matrix = np.asanyarray(matrix)
    if len(matrix.shape) < 2 or matrix.shape[-2:] != (3, 3):
        raise ValueError(
            'Input does not appear to be 3D rotation matrix, wrong size.')
    for i in np.ndindex(matrix.shape[:-2]):
        m = matrix[i + (Ellipsis,)]
        if not np.allclose(np.dot(m, m.transpose()),
                           np.identity(3), atol=0.2):
            raise ValueError('Input rotation matrix{} not orthogonal.'.format(
                '' if len(matrix.shape) == 2 else ' at ' + str(i)))
        det = np.linalg.det(m)
        if det < 0:
            raise ValueError('Input rotation matrix at {} not proper.'
                             .format(str(i)))
    inshape = matrix.shape
    # Flatten out most dimensions, easier to index
    matrix = np.reshape(matrix, (-1, 3, 3))
    # Indexing in Bar-Itzhack is reversed relative to numpy
    k = (matrix[..., 0, 0] - matrix[..., 1, 1] - matrix[..., 2, 2],
         matrix[..., 0, 1] + matrix[..., 1, 0],
         matrix[..., 0, 2] + matrix[..., 2, 0],
         matrix[..., 2, 1] - matrix[..., 1, 2],  # row 0
         matrix[..., 0, 1] + matrix[..., 1, 0],
         matrix[..., 1, 1] - matrix[..., 0, 0] - matrix[..., 2, 2],
         matrix[..., 1, 2] + matrix[..., 2, 1],
         matrix[..., 0, 2] - matrix[..., 2, 0],  # row 1
         matrix[..., 0, 2] + matrix[..., 2, 0],
         matrix[..., 1, 2] + matrix[..., 2, 1],
         matrix[..., 2, 2] - matrix[..., 0, 0] - matrix[..., 1, 1],
         matrix[..., 1, 0] - matrix[..., 0, 1],  # row 2
         matrix[..., 2, 1] - matrix[..., 1, 2],
         matrix[..., 0, 2] - matrix[..., 2, 0],
         matrix[..., 1, 0] - matrix[..., 0, 1],
         matrix[..., 0, 0] + matrix[..., 1, 1] + matrix[..., 2, 2]  # row 3
         )
    # Stack together on a new axis after the input sample, then
    # split that axis into the 4x4 k arrays
    k = np.reshape(np.stack(k, axis=1), (-1, 4, 4))
    k = k / 3.
    evals, evects = np.linalg.eig(k)
    evects = np.real(evects)
    idx = np.argmax(evals, axis=-1)
    Qout = evects[np.mgrid[0:len(idx)], :, idx]
    if scalarPos.lower() == 'first':
        Qout = np.roll(Qout, 1, axis=-1)
    # Recover original shape
    Qout = np.reshape(Qout, inshape[:-2] + (4,))
    return Qout


def quaternionToMatrix(Qin, scalarPos='last', normalize=True):
    '''
    Given an input quaternion, return the equivalent rotation matrix.

    The output has one more axis than the input (the last axis) and the
    shape is otherwise unchanged, allowing multi-dimensional quaternion input.

    Parameters
    ----------
    Qin : array_like
        input quaternion or array of quaternions, must be normalized.

    Returns
    -------
    out : array_like
        Rotation matrix

    Other Parameters
    ----------------
    scalarPos : str
        Location of the scalar component of the input quaternion, either
        'last' (default) or 'first'.

    normalize : True
        Normalize input quaternions before conversion (default). If False,
        raises error for non-normalized.

    Raises
    ------
    NotImplementedError
        for invalid values of ``scalarPos``.

    ValueError
        for inputs which are not valid normalized quaternions or arrays
        thereof: if the size doesn't end in (4), if the quaternion is not
        normalized and ``normalize`` is False.

    See Also
    --------
    quaternionFromMatrix

    Notes
    -----
    .. versionadded:: 0.2.2

    Implementation of the Euler–Rodrigues formula.

    Examples
    --------
    >>> import spacepy.coordinates
    >>> spacepy.coordinates.quaternionToMatrix([0.5, 0.5, 0.5, 0.5])
    array([[ 0.,  0.,  1.],
           [ 1.,  0.,  0.],
           [ 0.,  1.,  0.]])
    '''
    if scalarPos.lower() not in ('last', 'first'):
        raise NotImplementedError(
            'quaternionToMatrix: scalarPos must be set to "First" or "Last"')
    Qin = np.asanyarray(Qin)
    if Qin.shape[-1] != 4:
        raise ValueError('Input does not appear to be quaternion, wrong size.')
    if normalize:
        Qin = quaternionNormalize(Qin, scalarPos=scalarPos)
    if scalarPos.lower() == 'first':
        Qin = np.roll(Qin, -1, axis=-1)
    if not np.allclose(np.sum(Qin ** 2, axis=-1), 1):
        raise ValueError('Input quaternion not normalized.')
    # Maintain dimensions at end for stacking
    a, b, c, d = Qin[..., -1:], Qin[..., 0:1], Qin[..., 1:2], Qin[..., 2:3]
    # Output array, "flattened"
    # This can probably be written with a clever vector notation...
    out = np.concatenate((a ** 2 + b ** 2 - c ** 2 - d ** 2,
                          2 * (b * c - a * d),
                          2 * (b * d + a * c),
                          2 * (b * c + a * d),
                          a ** 2 + c ** 2 - b ** 2 - d ** 2,
                          2 * (c * d - a * b),
                          2 * (b * d - a * c),
                          2 * (c * d + a * b),
                          a ** 2 + d ** 2 - b ** 2 - c ** 2,
                          ), axis=-1)
    # And last two dims are 3x3 array for matrix
    out = np.reshape(out, out.shape[:-1] + (3, 3))
    return out
