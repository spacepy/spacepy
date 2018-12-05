#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of Coords class functions for coordinate transformations

Authors: Josef Koller and Steven Morley
Institution: Los ALamos National Laboratory
Contact: smorley@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.
"""

import numpy as np
from spacepy import help
import spacepy

__contact__ = 'Steven Morley, smorley@lanl.gov'

# -----------------------------------------------
# space coordinate class
# -----------------------------------------------
class Coords(object):
    '''
    a = Coords( data, dtype, carsph, [units, ticks] )

    A class holding spatial coordinates in Cartesian/spherical
    in units of Re and degrees

    Coordinate transforms are based on the IRBEM library; `its manual
    <http://svn.code.sf.net/p/irbem/code/tags/IRBEM-4.4.0/manual/user_guide.html>`_
    may prove useful. For a good reference on heliospheric and magnetospheric
    coordinate systems, see Franz & Harper, "Heliospheric Coordinate Systems",
    Planet. Space Sci., 50, pp 217-233, 2002.

    Parameters
    ==========
    data : list or ndarray, dim = (n,3)
        coordinate points [X,Y,Z] or [rad, lat, lon]
    dtype : string
        coordinate system; possible values are:

        * **GDZ** (Geodetic – WGS84),
        * **GEO** (Geographic Coordinate System),
        * **GSM** (Geocentric Solar Magnetospheric),
        * **GSE** (Geocentric Solar Ecliptic),
        * **SM** (Solar Magnetic),
        * **GEI** (Geocentric Equatorial Inertial),
        * **MAG** (Geomagnetic Coordinate System),
        * **SPH** (Spherical Coordinate System),
        * **RLL** (Radius, Latitude, Longitude)

    carsph : string
        Cartesian or spherical, 'car' or 'sph'
    units : list of strings, optional
        standard are  ['Re', 'Re', 'Re'] or ['Re', 'deg', 'deg'] depending on the carsph content
    ticks : Ticktock instance, optional
        used for coordinate transformations (see a.convert)

    Returns
    =======
    out : Coords instance
        instance with a.data, a.carsph, etc.

    See Also
    ========
    spacepy.time.Ticktock

    Examples
    ========
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
    .. automethod:: append
    .. automethod:: convert
    '''
    def __init__(self, data, dtype, carsph, units=None, ticks=None):

        from . import irbempy as op
        from spacepy.irbempy import SYSAXES_TYPES as typedict

        if isinstance(data[0], (float, int)):
            self.data = np.array([data])
        else:
            self.data = np.array(data)

        assert dtype in list(typedict.keys()), 'This dtype='+dtype+' is not supported. Only '+str(list(typedict.keys()))
        assert carsph in ['car','sph'], 'This carsph='+str(carsph)+' is not supported. Only "car" or "sph"'
        onerawarn = """Coordinate conversion to an ONERA-compatible system is required for any ONERA calls."""

        # add ticks
        if ticks:
            assert len(ticks) == len(self.data), 'Ticktock dimensions seem off'
        self.ticks = ticks

        # GEO,sph and SPH,sph are the same
        if dtype == 'GEO' and carsph == 'sph':
            dtype = 'SPH'
        self.sysaxes = typedict[dtype][carsph]

        #if self.sysaxes >= 10 and self.sysaxes < 20: #need sph2car
        #    try:
        #        self.data = op.sph2car(self.data)
        #        self.sysaxes -= 10
        #    except:
        #        print onerawarn
        #        self.sysaxes = None
        #if self.sysaxes >= 20: #need car2sph
        #    try:
        #        self.data = op.car2sph(self.data)
        #        self.sysaxes -= 20
        #    except:
        #        print onerawarn
        #        self.sysaxes = None

        self.dtype = dtype
        self.carsph = carsph
        # setup units
        self.Re = 6371000.0 #metres
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
        ## setup list for onera compatibility
        #self.sysaxes = op.get_sysaxes(dtype, carsph)
        self.shape = np.shape(self.data)
        return None

    # -----------------------------------------------
    def __str__(self):
        '''
        a.__str__() or a

        Will be called when printing Coords instance a

        Returns
        ========
        out : string
            string represenation of the instance

        Examples
        ========
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
        ==========
        idx : int
            integer numbers as index

        Returns
        =======
        out : numpy array
            new values

        Examples
        ========
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[0]
        Coords( [[1 2 4]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']
        '''
        arr = np.array(self.data)
        t_select = self.ticks[idx] if self.ticks else self.ticks

        return Coords(arr[idx].tolist(), self.dtype, self.carsph, self.units, t_select)

    # -----------------------------------------------
    def __setitem__(self, idx, vals):
        '''
        a.__setitem__(idx, vals) or a[idx] = vals

        Will be called setting items in this instance

        Parameters
        ==========
        idx : int
            integer numbers as index
        vals : numpy array or list
            new values

        Examples
        ========
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
        ========
        out : int
            length

        Examples
        ========
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
        ==========
        returntype : string
            coordinate system, possible are GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
        returncarsph : string
            coordinate type, possible 'car' for Cartesian and 'sph' for spherical

        Returns
        =======
        out : Coords object
            Coords object in the new coordinate system

        Examples
        ========
        >>> from spacepy.coordinates import Coords
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> from spacepy.time import Ticktock
        >>> y.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> x = y.convert('SM','car')
        >>> x
        Coords( [[ 0.81134097  2.6493305   3.6500375 ]
         [ 0.92060408  2.30678864  1.68262126]] ), dtype=SM,car, units=['Re', 'Re', 'Re']
        '''
        from . import irbempy as op
        from spacepy.irbempy import SYSAXES_TYPES as typedict

        #check return type/system is supported
        #if not (typedict[returntype][returncarsph]):
        #    raise NotImplementedError('System {0} is not supported in {1} coordinates'.format(returntype, returncarsph))

        # no change necessary
        if (self.dtype == returntype) and (self.carsph == returncarsph):
            return self

        # only car2sph or sph2car is needed
        if (self.dtype == returntype) and (self.carsph != returncarsph):
            if returncarsph == 'car':
                carsph = 'car'
                units = [self.units[0]]*3
                data = op.sph2car(self.data)
            else:
                carsph = 'sph'
                units =  [self.units[0], 'deg','deg']
                data = op.car2sph(self.data)
            return Coords(data, self.dtype, carsph, units, self.ticks)

        # check the length of ticks and do the more complex conversions
        if self.ticks:
            assert len(self.ticks) == len(self), 'Ticktock dimension does not match Coords dimensions'

        # check if car2sph is needed first for oneralib compatibility
        if (self.sysaxes is None) : # a car2sph or sph2car is needed
            if self.carsph == 'sph':
                carsph = 'car'
                units = [self.units[0]]*3
                data = op.sph2car(self.data)
            else:
                carsph = 'sph'
                units =  [self.units[0], 'deg','deg']
                data = op.car2sph(self.data)
        else:
            data = self.data
            units = self.units
            carsph = self.carsph

        NewCoords = Coords(data, self.dtype, carsph, units, self.ticks)

        # now convert to other coordinate system
        if (self.dtype != returntype) :
            assert NewCoords.ticks, "Time information required; add a.ticks attribute"
            NewCoords.data = op.coord_trans( NewCoords, returntype, returncarsph)
            NewCoords.dtype = returntype
            NewCoords.carsph = returncarsph
            NewCoords.sysaxes = op.get_sysaxes(returntype, returncarsph)

        # fix corresponding attributes
        if returncarsph == 'sph':
            NewCoords.units = [units[0], 'deg','deg']
            for k in ('x', 'y', 'z'):
                if hasattr(NewCoords, k):
                    delattr(NewCoords, k)
            NewCoords.radi = NewCoords.data[:,0]
            NewCoords.lati = NewCoords.data[:,1]
            NewCoords.long = NewCoords.data[:,2]
        else: # 'car'
            NewCoords.units =  [units[0]]*3
            for k in ('radi', 'lati', 'long'):
                if hasattr(NewCoords, k):
                    delattr(NewCoords, k)
            NewCoords.x = NewCoords.data[:,0]
            NewCoords.y = NewCoords.data[:,1]
            NewCoords.z = NewCoords.data[:,2]

        return NewCoords

    # -----------------------------------------------
    def __getstate__(self):
        '''
        Is called when pickling
        See Also http://docs.python.org/library/pickle.html
        '''
        odict = self.__dict__.copy() # copy the dict since we change it
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
        ==========
        other : Coords instance
            Coords instance to append

        Examples
        ========
        '''
        data = list(self.data)
        otherdata = other.convert(self.dtype, self.carsph)
        data.extend(list(otherdata.data))
        newobj = Coords(data, dtype=self.dtype, carsph=self.carsph)
        return newobj
