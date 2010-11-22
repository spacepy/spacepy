#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of Coords class functions

"""

from spacepy import help
__version__ = "$Revision: 1.9 $, $Date: 2010/11/22 22:01:31 $"
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'


# -----------------------------------------------
# space coordinate class
# -----------------------------------------------    
class Coords(object):
    """
    a = Coords( data, dtype, carsph, [units, ticks] )
    
    A class holding spatial coordinates in Cartesian/spherical
    in units of Re and degrees
        
    Input:
    ======
        - data (list or ndarray, dim = (n,3) ) : coordinate points 
        - dtype (string) :coordinate system, possible are GDZ, GEO, GSM, GSE, SM, GEI
                MAG, SPH, RLL
        - carsph (string) : Cartesian or spherical, 'car' or 'sph'
        - optional units (list of strings) : standard are  ['Re', 'Re', 'Re'] or 
            ['Re', 'deg', 'deg'] depending on the carsph content
        - optional ticks (Ticktock instance) : used for coordinate transformations (see a.convert)    
        
    Returns:
    ========
        - instance with a.data, a.carsph, etc.
    
    Example:
    ========  
    >>> cvals = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
    >>> cvals.x  # returns all x coordinates
    array([1, 1])
    >>> cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticks
    >>> newcoord = cvals.convert('GSM', 'sph')
    >>> newcoord
    
       
    See Also:
    =========
    time.Ticktock class
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 05-Mar-2010 (JK)
    
    """
    
    def __init__(self, data, dtype, carsph, units=None, ticks=None):
        
        import numpy as n
        from . import irbempy as op
        from spacepy.irbempy import SYSAXES_TYPES as typedict
        
        if isinstance(data[0], (float, int)):
            self.data = n.array([data])
        else:
            self.data = n.array(data)
        
        assert dtype in typedict.keys(), 'This dtype='+dtype+' is not supported. Only '+str(typedict.keys())
        assert carsph in ['car','sph'], 'This carsph='+str(carsph)+' is not supported. Only "car" or "sph"'
        onerawarn = """Coordinate conversion to an ONERA-compatible system is required for any ONERA calls."""
        
        # add ticks
        if ticks: assert len(ticks) == len(data), 'Ticktock dimensions seem off'
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
        if units==None and carsph=='car':
            # use standard units
            self.units = ['Re', 'Re', 'Re']
        elif units==None and carsph=='sph':
            self.units = ['Re', 'deg', 'deg']
        else:
            self.units = units
        if dtype == 'GDZ' and carsph == 'sph':
            self.units = ['km', 'deg', 'deg']
        # setup x,y,z etc
        if carsph == 'car':
            self.x = self.data[:,0]
            self.y = self.data[:,1]
            self.z = self.data[:,2]
        else:
            self.radi = self.data[:,0]
            self.lati = self.data[:,1]
            self.long = self.data[:,2]
        ## setup list for onera compatibility
        #self.sysaxes = op.get_sysaxes(dtype, carsph)
        self.shape = n.shape(self.data)
        return None
        
    # -----------------------------------------------    
    def __str__(self):
        """
        a.__str__() or a
        
        Will be called when printing Coords instance a
        
        Input:
        ======
            - a Coords class instance
 
        Returns:
        ========
            - output (string)          

        Example:
        ========
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y
        Coords( [[1 2 4]
         [1 2 2]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        """
 
        return 'Coords( '+str(self.data) + ' ), dtype='+self.dtype+','+self.carsph+', units='+\
            str(self.units)
    __repr__ = __str__
    
    # -----------------------------------------------    
    def __getitem__(self, idx):
        """
        a.__getitem__(idx) or a[idx]
        
        Will be called when requesting items in this instance 
        
        Input:
        ======
            - a Coords class instance
            - idx (int) : integer numbers as index

        Returns:
        ========
            - vals (numpy array) : new values         

        Example:
        ========
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[0] 
        array([1, 2, 4])
 
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        V1.1: 19-Apr-2010: now returns a Coords instance
 
        """
        import numpy as n
        
        arr = n.array(self.data)
        
        return Coords(arr[idx].tolist(), self.dtype, self.carsph, self.units, self.ticks)   
        
    # -----------------------------------------------    
    def __setitem__(self, idx, vals):
        """
        a.__setitem__(idx, vals) or a[idx] = vals
        
        Will be called setting items in this instance 
        
        Input:
        ======
            - a Coords class instance
            - idx (int) : integer numbers as index
            - vals (numpy array or list) : new values
            
        Example:
        ========
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[1] = [9,9,9]
        >>> y
        Coords( [[1 2 4]
         [9 9 9]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
 
        """
        self.data[idx] = vals
        return
        
    # -----------------------------------------------    
    def __len__(self):
        """
        a.__len__() or len(a)
        
        Will be called when requesting the length, i.e. number of items 
        
        Input:
        ======
            - a Coords class instance
            
        Returns:
        ========
            - length (int number)

        Example:
        ========
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> len(y)
        2

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
 
        """
         
        from numpy import ndarray 
        if isinstance(self.data, (list, ndarray)):
            return len(self.data)
        else:
            return 1
        return  
               
    # -----------------------------------------------    
    def convert(self, returntype, returncarsph):
        """
        a.convert( returntype, returncarsph )
        
        Can be used to create a new Coords instance with new coordinate types
        
        Input:
        ======
            - a Coords class instance
            - returntype (string) : coordinate system, possible are GDZ, GEO, GSM, GSE, SM, GEI
                MAG, SPH, RLL
            - returncarsph (string) : coordinate type, possible 'car' for Cartesian and 
                'sph' for spherical
 
        Returns:
        ========
            - a Coords object          

        Example:
        ========
        >>> y = Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> x = y.convert('SM','car')
        >>> x
        Coords( [[ 0.81134097  2.6493305   3.6500375 ]
         [ 0.92060408  2.30678864  1.68262126]] ), dtype=SM,car, units=['Re', 'Re', 'Re']


        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        
        """
        
        import numpy as n
        from . import irbempy as op
        import spacepy, spacepy.coordinates

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
        if (self.sysaxes == None) : # a car2sph or sph2car is needed            
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
        
        Coords = spacepy.coordinates.Coords(data, self.dtype, carsph, units, self.ticks)
        
        # now convert to other coordinate system
        if (self.dtype != returntype) : 
            assert Coords.ticks, "Time information required; add a.ticks attribute"
            Coords.data = op.coord_trans( Coords, returntype, returncarsph)
            Coords.dtype = returntype
            Coords.carsph = returncarsph
            Coords.sysaxes = op.get_sysaxes(returntype, returncarsph)

        # fix corresponding attributes         
        if returncarsph == 'sph':
            Coords.units = [units[0], 'deg','deg']
            if Coords.__dict__.has_key('x'): Coords.__delattr__('x')
            if Coords.__dict__.has_key('y'): Coords.__delattr__('y')
            if Coords.__dict__.has_key('z'): Coords.__delattr__('z')
            Coords.radi = Coords.data[:,0]
            Coords.lati = Coords.data[:,1]
            Coords.long = Coords.data[:,2]
        else: # 'car'
            Coords.units =  [units[0]]*3
            if Coords.__dict__.has_key('radi'): Coords.__delattr__('radi')
            if Coords.__dict__.has_key('lati'): Coords.__delattr__('lati')
            if Coords.__dict__.has_key('long'): Coords.__delattr__('long')
            Coords.x = Coords.data[:,0]
            Coords.y = Coords.data[:,1]
            Coords.z = Coords.data[:,2]
       
        return Coords

    # -----------------------------------------------    
    def __getstate__(self):
        """
        Is called when pickling
        Author: J. Koller
        See Also: http://docs.python.org/library/pickle.html
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        return odict

    def __setstate__(self, dict):
        """
        Is called when unpickling
        Author: J. Koller
        See Also: http://docs.python.org/library/pickle.html
        """
        self.__dict__.update(dict)
        return
    
    # -----------------------------------------------    
    def append(self, other):
        """
        a.append(other)
        
        Will be called when another Coords instance has to be appended to the current one

        Input:
        ======
            - a Coords class instance
            - other (Coords instance) 
      
        Example:
        ========
      
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 15-Sep-2010 (JK)
       
        """
        
        import numpy as n
        data = list(self.data)
        otherdata = other.convert(self.dtype, self.carsph)               
        data.extend(list(otherdata.data))
        newobj = Coords(data, dtype=self.dtype, carsph=self.carsph)
        return newobj
    
