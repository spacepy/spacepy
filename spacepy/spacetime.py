#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Implementation of TickTock class and SpaCo class functions

"""

from spacepy import help
__version__ = "$Revision: 1.1.1.1 $, $Date: 2010/05/20 17:19:44 $"
__log__ = """
$Log: spacetime.py,v $
Revision 1.1.1.1  2010/05/20 17:19:44  smorley
Pre-release repo import for SpacePy

Revision 1.19  2010/05/20 16:27:31  jkoller
moved testing routines into the specific modules

Revision 1.18  2010/05/19 22:30:50  smorley
Regenerated documentation

Revision 1.16  2010/05/19 00:27:28  smorley
Further mods to fix sub-second support - not yet complete

Revision 1.15  2010/05/18 18:44:15  smorley
Added sub-second support for spacetime class

Revision 1.14  2010/05/17 17:38:51  smorley
General updates

Revision 1.13  2010/05/11 15:54:20  smorley
Fixed cvs markup from previous version conflict

Revision 1.12  2010/05/11 15:53:12  smorley
Updated SpaCo class for smarter coordinate system handling

Revision 1.11  2010/05/10 22:03:58  balarsen
Missed a leap lear conversion, fixed and should be correct.

Revision 1.10  2010/05/10 21:57:33  balarsen
added GPS epoch to spacetime, also changed a few numpy.arange to the much faster and built-in xrange.

Revision 1.9  2010/04/28 20:55:26  smorley
Updates to radbelt documentation, spacetime, __init__ and setup

Revision 1.8  2010/04/20 23:17:06  jkoller
fixed a little bug in eDOY

Revision 1.7  2010/04/20 18:06:55  jkoller
added routine for elapsed DOY, called eDOY; DOY itself will return integer number according to convention

Revision 1.6  2010/04/12 17:54:55  smorley
Updated doy2date for datetime output and readepoch in seapy to accept any format

Revision 1.5  2010/04/02 23:55:41  smorley
Added datetime object output to doy2date in spacetime [SKM]

Revision 1.4  2010/03/16 15:14:09  jkoller
fixed floating point precision bug in tickrange and added __mul__ to TickDelta

Revision 1.3  2010/03/11 20:58:53  jkoller
fixed problem with timedelta; won't allow numpy.dtype('int64') numbers

Revision 1.2  2010/03/10 22:01:41  jkoller
added a tickrange function

Revision 1.1  2010/03/05 23:46:24  jkoller
First version of spacetime.py providing classes: SpaCo, TickTock, etc

"""
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'


# -----------------------------------------------
# space coordinate class
# -----------------------------------------------    
class SpaCo(object):
    """
    a = Spaco( data, dtype, carsph, [units, ticktock] )
    
    A class holding spatial coordinates in cartesian/spherical
    in units of Re and degrees
        
    Input:
    ======
        - data (list or ndarray, dim = (n,3) ) : coordinate points 
        - dtype (string) :coordinate system, possible are GDZ, GEO, GSM, GSE, SM, GEI
                MAG, SPH, RLL
        - carsph (string) : cartesian or spherical, 'car' or 'sph'
        - optional units (list of strings) : standard are  ['Re', 'Re', 'Re'] or 
            ['Re', 'deg', 'deg'] depending on the carsph content
        - optional ticktock (TickTock instance) : used for coordinate transformations (see a.convert)    
        
    Returns:
    ========
        - instance with a.data, a.carsph, etc.
    
    Example:
    ========  
    >>> coord = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
    >>> coord.x  # returns all x coordinates
    array([1, 1])
    >>> coord.ticktock = TickTock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
    >>> newcoord = coord.convert('GSM', 'sph')
    >>> newcoord
    
       
    See Also:
    =========
    spacetime.TickTock class
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 05-Mar-2010 (JK)
    
    """
    
    def __init__(self, data, dtype, carsph, units=None, ticktock=None):
        import numpy as n
        import onerapy as op
        if isinstance(data[0], (float, int)):
            self.data = n.array([data])
        else:
            self.data = n.array(data)
        
        typedict = {'GDZ': {'sph': 0, 'car': 10},
            'GEO': {'sph': 11, 'car': 1}, 'GSM': {'sph': 22, 'car': 2},
            'GSE': {'sph': 23, 'car': 3}, 'SM': {'sph': 24, 'car': 4},
            'GEI': {'sph': 25, 'car': 5}, 'MAG': {'sph': 26, 'car': 6},
            'SPH': {'sph': 7, 'car': 17}, 'RLL': {'sph': 8, 'car': 18}}
        assert dtype in typedict.keys(), 'This dtype='+dtype+' is not supported. Only '+str(typedict.keys())
        assert carsph in ['car','sph'], 'This carsph='+str(carsph)+' is not supported. Only "car" or "sph"'
        onerawarn = """Coordinate conversion to an ONERA-compatible system is required for any ONERA calls."""
        
        # add ticktock
        if ticktock: assert len(ticktock) == len(data), 'TickTock dimensions seem off'
        self.ticktock = ticktock
        
        self.sysaxes = typedict[dtype][carsph]
        if self.sysaxes >= 10 and self.sysaxes < 20: #need sph2car
            try:
                self.data = op.sph2car(self.data)
                self.sysaxes -= 10
            except:
                print onerawarn
                self.sysaxes = None
        if self.sysaxes >= 20: #need car2sph
            try:
                self.data = op.car2sph(self.data)
                self.sysaxes -= 20
            except:
                print onerawarn
                self.sysaxes = None
                
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
        
        return None
        
    # -----------------------------------------------    
    def __str__(self):
        """
        a.__str__() or a
        
        Will be called when printing SpaCo instance a
        
        Input:
        ======
            - a SpaCo class instance
 
        Returns:
        ========
            - output (string)          

        Example:
        ========
        >>> y = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y
        SpaCo( [[1 2 4]
         [1 2 2]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        """
 
        return 'SpaCo( '+str(self.data) + ' ), dtype='+self.dtype+','+self.carsph+', units='+\
            str(self.units)
    __repr__ = __str__
    
    # -----------------------------------------------    
    def __getitem__(self, idx):
        """
        a.__getitem__(idx) or a[idx]
        
        Will be called when requesting items in this instance 
        
        Input:
        ======
            - a SpaCo class instance
            - idx (int) : interger numbers as index

        Returns:
        ========
            - vals (numpy array) : new values         

        Example:
        ========
        >>> y = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[0] 
        array([1, 2, 4])
 
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        V1.1: 19-Apr-2010: now returns a SpaCo instance
 
        """
        import numpy as n
        
        arr = n.array(self.data)
        
        return SpaCo(arr[idx].tolist(), self.dtype, self.carsph, self.units, self.ticktock)   
        
    # -----------------------------------------------    
    def __setitem__(self, idx, vals):
        """
        a.__setitem__(idx, vals) or a[idx] = vals
        
        Will be called setting items in this instance 
        
        Input:
        ======
            - a SpaCo class instance
            - idx (int) : interger numbers as index
            - vals (numpy array or list) : new values
            
        Example:
        ========
        >>> y = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y[1] = [9,9,9]
        >>> y
        SpaCo( [[1 2 4]
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
            - a SpaCo class instance
            
        Returns:
        ========
            - length (int number)

        Example:
        ========
        >>> y = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
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
        
        Can be used to create a new SpaCo instance with new coordinate types
        
        Input:
        ======
            - a SpaCo class instance
            - returntype (string) : coordinate system, possible are GDZ, GEO, GSM, GSE, SM, GEI
                MAG, SPH, RLL
            - returncarsph (string) : coordinate type, possible 'car' for cartesian and 
                'sph' for spherical
 
        Returns:
        ========
            - a SpaCo object          

        Example:
        ========
        >>> y = SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
        >>> y.ticktock = TickTock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> x = y.convert('SM','car')
        >>> x
        SpaCo( [[ 0.81134097  2.6493305   3.6500375 ]
         [ 0.92060408  2.30678864  1.68262126]] ), dtype=SM,car, units=['Re', 'Re', 'Re']


        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        
        """
        
        import numpy as n
        import onerapy as op

        # no change necessary
        if (self.dtype is returntype) and (self.carsph is returncarsph):
            return self
            
        # only car2sph or sph2car is needed
        if (self.dtype is returntype) and (self.carsph is not returncarsph):
            if returncarsph == 'car':
                carsph = 'car'
                units = [self.units[0]]*3
                data = op.sph2car(self.data)
            else:
                carsph = 'sph'
                units =  [self.units[0], 'deg','deg']
                data = op.car2sph(self.data)
            return SpaCo(data, self.dtype, carsph, units, self.ticktock)

        # check the length of ticktock and do the more complex convertions
        if self.ticktock: 
            assert len(self.ticktock) == len(self), 'ticktock dimension does not match SpaCo dimensions'
        
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
        
        spaco = SpaCo(data, self.dtype, carsph, units, self.ticktock)
        
        # now convert to other coordinate system
        if (self.dtype is not returntype) : 
            assert spaco.ticktock, "Time information required; add a.ticktock attribute"
            spaco.data = op.coord_trans( spaco, returntype, returncarsph)
            spaco.dtype = returntype
            spaco.carsph = returncarsph
            spaco.sysaxes = op.get_sysaxes(returntype, returncarsph)

        # fix corresponding attributes         
        if returncarsph == 'sph':
            spaco.units = [units[0], 'deg','deg']
            if spaco.__dict__.has_key('x'): spaco.__delattr__('x')
            if spaco.__dict__.has_key('y'): spaco.__delattr__('y')
            if spaco.__dict__.has_key('z'): spaco.__delattr__('z')
            spaco.radi = spaco.data[:,0]
            spaco.lati = spaco.data[:,1]
            spaco.long = spaco.data[:,2]
        else: # 'car'
            spaco.units =  [units[0]]*3
            if spaco.__dict__.has_key('radi'): spaco.__delattr__('radi')
            if spaco.__dict__.has_key('lati'): spaco.__delattr__('lati')
            if spaco.__dict__.has_key('long'): spaco.__delattr__('long')
            spaco.x = spaco.data[:,0]
            spaco.y = spaco.data[:,1]
            spaco.z = spaco.data[:,2]
       
        return spaco

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
# TickDelta class
# -----------------------------------------------    
class TickDelta(object):
    """
    TickDelta( **kwargs )
    
    TickDelta class holding timedelta similar to datetime.timedelta
    This can be used to add/substract from TickTock objects
    
    Input:
    ======
        - days, hours, minutes and/or seconds (int or float) : time step
    
    Returns:
    ========
        - instance with self.days, self.secs, self.timedelta
    
    Example:
    ========  
    >>> dt = TickDelta(days=3.5, hours=12)
    >>> dt
    TickDelta( days=4.0 )
       
    See Also:
    =========
    spacetime.TickTock class
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 03-Mar-2010 (JK)
    """

    def __init__(self, **kwargs):
        import datetime as dt
        days, hours, minutes, seconds = [0,0,0,0]
        if kwargs.has_key('days'):  days = kwargs['days'] 
        if kwargs.has_key('hours'): hours = kwargs['hours']
        if kwargs.has_key('minutes'): minutes = kwargs['minutes']
        if kwargs.has_key('seconds'): seconds = kwargs['seconds']
        self.days = days + hours/24. + minutes/1440. + seconds/86400.
        self.hours = self.days*24.
        self.minutes = self.hours*60.
        self.seconds = self.minutes*60.
        self.timedelta = dt.timedelta(days=float(self.days))
        return  
        
    # -----------------------------------------------    
    def __str__(self):
        """
        dt.__str__() or dt
        
        Will be called when printing TickDelta instance dt
        
        Input:
        ======
            - a TickDelta class instance
 
        Returns:
        ========
            - output (string)          

        Example:
        ========
        >>> dt = TickDelta(3)
        >>> dt
        TickDelta( days=3 )
        
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
 
        """

        return 'TickDelta( days='+str(self.days) + ' )'  
    __repr__ = __str__
    
    # -----------------------------------------------    
    def __add__(self, other):
        """
        see TickTock.__add__
        
        """
        # call add routine from TickTock
        newobj = TickTock.__add__(other, self)
        return newobj

    # -----------------------------------------------    
    def __sub__(self, other):
        """
        see TickTock.__sub__
        
        """
        # call add routine from TickTock
        newobj = TickTock.__sub__(other, self)
        return newobj
     
    # -----------------------------------------------    
    def __mul__(self, other):
        """
        see TickTock.__sub__
        
        """
        # call add routine from TickTock
        newobj = self.timedelta.__mul__(other)
        days = newobj.days + newobj.seconds/86400.
        return TickDelta(days)
      
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
# TickTock class
# -----------------------------------------------    
class TickTock(object):
    """
    TickTock( data, dtype )
    
    TickTock class holding various time coordinate systems 
    (TAI, UTC, ISO, JD, MJD, UNX, RDT, CDF, DOY, eDOY)
    
    Possible data types:
    ISO: ISO standard format like '2002-02-25T12:20:30'
    UTC: datetime object with UTC time
    TAI: elapsed seconds since 1958/1/1 (includes leap seconds)
    UNX: elapsed seconds since 1970/1/1 (all days have 86400 secs sometimes unequal lenghts)
    JD: Julian days elapsed
    MJD: Modified Julian days
    RDT: Rata Die days elapsed since 1/1/1
    CDF: CDF epoch: milliseconds since 1/1/0000
    
    Input:
    ======
        - data (int, datetime, float, string) : time stamp
        - dtype (string) : data type for data, possible values: 
            CDF, ISO, UTC, TAI, UNX, JD, MJD, RDT
    
    Returns:
    ========
        - instance with self.data, self.dtype, self.UTC etc
    
    Example:
    ========
    
    >>> x=TickTock(2452331.0142361112, 'JD')
    >>> x.ISO
    '2002-02-25T12:20:30'
    >>> x.DOY # Day of year
    56
    
    See Also:
    =========
    a.getCDF, a.getISO, a.getUTC, etc.
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 20-Jan-2010
    V2: 25-Jan-2010: includes array support (JK)
    V3: 25-feb-2010: pulled functions into class (JK)
    V4: 19-May-2010: ISO format support (SM)
    """
    def __init__(self, data, dtype):
        from numpy import ndarray
        keylist = ['UTC','TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'GPS']
        assert dtype.upper() in keylist, "data type "+self.dtype+" not provided, only "+str(keylist)
        self.dtype = dtype.upper()
        self.__isoformatstr = {'seconds': '%Y-%m-%dT%H:%M:%S', 'microseconds': '%Y-%m-%dT%H:%M:%S.%f'}
        self.__isofmt = self.__isoformatstr['seconds']
        if isinstance(data, (list, ndarray)):
           self.data = data
        else:
            self.data = [data]
        if dtype.upper() == 'TAI': self.TAI = self.data
        if dtype.upper() == 'ISO':
            self.ISO = self.data
            self.update_items(self, 'data')
        if dtype.upper() == 'JD': self.JD = self.data
        if dtype.upper() == 'MJD': self.MJD = self.data
        if dtype.upper() == 'UNX': self.UNX = self.data
        if dtype.upper() == 'RDT': self.RDT = self.data
        if dtype.upper() == 'CDF': self.CDF = self.data
        self.UTC = self.getUTC()
        return
        
    # -----------------------------------------------    
    def __str__(self):
        """
        a.__str__() or a
        
        Will be called when printing TickTock instance a
        
        Input:
        ======
            - a TickTock class instance
 
        Returns:
        ========
            - output (string)          

        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:03', 'ISO')
        >>> a
        TickTock( ['2002-02-02T12:00:00'] ), dtype=ISO

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
 
        """

        return 'TickTock( '+str(self.data) + ' ), dtype='+str(self.dtype)  
    __repr__ = __str__
    
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
    def __getitem__(self, idx):
        """
        a.__getitem__(idx) or a[idx]
        
        Will be called when requesting items in this instance 
        
        Input:
        ======
            - a TickTock class instance
            - idx (int) : interger numbers as index

        Returns:
        ========
            - TickTock instance with requested values

        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:03', 'ISO')
        >>> a[0]
        '2002-02-02T12:00:00'
        
        See Also:
        =========
        a.__setitem__

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1.0: 03-Mar-2010 (JK)
        V1.1: 23-Mar-2010: now returns TickTock instance and can be indexed with arrays (JK)
 
        """
        import numpy as n
        
        arr = n.array(self.data)
        return TickTock(arr[idx].tolist(), self.dtype)   
    
    # -----------------------------------------------    
    def __setitem__(self, idx, vals):
        """
        a.__setitem__(idx, vals) or a[idx] = vals
        
        Will be called setting items in this instance 
        
        Input:
        ======
            - a TickTock class instance
            - idx (int) : integer numbers as index
            - vals (float, string, datetime) : new values
            

        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:03', 'ISO')
        >>> a[0] = '2003-03-03T00:00:00'
        
        See Also:
        =========
        a.__getitem__

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
 
        """
        self.data[idx] = vals
        self.update_items(self, 'data')
        
        return
        
    # -----------------------------------------------    
    def __len__(self):
        """
        a.__len__() or len(a)
        
        Will be called when requesting the length, i.e. number of items 
        
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - length (int number)

        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:03', 'ISO')
        >>> a.len
        1

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
 
        """
         
        from numpy import ndarray 
        if isinstance(self.data, (list, ndarray)):
            return len(self.data)
        else:
            return 1
        return  
    
    # -----------------------------------------------    
    def __cmp__(self, other):
        """
        a.__cmp__(other)
        
        Will be called when two TickTock instances are compared
        
        Input:
        ======
            - a TickTock class instance
            - other (TickTock instance) 
            
        Returns:
        ========
            - True or False

        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:03', 'ISO')
        >>> b = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a > b
        True
        
        See Also:
        =========
        a.__add__, a.__sub__

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
 
        """
        return cmp(self.UNX, other.UNX)  

    # -----------------------------------------------    
    def __gt__(self, other):
        return self.UNX > other.UNX
        
    def __lt__(self, other):
        return self.UNX < other.UNX
        
    def __ge__(self, other):
        return self.UNX >= other.UNX
        
    def __le__(self, other):
        return self.UNX <= other.UNX
        
    def __eq__(self, other):
        return self.UNX == other.UNX
        
    def __ne__(self, other):
        return self.UNX != other.UNX
    
    # -----------------------------------------------    
    def __sub__(self, other):
        """
        a.__sub__(other)
        
        Will be called if a TickDelta object is substracted to this instance and
        returns a new TickTock instance

        Input:
        ======
            - a TickTock class instance
            - other (TickDelta instance) 
      
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> dt = TickDelta(3)
        >>> a - dt
        TickTock( ['2002-02-05T12:00:00'] ), dtype=ISO


        See Also:
        =========
        __sub__
      
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
       
        """
        nTAI = len(self.data)
        newUTC = ['']*nTAI
        for i in range(nTAI):
            newUTC[i] = self.UTC[i] - other.timedelta
        newobj = TickTock(newUTC, 'UTC')
        newobj.data = eval('newobj.get'+self.dtype+'()')
        newobj.dtype = self.dtype
        newobj.update_items(self, 'data')
        return newobj

    
    # -----------------------------------------------    
    def __add__(self, other):
        """
        a.__add__(other)
        
        Will be called if a TickDelta object is added to this instance and
        returns a new TickTock instance

        Input:
        ======
            - a TickTock class instance
            - other (TickDelta instance) 
      
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> dt = TickDelta(3)
        >>> a + dt
        TickTock( ['2002-02-05T12:00:00'] ), dtype=ISO


        See Also:
        =========
        __sub__
      
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
       
        """

        nTAI = len(self.data)
        newUTC = ['']*nTAI
        for i in range(nTAI):
            newUTC[i] = self.UTC[i] + other.timedelta
        newobj = TickTock(newUTC, 'UTC')
        newobj.data = eval('newobj.get'+self.dtype+'()')
        newobj.dtype = self.dtype
        newobj.update_items(self, 'data')
        return newobj
       
    # -----------------------------------------------    
    def __getattr__(self, name):
        """
        a.__getattr__(name)
        
        Will be called if attribute "name" is not found in TickTock class instance.
        It will add TAI, RDT, etc

        Input:
        ======
            - a TickTock class instance
            - name (string) : a string from the list of time systems
                    'UTC', 'TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'DOY', 'eDOY', 'leaps'
        Returns:
        ========
            - requested values as either list/numpy array

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
       
        """

        keylist = ['TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'DOY', 'eDOY', 'leaps', 'GPS']
        assert name in keylist, "data type "+str(name)+" not provided, only "+str(keylist)
        if name.upper() == 'TAI': self.TAI = self.getTAI()
        if name.upper() == 'ISO': self.ISO = self.getISO()
        if name.upper() == 'JD': self.JD = self.getJD()
        if name.upper() == 'MJD': self.MJD = self.getMJD()
        if name.upper() == 'UNX': self.UNX = self.getUNX()
        if name.upper() == 'RDT': self.RDT = self.getRDT()
        if name.upper() == 'CDF': self.CDF = self.getCDF()
        if name.upper() == 'DOY': self.DOY = self.getDOY()
        if name.upper() == 'EDOY': self.eDOY = self.geteDOY()
        if name.upper() == 'GPS' : self.GPS = self.getGPS()
        #if name == 'isoformat': self.__isofmt = self.isoformat()
        if name == 'leaps': self.leaps = self.getleapsecs()
        return eval('self.'+name)
    
    # -----------------------------------------------
    def isoformat(self, fmt=None):
        """
        a.update_items(b, attrib)
        
        This changes the self.__isofmt attribute by and subsequently this 
        function will update the ISO attribute.

        Input:
        ======
            - a TickTock class instance
    
        Author:
        =======
        Steve Morley, Los Alamos National Lab (smorley@lanl.gov)
    
        Version:
        ========
        V1: 19-May-2010 (SM)
       
        """
        
        if not fmt:
            print 'Current ISO output format is %s' % self.__isofmt
            print 'Options are: %s' % [(k, self.__isoformatstr[k]) for k in self.__isoformatstr.keys()]
        else:
            try:
                self.__isofmt = self.__isoformatstr[fmt]
                self.update_items(self, 'data')
            except KeyError:
                print 'Not a valid option: Use %s' % self.__isoformatstr.keys()

        return
    
    # -----------------------------------------------
    def update_items(self, cls, attrib):
        """
        a.update_items(b, attrib)
        
        After changing the self.data attribute by either __setitem__ or __add__ etc
        this function will update all other attributes. This function is
        called automatically in __add__ and __setitem__

        Input:
        ======
            - a TickTock class instance
      
        See Also:
        =========
        __setitem__, __add__, __sub__
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 03-Mar-2010 (JK)
       
        """
        
        keylist = cls.__dict__.keys()
        keylist.remove('dtype')
        keylist.remove('data')
        if attrib is not 'data': keylist.remove(attrib)
        
        self.UTC = self.getUTC()
        for key in keylist:
            if key.upper() == 'TAI': self.TAI = self.getTAI()
            if key.upper() == 'ISO': self.ISO = self.getISO()
            if key.upper() == 'JD': self.JD = self.getJD()
            if key.upper() == 'MJD': self.MJD = self.getMJD()
            if key.upper() == 'UNX': self.UNX = self.getUNX()
            if key.upper() == 'RDT': self.RDT = self.getRDT()
            if key.upper() == 'CDF': self.CDF = self.getCDF()
            if key.upper() == 'DOY': self.DOY = self.getDOY()
            if key.upper() == 'eDOY': self.eDOY = self.geteDOY()
            if key.upper() == 'GPS' : self.GPS = self.getGPS()
            if key == 'leaps': self.leaps = self.getleapsecs()
            
        return

    # -----------------------------------------------    
    def convert(self, dtype):
        """
        a.convert(dtype)
        
        convert a TickTock instance into a new time coordinate system provided in dtype
        
        Input:
        ======
            - a TickTock class instance
            - dtype (string) : data type for new system, possible values are
                CDF, ISO, UTC, TAI, UNX, JD, MJD, RDT

        Returns:
        ========
            - newobj (TickTock instance) with new time coordinates
            
        Example:
        ========
        >>> a = TickTock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> s = a.convert('TAI')
        >>> type(s)
        <class 'spacetime.TickTock'>
        >>> s
        TickTock( [1391342432 1391342432] ), dtype=TAI

        See Also:
        =========
        a.CDF, a.ISO, a.UTC, etc. 

        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 05-Mar-2010 (JK)
        
        """
        newdat = eval('self.'+dtype)        
        return TickTock(newdat, dtype)

    # -----------------------------------------------    
    def append(self, other):
        """
        a.append(other)
        
        Will be called when another TickTock instance has to be appended to the current one

        Input:
        ======
            - a TickTock class instance
            - other (TickTock instance) 
      
        Example:
        ========
      
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 23-Mar-2010 (JK)
       
        """
        import numpy as n
        
        otherdata = eval('other.'+self.dtype)       
        newobj = TickTock(n.append(self.data, otherdata), dtype=self.dtype)
        return newobj
    
    # -----------------------------------------------
    def getCDF(self):
        """
        a.getCDF() or a.CDF
        
        Return CDF time which is milliseconds since 01-Jan-0000 00:00:00.000. 
        Year zero" is a convention chosen by NSSDC to measure epoch values. 
        This date is more commonly referred to as 1 BC. Remember that 1 BC was a leap year. 
        The CDF date/time calculations do not take into account the changes to the Gregorian 
        calendar, and cannot be directly converted into Julian date/times.

        Input:
        ======
            - a TickTock class instance
      
        Returns:
        ========
            - CDF (numpy array) : days elapsed since Jan. 1st
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.CDF
        array([  6.31798704e+13])
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getISO, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 02-Feb-2010 (JK)
        """
    
        RDTdata = self.getRDT()
        CDF = RDTdata*86400000.0 + 86400000.0*365.0
    
        self.CDF = CDF
        return CDF
    

    # -----------------------------------------------
    def getDOY(self):
        """
        a.DOY or a.getDOY()
        
        extract DOY (days since January 1st of given year)
        
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - DOY (numpy array int) : day of the year
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.DOY
        array([ 33])
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getISO, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 25-Jan-2010 (JK)
        V1.1: 20-Apr-2010: returns true DOY per definition as integer (JK)
        """ 
    
        import datetime
        import numpy as n
        
        nTAI = len(self.data)
        DOY = n.zeros(nTAI)
        
        for i in n.arange(nTAI):
            DOY[i] = self.UTC[i].toordinal() - datetime.date(self.UTC[i].year, 1, 1).toordinal() + 1
 
        self.DOY = DOY.astype(int)
        return DOY
        
    # -----------------------------------------------
    def geteDOY(self):
        """
        a.eDOY or a.geteDOY()
        
        extract eDOY (elapsed days since midnight January 1st of given year)
        
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - eDOY (numpy array) : days elapsed since midnight bbedJan. 1st
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.eDOY
        array([ 32.5])
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getISO, getDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Apr-2010 (JK)
        V2: 18-May-2010: Added microseconds (SM)
        """ 
    
        import datetime
        import numpy as n
        
        nTAI = len(self.data)
        eDOY = n.zeros(nTAI)
        
        for i in n.arange(nTAI):
            eDOY[i] = self.UTC[i].toordinal() - datetime.date(self.UTC[i].year, 1, 1).toordinal()
            eDOY[i] = eDOY[i] + self.UTC[i].hour/24. + self.UTC[i].minute/1440. + \
                self.UTC[i].second/86400. + self.UTC[i].microsecond/86400000000.

        self.eDOY = eDOY
        return eDOY


    # -----------------------------------------------
    def getJD(self):
        """
        a.JD or a.getJD()
        
        convert dtype data into Julian Date (JD) 
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - JD (numpy array) : elapsed days since 12:00 January 1, 4713 BC 
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.JD
        array([ 2452308.])

        See Also:
        =========
        getUTC, getUNX, getRDT, getISO, getMJD, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: added array support (JK)
        V3: 18-May-2010: added microseconds (SM)
        """ 
        import datetime, decimal
        import numpy as n
        
        nTAI = len(self.data)

        # convert all types in to UTC first and call again
        UTCdata = self.UTC

        if UTCdata[0] < datetime.datetime(1582,10,15):
            print "WARNING: Calendar date before the switch from Julian to Gregorian"
            print "    Calendar 1582-Oct-15: Use Julian Calendar dates as input"

        # include offset if given
        JD = n.zeros(nTAI)
        for i in n.arange(nTAI):
            offset = UTCdata[i].utcoffset()
            if offset:
                UTCdata[i] = UTCdata[i] - offset
        
            # extract year, month, day
            Y = int(UTCdata[i].year)
            M = int(UTCdata[i].month)
            D = int(UTCdata[i].day)

            # the following is from Wikipedia (but is wrong by 2 days)
            # JDN = D-32075+1461*(Y+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12-3*((Y+4900+(M-14)/12)/100)/4
            # JD = JDN + (data.hour-12)/24. + data.minute/1440. + data.second/86400.
        
            # following Press, "Numerical Recipes", Fct: JULDAY, p. 10
            igreg = 15+31*(10+12*1582)
            if M > 2:
                JY = Y
                JM = M+1
            else:
                JY = Y-1
                JM = M+13
            JD[i] = long(365.25*JY) + long(30.6001*JM) + D + 1720995
            c_val = (D+31*(M+12*Y))
            if c_val >= igreg: # yes if date after the Gregorian Switch in 1582-Oct-15
                JA = long(0.01*JY)
                JD[i] = JD[i]+2-JA+long(0.25*JA)

            # add this to num.recipes to get fractional days
            twelve, twofour, mind = decimal.Decimal('12.0'), decimal.Decimal('24.0'), decimal.Decimal('1440.0')
            sind, usind = decimal.Decimal('86400.0'), decimal.Decimal('86400000000.0')
            JD[i] = decimal.Decimal(str(JD[i])) + (decimal.Decimal(str(UTCdata[i].hour))-twelve)/twofour + \
                decimal.Decimal(str(UTCdata[i].minute/1440.)) + (decimal.Decimal(str(UTCdata[i].second))/sind) + \
                (decimal.Decimal(str(UTCdata[i].microsecond))/usind)
            JD[i] = float(JD[i])
            #JD[i] = JD[i] + (UTCdata[i].hour-12)/24. + UTCdata[i].minute/1440. + \
                #UTCdata[i].second/86400. + UTCdata[i].microsecond/86400000000.
    
        self.JD = JD
        return JD

    # -----------------------------------------------
    def getMJD(self):
        """
        a.MJD or a.getMJD()
        
        convert dtype data into MJD (modified Julian date)
        
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - MJD (numpy array) : elapsed days since November 17, 1858 
                    (Julian date was 2,400 000) 
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.MJD
        array([ 52307.5])

        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getISO, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: added support for arrays (JK)
        """

        import datetime

        if self.UTC[0] < datetime.datetime(1582,10,15):
            print "WARNING: Calendar date before the switch from Julian to Gregorian"
            print "Calendar 1582-Oct-15: Use Julian Calendar dates as input"
            
        MJD = self.JD - 2400000.5
    
        self.MJD = MJD
        return MJD

    # -----------------------------------------------
    def getUNX(self):
        """
        a.UNX or a.getUNX()
        
        convert dtype data into Unix Time (Posix Time)
        seconds since 1970-Jan-1 (not counting leap seconds)
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - UNX (numpy array) : elapsed secs since 1970/1/1 (not counting leap secs) 
    
        Example:
        ========
    
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.UNX
        array([  1.01265120e+09])
    
        See Also:
        =========
        getUTC, getISO, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: added array support (JK)
        V3: 18-May-2010: added sub-second support (SM)
        """
        
        import datetime
        import numpy as n
              
        nTAI = len(self.data)

        UNX0 = datetime.datetime(1970,1,1)
        d = ['']*nTAI
        UNX = n.zeros(nTAI)
        for i in n.arange(nTAI):
            d[i] = self.UTC[i] - UNX0 # timedelta object (only days, seconds, microsecs are stored)
            UNX[i] = (d[i].days)*86400 + d[i].seconds + d[i].microseconds/1.e6
    
        self.UNX = UNX
        return UNX
                
    # -----------------------------------------------
    def getRDT(self):
        """
        a.RDT or a.RDT()
        
        convert dtype data into Rata Die (lat.) Time (days since 1/1/0001)
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - RDT (numpy array) : elapsed days since 1/1/1 
    
        Example:
        ========
    
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.RDT
        array([ 730883.5])

        See Also:
        =========
        getUTC, getUNX, getISO, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: added array support (JK)
        V3: 17-May-2010: added microseconds (SM)
        """
        
        import datetime
        import numpy as n
        import matplotlib.dates as mpd
        
        nTAI = len(self.data)
        UTC = self.UTC
        #RDT = n.zeros(nTAI)
        RDT = mpd.date2num(UTC)
        #for i in n.arange(nTAI):
            #RDT[i] = UTC[i].toordinal() + UTC[i].hour/24. + UTC[i].minute/1440. + \
                #UTC[i].second/86400. + UTC[i].microsecond/86400000000.
    
        self.RDT = RDT
        return RDT
            
    # -----------------------------------------------
    def getUTC(self):
        """
        a.UTC or a.getUTC()
        
        convert dtype data into UTC object a la datetime()
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - UTC (list of datetime objects) : datetime object in UTC time
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.UTC
        [datetime.datetime(2002, 2, 2, 12, 0)]
    
        See Also:
        =========
        getISO, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: added array support (JK)
        V3: 17-May-2010: added microsecond ISO parsing (SM)
        """
    
        import datetime
        import numpy as n
    
        fmt,fmt2 = '%Y-%m-%dT%H:%M:%S','%Y-%m-%dT%H:%M:%S.%f'
        
        nTAI = len(self.data)
                       
        UTC = ['']*nTAI
        if self.dtype.upper() == 'UTC':
            UTC = self.data # return
            
        elif self.dtype.upper() == 'ISO':
            for i in n.arange(nTAI):
                if len(self.data[i])==19:
                    UTC[i] = datetime.datetime.strptime(self.data[i], fmt)
                else:
                    UTC[i] = datetime.datetime.strptime(self.data[i], fmt2)
            
        elif self.dtype.upper() == 'TAI':   
            for i in n.arange(nTAI):
                TAI0 = datetime.datetime(1958,1,1,0,0,0,0)
                for i in n.arange(nTAI):
                    UTC[i] = datetime.timedelta(seconds=float(self.data[i])) + TAI0                 
             # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in n.arange(nTAI):
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i]))

        elif self.dtype.upper() == 'GPS':   
            for i in n.arange(nTAI):
                GPS0 = datetime.datetime(1980,1,6,0,0,0,0)
                for i in xrange(nTAI):
                    UTC[i] = datetime.timedelta(seconds=float(self.data[i])) + GPS0                 
             # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in n.arange(nTAI):
                # there were 18 leap secinds before gps zero, need the -18 for that
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i])) + \
                    datetime.timedelta(seconds=19)           
 
        elif self.dtype.upper() == 'UNX':
            UNX0 = datetime.datetime(1970,1,1)     
            for i in n.arange(nTAI):
                UTC[i] = datetime.timedelta(seconds=self.data[i]) + UNX0 # timedelta object
    
        elif self.dtype.upper() == 'RDT':
            import matplotlib.dates as mpd
            UTC = mpd.num2date(self.data)
            UTC = [t.replace(tzinfo=None) for t in UTC]
            #for i in n.arange(nTAI):
                #UTC[i] = datetime.datetime(1,1,1) + \
                    #datetime.timedelta(days=n.floor(self.data[i])-1) +  \
                    #datetime.timedelta(microseconds=(self.data[i] - \
                        #self.data[i])*86400000.)
                # roundoff the microseconds
                #UTC[i] = UTC[i] - datetime.timedelta(microseconds=UTC[i].microsecond)
        
        elif self.dtype.upper() == 'CDF':
            for i in n.arange(nTAI):
                UTC[i] = datetime.timedelta(days=self.data[i]/86400000.) + \
                        datetime.datetime(1,1,1) - datetime.timedelta(days=366)
                #UTC[i] = datetime.timedelta(days=n.floor(self.data[i]/86400000.), \
                    #milliseconds=n.mod(self.data[i],86400000)) + \
                        #datetime.datetime(1,1,1) - datetime.timedelta(days=366)
                # the following has round off errors
                # UTC[i] = datetime.timedelta(data[i]/86400000.-366) + datetime.datetime(1,1,1)
        
        elif self.dtype.upper() in ['JD', 'MJD']:
            if self.dtype.upper() == 'MJD': 
                self.JD = n.array(self.data) + 2400000.5
            for i in n.arange(nTAI):
                # extract partial days
                ja = int(n.floor(self.JD[i]))
                p = self.JD[i] - n.floor(self.JD[i])
                # after Press: "Numerical Recipes"
                # http://www.rgagnon.com/javadetails/java-0506.html
                # only good for after 15-Oct-1582
                igreg = 15+31*(10+12*1582)
                if ja >= igreg: # after switching to Gregorian Calendar
                    jalpha = int(((ja-1867216)-0.25)/36524.25)
                    ja = ja + 1 + jalpha -jalpha/4
    
                jb = ja + 1524
                jc = int(6680.0 + ((jb - 2439870) - 122.1) / 365.25)
                jd = 365 * jc + jc/4
                je = int((jb - jd)/30.6001)
                day = jb - jd - int(30.6001 * je)
                month = je - 1
                if (month > 12): month = month - 12
                year = jc - 4715
                if (month > 2): year = year-1
                if (year <= 0): year = year-1
    
                # after http://aa.usno.navy.mil/faq/docs/JD_Formula.php
                # also good only for after 1582-Oct-15
                # L= JD+68569
                # N= 4*L/146097
                # = L-(146097*N+3)/4
                # I= 4000*(L+1)/1461001
                # L= L-1461*I/4+31
                # J= 80*L/2447
                # K= L-2447*J/80
                # L= J/11
                # J= J+2-12*L
                # I= 100*(N-49)+I+L   
            
                UTC[i] = datetime.datetime(year, month, day) + datetime.timedelta(hours=12) + \
                    datetime.timedelta(seconds = p*86400)
                if UTC[i] < datetime.datetime(1582,10,15):
                    print "WARNING: Calendar date before the switch from Julian to Gregorian"
                    print "Calendar 1582-Oct-15: Use Julian Calendar dates as input"
    
        else:
            print "ERROR: Data type ", self.dtype, ' in getUTC() not supported'
            return
    
        self.UTC = UTC
        return UTC

    # -----------------------------------------------
    def getGPS(self):
        """
        a.GPS or a.getGPS()
        
        return GPS epoch (0000 UT (midnight) on January 6, 1980)
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - GPS (numpy array) : elapsed secs since 6Jan1980 (excludes leap secs)
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.GPS
        array([])
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY
    
        Author:
        =======
        Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (BAL)
        V2: 17-May-2010: Added sub-second support (SM)
        """
    
        import datetime, time
        import numpy as np
    
        fmt = '%Y-%m-%dT%H:%M:%S'
        GPS0 = datetime.datetime(1980,1,6,0,0,0,0)
        
        nGPS = len(self.data)
        GPS = np.zeros(nGPS)
        UTC = self.UTC
        leapsec = self.getleapsecs()
        GPStup = ['']*nGPS
        for i in np.arange(nGPS):
            # get the leap seconds
            GPStup[i] = UTC[i] - GPS0 + datetime.timedelta(seconds=int(leapsec[i])) - datetime.timedelta(seconds=19)
            GPS[i] = GPStup[i].days*86400 + GPStup[i].seconds + GPStup[i].microseconds/1.e6

        self.GPS = np.array(GPS)#.astype(int)
        return self.GPS



    # -----------------------------------------------
    def getTAI(self):
        """
        a.TAI or a.getTAI()
        
        return TAI (International Atomic Time)
    
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - TAI (numpy array) : elapsed secs since 1958/1/1 (includes leap secs, 
                    i.e. all secs have equal lengths)
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.TAI
        array([1391342432])
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: include array support (JK)
        """
    
        import datetime, time
        import numpy as n
    
        fmt = '%Y-%m-%dT%H:%M:%S'
        TAI0 = datetime.datetime(1958,1,1,0,0,0,0)
        
        nTAI = len(self.data)
        TAI = n.zeros(nTAI)
        UTC = self.UTC
        leapsec = self.getleapsecs()
        TAItup = ['']*nTAI
        for i in n.arange(nTAI):
            #t = time.strptime(data[i], fmt)
            #dtimetup = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5])    
            # get the leap seconds
            TAItup[i] = UTC[i] - TAI0 + datetime.timedelta(seconds=int(leapsec[i]))
            TAI[i] = TAItup[i].days*86400 + TAItup[i].seconds + TAItup[i].microseconds/1.e6

        self.TAI = n.array(TAI)
        return self.TAI
    
    # -----------------------------------------------
    def getISO(self):
        """
        a.ISO or a.getISO()
        
        convert dtype data into ISO string
        
        Input:
        ======
            - a TickTock class instance
            
        Returns:
        ========
            - ISO (list of strings) : date in ISO format
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.ISO
        ['2002-02-02T12:00:00']
    
        See Also:
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010 (JK)
        V2: 25-Jan-2010: included arary support (JK)
        V3: 10-May-2010: speedup, arange to xrange (BAL)
        V4: 17-May-2010: switched to native formatting so sub-second is displayed
        """
    
        import datetime, time
    
        nTAI = len(self.data)
        ISO = ['']*nTAI
        
        for i in xrange(nTAI):
            ISO[i] = self.UTC[i].strftime(self.__isofmt)
    
        self.ISO = ISO
        return ISO

    # -----------------------------------------------
    def getleapsecs(self):
        """
        a.leaps or a.getleapsecs()
        
        retrieve leapseconds from lookup table, used in getTAI
    
        Input:
        ======
            - a TickTock class instance
    
        Returns:
        ========
            - secs (numpy array) : leap seconds
    
        Example:
        ========
        >>> a = TickTock('2002-02-02T12:00:00', 'ISO')
        >>> a.leaps
        array([32])
    
        See Also:
        =========
        getTAI
    
        Author:
        =======
        Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
        Version:
        ========
        V1: 20-Jan-2010: includes array support (JK)
        """
    
        from spacepy import __path__
        import numpy as n
        import datetime
    
        tup = self.UTC
        
        # so you don't have to read the file every single time
        global secs, year, mon, day
        
        try:
           leaps = secs[0]
    
        except:  # then we are calling this routine the 1st time
           # load current file
           fname = __path__[0]+'/data/tai-utc.dat'
           fh = open(fname)
           text = fh.readlines()
    
           secs = n.zeros(len(text))
           year = n.zeros(len(text))
           mon = n.zeros(len(text))
           day = n.zeros(len(text))
    
           months = n.array(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', \
                  'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])
        
           for line, i in zip(text, n.arange(len(secs))):
              secs[i] = int(float(line.split()[6]))  # truncate float seconds
              year[i] = int(line.split()[0])
              mon[i] = int(n.where(months == line.split()[1])[0][0] + 1)
              day[i] = int(line.split()[2])
    
        # check if array:
        if type(tup) == type(datetime.datetime(1,1,1)): # not an array of objects
            tup = [tup]
            nTAI = 1
            aflag = False
        else:
            nTAI = len(tup)
            aflag = True
    
        # convert them into a time tuple and find the correct leap seconds
        leaps = [secs[0]]*nTAI
        for i, itup in enumerate(tup):
            for y,m,d,s in zip(year, mon, day, secs):
                if tup[i] > datetime.datetime(int(y),int(m),int(d)):
                    leaps[i] = s
                else:
                    break
              
        #if datetime.datetime(1971,12,31) > tup[0]:
        #   print "WARNING: date before 1972/1/1; leap seconds are by fractions off"
           
        if aflag == False:
            self.leaps = int(leaps[0])
            return int(leaps[0])   # if you want to allow fractional leap seconds, remove 'int' here
        else:
            self.leaps = n.array(leaps, dtype=int)
            return self.leaps
        
# -----------------------------------------------
# End of TickTock class
# -----------------------------------------------


def doy2date(year,doy, dtobj=False):
    """
    convert day-of-year doy into a month and day
    after http://pleac.sourceforge.net/pleac_python/datesandtimes.html

    Input:
    ======
        - year (int or array of int) : year
        - doy (int or array of int) : day of year

    Returns:
    ========
        - month (int or array of int) : month as integer number
        - day (int or array of int) : as integer number

    Example:
    ========
    >>> month, day = doy2date(2002, 186)
    >>> dts = doy2date([2002]*4, range(186,190), dtobj=True)

    See Also:
    =========
    getDOY

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov 

    Version:
    ========
    V1: 24-Jan-2010: can handle arrays as input (JK)
    V2: 02-Apr-2010: option to return date objects (SM)
    V3: 07-Apr-2010: modified to return datetime objects (SM)
    """

    import datetime
    import numpy as n

    # check if array:
    try:
        nTAI = len(year)
        year = n.array(year, dtype=int)
    except:
        year = n.array([year], dtype=int)
        doy = n.array([doy], dtype=int)
        nTAI = 1 # flag is the values need to be returned as singles
        
    month = n.zeros(nTAI, dtype=int)
    day = n.zeros(nTAI, dtype=int)
    dateobj = ['']*nTAI
    for i, iyear, idoy in zip( n.arange(nTAI), year, doy): 
        dateobj[i] = datetime.datetime(int(year[i]),1,1) + datetime.timedelta(days=int(doy[i]-1))    
        month[i] = dateobj[i].month
        day[i] = dateobj[i].day
    
    if nTAI == 1:
        if dtobj:
            return dateobj[0]
        else:
            return month[0], day[0]
    else:
        if dtobj:
            return dateobj
        else:
            return month, day

# -----------------------------------------------
def tickrange(start, end, deltadays, dtype='ISO'):
    """
    return a TickTock range given the start, end, and delta

    Input:
    ======
        - start (string or number) : start time
        - end (string or number) : end time (inclusive)
        - deltadays: step in units of days (float); or datetime timedelta object
        - (optional) dtype (string) : data type for start, end; e.g. ISO, UTC, RTD, etc.
            see TickTock for all options

    Returns:
    ========
        - ticks (TickTock instance)
        
    Example:
    ========
    >>> ticks = st.tickrange('2002-02-01T00:00:00', '2002-02-10T00:00:00', deltadays = 1)
    >>> ticks
    TickTock( ['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00', 
    '2002-02-04T00:00:00'] ), dtype=ISO

    See Also:
    =========
    TickTock

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 10-Mar-2010: (JK)
    V1.1: 16-Mar-2010: fixed bug with floating point precision (JK)
    V1.2: 28-Apr-2010: added timedelta support for increment (SM)
    """
    
    import datetime as dt

    Tstart = TickTock(start, dtype)
    Tend = TickTock(end, dtype)
    diff = Tend.UTC[0] - Tstart.UTC[0]
    dmusec, dsec = diff.microseconds/86400000000., diff.seconds/86400.
    try:
        assert type(deltadays)==dt.timedelta
        musec, sec = deltadays.microseconds/86400000000., deltadays.seconds/86400.
        deltat = musec + sec + deltadays.days
        nticks = int((dmusec + dsec + diff.days)/deltat + 1)
        trange = [Tstart.UTC[0] + deltadays*n for n in range(nticks)]
    except:
        nticks = int((dmusec + dsec + diff.days)/float(deltadays) + 1)
        trange = [Tstart.UTC[0] + dt.timedelta(days=deltadays)*n for n in range(nticks)]
    ticks = TickTock(trange, 'UTC')
    ticks = eval('TickTock(ticks.'+dtype+',"'+dtype+'")')
    return ticks
    
# -----------------------------------------------
def test():
    """
    test all time conversions

    Returns:
    ========
        - nFAIL (int) : number of failures

    Example:
    ========

    >>> test()
    testing ticktock: PASSED TEST 1
    testing ticktock: PASSED TEST 2
    0

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)

    Version:
    ========
    V1: 20-Jan-2010
    """
    import spacetime as st
    import toolbox as tb
    import numpy as n
    import sys, datetime
    
    def timecomp(x,yarr):
        x=x[0]
        delta = datetime.timedelta(microseconds=1000000)
        truth = [(y[0] <= x + delta) & (y[0] >= x - delta) for y in yarr]
        return truth
        
    # number of failures
    nFAIL = 0

    # time types
    alldtypes = ['UNX', 'TAI', 'JD', 'MJD', 'RDT', 'UTC', 'ISO', 'CDF']

    # pick a start time
    data = '2002-02-25T12:20:30.1'
    dtype = 'ISO'
    t = st.TickTock(data,dtype)
    
    out = ['']*(len(alldtypes))
    back = ['']*(len(alldtypes))
    # cycle (keep the sequence the same as in alldtypes!!)
    out[0] = t.getUNX()[0]
    out[1] = t.getTAI()[0]
    out[2] = t.getJD()[0]
    out[3] = t.getMJD()[0]
    out[4] = t.getRDT()[0]
    out[5] = t.getUTC()[0]
    out[6] = t.getISO()[0]
    out[7] = t.getCDF()[0]
    
    # compare all to UTC (datetime format)
    for i, datatype in enumerate(alldtypes):
        back[i] = t.getUTC()[0]
    comp = [st.TickTock(data,'ISO').getUTC()[0]]*len(alldtypes)
    #
    if back == comp and tb.feq(round(out[2],5),2452331.01424) and \
    tb.feq(out[0],1014639630.1):  
        print "testing TickTock: PASSED TEST simple"
    else:
        print "testing TickTock: FAILED TEST simple"
        nFAIL =+ 1

    # now test the class
    foo = ['']*(len(alldtypes))
    bar = ['']*(len(alldtypes))
    UNX = ['']*(len(alldtypes))    
    TAI = ['']*(len(alldtypes))    
    JD  = ['']*(len(alldtypes))    
    MJD = ['']*(len(alldtypes))    
    RDT = ['']*(len(alldtypes))    
    UTC = ['']*(len(alldtypes))    
    ISO = ['']*(len(alldtypes))   
    CDF = ['']*(len(alldtypes))
    for i, dtype in enumerate(alldtypes):
        foo[i] = st.TickTock(out[i], dtype)
        foo[i].isoformat('microseconds')

    # now compare and make sure they are all the same
    for i, dtype in enumerate(alldtypes):
        UNX[i] = foo[i].getUNX()
        TAI[i] = foo[i].getTAI()
        JD[i] = foo[i].getJD()
        MJD[i] = foo[i].getMJD()
        RDT[i]  = foo[i].getRDT()
        UTC[i] = foo[i].getUTC()
        ISO[i] = foo[i].getISO()
        CDF[i] = foo[i].getCDF()
    
        
    prec = 5./86400000000.
    testTAI = n.where(tb.feq(TAI[0],TAI, precision=prec))
    testUNX = n.where(tb.feq(UNX[0],UNX, precision=prec))
    testJD = n.where(tb.feq(JD[0],JD, precision=prec))
    testMJD = n.where(tb.feq(MJD[0],MJD, precision=prec))
    testRDT = n.where(tb.feq(RDT[0],RDT, precision=prec))
    testUTC = timecomp(UTC[0],UTC)
    testCDF = n.where(tb.feq(CDF[0],CDF, precision=prec))
    try:
        assert len(testUNX[0]) == len(alldtypes)
        assert len(testTAI[0]) == len(alldtypes)
        assert len(testJD[0]) == len(alldtypes)
        assert len(testMJD[0]) == len(alldtypes)
        assert len(testRDT[0]) == len(alldtypes)
        assert False not in testUTC
        assert len(testCDF[0]) == len(alldtypes)
        assert len(n.unique1d(ISO)) == 1

        print "testing TickTock: PASSED TEST all combinations"
    except AssertionError:
        print "testing TickTock: FAILED TEST all combinations"
        nFAIL =+ 1

    # test with arrays 
    import numpy as n
    try:
        for i, dtype in enumerate(alldtypes):
            foo = st.TickTock( n.array([out[i], out[i], out[i]]), dtype)
        print "testing TickTock: PASSED TEST arrays"
    except:
        print "testing TickTock: FAILED TEST arrays"
        nFAIL =+ 1
        
    # test DOY
    try:
        foo.DOY
        print "testing TickTock: PASSED TEST DOY"
    except:
        print "testing TickTock: FAILED TEST DOY"
        nFAIL =+ 1
        
    return nFAIL
    
    
