#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Time conversion, manipulation and implementation of Ticktock class


Examples:
=========

>>> import spacepy.time as spt
>>> import datetime as dt

Day of year calculations

>>> dts = spt.doy2date([2002]*4, range(186,190), dtobj=True)
>>> dts
[datetime.datetime(2002, 7, 5, 0, 0),
datetime.datetime(2002, 7, 6, 0, 0),
datetime.datetime(2002, 7, 7, 0, 0),
datetime.datetime(2002, 7, 8, 0, 0)]

>>> dts = spt.Ticktock(dts,'UTC')
>>> dts.DOY
array([ 186.,  187.,  188.,  189.])

Ticktock object creation

>>> isodates = ['2009-12-01T12:00:00', '2009-12-04T00:00:00', '2009-12-06T12:00:00']
>>> dts = spt.Ticktock(isodates, 'ISO')

OR

>>> dtdates = [dt.datetime(2009,12,1,12), dt.datetime(2009,12,4), dt.datetime(2009,12,6,12)]
>>> dts = spt.Ticktock(dtdates, 'UTC')

ISO time formatting

>>> dts = spt.tickrange('2009-12-01T12:00:00','2009-12-06T12:00:00',2.5)

OR

>>> dts = spt.tickrange(dt.datetime(2009,12,1,12),dt.datetime(2009,12,6,12), \
    dt.timedelta(days=2, hours=12))

>>> dts
Ticktock( ['2009-12-01T12:00:00', '2009-12-04T00:00:00', '2009-12-06T12:00:00'] ), dtype=ISO

>>> dts.isoformat()
Current ISO output format is %Y-%m-%dT%H:%M:%S
Options are: [('seconds', '%Y-%m-%dT%H:%M:%S'), ('microseconds', '%Y-%m-%dT%H:%M:%S.%f')]

>>> dts.isoformat('microseconds')
>>> dts.ISO
['2009-12-01T12:00:00.000000',
 '2009-12-04T00:00:00.000000',
 '2009-12-06T12:00:00.000000']

Time manipulation

>>> tdelt  = spt.Tickdelta(days=1, hours=6)
>>> tdelt
Tickdelta( days=1.25 )

>>> new_dts = dts + tdelt
>>> new_dts.UTC
[datetime.datetime(2009, 12, 2, 18, 0),
 datetime.datetime(2009, 12, 5, 6, 0),
 datetime.datetime(2009, 12, 7, 18, 0)]

Other time formats

>>> dts.RDT  # Gregorian ordinal time
array([ 733742.5,  733745. ,  733747.5])

>>> dts.GPS # GPS time
array([  9.43704015e+08,   9.43920015e+08,   9.44136015e+08])

>>> dts.JD # Julian day
array([ 2455167. ,  2455169.5,  2455172. ])

And so on.

Authors:
--------
Josef Koller, Brian Larsen, Steve Morley, Jon Niehof

jkoller@lanl.gov,
Los Alamos National Laboratory

Copyright ©2010 Los Alamos National Security, LLC.
"""

from spacepy import help
import datetime
import datetime as dt
import numpy as np
__version__ = "$Revision: 1.37 $, $Date: 2011/05/31 20:25:32 $"


# -----------------------------------------------
# Tickdelta class
# -----------------------------------------------
class Tickdelta(object):
    """
    Tickdelta( **kwargs )

    Tickdelta class holding timedelta similar to datetime.timedelta
    This can be used to add/substract from Ticktock objects

    Parameters
    ==========
    days : float
        number of days in for the delta
    hours : float
        number of hours for the delta
    minutes : float
        number of minutes for the delta
    seconds : float
        number of secondes for the delta

    Returns
    =======
    out : Tickdelta
        instance with self.days, self.secs, self.timedelta

    Examples
    ========
    >>> dt = Tickdelta(days=3.5, hours=12)
    >>> dt
    Tickdelta( days=4.0 )

    See Also
    ========
    spacepy.time.Ticktock class
    """
    def __init__(self, **kwargs):
        days, hours, minutes, seconds = [0,0,0,0]
        if 'days' in kwargs:  days = kwargs['days']
        if 'hours' in kwargs: hours = kwargs['hours']
        if 'minutes' in kwargs: minutes = kwargs['minutes']
        if 'seconds' in kwargs: seconds = kwargs['seconds']
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

        Will be called when printing Tickdelta instance dt

        Returns
        =======
        out : string
            string representation of the time

        Examples
        ========
        >>> dt = Tickdelta(3)
        >>> dt
        Tickdelta( days=3 )

        """
        return 'Tickdelta( days='+str(self.days) + ' )'
    __repr__ = __str__

    # -----------------------------------------------
    def __add__(self, other):
        """
        See Also
        ========
        Ticktock.__add__

        """
        # call add routine from Ticktock
        newobj = Ticktock.__add__(other, self)
        return newobj

    # -----------------------------------------------
    def __sub__(self, other):
        """
        See Also
        ========
        Ticktock.__sub__

        """
        # call add routine from Ticktock
        newobj = Ticktock.__sub__(other, self)
        return newobj

    # -----------------------------------------------
    def __mul__(self, other):
        """
        See Also
        ========
        Ticktock.__sub__

        """
        # call add routine from Ticktock
        newobj = self.timedelta.__mul__(other)
        days = newobj.days + newobj.seconds/86400.
        return Tickdelta(days)

    # -----------------------------------------------
    def __getstate__(self):
        """
        Is called when pickling

        See Also
        ========
        http://docs.python.org/library/pickle.html
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        return odict

    def __setstate__(self, dict):
        """
        Is called when unpickling

        See Also
        ========
        http://docs.python.org/library/pickle.html
        """
        self.__dict__.update(dict)
        return


# -----------------------------------------------
# Ticktock class
# -----------------------------------------------
class Ticktock(object):
    """
    Ticktock( data, dtype )

    Ticktock class holding various time coordinate systems
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

    Parameters
    ==========
    data : array_like (int, datetime, float, string)
        time stamp
    dtype : string {`CDF`, `ISO`, `UTC`, `TAI`, `UNX`, `JD`, `MJD`, `RDT`}
        data type for data

    Returns
    =======
    out : Ticktock
        instance with self.data, self.dtype, self.UTC etc

    Examples
    ========
    >>> x=Ticktock([2452331.0142361112, 2452332.0142361112], 'JD')
    >>> x.ISO
    ['2002-02-25T12:20:30', '2002-02-26T12:20:30']
    >>> x.DOY # Day of year
    array([ 56.,  57.])

    See Also
    ========
    a.getCDF
    a.getISO
    a.getUTC
    """
    def __init__(self, data, dtype=None):
        from numpy import ndarray
        from datetime import datetime
        if isinstance(data, ndarray):
            self.data = data
        elif isinstance(data, (bytes, str)):
            self.data = [data]
        else:
            try:
                self.data = list(data)
            except TypeError:
                self.data = [data]
        # make some educated guess on the format if dtype not provided
        if isinstance(self.data[0], str):
            dtype = 'ISO'
        elif isinstance(self.data[0], datetime):
            dtype = 'UTC'
        elif self.data[0] > 1e13:
            dtype = 'CDF'
        keylist = ['UTC','TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'GPS']
        assert dtype.upper() in keylist, "data type "+self.dtype+" not provided, only "+str(keylist)
        self.dtype = dtype.upper()
        self.__isoformatstr = {'seconds': '%Y-%m-%dT%H:%M:%S', 'microseconds': '%Y-%m-%dT%H:%M:%S.%f'}
        self.__isofmt = self.__isoformatstr['seconds']

        if dtype.upper() == 'ISO':
            if self.data[0].find('Z'):
                for i in range(len(self.data)):
                    self.data[i] = self.data[i].split('Z')[0]
            if self.data[0].find('T') == -1: # then assume midnight
                for i in range(len(self.data)):
                    self.data[i] = self.data[i]+'T00:00:00'
            self.ISO = self.data
            self.update_items(self, 'data')
        if dtype.upper() == 'TAI': self.TAI = self.data
        if dtype.upper() == 'JD': self.JD = self.data
        if dtype.upper() == 'MJD': self.MJD = self.data
        if dtype.upper() == 'UNX': self.UNX = self.data
        if dtype.upper() == 'RDT': self.RDT = self.data
        if dtype.upper() == 'CDF': self.CDF = self.data
        self.UTC = self.getUTC()
        self.shape = (len(self.data), )
        return

    # -----------------------------------------------
    def __str__(self):
        """
        a.__str__() or a

        Will be called when printing Ticktock instance a

        Returns
        =======
        out : string
            string representaion of the class

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:03', 'ISO')
        >>> a
        Ticktock( ['2002-02-02T12:00:00'] ), dtype=ISO
        """
        return 'Ticktock( '+str(self.data) + ' ), dtype='+str(self.dtype)
    __repr__ = __str__

    # -----------------------------------------------
    def __getstate__(self):
        """
        Is called when pickling
        See Also http://docs.python.org/library/pickle.html
        """
        odict = self.__dict__.copy() # copy the dict since we change it
        return odict

    def __setstate__(self, dict):
        """
        Is called when unpickling
        See Also http://docs.python.org/library/pickle.html
        """
        self.__dict__.update(dict)
        return

    # -----------------------------------------------
    def __getitem__(self, idx):
        """
        a.__getitem__(idx) or a[idx]

        Will be called when requesting items in this instance

        Parameters
        ==========
        idx : int
            the item index to get

        Returns
        =======
        out : Ticktock
            Ticktock instance with requested values

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:03', 'ISO')
        >>> a[0]
        '2002-02-02T12:00:00'

        See Also
        ========
        a.__setitem__

        """
        arr = np.array(self.data)

        if isinstance(idx, int):
            return Ticktock(arr[idx], self.dtype)
        else:
            return Ticktock(arr[idx].tolist(), self.dtype)

    # -----------------------------------------------
    def __setitem__(self, idx, vals):
        """
        a.__setitem__(idx, vals) or a[idx] = vals

        Will be called setting items in this instance

        Parameters
        ==========
        idx : int
            integer numbers as index
        vals: {float, string, datetime}
            new values

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:03', 'ISO')
        >>> a[0] = '2003-03-03T00:00:00'

        See Also
        ========
        a.__getitem__
        """
        self.data[idx] = vals
        self.update_items(self, 'data')

        return

    # -----------------------------------------------
    def __len__(self):
        """
        a.__len__() or len(a)

        Will be called when requesting the length, i.e. number of items

        Returns
        =======
        out : int
            length

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:03', 'ISO')
        >>> a.len
        1
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

        Will be called when two Ticktock instances are compared

        Paramters
        =========
        other : Ticktock
            instance for comparison

        Returns
        ========
        out : boolean
            True or False

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:03', 'ISO')
        >>> b = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a > b
        True

        See Also
        ========
        a.__add__, a.__sub__
        """
        return cmp(self.UNX, other.UNX)

    # -----------------------------------------------
    def __gt__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX > other.UNX

    def __lt__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX < other.UNX

    def __ge__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX >= other.UNX

    def __le__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX <= other.UNX

    def __eq__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX == other.UNX

    def __ne__(self, other):
        if isinstance(other, datetime.datetime):
            other = Ticktock(other, 'UTC')
        return self.UNX != other.UNX

    # -----------------------------------------------
    def __sub__(self, other):
        """
        a.__sub__(other)

        Will be called if a Tickdelta object is substracted to this instance and
        returns a new Ticktock instance

        Paramters
        =========
        other : Ticktock or Tickdelta
            instance for comparison

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> dt = Tickdelta(3)
        >>> a - dt
        Ticktock( ['2002-02-05T12:00:00'] ), dtype=ISO

        See Also
        ========
        __sub__
        """
        nTAI = len(self.data)
        if isinstance(other, Tickdelta):
            newUTC = ['']*nTAI
            for i in range(nTAI):
                newUTC[i] = self.UTC[i] - other.timedelta
            newobj = Ticktock(newUTC, 'UTC')
            newobj.data = eval('newobj.get'+self.dtype+'()')
            newobj.dtype = self.dtype
            newobj.update_items(self, 'data')
            return newobj

        elif isinstance(other, Ticktock):
            newTAI = ['']*nTAI
            for i in range(nTAI):
                newTAI[i] = self.TAI[i] - other.TAI
            deltas = [Tickdelta(seconds=val) for val in newTAI ]
            return deltas

    # -----------------------------------------------
    def __add__(self, other):
        """
        a.__add__(other)

        Will be called if a Tickdelta object is substracted to this instance and
        returns a new Ticktock instance

        Paramters
        =========
        other : Ticktock or Tickdelta
            instance for comparison

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> dt = Tickdelta(3)
        >>> a + dt
        Ticktock( ['2002-02-05T12:00:00'] ), dtype=ISO


        See Also
        ========
        __sub__

        """
        nTAI = len(self.data)
        if isinstance(other, Tickdelta):
            newUTC = ['']*nTAI
            for i in range(nTAI):
                newUTC[i] = self.UTC[i] + other.timedelta
            newobj = Ticktock(newUTC, 'UTC')
            newobj.data = eval('newobj.get'+self.dtype+'()')
            newobj.dtype = self.dtype
            newobj.update_items(self, 'data')
            return newobj

        elif isinstance(other, Ticktock):
            newTAI = ['']*nTAI
            for i in range(nTAI):
                newTAI[i] = self.TAI[i] + other.TAI
            deltas = [Tickdelta(seconds=val) for val in newTAI ]
            return deltas

    # -----------------------------------------------
    def __getattr__(self, name):
        """
        a.__getattr__(name)

        Will be called if attribute "name" is not found in Ticktock class instance.
        It will add TAI, RDT, etc

        Parameters
        ==========
            name : string
                a string from the list of time systems
                    'UTC', 'TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'DOY', 'eDOY', 'leaps'
        Returns
        ========
            out: list, array
                requested values as either list/numpy array


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
    def sort(self):
        """
        a.sort()

        This will sort the Ticktock values in place
        """
        RDT = self.RDT
        RDTsorted = np.sort(RDT)
        tmp = Ticktock(RDTsorted, 'RDT').convert(self.dtype)
        self.data = tmp.data
        self.update_items(self, 'data')
        return

    # -----------------------------------------------
    def argsort(self):
        """
        idx = a.argsort()

        This will return the indices that would sort the Ticktock values

        Returns
        =======
        out : list
            indices that would sort the Ticktock values

        """
        RDT = self.RDT
        idx = np.argsort(RDT)
        return idx

    # -----------------------------------------------
    def isoformat(self, fmt=None):
        """
        a.update_items(b, attrib)

        This changes the self.__isofmt attribute by and subsequently this
        function will update the ISO attribute.

        Parameters
        ==========
        fmt : string, optional
        """
        if not fmt:
            print('Current ISO output format is %s' % self.__isofmt)
            print('Options are: %s' % [(k, self.__isoformatstr[k]) for k in list(self.__isoformatstr.keys())])
        else:
            try:
                self.__isofmt = self.__isoformatstr[fmt]
                self.update_items(self, 'data')
            except KeyError:
                print('Not a valid option: Use %s' % list(self.__isoformatstr.keys()))

        return

    # -----------------------------------------------
    def update_items(self, cls, attrib):
        """
        a.update_items(b, attrib)

        After changing the self.data attribute by either __setitem__ or __add__ etc
        this function will update all other attributes. This function is
        called automatically in __add__ and __setitem__

        Parameters
        ==========
        cls : Ticktock
        attrib : string
            attribute to update

        See Also
        ========
        __setitem__
        __add__
        __sub__
        """
        keylist = list(cls.__dict__.keys())
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

        convert a Ticktock instance into a new time coordinate system provided in dtype

        Parameters
        ==========
        dtype : string
            data type for new system, possible values are {`CDF`, `ISO`, `UTC`, `TAI`, `UNX`, `JD`, `MJD`, `RDT`}

        Returns
        =======
        out : Ticktock
            Ticktock instance with new time coordinates

        Examples
        ========
        >>> a = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO')
        >>> s = a.convert('TAI')
        >>> type(s)
        <class 'time.Ticktock'>
        >>> s
        Ticktock( [1391342432 1391342432] ), dtype=TAI

        See Also
        ========
        a.CDF
        a.ISO
        a.UTC
        """
        newdat = eval('self.'+dtype)
        return Ticktock(newdat, dtype)

    # -----------------------------------------------
    def append(self, other):
        """
        a.append(other)

        Will be called when another Ticktock instance has to be appended to the current one

        Parameters
        ==========
        other : Ticktock
            other (Ticktock instance)
        """
        otherdata = eval('other.'+self.dtype)
        newobj = Ticktock(np.append(self.data, otherdata), dtype=self.dtype)
        return newobj

    # -----------------------------------------------
    def getCDF(self):
        """
        a.getCDF() or a.CDF

        Return CDF time which is milliseconds since 01-Jan-0000 00:00:00.000.
        "Year zero" is a convention chosen by NSSDC to measure epoch values.
        This date is more commonly referred to as 1 BC. Remember that 1 BC was a leap year.
        The CDF date/time calculations do not take into account the changes to the Gregorian
        calendar, and cannot be directly converted into Julian date/times.

        Returns
        =======
        out : numpy array
            days elapsed since Jan. 1st

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.CDF
        array([  6.31798704e+13])

        See Also
        ========
        getUTC
        getUNX
        getRDT
        getJD
        getMJD
        getISO
        getTAI
        getDOY
        geteDOY
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

        Returns
        =======
        out : numpy array
            day of the year

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.DOY
        array([ 33])

        See Also
        ========
        getUTC
        getUNX
        getRDT
        getJD
        getMJD
        getISO
        getTAI
        getDOY
        geteDOY
        """
        nTAI = len(self.data)
        DOY = np.zeros(nTAI)

        for i in np.arange(nTAI):
            DOY[i] = self.UTC[i].toordinal() - datetime.date(self.UTC[i].year, 1, 1).toordinal() + 1

        self.DOY = DOY.astype(int)
        return DOY

    # -----------------------------------------------
    def geteDOY(self):
        """
        a.eDOY or a.geteDOY()

        extract eDOY (elapsed days since midnight January 1st of given year)

        Returns
        =======
        out : numpy array
            days elapsed since midnight bbedJan. 1st

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.eDOY
        array([ 32.5])

        See Also
        ========
        getUTC
        getUNX
        getRDT
        getJD
        getMJD
        getISO
        getTAI
        getDOY
        geteDOY
        """

        nTAI = len(self.data)
        eDOY = np.zeros(nTAI)

        for i in np.arange(nTAI):
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

        Returns
        =======
        out : numpy array
            elapsed days since 12:00 January 1, 4713 BC

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.JD
        array([ 2452308.])

        See Also
        ========
        getUTC
        getUNX
        getRDT
        getJD
        getMJD
        getISO
        getTAI
        getDOY
        geteDOY
        """
        import decimal

        nTAI = len(self.data)

        # convert all types in to UTC first and call again
        UTCdata = self.UTC

        if UTCdata[0] < datetime.datetime(1582,10,15):
            print("WARNING: Calendar date before the switch from Julian to Gregorian")
            print("    Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        # include offset if given
        JD = np.zeros(nTAI)
        for i in np.arange(nTAI):
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
            JD[i] = int(365.25*JY) + int(30.6001*JM) + D + 1720995
            c_val = (D+31*(M+12*Y))
            if c_val >= igreg: # yes if date after the Gregorian Switch in 1582-Oct-15
                JA = int(0.01*JY)
                JD[i] = JD[i]+2-JA+int(0.25*JA)

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

        Returns
        ========
            out : numpy array
                elapsed days since November 17, 1858
                    (Julian date was 2,400 000)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.MJD
        array([ 52307.5])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getISO, getCDF, getTAI, getDOY, geteDOY

        """

        if self.UTC[0] < datetime.datetime(1582,10,15):
            print("WARNING: Calendar date before the switch from Julian to Gregorian")
            print("Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        MJD = self.JD - 2400000.5

        self.MJD = MJD
        return MJD

    # -----------------------------------------------
    def getUNX(self):
        """
        a.UNX or a.getUNX()

        convert dtype data into Unix Time (Posix Time)
        seconds since 1970-Jan-1 (not counting leap seconds)

        Returns
        ========
            out : numpy array
                elapsed secs since 1970/1/1 (not counting leap secs)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.UNX
        array([  1.01265120e+09])

        See Also
        =========
        getUTC, getISO, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """

        nTAI = len(self.data)

        UNX0 = datetime.datetime(1970,1,1)
        d = ['']*nTAI
        UNX = np.zeros(nTAI)
        for i in np.arange(nTAI):
            d[i] = self.UTC[i] - UNX0 # timedelta object (only days, seconds, microsecs are stored)
            UNX[i] = (d[i].days)*86400 + d[i].seconds + d[i].microseconds/1.e6

        self.UNX = UNX
        return UNX

    # -----------------------------------------------
    def getRDT(self):
        """
        a.RDT or a.RDT()

        convert dtype data into Rata Die (lat.) Time (days since 1/1/0001)

        Returns
        ========
            out : numpy array
                elapsed days since 1/1/1

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.RDT
        array([ 730883.5])

        See Also
        =========
        getUTC, getUNX, getISO, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """

        import matplotlib.dates as mpd

        nTAI = len(self.data)
        UTC = self.UTC
        #RDT = np.zeros(nTAI)
        RDT = mpd.date2num(UTC)
        #for i in np.arange(nTAI):
            #RDT[i] = UTC[i].toordinal() + UTC[i].hour/24. + UTC[i].minute/1440. + \
                #UTC[i].second/86400. + UTC[i].microsecond/86400000000.

        self.RDT = RDT
        return RDT

    # -----------------------------------------------
    def getUTC(self):
        """
        a.UTC or a.getUTC()

        convert dtype data into UTC object a la datetime()

        Returns
        ========
            out : list of datetime objects
                datetime object in UTC time

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.UTC
        [datetime.datetime(2002, 2, 2, 12, 0)]

        See Also
        =========
        getISO, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """

        fmt,fmt2 = '%Y-%m-%dT%H:%M:%S','%Y-%m-%dT%H:%M:%S.%f'

        nTAI = len(self.data)

        UTC = ['']*nTAI
        if self.dtype.upper() == 'UTC':
            UTC = self.data # return

        elif self.dtype.upper() == 'ISO':
            for i in np.arange(nTAI):
                if len(self.data[i])==19:
                    UTC[i] = datetime.datetime.strptime(self.data[i], fmt)
                else:
                    UTC[i] = datetime.datetime.strptime(self.data[i], fmt2)

        elif self.dtype.upper() == 'TAI':
            TAI0 = datetime.datetime(1958,1,1,0,0,0,0)
            for i in np.arange(nTAI):
                UTC[i] = datetime.timedelta(seconds=float(self.data[i])) + TAI0
             # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i]))
                tmpleaps = Ticktock(self.UTC[i]).leaps
                if tmpleaps == leapsecs[i]-1: self.UTC[i] = self.UTC[i]+datetime.timedelta(seconds=1)

        elif self.dtype.upper() == 'GPS':
            GPS0 = datetime.datetime(1980,1,6,0,0,0,0)
            for i in np.arange(nTAI):
                UTC[i] = datetime.timedelta(seconds=float(self.data[i])) + GPS0
             # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                # there were 18 leap secinds before gps zero, need the -18 for that
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i])) + \
                    datetime.timedelta(seconds=19)

        elif self.dtype.upper() == 'UNX':
            UNX0 = datetime.datetime(1970,1,1)
            for i in np.arange(nTAI):
                UTC[i] = datetime.timedelta(seconds=self.data[i]) + UNX0 # timedelta object

        elif self.dtype.upper() == 'RDT':
            import matplotlib.dates as mpd
            UTC = mpd.num2date(self.data)
            UTC = [t.replace(tzinfo=None) for t in UTC]
            #for i in np.arange(nTAI):
                #UTC[i] = datetime.datetime(1,1,1) + \
                    #datetime.timedelta(days=np.floor(self.data[i])-1) +  \
                    #datetime.timedelta(microseconds=(self.data[i] - \
                        #self.data[i])*86400000.)
                # roundoff the microseconds
                #UTC[i] = UTC[i] - datetime.timedelta(microseconds=UTC[i].microsecond)

        elif self.dtype.upper() == 'CDF':
            for i in np.arange(nTAI):
                UTC[i] = datetime.timedelta(days=self.data[i]/86400000.) + \
                        datetime.datetime(1,1,1) - datetime.timedelta(days=366)
                #UTC[i] = datetime.timedelta(days=np.floor(self.data[i]/86400000.), \
                    #milliseconds=np.mod(self.data[i],86400000)) + \
                        #datetime.datetime(1,1,1) - datetime.timedelta(days=366)
                # the following has round off errors
                # UTC[i] = datetime.timedelta(data[i]/86400000.-366) + datetime.datetime(1,1,1)

        elif self.dtype.upper() in ['JD', 'MJD']:
            if self.dtype.upper() == 'MJD':
                self.JD = np.array(self.data) + 2400000.5
            for i in np.arange(nTAI):
                # extract partial days
                ja = int(np.floor(self.JD[i]))
                p = self.JD[i] - np.floor(self.JD[i])
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
                    print("WARNING: Calendar date before the switch from Julian to Gregorian")
                    print("Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        else:
            print("ERROR: Data type ", self.dtype, ' in getUTC() not supported')
            return

        self.UTC = UTC
        return UTC

    # -----------------------------------------------
    def getGPS(self):
        """
        a.GPS or a.getGPS()

        return GPS epoch (0000 UT (midnight) on January 6, 1980)

        Returns
        ========
            out : numpy array
                elapsed secs since 6Jan1980 (excludes leap secs)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.GPS
        array([])

        See Also
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY

        """

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

        Returns
        =======
            out : numpy array
                elapsed secs since 1958/1/1 (includes leap secs,
                    i.e. all secs have equal lengths)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.TAI
        array([1391342432])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY

        """

        fmt = '%Y-%m-%dT%H:%M:%S'
        TAI0 = datetime.datetime(1958,1,1,0,0,0,0)

        nTAI = len(self.data)
        TAI = np.zeros(nTAI)
        UTC = self.UTC
        leapsec = self.getleapsecs()
        TAItup = ['']*nTAI
        for i in np.arange(nTAI):
            #t = time.strptime(data[i], fmt)
            #dtimetup = datetime.datetime(t[0], t[1], t[2], t[3], t[4], t[5])
            # get the leap seconds
            TAItup[i] = UTC[i] - TAI0 + datetime.timedelta(seconds=int(leapsec[i]))
            TAI[i] = TAItup[i].days*86400 + TAItup[i].seconds + TAItup[i].microseconds/1.e6

        self.TAI = np.array(TAI)
        return self.TAI

    # -----------------------------------------------
    def getISO(self):
        """
        a.ISO or a.getISO()

        convert dtype data into ISO string

        Returns
        =======
            out : list of strings
                date in ISO format

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.ISO
        ['2002-02-02T12:00:00']

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """

        nTAI = len(self.data)
        ISO = ['']*nTAI
        for i in range(nTAI):
            ISO[i] = self.UTC[i].strftime(self.__isofmt)

            if self.TAI[i] in self.TAIleaps:
                tmptick = Ticktock(self.UTC[i] - datetime.timedelta(seconds=1), 'UTC')
                a,b,c = tmptick.ISO[0].split(':')
                cnew = c.replace('59','60')
                ISO[i] = a+':'+b+':'+cnew

        self.ISO = ISO
        return ISO

    # -----------------------------------------------
    def getleapsecs(self):
        """
        a.leaps or a.getleapsecs()

        retrieve leapseconds from lookup table, used in getTAI

        Returns
        =======
            out : numpy array
                leap seconds

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.leaps
        array([32])

        See Also
        =========
        getTAI

        """

        import os
        from spacepy import DOT_FLN


        tup = self.UTC
        # so you don't have to read the file every single time
        global secs, year, mon, day, TAIleaps

        try:
           leaps = secs[0]

        except:  # then we are calling this routine the 1st time
           # load current file
           fname = DOT_FLN+'/data/tai-utc.dat'
           fh = open(fname)
           text = fh.readlines()

           secs = np.zeros(len(text))
           year = np.zeros(len(text))
           mon = np.zeros(len(text))
           day = np.zeros(len(text))

           months = np.array(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', \
                  'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])

           for line, i in zip(text, np.arange(len(secs))):
              secs[i] = int(float(line.split()[6]))  # truncate float seconds
              year[i] = int(line.split()[0])
              mon[i] = int(np.where(months == line.split()[1])[0][0] + 1)
              day[i] = int(line.split()[2])

           TAIleaps = np.zeros(len(secs))
           TAItup = ['']*len(secs)
           TAI0 = datetime.datetime(1958,1,1,0,0,0,0)
           for i in np.arange(len(secs)):
                TAItup[i] = datetime.datetime(int(year[i]), int(mon[i]), int(day[i])) - TAI0 + datetime.timedelta(seconds=int(secs[i])-1)
                TAIleaps[i] = TAItup[i].days*86400 + TAItup[i].seconds + TAItup[i].microseconds/1.e6

        # check if array:
        if type(tup) == type(datetime.datetime(1,1,1)): # not an array of objects
            tup = [tup]
            nTAI = 1
            aflag = False
        else:
            nTAI = len(tup)
            aflag = True

        # convert them into a time tuple and find the correct leap seconds
        self.TAIleaps = TAIleaps
        leaps = [secs[0]]*nTAI
        for i, itup in enumerate(tup):
            for y,m,d,s in zip(year, mon, day, secs):
                if tup[i] >= datetime.datetime(int(y),int(m),int(d)):
                    leaps[i] = s
                else:
                    break

        #if datetime.datetime(1971,12,31) > tup[0]:
        #   print "WARNING: date before 1972/1/1; leap seconds are by fractions off"

        if aflag == False:
            self.leaps = int(leaps[0])
            return int(leaps[0])   # if you want to allow fractional leap seconds, remove 'int' here
        else:
            self.leaps = np.array(leaps, dtype=int)
            return self.leaps

    # -----------------------------------------------
    @classmethod
    def now(self):
        """
        Creates a Ticktock object with the current time, equivalent to datetime.now()

        Returns
        =======
            out : ticktock
                Ticktock object with the current time, equivalent to datetime.now()

        See Also
        ========
        datetime.datetime.now()

        """
        dt = datetime.datetime.now()
        return Ticktock(dt, 'utc')



# -----------------------------------------------
# End of Ticktock class
# -----------------------------------------------


def doy2date(year, doy, dtobj=False, flAns=False):
    """
    convert integer day-of-year doy into a month and day
    after http://pleac.sourceforge.net/pleac_python/datesandtimes.html

    Parameters
    ==========
        year : int or array of int
            year
        doy : int or array of int
            day of year

    Returns
    =======
        month : int or array of int
            month as integer number
        day : int or array of int
            as integer number

    Examples
    ========
    >>> month, day = doy2date(2002, 186)
    >>> dts = doy2date([2002]*4, range(186,190), dtobj=True)

    See Also
    ========
    getDOY

    """
    try:
        n_year = len(year)
    except TypeError:
        n_year = -1 #Special case: this is a scalar
    try:
        n_doy = len(doy)
    except TypeError:
        n_doy = -1
    if n_doy != n_year:
        raise ValueError('Day of year and year must have same length')

    if n_year == -1:
        if doy < 1:
            raise ValueError(
                'Day-of-Year less than 1 detected: DOY starts from 1')
        if not flAns:
            dt = datetime.datetime(int(year), 1, 1) + \
                 datetime.timedelta(days=int(doy) - 1)
        else:
            dt = datetime.datetime(year, 1, 1) + \
                 datetime.timedelta(days=float(doy) - 1)
        if dtobj:
            return dt
        else:
            return dt.month, dt.day

    if min(doy) < 1:
        raise ValueError('Day-of-Year less than 1 detected: DOY starts from 1')

    if flAns:
        dateobj = [datetime.datetime(year[i], 1, 1) +
                   datetime.timedelta(days=float(doy[i]) - 1)
                   for i in range(n_year)]
    else:
        dateobj = [datetime.datetime(int(year[i]), 1, 1) +
                   datetime.timedelta(days=int(doy[i]) - 1)
                   for i in range(n_year)]
    if dtobj:
        return dateobj
    else:
        return [dt.month for dt in dateobj], [dt.day for dt in dateobj]



# -----------------------------------------------
def tickrange(start, end, deltadays, dtype='ISO'):
    """
    return a Ticktock range given the start, end, and delta

    Parameters
    ==========
        start : string or number
            start time
        end : string or number
            end time (inclusive)
        deltadays : float or timedelta
            step in units of days (float); or datetime timedelta object
        dtype : string (optional)
            data type for start, end; e.g. ISO, UTC, RTD, etc. see Ticktock for all options

    Returns
    =======
        out : Ticktock instance
            ticks

    Examples
    ========
    >>> ticks = st.tickrange('2002-02-01T00:00:00', '2002-02-10T00:00:00', deltadays = 1)
    >>> ticks
    Ticktock( ['2002-02-01T00:00:00', '2002-02-02T00:00:00', '2002-02-03T00:00:00',
    '2002-02-04T00:00:00'] ), dtype=ISO

    See Also
    ========
    Ticktock

    """
    Tstart = Ticktock(start, dtype)
    Tend = Ticktock(end, dtype)
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
    ticks = Ticktock(trange, 'UTC')
    ticks = eval('Ticktock(ticks.'+dtype+',"'+dtype+'")')
    return ticks


def sec2hms(sec, rounding=True, days=False, dtobj=False):
    """Convert seconds of day to hours, minutes, seconds

    Parameters
    ==========
        sec : float
            Seconds of day

    Other Parameters
    ================
        rounding : boolean
            set for integer seconds
        days : boolean
            set to wrap around day (i.e. modulo 86400)
        dtobj : boolean
            set to return a timedelta object

    Returns
    =======
        out : [hours, minutes, seconds] or datetime.timedelta

    """
    if not days:
        try:
            assert sec <= 86400
        except:
            print("Warning: Number of seconds > seconds in day. Try days keyword.")
    else:
        sec %= 86400

    hours = int(sec)//3600
    try:
        minutes = int((sec - hours*3600) // 60) % 60
    except ZeroDivisionError:
        minutes = 0
    seconds = sec % 60
    if rounding:
        seconds = int(round(seconds))

    if dtobj:
        return dt.timedelta(hours=hours, minutes=minutes, seconds=seconds)
    else:
        return [hours, minutes, seconds]


# -----------------------------------------------
def test():
    """
    test all time conversions

    Returns
    ========
        out : int
            number of failures

    Examples
    ========
    >>> test()
    testing ticks: PASSED TEST 1
    testing ticks: PASSED TEST 2
    0

    """
    from . import time as st
    from . import toolbox as tb
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
    t = st.Ticktock(data,dtype)

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
    comp = [st.Ticktock(data,'ISO').getUTC()[0]]*len(alldtypes)
    #
    if back == comp and tb.feq(round(out[2],5),2452331.01424) and \
    tb.feq(out[0],1014639630.1):
        print("testing Ticktock: PASSED TEST simple")
    else:
        print("testing Ticktock: FAILED TEST simple")
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
        foo[i] = st.Ticktock(out[i], dtype)
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
    testTAI = np.where(tb.feq(TAI[0],TAI, precision=prec))
    testUNX = np.where(tb.feq(UNX[0],UNX, precision=prec))
    testJD = np.where(tb.feq(JD[0],JD, precision=prec))
    testMJD = np.where(tb.feq(MJD[0],MJD, precision=prec))
    testRDT = np.where(tb.feq(RDT[0],RDT, precision=prec))
    testUTC = timecomp(UTC[0],UTC)
    testCDF = np.where(tb.feq(CDF[0],CDF, precision=prec))
    try:
        assert len(testUNX[0]) == len(alldtypes)
        assert len(testTAI[0]) == len(alldtypes)
        assert len(testJD[0]) == len(alldtypes)
        assert len(testMJD[0]) == len(alldtypes)
        assert len(testRDT[0]) == len(alldtypes)
        assert False not in testUTC
        assert len(testCDF[0]) == len(alldtypes)
        #assert len(np.unique(ISO)) == 1

        print("testing Ticktock: PASSED TEST all combinations")
    except AssertionError:
        print("testing Ticktock: FAILED TEST all combinations")
        nFAIL =+ 1

    # test with arrays
    try:
        for i, dtype in enumerate(alldtypes):
            foo = st.Ticktock( np.array([out[i], out[i], out[i]]), dtype)
        print("testing Ticktock: PASSED TEST arrays")
    except:
        print("testing Ticktock: FAILED TEST arrays")
        nFAIL =+ 1

    # test DOY
    try:
        foo.DOY
        print("testing Ticktock: PASSED TEST DOY")
    except:
        print("testing Ticktock: FAILED TEST DOY")
        nFAIL =+ 1

    return nFAIL
