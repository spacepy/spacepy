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

.. currentmodule:: spacepy.time

.. NOTE... there is an error with this reference

Authors: Josef Koller, Brian Larsen, Steve Morley, Jon Niehof
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov,


Copyright 2010 Los Alamos National Security, LLC.

"""

from spacepy import help
import spacepy.datamodel

import datetime, collections
import dateutil.parser as dup
import warnings

import numpy as np

__contact__ = 'Josef Koller, jkoller@lanl.gov'

# -----------------------------------------------
# Tickdelta class
# -----------------------------------------------
class Tickdelta(object):
    """
    Tickdelta( **kwargs )

    Tickdelta class holding timedelta similar to datetime.timedelta
    This can be used to add/subtract from Ticktock objects

    .. deprecated:: 0.1.3
       Use :class:`datetime.timedelta` instead.

    Parameters
    ==========
    days : float
        number of days in for the delta
    hours : float
        number of hours for the delta
    minutes : float
        number of minutes for the delta
    seconds : float
        number of seconds for the delta

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
    Ticktock
    """
    def __init__(self, **kwargs):
        warnings.warn('Tickdelta is deprecated, use datetime.timedelta',
                  DeprecationWarning)
        days, hours, minutes, seconds = [0,0,0,0]
        if 'days' in kwargs:  days = kwargs['days']
        if 'hours' in kwargs: hours = kwargs['hours']
        if 'minutes' in kwargs: minutes = kwargs['minutes']
        if 'seconds' in kwargs: seconds = kwargs['seconds']
        self.days = days + hours/24. + minutes/1440. + seconds/86400.
        self.hours = self.days*24.
        self.minutes = self.hours*60.
        self.seconds = self.minutes*60.
        self.timedelta = datetime.timedelta(days=float(self.days))
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
class Ticktock(collections.MutableSequence):
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

    .. currentmodule:: spacepy.time
    .. autosummary::
        ~Ticktock.append
        ~Ticktock.argsort
        ~Ticktock.convert
        ~Ticktock.getCDF
        ~Ticktock.getDOY
        ~Ticktock.getGPS
        ~Ticktock.getISO
        ~Ticktock.getJD
        ~Ticktock.getMJD
        ~Ticktock.getRDT
        ~Ticktock.getTAI
        ~Ticktock.getUNX
        ~Ticktock.getUTC
        ~Ticktock.geteDOY
        ~Ticktock.getleapsecs
        ~Ticktock.isoformat
        ~Ticktock.now
        ~Ticktock.sort
        ~Ticktock.update_items
    .. automethod:: append
    .. automethod:: argsort
    .. automethod:: convert
    .. automethod:: getCDF
    .. automethod:: getDOY
    .. automethod:: getGPS
    .. automethod:: getISO
    .. automethod:: getJD
    .. automethod:: getMJD
    .. automethod:: getRDT
    .. automethod:: getTAI
    .. automethod:: getUNX
    .. automethod:: getUTC
    .. automethod:: geteDOY
    .. automethod:: getleapsecs
    .. automethod:: isoformat
    .. automethod:: now
    .. automethod:: sort
    .. automethod:: update_items
    """
    def __init__(self, data, dtype=None):
        self._keylist = ['UTC','TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'GPS', 'DOY', 'eDOY', 'leaps']

        if isinstance(data, Ticktock):
            dtype = data.data.attrs['dtype']
            self.data = data.data
        else:
            try:
                spacepy.datamodel.dmarray(data)[0]
            except IndexError:
                self.data = spacepy.datamodel.dmarray([data])
            else:
                self.data = spacepy.datamodel.dmarray(data)

            if isinstance(self.data[0], str):
                dtype = 'ISO'
            elif isinstance(self.data[0], datetime.datetime):
                dtype = 'UTC'
            elif self.data[0] > 1e13:
                dtype = 'CDF'
            assert dtype.upper() in self._keylist, "data type " + dtype +" not provided, only "+str(self._keylist)
        self.data.attrs['dtype'] = dtype.upper()
        self.__isoformatstr = {'seconds': '%Y-%m-%dT%H:%M:%S', 'microseconds': '%Y-%m-%dT%H:%M:%S.%f'}
        self.__isofmt = self.__isoformatstr['seconds']

        if dtype.upper() == 'ISO':
            if self.data[0].find('Z'):
                for i in range(len(self.data)):
                    self.data[i] = self.data[i].split('Z')[0]
            if self.data[0].find('T') == -1: # then assume midnight
                self.data = spacepy.datamodel.dmarray([el + 'T00:00:00' for el in self.data], attrs={'dtype': dtype.upper()})
            self.ISO = self.data
        self.update_items(self, 'data')
        if dtype.upper() == 'TAI': self.TAI = self.data
        if dtype.upper() == 'JD' : self.JD = self.data
        if dtype.upper() == 'MJD': self.MJD = self.data
        if dtype.upper() == 'UNX': self.UNX = self.data
        if dtype.upper() == 'RDT': self.RDT = self.data
        if dtype.upper() == 'CDF': self.CDF = self.data
        if dtype.upper() == 'UTC': self.UTC = self.data


## Brian and Steve were looking at this to see about making plot work directly on the object
## is also making iterate as an array of datetimes
#    def __iter__(self):
#        i = 0
#        try:
#            while True:
#                v = self[i].UTC[0]
#                yield v
#                i += 1
#        except IndexError:
#            return

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
        Ticktock( ['2002-02-02T12:00:00'], dtype=ISO)
        """
        return 'Ticktock( '+str(self.data) + ', dtype='+str(self.data.attrs['dtype'] +')')
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
        return Ticktock(self.UTC[idx])

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
        tmp = Ticktock(vals)
        if len(tmp) > 1:
            self.data[idx] = tmp.__getattribute__(self.data.attrs['dtype'])[:]
        else:
            self.data[idx] = tmp.__getattribute__(self.data.attrs['dtype'])[0]
        self.update_items(self, 'data')

    # -----------------------------------------------
    def __delitem__(self, idx):
        """
        a.__delitem(index)

        will be called when deleting items in the sequence
        """
        self.data = np.delete(self.data, idx)

        self.update_items(self, 'data')

        #del self.data[idx]
        #self.update_items(self, 'data')

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
        return len(self.data)

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

        Will be called if a Tickdelta object is subtracted from this instance and
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
        __add__
        """
        if isinstance(other, datetime.timedelta):
            newobj = Ticktock(self.UTC - other, 'UTC')
        elif isinstance(other, Tickdelta):
            newUTC = [t - other.timedelta for t in self.UTC]
            newobj = Ticktock(newUTC, 'UTC')
            newobj.data = getattr(newobj, 'get' + self.data.attrs['dtype'])()
            newobj.data.attrs['dtype'] = self.data.attrs['dtype']
            newobj.dtype = self.data.attrs['dtype']
            newobj.update_items(self, 'data')
        elif isinstance(other, Ticktock):
            return [datetime.timedelta(seconds=t - other.TAI[0])
                    for t in self.TAI]
        else:
            raise TypeError("unsupported operand type(s) for -: {0} and {1}".format(type(other),type(self)))
        return newobj

    # -----------------------------------------------
    def __add__(self, other):
        """
        a.__add__(other)

        Will be called if a Tickdelta object is added to this instance and
        returns a new Ticktock instance

        Parameters
        ==========
        other : datetime.timedelta
            instance for comparison

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> import datetime as dt
        >>> delt = dt.timedelta(minutes=1)
        >>> a + delt
        Ticktock( ['2002-02-02T12:01:00'] , dtype=ISO)


        See Also
        ========
        __sub__

        """
        if isinstance(other, datetime.timedelta):
            newobj = Ticktock(self.UTC + other, 'UTC')
        elif isinstance(other, Tickdelta):
            newUTC = ['']*len(self.data)
            for i in range(len(self.data)):
                newUTC[i] = self.UTC[i] + other.timedelta
            newobj = Ticktock(newUTC, 'UTC')
            newobj.data = eval('newobj.get'+self.data.attrs['dtype']+'()')
            newobj.dtype = self.data.attrs['dtype']
            newobj.update_items(self, 'data')
        else:
            raise TypeError("unsupported operand type(s) for +: {0} and {1}".format(type(other),type(self)))

        return newobj

    def __radd__(self, other):
        """
        a.__radd__(other)

        reverse add -- Will be called if this object is added to a datetime.timedelta and
        returns a new Ticktock instance

        Parameters
        ==========
        other : datetime.timedelta
            instance for comparison

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> import datetime as dt
        >>> delt = dt.timedelta(minutes=1)
        >>> a + delt
        Ticktock( ['2002-02-02T12:01:00'] , dtype=ISO)


        See Also
        ========
        __sub__

        """
        return self.__add__(other)

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
        assert name in self._keylist, "data type "+str(name)+" not provided, only "+str(self._keylist)
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
    def insert(self, idx, val, dtype=None):
        """
        insert values into the TickTock object

        .. note:: If more than one value to insert a slice must be specified
        as the index.  See numpy.insert

        Parameters
        ==========
            idx : int, slice or sequence of ints
                Object that defines the index or indices before which `val` is inserted.
            val : array_like
                values to insert
            dtype : str (optional)
                must be specified if not CDF, ISO, or UTC
        """
        fmt = self.data.attrs['dtype']
        if not dtype:
           dum = Ticktock(val)
        else:
           dum = Ticktock(val, dtype=dtype)
        ival = eval('dum.{0}'.format(fmt))
        self.data = np.insert(self.data, idx, ival)

        self.update_items(self, 'data')

    # -----------------------------------------------
    def remove(self, idx):
        """
        a.remove(idx)

        This will remove the Ticktock value at index idx
        """
        del self[idx]

    # -----------------------------------------------
    def sort(self):
        """
        a.sort()

        This will sort the Ticktock values in place
        """
        RDT = self.RDT
        RDTsorted = np.sort(RDT)
        tmp = Ticktock(RDTsorted, 'RDT').convert(self.data.attrs['dtype'])
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
                raise(ValueError('Not a valid option: Use %s' % list(self.__isoformatstr.keys())))

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
        spacepy.Ticktock.__setitem__
        spacepy.Ticktock.__add__
        spacepy.Ticktock.__sub__
        """
        keylist = list(cls.__dict__.keys())
        #keylist.remove('dtype')
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
        CDF
        ISO
        UTC
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
        otherdata = eval('other.'+self.data.attrs['dtype'])
        newobj = Ticktock(np.append(self.data, otherdata), dtype=self.data.attrs['dtype'])
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
        DOY = spacepy.datamodel.dmarray(np.zeros(nTAI))

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
        eDOY = spacepy.datamodel.dmarray(np.zeros(nTAI))

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
            warnings.warn("Calendar date before the switch from Julian to Gregorian\n" +
                "    Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        # include offset if given
        JD = spacepy.datamodel.dmarray(np.zeros(nTAI))

        twelve, twofour, mind = decimal.Decimal('12.0'), decimal.Decimal('24.0'), decimal.Decimal('1440.0')
        sind, usind = decimal.Decimal('86400.0'), decimal.Decimal('86400000000.0')

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
            # twelve, twofour, mind = decimal.Decimal('12.0'), decimal.Decimal('24.0'), decimal.Decimal('1440.0')
            # sind, usind = decimal.Decimal('86400.0'), decimal.Decimal('86400000000.0')
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
            elapsed days since November 17, 1858 (Julian date was 2,400 000)

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
            warnings.warn("WARNING: Calendar date before the switch from Julian to Gregorian\n" +
                "Calendar 1582-Oct-15: Use Julian Calendar dates as input")

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
        UNX = spacepy.datamodel.dmarray(np.zeros(nTAI))
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
        from matplotlib.dates import date2num, num2date
        #        import matplotlib.dates as mpd

        # nTAI = len(self.data)
        UTC = self.UTC
        #RDT = np.zeros(nTAI)
        RDT = date2num(UTC)
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
        from matplotlib.dates import date2num, num2date

        nTAI = len(self.data)

        if self.data.attrs['dtype'].upper() == 'UTC':
            UTC = self.data # return

        elif self.data.attrs['dtype'].upper() == 'ISO':
            self.ISO = self.data
            UTC = [dup.parse(isot) for isot in self.data]

        elif self.data.attrs['dtype'].upper() == 'TAI':
            self.TAI = self.data
            TAI0 = datetime.datetime(1958,1,1,0,0,0,0)
            UTC = [datetime.timedelta(seconds=float(tait)) + TAI0 for tait in self.data]
            # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i]))
                tmpleaps = Ticktock(self.UTC[i]).leaps
                if tmpleaps == leapsecs[i]-1: self.UTC[i] = self.UTC[i]+datetime.timedelta(seconds=1)

        elif self.data.attrs['dtype'].upper() == 'GPS':
            self.GPS = self.data
            GPS0 = datetime.datetime(1980,1,6,0,0,0,0)
            UTC = [datetime.timedelta(seconds=float(gpst)) + GPS0 for gpst in self.data]
            # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                # there were 18 leap secinds before gps zero, need the -18 for that
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i])) + \
                    datetime.timedelta(seconds=19)

        elif self.data.attrs['dtype'].upper() == 'UNX':
            self.UNX = self.data
            UNX0 = datetime.datetime(1970,1,1)
            UTC = [datetime.timedelta(seconds=unxt) + UNX0 for unxt in self.data]

        elif self.data.attrs['dtype'].upper() == 'RDT':
            self.RDT = self.data
            # import matplotlib.dates as mpd
            UTC = num2date(self.data)
            UTC = no_tzinfo(UTC)
            #for i in np.arange(nTAI):
                #UTC[i] = datetime.datetime(1,1,1) + \
                    #datetime.timedelta(days=np.floor(self.data[i])-1) +  \
                    #datetime.timedelta(microseconds=(self.data[i] - \
                        #self.data[i])*86400000.)
                # roundoff the microseconds
                #UTC[i] = UTC[i] - datetime.timedelta(microseconds=UTC[i].microsecond)

        elif self.data.attrs['dtype'].upper() == 'CDF':
            self.CDF = self.data
            UTC = [datetime.timedelta(days=cdft/86400000.) +
                        datetime.datetime(1,1,1) - datetime.timedelta(days=366) for cdft in self.data]
                #UTC[i] = datetime.timedelta(days=np.floor(self.data[i]/86400000.), \
                    #milliseconds=np.mod(self.data[i],86400000)) + \
                        #datetime.datetime(1,1,1) - datetime.timedelta(days=366)
                # the following has round off errors
                # UTC[i] = datetime.timedelta(data[i]/86400000.-366) + datetime.datetime(1,1,1)

        elif self.data.attrs['dtype'].upper() in ['JD', 'MJD']:
            if self.data.attrs['dtype'].upper() == 'MJD':
                self.JD = self.data + 2400000.5
                self.MJD = self.data
            else:
                self.JD = self.data
            UTC = ['']*nTAI
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
                    warnings.warn("WARNING: Calendar date before the switch from Julian to Gregorian\n" +
                       "Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        else:
            print("ERROR: Data type ", self.data.attrs['dtype'], ' in getUTC() not supported')
            return

        UTC = spacepy.datamodel.dmarray(UTC, attrs={'dtype': 'UTC'})
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
        # fmt = '%Y-%m-%dT%H:%M:%S'
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

        self.GPS = spacepy.datamodel.dmarray(GPS)#.astype(int)
        return self.GPS


    # -----------------------------------------------
    def getTAI(self):
        """
        a.TAI or a.getTAI()

        return TAI (International Atomic Time)

        Returns
        =======
        out : numpy array
            elapsed secs since 1958/1/1 (includes leap secs, i.e. all secs have equal lengths)

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

        self.TAI = spacepy.datamodel.dmarray(TAI)
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
        dmarray(['2002-02-02T12:00:00'])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """

        nTAI = len(self.data)
        ISO = ['']*nTAI
        self.TAI = self.getTAI()
        for i in range(nTAI):
            ISO[i] = self.UTC[i].strftime(self.__isofmt)

            if self.TAI[i] in self.TAIleaps:
                tmptick = Ticktock(self.UTC[i] - datetime.timedelta(seconds=1), 'UTC')
                a,b,c = tmptick.ISO[0].split(':')
                cnew = c.replace('59','60')
                ISO[i] = a+':'+b+':'+cnew

        self.ISO = spacepy.datamodel.dmarray(ISO)
        return self.ISO

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

#TODO add the deprication decorator when we have it
def num2date(*args, **kwargs):
    """
    Convert matplotlib epoch to datetime fast using an extension module

    .. deprecated:: 0.1.3

    Equivalent functionality to matplotlib.dates.num2date

    Parameters
    ==========
    mplnum : float or iterable of floats
        matplotlib epoch or iterable of epochs

    Returns
    =======
    out : np.array
        Array of datetime objects accosicated with the matplotlib epochs

    See Also
    ========
    matplotlib.dates.num2date
    """
    warnings.warn('num2date has been deprecated, see matplotlib.dates.num2date', DeprecationWarning)
    from matplotlib.dates import num2date
    return num2date(*args, **kwargs)

def date2num(*args, **kwargs):
    """
    Convert datetimes to matplotlib fast using an extension module

    .. deprecated:: 0.1.3

    Equivalent functionality to matplotlib.dates.date2num

    Parameters
    ==========
    dates : datetime.datetime or iterable of datetime.datetime
        datetime objects to convert to matplotlib epochs

    Returns
    =======
    out : np.array
        Array of floats accosicated with the datetimes

    See Also
    ========
    matplotlib.dates.date2num
    """
    warnings.warn('date2num has been deprecated, see matplotlib.dates.date2num', DeprecationWarning)
    from matplotlib.dates import date2num
    return date2num(*args, **kwargs)

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
    Ticktock.getDOY

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
def tickrange(start, end, deltadays, dtype='UTC'):
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
    '2002-02-04T00:00:00'] , dtype=ISO)

    See Also
    ========
    Ticktock

    """
    Tstart = Ticktock(start, dtype)
    Tend = Ticktock(end, dtype)
    diff = Tend.UTC[0] - Tstart.UTC[0]
    dmusec, dsec = diff.microseconds/86400000000., diff.seconds/86400.
    try:
        assert type(deltadays)==datetime.timedelta
        musec, sec = deltadays.microseconds/86400000000., deltadays.seconds/86400.
        deltat = musec + sec + deltadays.days
        nticks = int((dmusec + dsec + diff.days)/deltat + 1)
        trange = [Tstart.UTC[0] + deltadays*n for n in range(nticks)]
    except:
        nticks = int((dmusec + dsec + diff.days)/float(deltadays) + 1)
        trange = [Tstart.UTC[0] + datetime.timedelta(days=deltadays)*n for n in range(nticks)]
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
        if sec > 86400:
            warnings.warn("Number of seconds > seconds in day. "
                          "Try days keyword.")
    else:
        sec %= 86400

    hours = int(sec)//3600
    minutes = int((sec - hours*3600) // 60) % 60
    seconds = sec % 60
    if rounding:
        seconds = int(round(seconds))

    if dtobj:
        return datetime.timedelta(hours=hours, minutes=minutes, seconds=seconds)
    else:
        return [hours, minutes, seconds]


def no_tzinfo(dt):
    """
    take in an arraylike of datetime objects and return them without any tzinfo

    Parameters
    ==========
    dt : iterable
        iterable of datetime.datetime objects

    Returns
    =======
    out : list
        list of datetime.datetime without tzinfo
    """
    try:
        return [val.replace(tzinfo=None) for val in dt]
    except TypeError: # was not an iterable
        return dt.replace(tzinfo=None)


def leapyear(year, numdays=False):
    """
    return an array of boolean leap year,
    a lot faster than the mod method that is normally seen

    Parameters
    ==========
    year : array_like
        array of years
    numdays : boolean (optional)
        optionally return the number of days in the year

    Returns
    =======
    out : numpy array
        an array of boolean leap year, or array of number of days

    Examples
    ========
    >>> import numpy
    >>> import spacepy.time
    >>> spacepy.time.leapyear(numpy.arange(15)+1998)
    [False, False,  True, False, False, False,  True, False, False,
          False,  True, False, False, False,  True]
    """
    if not isinstance(year, (tuple, np.ndarray, list)):
        year = [year]
    mask400 = [(val % 400) == 0 for val in year]   # this is a leap year
    mask100 = [(val % 100) == 0 for val in year ]   # these are not leap years
    mask4   = [(val % 4) == 0 for val in year ]   # this is a leap year
    if numdays:
        numdays=365
        ans = [numdays + ((val[0] | val[2]) & (~val[1] | val[0])) for val in zip(mask400, mask100, mask4)]
        if len(ans) == 1:
            return ans[0]
        else:
            return ans
    else:
        ans = [bool(((val[0] | val[2]) & (~val[1] | val[0]))) for val in zip(mask400, mask100, mask4)]
        if len(ans) == 1:
            return ans[0]
        else:
            return ans


def randomDate(dt1, dt2, N=1, tzinfo=False, sorted=False):
    """
    Return a (or many) random datetimes between two given dates, this is done under the convention dt <=1 rand < dt2

    Parameters
    ==========
    dt1 : datetime.datetime
        start date for the the random date
    dt2 : datetime.datetime
        stop date for the the random date

    Other Parameters
    ================
    N : int (optional)
        the number of random dates to generate (defualt=1)
    tzinfo : bool (optional)
        maintain the tzinfo of the input datetimes (default=False)
    sorted : bool (optional)
        return the times sorted (default=False)

    Returns
    =======
    out : datetime.datetime or numpy.ndarray of datetime.datetime
        the new time for the next call to EventTimer

    Examples
    ========
    """
    from matplotlib.dates import date2num, num2date

    if dt1.tzinfo != dt2.tzinfo:
        raise(ValueError('tzinfo for the input and output datetimes must match'))
    dt1n = date2num(dt1)
    dt2n = date2num(dt2)
    rnd_tn = np.random.uniform(dt1n, dt2n, size=N)
    rnd_t = num2date(rnd_tn)
    if not tzinfo:
        tzinfo = None
    else:
        tzinfo = dt1.tzinfo
    rnd_t = np.asarray([val.replace(tzinfo=tzinfo) for val in rnd_t])
    if sorted:
        rnd_t.sort()
    return rnd_t
