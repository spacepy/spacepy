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

Authors: Steve Morley, Josef Koller, Brian Larsen, Jon Niehof
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov,


Copyright 2010 Los Alamos National Security, LLC.
"""
import bisect
try:
    from collections.abc import Callable, MutableSequence
except ImportError:
    from collections import Callable, MutableSequence
import datetime
import decimal

try:
    from itertools import izip as zip
except ImportError:
    pass  # just use system zip. In python3 itertools.izip is just python zip
import os.path
import re
import warnings

import dateutil.parser as dup
import numpy as np

from spacepy import help
import spacepy.datamodel

__contact__ = 'Steve Morley, smorley@lanl.gov'

"""
Notes:
Python's assert statement is a good way to catch situations that should never happen. And it can be removed by Python
optimization when the code is trusted to be correct. Assert is not to be used in program flow.

In Python, on the other hand, there is no strict distinction between debug and release mode. The interpreter features
an "optimization flag" (-O), but currently this does not actually optimize the byte code, but only removes asserts.
"""


# -----------------------------------------------
# Ticktock class
# -----------------------------------------------
class Ticktock(MutableSequence):
    """
    Ticktock( data, dtype )

    Ticktock class holding various time coordinate systems
    (TAI, UTC, ISO, JD, MJD, UNX, RDT, CDF, DOY, eDOY)

    Possible input data types:

    ISO
        ISO standard format like '2002-02-25T12:20:30'
    UTC
        datetime object with UTC time
    TAI
        Elapsed seconds since 1958-1-1 (includes leap seconds)
    UNX
        Elapsed seconds since 1970-1-1 ignoring leapseconds (all days have
        86400 secs).
    JD
        Julian days elapsed
    MJD
        Modified Julian days
    RDT
        Rata Die days elapsed since 0001-1-1
    CDF
        CDF Epoch type: float milliseconds since 0000-1-1 ignoring leapseconds

    Possible output data types: All those listed above, plus:

    DOY
        Integer day of year, starts with day 1
    eDOY
        Fractional day of year, starts at day 0

    It is strongly recommended to access various time systems via the
    attributes listed above, as in the examples. They will be calculated
    automatically if necessary. Using the ``get`` methods will force a
    recalculation.

    The original input data will always be available as the ``data``
    attribute.

    Parameters
    ==========
    data : array_like (int, datetime, float, string)
        time stamp
    dtype : string {`CDF`, `ISO`, `UTC`, `TAI`, `UNX`, `JD`, `MJD`, `RDT`} or function
        data type for data, if a function it must convert input time format to Python datetime

    Returns
    =======
    out : Ticktock
        instance with self.data, self.dtype, self.UTC etc

    Notes
    =====
    UTC data type is implemented as Python datetime, which cannot represent
    leap seconds. The time within a leap second is regarded as not happening.

    The CDF data type is the older CDF_EPOCH time type, not the newer
    CDF_TIME_TT2000. It similarly cannot represent leap seconds. Year
    0 is considered a leap year.

    Examples
    ========
    >>> x = Ticktock([2452331.0142361112, 2452332.0142361112], 'JD')
    >>> x.ISO
    dmarray(['2002-02-25T12:20:30', '2002-02-26T12:20:30'], dtype='|S19')
    >>> x.DOY # Day of year
    dmarray([ 56.,  57.])
    >>> y = Ticktock(['01-01-2013', '20-03-2013'], lambda x: datetime.datetime.strptime(x, '%d-%m-%Y'))
    >>> y.UTC
    dmarray([2013-01-01 00:00:00, 2013-03-20 00:00:00], dtype=object)
    >>> y.DOY # Day of year
    dmarray([  1.,  79.])

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
    _keylist = ['UTC', 'TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'GPS', 'DOY', 'eDOY', 'leaps']
    _keylist_upper = [key.upper() for key in _keylist]
    _isoformatstr = {'seconds': '%Y-%m-%dT%H:%M:%S', 'microseconds': '%Y-%m-%dT%H:%M:%S.%f'}

    def __init__(self, data, dtype=None):
        self._isofmt = Ticktock._isoformatstr['seconds']

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

            if not isinstance(dtype, Callable):
                if isinstance(self.data[0], (str, bytes)):
                    dtype = 'ISO'
                elif str is bytes and isinstance(self.data[0], unicode): #Py2k
                    dtype = 'ISO'
                elif isinstance(self.data[0], datetime.datetime):
                    dtype = 'UTC'
                elif self.data[0] > 1e13:
                    dtype = 'CDF'
                elif dtype is None:
                    raise ValueError('Unable to guess dtype from data; '
                                     'please specify dtype.')
                if dtype.upper() not in Ticktock._keylist_upper:
                    raise ValueError("data type " + dtype + " not provided, only " + str(Ticktock._keylist))
            else:
                # process input data using callable dtype to convert to datetime/UTC
                dtype_func = np.vectorize(dtype)
                self.data = dtype_func(self.data)
                self.UTC = no_tzinfo(self.data)

        try:
            self.data.attrs['dtype'] = dtype.upper()
        except AttributeError:
            self.data.attrs['dtype'] = str(dtype_func)
        else:
            if dtype.upper() == 'ISO': self.ISO = self.data
            self.update_items('data')
            if dtype.upper() == 'TAI':
                self.TAI = self.data
            elif dtype.upper() == 'JD':
                self.JD = self.data
            elif dtype.upper() == 'MJD':
                self.MJD = self.data
            elif dtype.upper() == 'UNX':
                self.UNX = self.data
            elif dtype.upper() == 'RDT':
                self.RDT = self.data
            elif dtype.upper() == 'CDF':
                self.CDF = self.data
            elif dtype.upper() == 'UTC':
                self.UTC = no_tzinfo(self.data)

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
        return 'Ticktock( ' + str(self.data) + ', dtype=' + str(self.data.attrs['dtype'] + ')')

    __repr__ = __str__

    # -----------------------------------------------
    def __getstate__(self):
        """
        Is called when pickling
        See Also http://docs.python.org/library/pickle.html
        """
        odict = self.__dict__.copy()  # copy the dict since we change it
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
        self.update_items('data')

    # -----------------------------------------------
    def __delitem__(self, idx):
        """
        a.__delitem(index)

        will be called when deleting items in the sequence
        """
        self.data = np.delete(self.data, idx)
        self.update_items('data')

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

        Will be called if a timedelta object is subtracted from this instance and
        returns a new Ticktock instance. If a Ticktock is subtracted from another
        Ticktock then a list of timedeltas is returned.

        Parameters
        ==========
        other : Ticktock or datetime.timedelta
            instance for comparison

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> dt = datetime.timedelta(3)
        >>> a - dt
        Ticktock( ['2002-02-05T12:00:00'] ), dtype=ISO

        See Also
        ========
        __add__
        """
        if isinstance(other, datetime.timedelta):
            newobj = Ticktock(self.UTC - other, 'UTC')
        elif isinstance(other, Ticktock):
            if not (len(other) == len(self.data)) and not (len(other) == 1):
                raise ValueError('Ticktock lengths are mismatched, subtraction is not possible')
            same = True
            if len(other) == 1:
                same = False
            if same:
                return [datetime.timedelta(seconds=t - other.TAI[i])
                        for i, t in enumerate(self.TAI)]
            else:
                return [datetime.timedelta(seconds=t - other.TAI[0])
                        for t in self.TAI]
        elif hasattr(other, '__iter__'):
            if not isinstance(other[0], datetime.timedelta):
                raise TypeError("Data supplied for addition is of the wrong type")
            if not (len(other) == len(self.data)) or (len(other) == 1):
                raise TypeError("Data supplied for addition is of the wrong shape")
            same = True
            if len(other) == 1: same = False
            if same:
                newUTC = [utc - o for utc, o in zip(self.UTC, other)]
            else:
                newUTC = [utc - other for utc in self.UTC]
            newobj = Ticktock(newUTC, 'UTC')
        else:
            raise TypeError("unsupported operand type(s) for -: {0} and {1}".format(type(other), type(self)))
        return newobj

    # -----------------------------------------------
    def __add__(self, other):
        """
        a.__add__(other)

        Will be called if an iterable of datetime.timedeltas is added to this instance and
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
        elif hasattr(other, '__iter__'):
            if not isinstance(other[0], datetime.timedelta):
                raise TypeError("Data supplied for addition is of the wrong type")
            if not (len(other) == len(self.data)) or (len(other) == 1):
                raise TypeError("Data supplied for addition is of the wrong shape")
            same = True
            if len(other) == 1: same = False
            if same:
                newUTC = [utc + o for utc, o in zip(self.UTC, other)]
            else:
                newUTC = [utc + other for utc in self.UTC]
            newobj = Ticktock(newUTC, 'UTC')
        else:
            raise TypeError("unsupported operand type(s) for +: {0} and {1}".format(type(other), type(self)))

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
        assert name in Ticktock._keylist, "data type " + str(name) + " not provided, only " + str(Ticktock._keylist)
        if name.upper() == 'TAI': self.TAI = self.getTAI()
        if name.upper() == 'ISO': self.ISO = self.getISO()
        if name.upper() == 'JD': self.JD = self.getJD()
        if name.upper() == 'MJD': self.MJD = self.getMJD()
        if name.upper() == 'UNX': self.UNX = self.getUNX()
        if name.upper() == 'RDT': self.RDT = self.getRDT()
        if name.upper() == 'CDF': self.CDF = self.getCDF()
        if name.upper() == 'DOY': self.DOY = self.getDOY()
        if name.upper() == 'EDOY': self.eDOY = self.geteDOY()
        if name.upper() == 'GPS': self.GPS = self.getGPS()
        # if name == 'isoformat': self.__isofmt = self.isoformat()
        if name == 'leaps': self.leaps = self.getleapsecs()
        return getattr(self, name)

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
        ival = getattr(dum, fmt)
        self.data = np.insert(self.data, idx, ival)

        self.update_items('data')

    # -----------------------------------------------
    def remove(self, idx):
        """
        a.remove(idx)

        This will remove the Ticktock value at index idx
        """
        del self[idx]

    # -----------------------------------------------
    def sort(self, kind='quicksort'):
        """
        a.sort()

        This will sort the Ticktock values in place based on the values
        in `data`. If you need a stable sort use kind='mergesort'
        """
        idx = self.argsort(kind=kind)
        self.data = self.data[idx]
        self.update_items('data')

    # -----------------------------------------------
    def argsort(self, kind='quicksort'):
        """
        idx = a.argsort()

        This will return the indices that would sort the Ticktock values

        Returns
        =======
        out : list
            indices that would sort the Ticktock values

        """
        return np.argsort(self.TAI, kind=kind)

    # -----------------------------------------------
    def isoformat(self, fmt=None):
        """
        a.update_items(b, attrib)

        This changes the self._isofmt attribute by and subsequently this
        function will update the ISO attribute.

        Parameters
        ==========
        fmt : string, optional
        """
        if fmt is None:
            print('Current ISO output format is %s' % self._isofmt)
            print('Options are: {0}'.format([(k, Ticktock._isoformatstr[k]) for k in list(Ticktock._isoformatstr.keys())]))
        else:
            try:
                self._isofmt = Ticktock._isoformatstr[fmt]
                self.update_items('data')
            except KeyError:
                raise (ValueError('Not a valid option: Use {0}'.format(list(Ticktock._isoformatstr.keys()))))

    # -----------------------------------------------
    def update_items(self, attrib, cls=None):
        """
        a.update_items(attrib)

        After changing the self.data attribute by either __setitem__ or __add__ etc
        this function will update all other attributes. This function is
        called automatically in __add__, __init__, and __setitem__.

        ``UTC`` is always updated (even if it was not previously set.)

        Parameters
        ==========
        attrib : str
            attribute that was updated; update others from this

        Other Parameters
        ================
        cls : type
            .. deprecated:: 0.2.2
                Class is now taken from the ``self`` of the bound instance.

            Class to use for finding possible attributes; now ignored.

        See Also
        ========
        spacepy.Ticktock.__setitem__
        spacepy.Ticktock.__add__
        spacepy.Ticktock.__sub__
        """
        keylist = dir(self)
        # Formerly took arguments (cls, attrib) but there's nothing from
        # the class we can't get from the instance, so removed cls.
        # If we got two position arguments, though, that indicates the cls
        if cls is not None:
            attrib = cls
            warnings.warn(
                'cls argument of update_items was deprecated in 0.2.2'
                ' and will be ignored.',
                DeprecationWarning)
        keylist.remove('data')
        if attrib != 'data':
            keylist.remove(attrib)
            attrib = attrib.upper()
            if attrib != self.data.attrs['dtype']:
                # Repopulating based on a different dtype, so make a temp
                # Ticktock to do the conversion.
                cls = type(self)
                dt = self.data.attrs['dtype']
                self.data = getattr(
                    cls(getattr(self, attrib), dtype=attrib), dt)
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
            if key.upper() == 'EDOY': self.eDOY = self.geteDOY()
            if key.upper() == 'GPS': self.GPS = self.getGPS()
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
        newdat = getattr(self, dtype)
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
        otherdata = getattr(other, self.data.attrs['dtype'])
        return Ticktock(np.append(self.data, otherdata), dtype=self.data.attrs['dtype'])

    # -----------------------------------------------
    def getCDF(self):
        """
        a.getCDF() or a.CDF

        Return CDF Epoch time which is milliseconds since 0000-1-1 at
        00:00:00.000. "Year zero" is a convention chosen by NSSDC to measure
        epoch values. This date is more commonly referred to as 1 BC and is
        considered a leap year.

        The CDF date/time calculations do not take into account the change
        to the Gregorian calendar or leap seconds, and cannot be directly
        converted into Julian date/times.

        Returns ``data`` if it was provided in CDF; otherwise always
        recalculates from ``RDT`` and calls :meth:`getRDT` to do so,
        updating the ``RDT`` attribute.

        Updates the ``CDF`` attribute.

        Returns
        =======
        out : numpy array
            milliseconds since 0000-01-01T00:00:00 assuming no discontinuities.

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
        if self.data.attrs['dtype'] == 'CDF':
            # This should be the case from the constructor
            self.CDF = self.data
            return self.CDF
        RDTdata = self.getRDT()
        # RDT has 0001-01-01 as day 1, but this is day 3666
        # of CDF Epoch (since 0000-01-01 is day 0, and a leap year).
        CDF = (RDTdata + 365) * 86400000.0
        self.CDF = CDF
        return self.CDF

    # -----------------------------------------------
    def getDOY(self):
        """
        a.DOY or a.getDOY()

        extract DOY (days since January 1st of given year)

        Always recalculates from the current value of ``UTC``, which will
        be created if necessary.

        Updates the ``DOY`` attribute.

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
        DOY = [utc.toordinal() - datetime.date(utc.year, 1, 1).toordinal() + 1 for utc in self.UTC]
        self.DOY = spacepy.datamodel.dmarray(DOY, attrs={'dtype': 'DOY'}).astype(int)
        return self.DOY

    # -----------------------------------------------
    def geteDOY(self):
        """
        a.eDOY or a.geteDOY()

        extract eDOY (elapsed days since midnight January 1st of given year)

        Always recalculates from the current value of ``UTC``, which will
        be created if necessary.

        Updates the ``eDOY`` attribute.

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

        eDOY = [utc.toordinal() - datetime.date(utc.year, 1, 1).toordinal() for utc in self.UTC]
        eDOY = [edoy + utc.hour / 24. + utc.minute / 1440. + utc.second / 86400. + utc.microsecond / 86400000000.
                for edoy, utc in zip(eDOY, self.UTC)]

        self.eDOY = spacepy.datamodel.dmarray(eDOY, attrs={'dtype': 'eDOY'})
        return self.eDOY

    # -----------------------------------------------
    def getJD(self):
        """
        a.JD or a.getJD()

        convert dtype data into Julian Date (JD)


        Returns ``data`` if it was provided in JD; otherwise always
        recalculates from the current value of ``UTC`` (thus is
        not leap-second aware), which will be created if necessary.

        Updates the ``JD`` attribute.

        Returns
        =======
        out : numpy array
            elapsed days since 4713 BCE 01-01T12:00

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
        if self.data.attrs['dtype'] == 'JD':
            # This should be the case from the constructor
            self.JD = self.data
            return self.JD

        nTAI = len(self.data)

        # convert all types in to UTC first and call again
        UTCdata = self.UTC

        #In TAI terms, if self.TAI[0] < -11840601564.0:
        if UTCdata[0] < datetime.datetime(1582, 10, 15):
            warnings.warn("Calendar date before the switch from Julian to Gregorian\n" +
                          "    Calendar 1582-10-15: Use Julian Calendar dates as input")

        # include offset if given
        JD = spacepy.datamodel.dmarray(np.zeros(nTAI), attrs={'dtype': 'JD'})

        twelve, twofour, mind = decimal.Decimal('12.0'), decimal.Decimal('24.0'), decimal.Decimal('1440.0')
        sind, usind = decimal.Decimal('86400.0'), decimal.Decimal('86400000000.0')

        for i in np.arange(nTAI):
            offset = UTCdata[i].utcoffset()
            if offset:
                UTCdata[i] = UTCdata[i] - offset

            # extract year, month, day
            Y = UTCdata[i].year
            M = UTCdata[i].month
            D = UTCdata[i].day

            # the following is from Wikipedia (but is wrong by 2 days)
            # JDN = D-32075+1461*(Y+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12-3*((Y+4900+(M-14)/12)/100)/4
            # JD = JDN + (data.hour-12)/24. + data.minute/1440. + data.second/86400.

            # following Press, "Numerical Recipes", Fct: JULDAY, p. 10
            igreg = 15 + 31 * (10 + 12 * 1582)
            if M > 2:
                JY = Y
                JM = M + 1
            else:
                JY = Y - 1
                JM = M + 13
            JD[i] = int(365.25 * JY) + int(30.6001 * JM) + D + 1720995
            c_val = (D + 31 * (M + 12 * Y))
            if c_val >= igreg:  # yes if date after the Gregorian Switch in 1582-Oct-15
                JA = int(0.01 * JY)
                JD[i] = JD[i] + 2 - JA + int(0.25 * JA)

            # add this to num.recipes to get fractional days
            # twelve, twofour, mind = decimal.Decimal('12.0'), decimal.Decimal('24.0'), decimal.Decimal('1440.0')
            # sind, usind = decimal.Decimal('86400.0'), decimal.Decimal('86400000000.0')
            JD[i] = decimal.Decimal(str(JD[i])) + (decimal.Decimal(UTCdata[i].hour) - twelve) / twofour + \
                    decimal.Decimal(str(UTCdata[i].minute / 1440.)) + (decimal.Decimal(UTCdata[i].second) / sind) + \
                    (decimal.Decimal(UTCdata[i].microsecond) / usind)
            # JD[i] = JD[i] + (UTCdata[i].hour-12)/24. + UTCdata[i].minute/1440. + \
            # UTCdata[i].second/86400. + UTCdata[i].microsecond/86400000000.

        self.JD = JD
        return self.JD

    # -----------------------------------------------
    def getMJD(self):
        """
        a.MJD or a.getMJD()

        convert dtype data into MJD (modified Julian date)

        Returns ``data`` if it was provided in MJD; otherwise always
        recalculates from the current value of ``JD`` which will be
        created if necessary.

        Updates the ``MJD`` attribute.

        Returns
        ========
        out : numpy array
            elapsed days since 1858-11-17T00:00
            (Julian date of 1858-11-17T12:00 was 2 400 000)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.MJD
        array([ 52307.5])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getISO, getCDF, getTAI, getDOY, geteDOY
        """
        if self.data.attrs['dtype'] == 'MJD':
            # This should be the case from the constructor
            self.MJD = self.data
            return self.MJD

        self.MJD = self.JD - 2400000.5
        return self.MJD

    # -----------------------------------------------
    def getUNX(self):
        """
        a.UNX or a.getUNX()

        convert dtype data into Unix Time (Posix Time)
        seconds since 1970-1-1 (not counting leap seconds)

        Returns ``data`` if it was provided in UNX; otherwise always
        recalculates from the current value of ``UTC``, which will be
        created if necessary.

        Updates the ``UNX`` attribute.

        Returns
        ========
        out : numpy array
            elapsed secs since 1970-1-1 (not counting leap secs)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.UNX
        array([  1.01265120e+09])

        See Also
        =========
        getUTC, getISO, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """
        if self.data.attrs['dtype'] == 'UNX':
            # This should be the case from the constructor
            self.UNX = self.data
            return self.UNX

        UNX0 = datetime.datetime(1970, 1, 1)
        d = [utc - UNX0 for utc in self.UTC]
        UNX = [dd.days * 86400 + dd.seconds + dd.microseconds / 1.e6 for dd in d]

        self.UNX = spacepy.datamodel.dmarray(UNX, attrs={'dtype': 'UNX'})
        return self.UNX

    # -----------------------------------------------
    def getRDT(self):
        """
        a.RDT or a.RDT()

        convert dtype data into Rata Die (lat.) Time, or elapsed days
        counting 0001-01-01 as day 1. This is a naive conversion: it
        ignores the existence of leapseconds for fractional days and
        ignores the conversion from Julian to Gregorian calendar, i.e.
        it assumes Gregorian calendar infinitely into the past.

        Returns ``data`` if it was provided in RDT; otherwise always
        recalculates from the current value of ``UTC``, which will be
        created if necessary.

        Updates the ``RDT`` attribute.

        Returns
        ========
        out : numpy array
            elapsed days counting 1/1/1 as day 1.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.RDT
        array([ 730883.5])

        See Also
        =========
        getUTC, getUNX, getISO, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY

        """
        if self.data.attrs['dtype'] == 'RDT':
            # This should be the case from the constructor
            self.RDT = self.data
            return self.RDT

        from matplotlib.dates import date2num

        # This is how to do this without date2num
        # nTAI = len(self.data)
        # RDT = np.zeros(nTAI)
        # for i in np.arange(nTAI):
        # RDT[i] = UTC[i].toordinal() + UTC[i].hour/24. + UTC[i].minute/1440. + \
        # UTC[i].second/86400. + UTC[i].microsecond/86400000000.

        self.RDT = spacepy.datamodel.dmarray(date2num(self.UTC), attrs={'dtype': 'RDT'})
        return self.RDT

    # -----------------------------------------------
    def getUTC(self):
        """
        a.UTC or a.getUTC()

        convert dtype data into UTC object a la datetime()

        Always recalculates from ``data``, the provided input data.

        Updates the ``UTC`` attribute.

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
        from matplotlib.dates import num2date

        nTAI = len(self.data)

        # if already UTC, we are done, no conversion
        if self.data.attrs['dtype'].upper() == 'UTC':
            UTC = self.data

        elif self.data.attrs['dtype'].upper() == 'ISO':
            self.ISO = self.data
            data = self.data if isinstance(self.data[0], str) \
                   else self.data.astype('S' if str is bytes else 'U')
            # try a few special cases that are faster than dateutil.parser
            for strfmt in ('%Y-%m-%dT%H:%M:%S',
                           '%Y-%m-%dT%H:%M:%SZ',
                           '%Y-%m-%d',
                           '%Y%m%d',
                           '%Y%m%d %H:%M:%S'):
                try:
                    UTC = [datetime.datetime.strptime(isot, strfmt)
                           for isot in data]
                    break
                except ValueError:
                    continue
            else:
                UTC = [dup.parse(isot) for isot in data]
        elif self.data.attrs['dtype'].upper() == 'TAI':
            self.TAI = self.data
            TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)
            UTC = [datetime.timedelta(seconds=float(tait)) + TAI0 for tait in self.data]
            # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i]))
                tmpleaps = Ticktock(self.UTC[i]).leaps
                if tmpleaps == leapsecs[i] - 1: self.UTC[i] = self.UTC[i] + datetime.timedelta(seconds=1)

        elif self.data.attrs['dtype'].upper() == 'GPS':
            self.GPS = self.data
            GPS0 = datetime.datetime(1980, 1, 6, 0, 0, 0, 0)
            UTC = [datetime.timedelta(seconds=float(gpst)) + GPS0 for gpst in self.data]
            # add leap seconds after UTC is created
            self.UTC = UTC
            leapsecs = self.getleapsecs()
            for i in np.arange(nTAI):
                # there were 18 leap seconds before gps zero, need the -18 for that
                self.UTC[i] = UTC[i] - datetime.timedelta(seconds=float(leapsecs[i])) + \
                              datetime.timedelta(seconds=19)

        elif self.data.attrs['dtype'].upper() == 'UNX':
            self.UNX = self.data
            UNX0 = datetime.datetime(1970, 1, 1)
            if np.issubdtype(self.data.dtype, np.integer)\
                     and not issubclass(np.int64, int):
                # numpy integer will not be accepted by timedelta
                UTC = [datetime.timedelta(seconds=unxt.item()) + UNX0
                       for unxt in self.data]
            else:
                UTC = [datetime.timedelta(seconds=unxt) + UNX0
                       for unxt in self.data]

        elif self.data.attrs['dtype'].upper() == 'RDT':
            self.RDT = self.data
            UTC = num2date(self.data)
            UTC = no_tzinfo(UTC)
            # for i in np.arange(nTAI):
            # UTC[i] = datetime.datetime(1,1,1) + \
            # datetime.timedelta(days=np.floor(self.data[i])-1) +  \
            # datetime.timedelta(microseconds=(self.data[i] - \
            # self.data[i])*86400000.)
            # roundoff the microseconds
            # UTC[i] = UTC[i] - datetime.timedelta(microseconds=UTC[i].microsecond)

        elif self.data.attrs['dtype'].upper() == 'CDF':
            self.CDF = self.data
            UTC = [datetime.timedelta(days=cdft / 86400000.) +
                   datetime.datetime(1, 1, 1) - datetime.timedelta(days=366) for cdft in self.data]
            # UTC[i] = datetime.timedelta(days=np.floor(self.data[i]/86400000.), \
            # milliseconds=np.mod(self.data[i],86400000)) + \
            # datetime.datetime(1,1,1) - datetime.timedelta(days=366)
            # the following has round off errors
            # UTC[i] = datetime.timedelta(data[i]/86400000.-366) + datetime.datetime(1,1,1)

        elif self.data.attrs['dtype'].upper() in ['JD', 'MJD']:
            if self.data.attrs['dtype'].upper() == 'MJD':
                self.JD = self.data + 2400000.5
                self.MJD = self.data
            else:
                self.JD = self.data
            UTC = [''] * nTAI
            for i in np.arange(nTAI):
                # extract partial days
                ja = int(np.floor(self.JD[i]))
                p = self.JD[i] - np.floor(self.JD[i])
                # after Press: "Numerical Recipes"
                # http://www.rgagnon.com/javadetails/java-0506.html
                # only good for after 15-Oct-1582
                igreg = 15 + 31 * (10 + 12 * 1582)
                if ja >= igreg:  # after switching to Gregorian Calendar
                    jalpha = int(((ja - 1867216) - 0.25) / 36524.25)
                    ja = ja + 1 + jalpha - jalpha // 4

                jb = ja + 1524
                jc = int(6680.0 + ((jb - 2439870) - 122.1) / 365.25)
                jd = 365 * jc + jc // 4
                je = int((jb - jd) / 30.6001)
                day = jb - jd - int(30.6001 * je)
                month = je - 1
                if (month > 12): month = month - 12
                year = jc - 4715
                if (month > 2): year = year - 1
                if (year <= 0): year = year - 1

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

                UTC[i] = datetime.datetime(year, month, int(day)) + \
                         datetime.timedelta(hours=12) + \
                         datetime.timedelta(seconds=p * 86400)
                if UTC[i] < datetime.datetime(1582, 10, 15):
                    warnings.warn("WARNING: Calendar date before the switch from Julian to Gregorian\n" +
                                  "Calendar 1582-Oct-15: Use Julian Calendar dates as input")

        else:
            warnstr1 = 'Input data type {0} does not support calculation of UTC times'.format(self.data.attrs['dtype'])
            warnstr2 = 'Valid input dtypes are: {0}'.format(
                ', '.join([kk for kk in Ticktock._keylist if kk not in ['DOY', 'eDOY', 'leaps']]))
            raise TypeError('{0}\n{1}'.format(warnstr1, warnstr2))

        UTC = no_tzinfo(UTC)
        self.UTC = spacepy.datamodel.dmarray(UTC, attrs={'dtype': 'UTC'})
        return self.UTC

    # -----------------------------------------------
    def getGPS(self):
        """
        a.GPS or a.getGPS()

        return GPS epoch (1980-1-6T00:00 UT)

        Always recalculates from the current value of ``UTC``, which
        will be created if necessary.

        Updates the ``GPS`` attribute.

        Returns
        ========
            out : numpy array
                elapsed secs since 1980-1-6 (excludes leap secs)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.GPS
        array([])

        See Also
        =========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY

        """
        GPS0 = datetime.datetime(1980, 1, 6, 0, 0, 0, 0)
        leapsec = self.getleapsecs()

        GPStup = [utc - GPS0 + datetime.timedelta(seconds=int(ls)) - datetime.timedelta(seconds=19)
                  for utc, ls in zip(self.UTC, leapsec)]
        GPS = [gps.days * 86400 + gps.seconds + gps.microseconds / 1.e6 for gps in GPStup]
        self.GPS = spacepy.datamodel.dmarray(GPS, attrs={'dtype': 'GPS'})  # .astype(int)
        return self.GPS

    # -----------------------------------------------
    def getTAI(self):
        """
        a.TAI or a.getTAI()

        return TAI (International Atomic Time), elapsed secs since 1958-1-1
        (leap seconds are counted.) Ticktock's handling of TAI and UTC
        conversions treats the UTC second as always equal in length to the SI
        second (and thus TAI), ignoring frequency changes and fractional
        leap seconds from 1958 through 1972, i.e. the UTC to TAI offset
        is always treated as an integer, truncated (not rounded) from the
        value at the most recent leap second (or fraction thereof).

        Returns ``data`` if it was provided in TAI; otherwise always
        recalculates from the current value of ``UTC``, which will be
        created if necessary.

        Updates the ``TAI`` attribute.

        Returns
        =======
        out : numpy array
            TAI as seconds since 1958-1-1.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.TAI
        array([1391342432])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY

        """
        if self.data.attrs['dtype'] == 'TAI':
            # This should be the case from the constructor
            self.TAI = self.data
            return self.TAI

        TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)

        leapsec = self.getleapsecs()
        TAItup = [utc - TAI0 + datetime.timedelta(seconds=int(ls)) for utc, ls in zip(self.UTC, leapsec)]
        TAI = [tai.days * 86400 + tai.seconds + tai.microseconds / 1.e6 for tai in TAItup]

        self.TAI = spacepy.datamodel.dmarray(TAI, attrs={'dtype': 'TAI'})
        return self.TAI

    # -----------------------------------------------
    def getISO(self):
        """
        a.ISO or a.getISO()

        convert dtype data into ISO string

        Always recalculates from the current value of ``UTC``, which
        will be created if necessary.

        Updates the ``ISO`` attribute.

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
        self.TAI = self.getTAI()
        self.ISO = spacepy.datamodel.dmarray([utc.strftime(self._isofmt) for utc in self.UTC], attrs={'dtype': 'ISO'})
        for i in range(nTAI):
            if self.TAI[i] in self.TAIleaps:
                tmpdt = self.UTC[i] - datetime.timedelta(seconds=1)
                a, b, c = tmpdt.strftime(self._isofmt).split(':')
                cnew = c.replace('59', '60')
                self.ISO[i] = a + ':' + b + ':' + cnew

        return self.ISO

    # -----------------------------------------------
    def getleapsecs(self):
        """
        a.leaps or a.getleapsecs()

        retrieve leapseconds from lookup table, used in getTAI

        Always recalculates from the current value of ``UTC``, which
        will be created if necessary.

        Updates the ``leaps`` attribute.

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
            fname = os.path.join(DOT_FLN, 'data', 'tai-utc.dat')
            with open(fname) as fh:
                text = fh.readlines()
            # Some files have a "last checked" line at the top
            if text[0].startswith('Checked'):
                del text[0]

            secs = np.zeros(len(text))
            year = np.zeros(len(text))
            mon = np.zeros(len(text))
            day = np.zeros(len(text))

            months = np.array(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])

            for line, i in zip(text, np.arange(len(secs))):
                secs[i] = int(float(line.split()[6]))  # truncate float seconds
                year[i] = int(line.split()[0])
                mon[i] = int(np.where(months == line.split()[1])[0][0] + 1)
                day[i] = int(line.split()[2])

            TAIleaps = np.zeros(len(secs))
            TAItup = [''] * len(secs)
            TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)
            for i in np.arange(len(secs)):
                TAItup[i] = datetime.datetime(int(year[i]), int(mon[i]), int(day[i])) - TAI0 + datetime.timedelta(
                    seconds=int(secs[i]) - 1)
                TAIleaps[i] = TAItup[i].days * 86400 + TAItup[i].seconds + TAItup[i].microseconds / 1.e6

        # check if array:
        if isinstance(tup, datetime.datetime):  # not an array of objects
            tup = [tup]
            nTAI = 1
            aflag = False
        else:
            nTAI = len(tup)
            aflag = True

        # convert them into a time tuple and find the correct leap seconds
        self.TAIleaps = TAIleaps
        leaps = [secs[0]] * nTAI
        leap_dates = [datetime.datetime(int(y), int(m), int(d)) for
                      y, m, d, s in zip(year, mon, day, secs)]
        for i, itup in enumerate(tup):
            ind = bisect.bisect_right(leap_dates, tup[i])
            leaps[i] = secs[ind - 1] if ind > 0 else 0

        ## ldatetime = datetime.datetime # avoid an expensive lookup below
        ## for i, itup in enumerate(tup):
        ##     for y,m,d,s in zip(year, mon, day, secs):
        ##         if tup[i] >= ldatetime(int(y),int(m),int(d)):
        ##             leaps[i] = s
        ##         else:
        ##             break

        # if datetime.datetime(1971,12,31) > tup[0]:
        #   print "WARNING: date before 1972/1/1; leap seconds are by fractions off"

        if aflag == False:
            self.leaps = int(leaps[0])
            return int(leaps[0])  # if you want to allow fractional leap seconds, remove 'int' here
        else:
            self.leaps = np.array(leaps, dtype=int)
            return self.leaps

    # -----------------------------------------------
    @classmethod
    def now(cls):
        """
        Creates a Ticktock object with the current time, equivalent to datetime.now()

        Returns
        =======
            out : ticktock
                Ticktock object with the current time, equivalent to datetime.now()

        See Also
        ========
        datetime.datetime.now

        """
        dt = datetime.datetime.now()
        return Ticktock(dt, 'UTC')

    # -----------------------------------------------
    @classmethod
    def today(cls):
        """
        Creates a Ticktock object with the current date and time set to 00:00:00, equivalent to date.today() with time
        included

        Returns
        =======
            out : ticktock
                Ticktock object with the current time, equivalent to date.today() with time included

        See Also
        ========
        datetime.date.today()

        """
        dt = datetime.datetime.now()
        dt = dt.replace(hour=0, minute=0, second=0, microsecond=0)
        return Ticktock(dt, 'UTC')


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
    Ticktock.getDOY

    """
    try:
        n_year = len(year)
    except TypeError:
        n_year = -1  # Special case: this is a scalar
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
        dateobj = spacepy.datamodel.dmarray([datetime.datetime(year[i], 1, 1) +
                                             datetime.timedelta(days=float(doy[i]) - 1)
                                             for i in range(n_year)])
    else:
        dateobj = spacepy.datamodel.dmarray([datetime.datetime(int(year[i]), 1, 1) +
                                             datetime.timedelta(days=int(doy[i]) - 1)
                                             for i in range(n_year)])
    if dtobj:
        return dateobj
    else:
        return (spacepy.datamodel.dmarray([dt.month for dt in dateobj]),
                spacepy.datamodel.dmarray([dt.day for dt in dateobj]))


# -----------------------------------------------
def tickrange(start, end, deltadays, dtype=None):
    """
    return a Ticktock range given the start, end, and delta

    Parameters
    ==========
    start : string or number
        start time (ISO standard string and UTC/datetime do not require a dtype)
    end : string or number
        last possible time in series (excluded unless end=start+n*step for integer n)
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
    dmusec, dsec = diff.microseconds / 86400000000., diff.seconds / 86400.
    if isinstance(deltadays, datetime.timedelta):
        musec, sec = deltadays.microseconds / 86400000000., deltadays.seconds / 86400.
        deltat = musec + sec + deltadays.days
        nticks = int((dmusec + dsec + diff.days) / deltat + 1)
        trange = [Tstart.UTC[0] + deltadays * n for n in range(nticks)]
    else:
        nticks = int((dmusec + dsec + diff.days) / float(deltadays) + 1)
        trange = [Tstart.UTC[0] + datetime.timedelta(days=deltadays) * n for n in range(nticks)]
    ticks = Ticktock(trange, 'UTC')
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
    if rounding:
        sec = int(round(sec))
    if not days:
        if sec > 86400:
            warnings.warn("Number of seconds > seconds in day. "
                          "Try days keyword.")
    else:
        sec %= 86400
    if dtobj:  # no need to do the computation
        return datetime.timedelta(seconds=sec)
    else:
        hours = sec // 3600
        minutes = ((sec - hours * 3600) // 60) % 60
        seconds = sec % 60
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
    returnclass = type(dt)
    try:
        retval = [val.replace(tzinfo=None) for val in dt]
    except TypeError:  # was not an iterable, assume datetime
        return dt.replace(tzinfo=None)
    #special case: numpy ndarray - dmarray and masked array work, but not ndarray
    if returnclass is not np.ndarray:
        retval = returnclass(retval)
        if hasattr(dt, 'attrs') and hasattr(retval, 'attrs'):
            retval.attrs.update(dt.attrs)
        return retval
    else:
        return np.asarray(retval)


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
    year = np.asanyarray(year)
    mask400 = (year % 400) == 0
    mask100 = (year % 100) == 0
    mask4 = (year % 4) == 0
    isleap = (mask400 | mask4) & (~mask100 | mask400)
    if numdays:
        isleap = isleap.astype(int) + 365
    return isleap


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
        raise (ValueError('tzinfo for the input and output datetimes must match'))
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


def extract_YYYYMMDD(filename):
    """
    go through a string and extract the first valid YYYYMMDD as a datetime

    Parameters
    ==========
    filename : str
        string to parse for a YYYYMMDD format

    Returns
    =======
    out : (None, datetime.datetime)
        the datetime found in the string or None
    """
    # return a datetime if there is one from YYYYMMDD
    # be picky so don't match random numbers (1950-2049)
    m = re.search(r"(19[5-9]|20[0-4])\d(0\d|1[0-2])([0-2]\d|3[01])",
                  os.path.basename(filename))
    if not m:
        return None
    else:
        return datetime.datetime.strptime(m.group(), '%Y%m%d')


def valid_YYYYMMDD(inval):
    """
    if inval is valid YYYYMMDD return True, False otherwise
    """
    if re.search(r"(19[5-9]|20[0-4])\d(0\d|1[0-2])([0-2]\d|3[01])",
                 inval):
        return True
    else:
        return False
