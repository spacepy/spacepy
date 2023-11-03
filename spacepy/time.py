#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Time conversion, manipulation and implementation of Ticktock class

Notes
=====
The handling of time, in particular the conversions between representations,
can be more complicated than it seems on the surface. This can result in
some surprising behavior, particularly when requiring second-level accuracy and
converting between time systems outside of the period 1972 to present.
It is strongly recommended to use TAI if transferring times between SpacePy
and other libraries. TAI has a consistent, unambiguous definition and no
discontinuities.

Some time systems (e.g. the UTC representation via datetime) cannot represent
times during a leapsecond. SpacePy represents all these times as the latest
representable time in the day, e.g.::

    >>> spacepy.time.Ticktock('2008-12-31T23:59:60').UTC[0]
    datetime.datetime(2008, 12, 31, 23, 59, 59, 999999)

Conversions between continuous time representations (e.g. TAI), leap second
aware representations (e.g. ISO timestrings), and those that ignore leap
seconds (e.g. UTC datetime, Unix time) are well-defined between the
introduction of the leap second system to UTC in 1972 and the present.
For systems that cannot represent leap seconds, the leap second moment is
considered not to exist. For example, from 23:59:59 on 2008-12-31 to 00:00:00
on 2009-01-01 is two seconds, but only represents a one-second increment in
Unix time. Details are also discussed in the individual time representations.

UTC times more than six months in the future are not well-defined, since
the schedule of leap second insertion is not known in advance. SpacePy
performs conversions assuming there are no leapseconds after those which have
been announced by IERS.

Between 1960 and 1972, UTC was defined by means of fractional leap
seconds and a varying-length second. From 1958 (when UTC was set equal
to TAI) and 1972, SpacePy treats UTC time similar to after 1972, with
a consistent second the same length of the SI second, and applying a
full leap second before the beginning of January and July if UTC - UT1
exceeded 0.4s. The difference with other methods of calculating UTC is
less than half a second.

.. versionchanged:: 0.2.3
   The application of post-1972 rules to 1958-1927 is new in
   0.2.3. Before, SpacePy applied leap seconds wherever there was an
   entry in the USNO record of TAI-UTC, rounding fractional total leap
   second counts to the integer (0.5 rounds up). The UTC second was still
   treated as the same length as the SI second (i.e., rate changed were
   not applied.) This resulted in the application of six leap seconds at
   the beginning of 1972. The discrepancy with other means of calculating
   TAI-UTC was as much as five seconds by the end of this period.

.. versionchanged:: 0.2.2
   Before 0.2.2, SpacePy truncated fractional leapseconds rather than rounding.

Before 1958, UTC is not defined. SpacePy assumes days of constant length
86400 seconds, equal to the SI second. This is almost guaranteed to be wrong;
for times well out of the space era, it is strongly recommended to work
consistently in either a continuous time system (e.g. TAI) or a day-based
system (e.g. JD).

SpacePy assumes dates including and after 1582-10-15 to be in the Gregorian
calendar and dates including and before 1582-10-04 to be Julian. 10-05 through
10-14 do not exist. This change is ignored for continuously-running non leap
second aware timebases: CDF and RDT.

See the :class:`Ticktock` documentation and its various ``get`` functions for
more details on the exact definitions of time systems used by SpacePy.

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
from collections.abc import Callable, MutableSequence
import datetime

import os.path
import re
import time
import warnings

try:
    import astropy.time
    HAVE_ASTROPY = True
except ImportError:
    HAVE_ASTROPY = False
import dateutil.parser as dup
import numpy as np

import spacepy
from spacepy import help
import spacepy.datamodel

__contact__ = 'Steve Morley, smorley@lanl.gov'


# -----------------------------------------------
# Ticktock class
# -----------------------------------------------
class Ticktock(MutableSequence):
    """
    Ticktock( data, dtype )

    Ticktock class holding various time coordinate systems
    (TAI, UTC, ISO, JD, MJD, GPS, UNX, RDT, CDF, DOY, eDOY, APT)

    Possible input data types:

    ISO
        ISO standard format like '2002-02-25T12:20:30'
    UTC
        datetime object with UTC time
    TAI
        Elapsed seconds since 1958-1-1 (includes leap seconds)
    GPS
        Elapsed seconds since 1980-1-6 (includes leap seconds)
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
    APT
        AstroPy :class:`~astropy.time.Time`. Requires AstroPy 1.0.
        (New in version 0.2.2.)

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

    .. versionchanged:: 0.2.2
       In earlier versions of SpacePy, most values were derived from the
       ``datetime``-based ``UTC`` representation. This did not properly
       handle leap seconds in many cases. Now most are derived from ``TAI``
       (exceptions being ``DOY`` and ``eDOY``). In addition to differences
       around actual leap seconds, this may result in small differences
       between versions of SpacePy, with relative magnitude on the order of the
       resolution of a 64-bit float (2e-16). For times in the modern era, this
       is about 50 microseconds (us) for ``JD``, 15 us for ``CDF``, 1.5 us
       for ``RDT``, 1 us for ``MJD``, and 360 *nanoseconds* for ``TAI``.

    The relationships between parameters and how they are calculated are
    listed in the ``get`` methods and illustrated below.

    .. only:: latex

        .. image:: ../images/ticktock_relationships.*

    .. only:: html

        .. image:: ../images/ticktock_relationships.svg
            :target: ../_images/ticktock_relationships1.svg

    Parameters
    ==========
    data : array_like (int, datetime, float, string)
        time stamp
    dtype : string {`CDF`, `ISO`, `UTC`, `TAI`, 'GPS', `UNX`, `JD`, `MJD`, `RDT`, `APT`} or function
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

    Other Parameters
    ================
    isoformat : str, optional

        .. versionadded:: 0.2.2

        Format string used for parsing and outputting ISO format. Input is
        not forced to be in this format; it is tried first, and other
        common formats tried if parsing fails. Because of this, if ISO input
        is in a consistent format, specifying this can speed up the input
        parsing. Can be changed on existing ``Ticktock`` with :meth:`isoformat`
        method. Default ``'%Y-%m-%dT%H:%M:%S'``.

    See Also
    ========
    datetime.datetime.strptime, isoformat

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
        ~Ticktock.getAPT
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
        ~Ticktock.today
        ~Ticktock.update_items
    .. automethod:: append
    .. automethod:: argsort
    .. automethod:: convert
    .. automethod:: getAPT
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
    .. automethod:: today
    .. automethod:: update_items

    """
    _keylist = ['UTC', 'TAI', 'ISO', 'JD', 'MJD', 'UNX', 'RDT', 'CDF', 'GPS', 'DOY', 'eDOY', 'leaps']
    if HAVE_ASTROPY:
        _keylist.append('APT')
    _keylist_upper = [key.upper() for key in _keylist]
    _isoformatstr = {'seconds': '%Y-%m-%dT%H:%M:%S', 'microseconds': '%Y-%m-%dT%H:%M:%S.%f'}

    def __init__(self, data, dtype=None, isoformat=None):

        self._isofmt = isoformat or self._isoformatstr['seconds']

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
                elif isinstance(self.data[0], datetime.datetime):
                    dtype = 'UTC'
                elif HAVE_ASTROPY and isinstance(self.data[0],
                                                 astropy.time.Time):
                    dtype = 'APT'
                    # Recover original input (not dmarray), add axis if scalar
                    self.data = astropy.time.Time([data]) if data.shape == ()\
                                else data
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
        # AstroPy time objects don't have attrs, but will accept one.
        if not hasattr(self.data, 'attrs'):
            self.data.attrs = {}
        try:
            self.data.attrs['dtype'] = dtype.upper()
        except AttributeError:
            self.data.attrs['dtype'] = str(dtype_func)
        else:
            # ISO is populated by update_items
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
            elif dtype.upper() == 'APT':
                self.APT = self.data

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
        =======
            out: list, array
                requested values as either list/numpy array


        """
        if name not in Ticktock._keylist:
            raise AttributeError("data type {} not provided, only {}".format(
                str(name), str(Ticktock._keylist)))
        if name.upper() == 'TAI': self.TAI = self.getTAI()
        if name.upper() == 'UTC': self.UTC = self.getUTC()
        if name.upper() == 'ISO': self.ISO = self.getISO()
        if name.upper() == 'JD': self.JD = self.getJD()
        if name.upper() == 'MJD': self.MJD = self.getMJD()
        if name.upper() == 'UNX': self.UNX = self.getUNX()
        if name.upper() == 'RDT': self.RDT = self.getRDT()
        if name.upper() == 'CDF': self.CDF = self.getCDF()
        if name.upper() == 'DOY': self.DOY = self.getDOY()
        if name.upper() == 'EDOY': self.eDOY = self.geteDOY()
        if name.upper() == 'GPS': self.GPS = self.getGPS()
        if name.upper() == 'APT': self.APT = self.getAPT()
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

        Other Parameters
        ================
        kind : str
            Sort algorithm to use, default 'quicksort'.

        See Also
        ========
        argsort, numpy.argsort
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

        Other Parameters
        ================
        kind : str, optional
            Sort algorithm to use, default 'quicksort'.

            .. versionchanged:: 0.2.2
                Default is now 'quicksort' to match numpy default;
                previously was 'mergesort'.

        See Also
        ========
        argsort, numpy.argsort
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
    def update_items(self, attrib):
        """
        a.update_items(attrib)

        After changing the self.data attribute by either __setitem__ or __add__ etc
        this function will update all other attributes. This function is
        called automatically in __add__, __init__, and __setitem__.

        Parameters
        ==========
        attrib : str
            attribute that was updated; update others from this

        See Also
        ========
        spacepy.Ticktock.__setitem__
        spacepy.Ticktock.__add__
        spacepy.Ticktock.__sub__
        """
        keylist = dir(self)
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
        if self.data.attrs['dtype'] in (
                'TAI', 'GPS', 'JD', 'MJD', 'RDT', 'CDF', 'UNX', 'ISO', 'APT'):
            if self.data.attrs['dtype'] == 'ISO':
                if 'UTC' in keylist:
                    del self.UTC # Force recalc of UTC in TAI calc
                if 'ISO' in keylist:
                    del self.ISO # Also will recalc the ISO
                    del keylist[keylist.index('ISO')] # So no need to calc again
            self.TAI = self.getTAI()
            if 'UTC' in keylist and self.data.attrs['dtype'] != 'ISO':
                self.UTC = self.getUTC()
        else:
            self.UTC = self.getUTC()
            if 'TAI' in keylist:
                self.TAI = self.getTAI()
        for key in keylist:
            if key.upper() == 'ISO': self.ISO = self.getISO()
            if key.upper() == 'JD': self.JD = self.getJD()
            if key.upper() == 'MJD': self.MJD = self.getMJD()
            if key.upper() == 'UNX': self.UNX = self.getUNX()
            if key.upper() == 'RDT': self.RDT = self.getRDT()
            if key.upper() == 'CDF': self.CDF = self.getCDF()
            if key.upper() == 'DOY': self.DOY = self.getDOY()
            if key.upper() == 'EDOY': self.eDOY = self.geteDOY()
            if key.upper() == 'GPS': self.GPS = self.getGPS()
            if key.upper() == 'APT': self.APT = self.getAPT()
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
        recalculates from the current value of ``TAI``, which will be
        created if necessary.

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
        getAPT
        """
        if self.data.attrs['dtype'] == 'CDF':
            # This should be the case from the constructor
            self.CDF = self.data
            return self.CDF
        CDFofTAI0 = 61788528000000.
        naive_tai = _tai_real_to_naive(self.TAI)
        cdf = naive_tai * 1e3 + CDFofTAI0
        self.CDF = spacepy.datamodel.dmarray(cdf, attrs={'dtype': 'CDF'})
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
        getAPT
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
        getAPT
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
        recalculates from the current value of ``TAI``, which will be
        created if necessary.

        Updates the ``JD`` attribute.

        Returns
        =======
        out : numpy array
            elapsed days since 4713 BCE 01-01T12:00

        Notes
        =====
        This is based on the UTC day, defined as JD(UTC),
        per the recommendation
        of `IAU General Assembly XXIII resolution B1
        <https://www.iers.org/IERS/EN/Science/Recommendations/
        resolutionB1.html>`_.
        Julian days with leapseconds are 86401 seconds long and each second
        is a smaller fraction of the day. Note this "stretching" is across
        the *Julian* Day, noon to noon.

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
        getAPT
        """
        if self.data.attrs['dtype'] == 'JD':
            # This should be the case from the constructor
            self.JD = self.data
            return self.JD
        # 1958-01-01T12:00 is JD 2436205.0
        JDofTAI0 = 2436205.0 # Days-since-1958 is relative to noon
        self.JD = _days1958(self.TAI, leaps='rubber') + JDofTAI0
        return self.JD

    # -----------------------------------------------
    def getMJD(self):
        """
        a.MJD or a.getMJD()

        convert dtype data into MJD (modified Julian date)

        Returns ``data`` if it was provided in MJD; otherwise always
        recalculates from the current value of ``TAI`` which will be
        created if necessary.

        Updates the ``MJD`` attribute.

        Returns
        =======
        out : numpy array
            elapsed days since 1858-11-17T00:00
            (Julian date of 1858-11-17T12:00 was 2 400 000)

        Notes
        =====
        This is based on the UTC day, defined as JD(UTC) - 2 400 000.5,
        per the recommendation
        of `IAU General Assembly XXIII resolution B1
        <https://www.iers.org/IERS/EN/Science/Recommendations/
        resolutionB1.html>`_.
        Julian days with leapseconds are 86401 seconds long and each second
        is a smaller fraction of the day. Note this "stretching" is across
        the *Julian* Day not the MJD, so it will affect the last half of
        the MJD before the leap second and the first half of the following
        MJD, so that MJD is always JD - 2 400 000.5 This also means that
        the MJD following a leap second does not begin exactly at midnight.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.MJD
        array([ 52307.5])
        >>> a = Ticktock('2009-01-01T00:00:00', 'ISO')
        >>> a.MJD
        array([ 54832.00000579])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getISO, getCDF, getTAI, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'MJD':
            # This should be the case from the constructor
            self.MJD = self.data
            return self.MJD

        # 1958-01-01T12:00 is MJD 36204.5 (days since 1858-11-17T00:00)
        self.MJD = _days1958(self.TAI, leaps='rubber') + 36204.5
        return self.MJD

    # -----------------------------------------------
    def getUNX(self):
        """
        a.UNX or a.getUNX()

        convert dtype data into Unix Time (Posix Time)
        seconds since 1970-1-1 (not counting leap seconds)

        Returns ``data`` if it was provided in UNX; otherwise always
        recalculates from the current value of ``TAI``, which will be
        created if necessary.

        Updates the ``UNX`` attribute.

        Returns
        =======
        out : numpy array
            elapsed secs since 1970-1-1 (not counting leap secs)

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.UNX
        array([  1.01265120e+09])

        See Also
        ========
        getUTC, getISO, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'UNX':
            # This should be the case from the constructor
            self.UNX = self.data
            return self.UNX
        naive_tai = _tai_real_to_naive(self.TAI)
        UNXofTAI0 = -378691200.
        unx = naive_tai + UNXofTAI0

        self.UNX = spacepy.datamodel.dmarray(unx, attrs={'dtype': 'UNX'})
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
        recalculates from the current value of ``TAI``, which will be
        created if necessary.

        Updates the ``RDT`` attribute.

        Returns
        =======
        out : numpy array
            elapsed days counting 1/1/1 as day 1.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.RDT
        array([ 730883.5])

        See Also
        ========
        getUTC, getUNX, getISO, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'RDT':
            # This should be the case from the constructor
            self.RDT = self.data
            return self.RDT

        # RDT date at 1958-1-1T00
        RDTTAI0 = 714780.0
        RDT = _days1958(self.TAI, leaps='drop', midnight=True) + RDTTAI0
        # RDT can represent 1582-10-5 through 1582-10-14, which do not exist.
        # So everything before 1582-10-15 (RDT day 577736) is ten days earlier.
        RDT[RDT < 577736.0] -= 10

        self.RDT = spacepy.datamodel.dmarray(RDT, attrs={'dtype': 'RDT'})
        return self.RDT

    # -----------------------------------------------
    def getUTC(self):
        """
        a.UTC or a.getUTC()

        convert dtype data into UTC object a la datetime()

        Return value comes from (in priority order):

            1. If ``data`` was provided in UTC, returns ``data``.
            2. Else recalculates directly from ``data`` if it was
               provided in ISO.
            3. Else calculates from current value of ``TAI``, which
               will be created if necessary. (``data`` is TAI, GPS,
               JD, MJD, RDT, CDF, UNX)).

        Updates the ``UTC`` attribute.

        Returns
        =======
        out : list of datetime objects
            datetime object in UTC time

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.UTC
        [datetime.datetime(2002, 2, 2, 12, 0)]

        See Also
        ========
        getISO, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY,
        getAPT
        """
        nTAI = len(self.data)

        # if already UTC, we are done, no conversion
        if self.data.attrs['dtype'].upper() == 'UTC':
            UTC = self.data

        elif self.data.attrs['dtype'].upper() == 'ISO':
            self.ISO = self.data
            _, UTC, _ = dtstr2iso(self.data, fmt=self._isofmt)

        elif self.data.attrs['dtype'].upper() in (
                'TAI', 'GPS', 'JD', 'MJD', 'RDT', 'CDF', 'UNX', 'APT'):
            TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)
            UTC = [datetime.timedelta(
                # Before 1582-10-15, UTC 10 days earlier than naive conversion
                # since those dates are Julian not Gregorian.
                seconds=float(tait - (864000 if tait < -11840601600.0 else 0)))
                   + TAI0 for tait in self.TAI]
            for i in np.arange(nTAI):
                # This is the index of number of seconds to subtract from
                # "naive" UTC, not just TAI - UTC.
                # TAI of leap second does not have a new TAI - UTC, this
                # is because UTC seconds = 60, but need to subtract off
                # one more to make the UTC seconds = 59 in that case, thus
                # "flip" to next leap second count 1s earlier.
                idx = np.searchsorted(TAIleaps, self.TAI[i], side='right') - 1
                UTC[i] = UTC[i] - datetime.timedelta(seconds=secs[idx]
                                                          if idx > 0 else 0)
                if int(self.TAI[i]) == TAIleaps[idx]:
                    # TAI is in leap second
                    UTC[i] = UTC[i].replace(
                        second=59, microsecond=999999)

        else:
            warnstr1 = 'Input data type {0} does not support calculation of UTC times'.format(self.data.attrs['dtype'])
            warnstr2 = 'Valid input dtypes are: {0}'.format(
                ', '.join([kk for kk in Ticktock._keylist if kk not in ['DOY', 'eDOY', 'leaps']]))
            raise TypeError('{0}\n{1}'.format(warnstr1, warnstr2))

        self.UTC = spacepy.datamodel.dmarray(UTC, attrs={'dtype': 'UTC'})
        return self.UTC

    # -----------------------------------------------
    def getGPS(self):
        """
        a.GPS or a.getGPS()

        Return seconds since the GPS epoch (1980-1-6T00:00 UT)

        Always recalculates from the current value of ``TAI``, which
        will be created if necessary.

        Updates the ``GPS`` attribute.

        Returns
        =======
            out : numpy array
                elapsed secs since 1980-1-6. Leap seconds are counted;
                i.e. there are no discontinuities.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.GPS
        dmarray([6.96686413e+08])

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'GPS':
            # This should be the case from the constructor
            self.GPS = self.data
            return self.GPS
        GPS0 = 694656019
        self.GPS = spacepy.datamodel.dmarray(
            self.TAI - GPS0, attrs={'dtype': 'GPS'})
        return self.GPS

    # -----------------------------------------------
    def getAPT(self):
        """
        a.APT or a.getAPT()

        Return AstroPy time object.

        Always recalculates from the current value of ``TAI``, which
        will be created if necessary.

        Updates the ``APT`` attribute.

        Returns
        =======
            out : astropy.time.Time
                AstroPy Time object

        Notes
        =====
        .. versionadded:: 0.2.2

        Requires AstroPy 1.0. The returned value will be on the ``tai``
        scale in ``gps`` format (unless the :class:`Ticktock` was created
        from a :class:`~astropy.time.Time` object, in which case it will
        be the original input.) See the :mod:`astropy.time` docs for
        conversion to other scales and formats.

        Examples
        ========
        >>> a = Ticktock('2002-02-02T12:00:00', 'ISO')
        >>> a.APT
        <Time object: scale='tai' format='gps' value=696686413.0>

        See Also
        ========
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY,
        getGPS
        """
        if self.data.attrs['dtype'] == 'APT':
            # This should be the case from the constructor
            self.APT = self.data
            return self.APT
        if not HAVE_ASTROPY:
            raise RuntimeError('Import of astropy.time failed.')
        GPS0 = 694656019
        self.APT = astropy.time.Time(
            self.TAI - GPS0, scale='tai', format='gps')
        self.APT.attrs = {'dtype': 'APT'}
        return self.APT

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

        Return value comes from (in priority order):

            1. If ``data`` was provided in TAI, returns ``data``.
            2. Else recalculates directly from ``data`` if it was
               provided in APT, CDF, GPS, ISO, JD, MJD, RDT, or UNX.
            3. Else calculates from current value of ``UTC``, which
               will be created if necessary.

        Updates the ``TAI`` attribute; will also create the ``UTC``
        and ``ISO`` attributes from ``data`` if input is in ``ISO``
        (but will not overwrite an existing ``ISO`` or ``UTC``). This is
        for efficiency, as computation from ISO requires calculating UTC
        and makes creating a formatted ISO string easy.

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
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getISO, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'TAI':
            # This should be the case from the constructor
            self.TAI = self.data
            return self.TAI
        if self.data.attrs['dtype'] == 'GPS':
            GPS0 = 694656019
            self.TAI = spacepy.datamodel.dmarray(
                self.data + GPS0, attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'MJD':
            MJDofTAI0 = 36204.5 # Days-since-1958 is relative to noon
            self.TAI = spacepy.datamodel.dmarray(
                _days1958totai(
                    np.require(self.data, dtype=np.float64)
                    - MJDofTAI0, leaps='rubber', midnight=False),
                attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'JD':
            JDofTAI0 = 2436205.0 # Days-since-1958 is relative to noon
            self.TAI = spacepy.datamodel.dmarray(
                _days1958totai(
                    np.require(self.data, dtype=np.float64)
                    - JDofTAI0, leaps='rubber', midnight=False),
                attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'RDT':
            RDTofTAI0 = 714780. #RDT date of 1958-1-1T00
            tai = _days1958totai(
                np.require(self.data, dtype=np.float64)
                - RDTofTAI0, leaps='drop', midnight=True)
            # Anything before 1582-10-5 has TAI ten days later than the
            # naive conversion, because RDT has ten days that are not in TAI.
            tai[tai < -11840601600.0] += 864000
            self.TAI = spacepy.datamodel.dmarray(tai, attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'CDF':
            CDFofTAI0 = 61788528000000.0
            # Naive TAI conversion
            tai = (self.data - CDFofTAI0) / 1.e3
            tai = _tai_naive_to_real(tai)
            self.TAI = spacepy.datamodel.dmarray(tai, attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'UNX':
            UNXofTAI0 = -378691200.
            # Naive TAI
            tai = self.data - UNXofTAI0
            tai = _tai_naive_to_real(tai)
            self.TAI = spacepy.datamodel.dmarray(tai, attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'APT':
            # AstroPy time doesn't support TAI directly, have to go relative
            # to GPS. But since it's a simple offset, just do the math here
            # instead of via the GPS attribute.
            GPS0 = 694656019
            self.TAI = spacepy.datamodel.dmarray(
                self.data.gps + GPS0, attrs={'dtype': 'TAI'})
            return self.TAI
        if self.data.attrs['dtype'] == 'ISO':
            isoout, UTC, offset = dtstr2iso(self.data, self._isofmt)
            if 'UTC' not in dir(self):
                self.UTC = spacepy.datamodel.dmarray(
                    UTC, attrs={'dtype': 'UTC'})
            if 'ISO' not in dir(self):
                self.ISO = spacepy.datamodel.dmarray(
                    isoout, attrs={'dtype': 'ISO'})
        else:
            UTC = self.UTC
            offset = None

        TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)

        leapsec = self.getleapsecs()
        TAItup = [utc - TAI0 + datetime.timedelta(seconds=int(ls)) for utc, ls in zip(UTC, leapsec)]
        if offset is not None:
            TAI = [tai.days * 86400 + tai.seconds
                   + (tai.microseconds + offset[i])/ 1.e6
                   for i, tai in enumerate(TAItup)]
        else:
            TAI = [tai.days * 86400 + tai.seconds + tai.microseconds / 1.e6
                   for tai in TAItup]

        TAI = spacepy.datamodel.dmarray(TAI, attrs={'dtype': 'TAI'})
        # 1582-10-5 through 1582-10-14 do not exist, so anything
        # before 1582-10-15 is 10 TAI days later than the naive conversion.
        TAI[TAI < -11840601600.0] += (86400 * 10)

        self.TAI = TAI
        return self.TAI

    # -----------------------------------------------
    def getISO(self):
        """
        a.ISO or a.getISO()

        convert dtype data into ISO string

        Always recalculates from the current value of ``UTC``, which
        will be created if necessary. Applies leapsecond correction
        based on ``TAI``, also created as necessary.

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
        getUTC, getUNX, getRDT, getJD, getMJD, getCDF, getTAI, getDOY, geteDOY,
        getAPT
        """
        if self.data.attrs['dtype'] == 'ISO': # Convert the string directly.
            self.ISO = spacepy.datamodel.dmarray(
                dtstr2iso(self.data, fmt=self._isofmt)[0],
                attrs={'dtype': 'ISO'})
            return self.ISO
        nTAI = len(self.data)
        self.TAI = self.getTAI()
        try:
            iso = [utc.strftime(self._isofmt) for utc in self.UTC]
        except ValueError: # Python before 3.3 fails on strftime before 1900.
            iso = [utc.replace(year=1900).strftime(
                self._isofmt.replace('%Y', str(utc.year))) for utc in self.UTC]
        self.ISO = spacepy.datamodel.dmarray(iso, attrs={'dtype': 'ISO'})
        for i in range(nTAI):
            if int(self.TAI[i]) in TAIleaps:
                # UTC is 23:59:59.9999, get correct number of microseconds
                tmpdt = self.UTC[i].replace(
                    microsecond=int((self.TAI[i] % 1) * 1e6))
                # And fudge the second
                a, b, c = tmpdt.strftime(self._isofmt).split(':')
                cnew = c.replace('59', '60')
                self.ISO[i] = a + ':' + b + ':' + cnew

        return self.ISO

    # -----------------------------------------------
    def getleapsecs(self):
        """
        a.leaps or a.getleapsecs()

        retrieve leapseconds from lookup table, used in getTAI

        Always recalculates from current value of ``TAI`` if ``data``
        is dtype ``TAI``, otherwise from the current value of ``UTC``,
        which will be created if necessary.

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
        ========
        getTAI

        """
        if self.data.attrs['dtype'] == 'TAI':
            # TAIleaps contains the TAI which IS the leap second.
            # The leap second count increments in the NEXT second.
            # (TAI - UTC changes at end of leap second.)
            # So find the latest index where the TAI-after-leap-second
            # is less than current TAI (i.e. we are not after leap second yet).
            idx = np.searchsorted(TAIleaps + 1, self.data, side='right') - 1
            return secs[idx]
        tup = self.UTC

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
        Create Ticktock with the current UTC time.

        Equivalent to datetime.utcnow()

        .. versionchanged:: 0.2.2
            This now returns a UTC time; previously it returned a Ticktock
            UTC object, but in the local timezone, which made all conversions
            incorrect.

        Returns
        =======
        out : ticktock
            Ticktock object with the current time, equivalent to datetime.utcnow()

        See Also
        ========
        datetime.datetime.now, datetime.datetime.utcnow

        """
        try:
            dt = datetime.datetime.now(datetime.UTC).replace(tzinfo=None)
        except AttributeError:
            dt = datetime.datetime.utcnow()
        return Ticktock(dt, 'UTC')

    # -----------------------------------------------
    @classmethod
    def today(cls):
        """
        Create Ticktock with the current UTC date and time set to 00:00:00

        Similar to date.today() with time included but in UTC and with the
        time included (zero hours, minutes, seconds)

        .. versionchanged:: 0.2.2
            This now returns the UTC day; previously it returned a Ticktock
            UTC object, but in the local timezone, which made all conversions
            incorrect.

        Returns
        =======
            out : ticktock
                Ticktock object with the current UTC day

        See Also
        ========
        datetime.date.today

        """
        warnings.warn('today() returns UTC day as of 0.2.2.',
                      DeprecationWarning)
        try:
            dt = datetime.datetime.now(datetime.UTC).replace(tzinfo=None)
        except AttributeError:
            dt = datetime.datetime.utcnow()
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


def dtstr2iso(dtstr, fmt='%Y-%m-%dT%H:%M:%S'):
    """Convert a datetime string to a standard format

    Attempts to maintain leap second representation while converting
    time strings to the specified format (by default, ISO8601-like.)
    Only handles a single positive leap second; negative leap seconds
    require no special handling and policy is for UTC-UT1 not to
    exceed 0.9.

    Parameters
    ==========
    dtstr : sequence of str
        Date + time representation, format is fairly open.

    Returns
    =======
    isostr : array of str
        Representation of `dtstr` formatted according to `fmt`.
        Always a new sequence even if contents are identical to `dtstr`.
    UTC : array of datetime.datetime
        The closest-possible rendering of UTC time before or equal to `dtstr`.
    offset : array of int
        Amount (in microseconds) to add to `UTC` to get the real time.

    Other Parameters
    ================
    fmt : str, optional
        Format appropriate for :meth:`~datetime.datetime.strftime` for
        rendering the output time.
    """
    indtstr = dtstr
    # Will be editing this, so force it to own its data (but keep copy!)
    dtstr = np.require(indtstr, dtype='U', requirements='O')
    dtstr = dtstr.copy() if dtstr is indtstr else dtstr
    # Replace leapsecond with a valid "59"
    # Indices of every place that might be leap second
    # Leap second is sec==60, must come at LEAST after YYMMDDHHMM, if
    # being very terse, but this still will avoid YYYY-060
    lsidx = np.char.rfind(dtstr, '60') >= 10
    if lsidx.shape == ():
        lsidx = lsidx.reshape((1,)) #atleast1d is only in new numpy
    leapidx = np.transpose(np.nonzero(lsidx))
    # Add offset to datetime value to get actual UTC,
    # in integer microseconds.
    offset = np.zeros(shape=dtstr.shape, dtype=np.uint32)
    # Actual leap second indices (not just suspected)
    realleap = []
    for j, idx in enumerate(leapidx): # Should be short list, so loop it
        # Get the index to scalar if necessary
        i = tuple(idx) if dtstr.shape else ()
        # The leap second must be preceded by at least 10 char (above),
        # followed by either nothing or fractional seconds,
        # preceded by 59 minutes (and optional separator)
        # It must also not be preceded by a . (to avoid replacing
        # 60 milliseconds with 59 milliseconds).
        dtstr[i], count = re.subn(
            r'^([^\.]{8,}59[^\d]?)60((?:\.\d+)?)$', r'\g<1>59\g<2>', dtstr[i])
        # Doing this subtracted one second.
        if count:
            realleap.append(idx)
            offset[i] = 1e6
    # Cut index of leap seconds down to real ones.
    leapidx = np.array(realleap)
    # try a few special cases that are faster than dateutil.parser
    strfmts = ['%Y-%m-%dT%H:%M:%S',
               '%Y-%m-%dT%H:%M:%SZ',
               '%Y-%m-%d',
               '%Y%m%d',
               '%Y%m%d %H:%M:%S']
    if fmt not in strfmts:
        strfmts.insert(0, fmt)
    for strfmt in strfmts:
        try:
            UTC = np.frompyfunc(
                lambda x: datetime.datetime.strptime(x, strfmt), 1, 1)(dtstr)
            break
        except ValueError:
            continue
    else:
        UTC = np.frompyfunc(dup.parse, 1, 1)(dtstr)
    isostr = np.vectorize(lambda x: x.strftime(fmt), otypes=['U'])(UTC)
    # Check that leap seconds are actually valid
    if len(leapidx):
        # Day that ends in leap second *entry* (may not be leap second)
        # Unbelievably these are read as floats, and other places rely on that.
        leapsecday = np.array([
            datetime.date(int(y), int(m), int(d))- datetime.timedelta(days=1)
            for y, m , d in zip(year, mon, day)])
        # Find only those leap seconds that are really changes
        idx = np.nonzero(np.diff(np.concatenate(([0], secs))))[0]
        leapsecday = leapsecday[idx]
        if dtstr.shape == ():
            if UTC.date() not in leapsecday:
                raise ValueError('{} is not a valid leapsecond.'.format(
                    indtstr))
        else:
            utcday = np.frompyfunc(lambda x: x.date(), 1, 1)(UTC[leapidx])
            # Last day w/leapsecond that comes before supposed leapsec day
            closestday = np.clip(
                np.searchsorted(leapsecday, utcday), None, len(leapsecday) - 1)
            bad = utcday != leapsecday[closestday]
            if bad.any():
                raise ValueError('{} is not a valid leapsecond.'.format(
                    indtstr[leapidx[bad][0]]))
    #Peg all leap seconds to the end of the previous second
    for idx in leapidx:
        i = tuple(idx) if dtstr.shape else ()
        if not offset[i]:
            continue
        # The time string is off by one second.
        if i == (): # Scalar input.
            isostr = UTC.strftime(fmt.replace('%S', '60'))
            us = UTC.microsecond
            UTC = UTC + datetime.timedelta(microseconds=(1e6 - us - 1))
            offset = us + 1
        else:
            isostr[i] = UTC[i].strftime(fmt.replace('%S', '60'))
            us = UTC[i].microsecond
            UTC[i] = UTC[i] + datetime.timedelta(microseconds=(1e6 - us - 1))
            offset[i] = us + 1
    return isostr, UTC, offset


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
    Return a (or many) random datetimes between two given dates

    Convention used is dt1 <= rand < dt2. Leap second times will
    not be returned.

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
    if dt1.tzinfo != dt2.tzinfo:
        raise ValueError('tzinfo for the input and output datetimes must match')
    tt = Ticktock([dt1, dt2]).RDT
    rnd_tn = np.random.uniform(tt[0], tt[1], size=N)
    rnd_t = Ticktock(rnd_tn, dtype='RDT')
    if sorted:
        rnd_t.sort()
    rnd_t = np.asarray([val.replace(tzinfo=dt1.tzinfo if tzinfo else None)
                        for val in rnd_t.UTC])
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


def _leapsgood(now, filetime, lastleap):
    """Determines if leap seconds are up-to-date

    Parameters
    ----------
    now : datetime.datetime
        Current time (UTC assumed)
    filetime : datetime.datetime
        Timestamp (mtime) of the leapsecond file
    lastleap : datetime.datetime
        Last leapsecond in the file (really the moment after the leap)

    Returns
    -------
    bool
        True if leap second file is up-to-date, False if might not be.
    """
    # There cannot be a leap second until AFTER the ls after the
    # bulletin (e.g. once past 1 Jan, know the next possible is next 1 Jan)
    goodto_mtime = datetime.datetime(
        filetime.year + 1, 7 if filetime.month > 6 else 1, 1)
    # The next possible leapsecond is 6mo after the previous known leap,
    # which is (technically just before) 1 Jan or 1 July
    goodto_ls = datetime.datetime(
        lastleap.year + int(lastleap.month > 6),
        1 if lastleap.month > 6 else 7, 1)
    goodto = max(goodto_ls, goodto_mtime)
    return now < goodto


def _read_leaps(oldstyle=False):
    """Read leapseconds in from spacepy tai-utc file

    Populates module-global variables with leapsecond information:
    secs, year, mon, day, TAIleaps.
    Called on import of this module.

    Other Parameters
    ----------------
    oldstyle : bool

        .. versionadded:: 0.2.3

        Treat leapseconds as in SpacePy 0.2.2 (default False). Default is
        to ignore the file contents through 1 Jan 1972 and use SpacePy's own
        list of integral leapseconds, which uses the post-1972 standard
        (add a leapsecond on January 1/July 1 if UTC - UT1 > 0.4s).
    """
    global secs, year, mon, day, TAIleaps
    # load current file
    fname = os.path.join(spacepy.DOT_FLN, 'data', 'tai-utc.dat')
    try:
        with open(fname) as fh:
            text = fh.readlines()
        mtime = datetime.datetime(*time.gmtime(os.path.getmtime(fname))[:6])
    except IOError:
        warnings.warn('Cannot read leapsecond file. Use'
                      ' spacepy.toolbox.update(leapsecs=True).')
        text = [] # Use built-in pre-1972 leaps
        mtime = None
    # Some files have a "last checked" line at the top
    if text and text[0].startswith('Checked'):
        del text[0]

    months = np.array(['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                       'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'])

    if not oldstyle: # Use internal integral leapsecond count until 1972
        keep_idx = next((i for i, l in enumerate(text)
                         if int(l[:5]) > 1972
                         or int(l[:5]) == 1972 and l[6:9] == 'JUL'), len(text))
        leaps = [(1959, 1, 36569), (1961, 1, 37300), (1963, 7, 38211),
                 (1965, 1, 38761), (1966, 7, 39307), (1967, 7, 39672),
                 (1968, 7, 40038), (1969, 7, 40403), (1970, 7, 40768),
                 (1971, 7, 41133)]
        # 1972 Jan NOT a leapsecond in this formulation; TAI-UTC=10s at 71 Jul
        text = [
            ' {Y} {M}  1 =JD {JD}  TAI-UTC={L:12.7f} S '
            '+ (MJD - 00000.) X 0.0000000 S   ACTUAL\n'.format(
                Y=l[0], M=months[l[1] - 1], JD=l[2] + 2400000.5, L=i + 1)
            for i, l in enumerate(leaps)] \
                + text[keep_idx:]
    secs = np.zeros(len(text))
    year = np.zeros(len(text))
    mon = np.zeros(len(text))
    day = np.zeros(len(text))

    for line, i in zip(text, np.arange(len(secs))):
        # Round float seconds (0.5 always rounds up.)
        secs[i] = int(float(line.split()[6]) + 0.5)
        year[i] = int(line.split()[0])
        mon[i] = int(np.where(months == line.split()[1])[0][0] + 1)
        day[i] = int(line.split()[2])

    # Check for out of date. The leap second bulletin comes every
    # six months, and that contains information through the following
    # leap second (end of June/Dec)
    try:
        utcnow = datetime.datetime.now(datetime.UTC).replace(tzinfo=None)
    except AttributeError:
        utcnow = datetime.datetime.utcnow()
    if mtime is not None and spacepy.config['enable_old_data_warning'] \
       and not _leapsgood(
           utcnow, mtime,
           datetime.datetime(int(year[-1]), int(mon[-1]), int(day[-1]))):
            warnings.warn('Leapseconds may be out of date.'
                          ' Use spacepy.toolbox.update(leapsecs=True)')

    TAIleaps = np.zeros(len(secs))
    TAItup = [''] * len(secs)
    TAI0 = datetime.datetime(1958, 1, 1, 0, 0, 0, 0)
    for i in np.arange(len(secs)):
        TAItup[i] = datetime.datetime(int(year[i]), int(mon[i]), int(day[i])) - TAI0 + datetime.timedelta(
            seconds=int(secs[i]) - 1)
        TAIleaps[i] = TAItup[i].days * 86400 + TAItup[i].seconds + TAItup[i].microseconds / 1.e6


def _days1958(tai, leaps='rubber', midnight=False):
    """Calculate days and fractional days since 1958-01-01T12:00

    This is basically a Julian Date but baselined from the start of TAI.
    Since it is calculated from TAI, using this as day 0 maximizes the
    resolution of the resulting values.

    Parameters
    ==========
    tai : sequence of float
        TAI seconds (i.e. continuous SI seconds relative to 1958-01-01T00:00)

    leaps : str, optional
        How to treat days with leapseconds. Since the Julian Date runs
        noon (inclusive) to noon (exclusive), this affects the back
        half of the date with the leapsecond, and the first half of the
        next. Accepted values are:

        rubber
            Consider these days to be 86401 seconds long and thus treat
            each second as a slightly smaller fraction of the day,
            so that all times are represented and evently spaced within
            the day. I.e. the day is "stretched" across more seconds, by
            analogy with the "rubber second" of 1960-1972. (default)

        drop
            Treate time as a flow of SI seconds that is suspended during
            leapseconds, i.e. leapsecond time does not exist. Time during
            leap seconds is pinned to the last microsecond of the
            previous second.

        continuous
            Treat time as a continuous flow of SI seconds with days always
            86400 seconds long. This is the proper handling for true
            Julian Dates, i.e. JD(TAI), as recommended by IAU General Assembly
            XXIII, resolution B1.

    Returns
    =======
    sequence of float
        Days, including fraction, relative to 1958-01-01T12:00

    Other parameters
    ================
    midnight : bool, optional
        Start the day at midnight instead of noon. This affects the allocation
        of leap seconds, and of course days are relative to 1958-01-01T00:00
    """
    off = 0. if midnight else 43200. # Offset from midnight
    # Shift to time-since-noon, if desired (also makes copy), call delta-TAI
    dtai = np.require(tai, dtype=np.float64) - off
    leap_dtai, taiutc = _changed_leaps()
    leap_dtai -= off # delta-TAI of start of leap second
    if leaps in ('rubber', 'drop'):
        # Index of leapsecond equal to or before each time record
        lidx = np.searchsorted(leap_dtai, dtai, side='right') - 1
    elif leaps != 'continuous':
        raise ValueError('leaps handling {} not recognized.'.format(leaps))
    if leaps == 'rubber':
        # d-TAI of start-of-day (noon or midnight) before/after each leapsecond
        # because leapseconds always start at 23:59:60 (43200s after noon).
        if midnight:
            daystart_before_leap = leap_dtai - 86400
            daystart_after_leap = leap_dtai + 1
        else:
            daystart_before_leap = leap_dtai - 43200
            daystart_after_leap = leap_dtai + 43201
        # Closest leapsecond-day-start before record.
        # May be greater than lidx, if near but before leapsecond
        ldidx = np.searchsorted(daystart_before_leap, dtai, side='right') - 1
        # All records that happen on leap second days.
        leap_sec_day = (daystart_before_leap[ldidx] <= dtai) \
                       & (dtai < daystart_after_leap[ldidx])
        # Save SSD on leap second days, about to be destroyed
        if leap_sec_day.shape == (): # Scalar
            ssd_leap_sec_day = dtai - daystart_before_leap[ldidx]
        else:
            ssd_leap_sec_day = dtai[leap_sec_day] \
                               - daystart_before_leap[ldidx[leap_sec_day]]
        # Destroy leap seconds, so days always break at 86400s
        dtai -= taiutc[lidx]
    elif leaps == 'drop':
        # Where in a leapsecond, pin to previous second
        inleap = dtai - leap_dtai[lidx] < 1
        if inleap.shape == (): # Scalar
            if inleap:
                dtai = np.floor(dtai) - .000001
                lidx -= 1 # Associated with previous TAI - UTC
        else:
            dtai[inleap] = np.floor(dtai[inleap]) - .000001
            # Already subtracted ~1sec, now associated with previous TAI - UTC
            lidx[inleap] -= 1
        # Remove all of the TAI seconds that "disappeared".
        dtai -= taiutc[lidx]
    day = np.floor(dtai / 86400)
    ssd = dtai - day * 86400 # Mod does wrong thing if negative.
    if leaps == 'rubber':
        # Patch back in the SSD on leap-second days.
        if leap_sec_day.shape == (): # Scalar
            if leap_sec_day:
                ssd = ssd_leap_sec_day
        else:
            ssd[leap_sec_day] = ssd_leap_sec_day
        daylen = np.choose(leap_sec_day, (86400., 86401.))
    else:
        daylen = 86400
    return day + ssd / daylen


def _days1958totai(days, leaps='rubber', midnight=False):
    """Calculate TAI from days and fractional days since 1958-01-01T12:00

    Input is basically a Julian Date but baselined from the start of TAI.
    Since it is calculated from TAI, using this as day 0 maximizes the
    resolution of the resulting values.

    Parameters
    ==========
    days : sequence of float
        Days, including fraction, relative to 1958-01-01T12:00

    leaps : str, optional
        How to treat days with leapseconds. Since the Julian Date runs
        noon (inclusive) to noon (exclusive), this affects the back
        half of the date with the leapsecond, and the first half of the
        next. Accepted values are:

        rubber
            Consider these days to be 86401 seconds long and thus treat
            each second as a slightly smaller fraction of the day,
            so that all times are represented and evently spaced within
            the day. I.e. the day is "stretched" across more seconds, by
            analogy with the "rubber second" of 1960-1972. (default)

        drop
            Treate time as a flow of SI seconds that is suspended during
            leapseconds, i.e. leapsecond time does not exist. Time during
            leap seconds is not recovered (no output time during leaps.)

        continuous
            Treat time as a continuous flow of SI seconds with days always
            86400 seconds long. This is the proper handling for true
            Julian Dates, i.e. JD(TAI), as recommended by IAU General Assembly
            XXIII, resolution B1.

    Returns
    =======
    sequence of float
        TAI seconds (i.e. continuous SI seconds relative to 1958-01-01T00:00)

    Other parameters
    ================
    midnight : bool, optional
        Start the day at midnight instead of noon. This affects the allocation
        of leap seconds, and of course `days` are relative to 1958-01-01T00:00
    """
    days = np.require(days, dtype=np.float64)
    off = 0. if midnight else 43200. # Offset from midnight
    if leaps == 'continuous': # Very simple case.
        return days * 86400 + off
    elif leaps not in ('rubber', 'drop'):
        raise ValueError('leaps handling {} not recognized.'.format(leaps))
    leap_tai, taiutc = _changed_leaps()
    leap_dtai = leap_tai - off # delta-TAI of start of leap second
    # Days with leap second. Leap second is always late in day: end
    # of day if doing midnight-to-midnight, halfway if noon-to-noon.
    leap_day = np.floor((leap_dtai - taiutc) / 86400)
    # Closest leapsecond-day before record.
    ldidx = np.searchsorted(leap_day, days) - 1
    # All records that happen on leap second days.
    leap_sec_day = (ldidx > 0) \
                   & (leap_day[ldidx] <= days) \
                   & (days < leap_day[ldidx] + 1)
    if leaps == 'rubber':
        daylen = np.choose(leap_sec_day, (86400., 86401.))
    else:
        daylen = 86400.
    # Naive TAI at start of day (no leapsecond corrections).
    dayint = np.floor(days)
    daystart = dayint * 86400
    # Seconds from the start of naive day.
    ssd = (days - dayint) * daylen # Mod does wrong thing if negative
    # Index TAIUTC for every input. Start of day w/LS is still previous value.
    taiutcidx = ldidx - np.require(leap_sec_day, dtype=np.intp)
    # Keep small number adds together: ssd and leapseconds.
    taiout = daystart + (ssd + off + taiutc[taiutcidx])
    if leaps == 'drop':
        # Any times after leap-second skip need an extra second.
        taiutcidx += np.require(leap_sec_day & (taiout >= leap_tai[ldidx]),
                                dtype=np.intp)
        # Recalculate to keep small-number addition together
        taiout = daystart + (ssd + off + taiutc[taiutcidx])
    return taiout


def _changed_leaps():
    """Find only those times where leap seconds actually changed

    There are leap second records where there is no actual change
    in leap seconds, and there isn't a record for TAI-UTC == 0; this
    adds a record (where the TAI of the leap second time is -inf) and
    eliminated those with no actual change.

    Returns
    =======
    leap_tai : sequence of float
        TAI of the start of every leap second.
    taiutc : sequence of float
        TAI - UTC at the end of the leap second.
    """
    # Find only those leap seconds that are really changes
    idx = np.nonzero(np.diff(np.concatenate(([0], secs))))[0]
    leap_tai = TAIleaps[idx] # TAI of start of leap second
    taiutc = secs[idx] # TAI - UTC after leap second (same index as leap_tai)
    # Add fake TAI - UTC = 0 at the Big Bang
    leap_tai = np.concatenate(([-np.inf], leap_tai))
    taiutc = np.concatenate(([0], taiutc))
    return leap_tai, taiutc


def _tai_naive_to_real(tai):
    """Convert naive TAI to real TAI

    Convert a TAI on a continuous timescale that skips over leapseconds
    and is unaware of calendar conversion to actual TAI, which includes
    leapseconds and is continuous across the Gregorian conversion (naive
    has a gap, i.e. dates it can represent that don't exist.)

    Parameters
    ==========
    tai : sequence of float
        Naive TAI

    Returns
    =======
    tai : sequence of float
        TAI
    """
    # This is the ACTUAL TAI and TAI-UTC at the end of that TAI
    leap_tai, taiutc = _changed_leaps()
    # Naive TAI and index that corresponds to TAI - UTC
    naive_leap_tai = leap_tai - taiutc + 1
    taiutcidx = np.searchsorted(naive_leap_tai, tai, side='right') - 1
    realtai = tai + taiutc[taiutcidx]
    # Anything before 1582-10-5 has TAI ten days later than the
    # naive conversion, because naive has ten days that are not in TAI.
    realtai[realtai < -11840601600.0] += 864000
    return realtai


def _tai_real_to_naive(tai):
    """Convert naive TAI to actual TAI

    Convert actual TAI, which includes leapseconds and is continuous
    across the Gregorian calendar conversion, to naive TAI, which skips
    over leapseconds and has times in the calendar conversion that
    do not actually exist.

    Parameters
    ==========
    tai : sequence of float
        TAI

    Returns
    =======
    tai : sequence of float
        Naive TAI
    """
    # ACTUAL TAI and TAI-UTC at the end of that TAI
    leap_tai, taiutc = _changed_leaps()
    # Points to largest leap-TAI less-than input TAI, thus also TAI-UTC
    lidx = np.searchsorted(leap_tai, tai, side='right') - 1
    # Records in a leap second
    inleap = tai < leap_tai[lidx] + np.diff(taiutc)[lidx - 1]
    naive_tai = np.choose(inleap, (
        tai - taiutc[lidx], # Just subtract off LS
        np.floor(tai) + (.999 - taiutc[lidx]) # Peg to end of sec
        ))
    # Anything before 1582-10-15 needs to skip 10 CDF days backward,
    # since naive has ten days that are not in TAI
    naive_tai[tai < -11839737600.0] -= 864000
    return naive_tai


_read_leaps()
