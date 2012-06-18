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

from time import *

__contact__ = 'Josef Koller, jkoller@lanl.gov'
