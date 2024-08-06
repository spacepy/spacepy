#################################################################################
time - Time conversion, manipulation and implementation of Ticktock class
#################################################################################

Contents
--------

.. contents::
   :depth: 2
   :local:

.. currentmodule:: spacepy.time

See also the `full API documentation <spacepy.time>`.

Introduction
============
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

Examples
========

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