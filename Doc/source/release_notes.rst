=============
Release Notes
=============

This document presents user-visible changes in each release of SpacePy.

.. contents::
   :depth: 2
   :local:


0.3 Series
==========
0.3.0 (2021-xx-xx)
------------------

Deprecations and removals
*************************
Colourmaps have been removed from :class:`~spacepy.plot`. The same
colourmaps (``plasma`` and ``viridis``) have been available in
matplotlib since at least 1.5. (Deprecated in 0.2.3.)

The old name ``spectrogram`` for :class:`~spacepy.plot.Spectrogram`
has been removed. (Deprecated in 0.2.2.)

The ``read_ram_dst`` function has been removed from
:mod:`~spacepy.pybats.ram`, as it operates on files that are no longer
written by RAM-SCB. (Deprecated in 0.1.6.)

The ``fix_format`` function has been removed from
:mod:`~spacepy.pybats.rim`; :class:`~spacepy.pybats.rim.Iono` can now
read these files directly. (Deprecated in 0.2.2.)

The ``from_dict`` method of CDF attribute lists
(:meth:`~spacepy.pycdf.gAttrList`, :meth:`~spacepy.pycdf.zAttrList`)
has been removed. Use :meth:`~spacepy.pycdf.AttrList.clone`, which
supports cloning from dictionaries. (Deprecated in 0.1.5.)

The ``feq`` function has been removed from :mod:`~spacepy.toolbox`;
use :func:`numpy.isclose`. (Deprecated in 0.2.2.)

Quaternion math functions have been removed from
:mod:`~spacepy.toolbox`; they are available in
:mod:`~spacepy.coordinates`. (Deprecated in 0.2.2.)

Other changes
*************
:mod:`~spacepy.pycdf` now defaults to creating version 3 (not
backward-compatible) CDFs if the backward compatible mode is not
explicitly set (:meth:`~spacepy.pycdf.Library.set_backward`). It still
issues a warning when creating a CDF if this is not set; this warning
will be removed in 0.4.0. (Warning added in 0.2.2.)

Similarly, :mod:`~spacepy.pycdf` defaults to TIME_TT2000 when creating
a time variable or attribute without specifying a type (EPOCH or
EPOCH16 are used if TT2000 isn't available). A warning is issued when
doing so; this warning will be removed in 0.4.0. (Warning added in 0.2.2.)

0.2 Series
==========

0.2.3 (2021-10-30)
------------------
This is the last release of the 0.2 series and the last with full
support for :doc:`Python 2 <py2k_eol>`. Binary installers (including
wheels) for :doc:`32-bit Windows <install_windows>` will also end
after the 0.2 series, as will Windows installers. The only binaries
for Windows will be 64-bit wheels, installable with ``pip``.

New features
************
:mod:`~spacepy.pycdf` now supports variables with sparse records, including
enabling/disabling sparse records (:meth:`~spacepy.pycdf.Var.sparse`) and
setting the pad value (:meth:`~spacepy.pycdf.Var.pad`). Thanks Antoine Brunet.

Deprecations and removals
*************************
The colourmaps provided in the :mod:`~spacepy.plot` module have been
deprecated. The same colourmaps have been available in matplotlib since
at least 1.5, and users who do not directly import the colourmaps should
see no impact.

Major bugfixes
**************
The passing of keyword arguments from :func:`~spacepy.toolbox.bootHisto`
to :func:`numpy.histogram` and :func:`matplotlib.pyplot.bar` has been fixed.

The check for out-of-date leapseconds in :mod:`~spacepy.time` has been
fixed (previously warned even when the file was up to date.)

Fixed installation on new versions of setuptools, which removed
``bdist_wininst`` support (`#530
<https://github.com/spacepy/spacepy/issues/530>`_).

The handling of library paths on Windows has been updated. This should
fix situations where :mod:`~spacepy.irbempy` would not import on
Windows with Python 3.8 or later. This did not seem to be a problem
with Anaconda, but would sometimes manifest with Python from the app
store or from `<http://python.org/>`_ (`#507
<https://github.com/spacepy/spacepy/issues/507>`_)

Other changes
*************
Modern leapsecond rules are applied from 1958-1972 rather than
rounding fractional leapseconds. See :mod:`~spacepy.time` for full
discussion of leap seconds and other conversion considerations.

The handling of the ``.spacepy`` directory (see :doc:`configuration`)
has been improved. If the ``SPACEPY`` environment variable is used,
the directory will be created. The import process also is less fragile
in the case of a partially-created ``.spacepy`` directory or an
invalid (e.g. empty) ``spacepy.rc``.

0.2.2 (2020-12-29)
------------------

The 0.2 series will be the last with full support for :doc:`Python 2
<py2k_eol>`; 0.2.3 will likely be the last release. Binary installers
for :doc:`32-bit Windows <install_windows>` will also end after the 0.2
series.

New features
************
:mod:`~spacepy.irbempy` incorporates upstream IRBEMlib rev620. This
adds IGRF13 coefficients. :mod:`~spacepy.coordinates` and
:mod:`~spacepy.irbempy` now also support using all supported
coordinate systems as inputs to routines; if a routine does not
support an input system, it will be automatically converted.

:class:`~spacepy.time.Ticktock` supports conversions to and from
:class:`astropy.time.Time`.

The following classes, functions, and methods are new:

.. autosummary::
   ~spacepy.coordinates.quaternionFromMatrix
   ~spacepy.coordinates.quaternionToMatrix
   ~spacepy.datamanager.rebin
   ~spacepy.plot.utils.add_arrows
   ~spacepy.pycdf.concatCDF
   ~spacepy.pycdf.istp.nanfill
   ~spacepy.pycdf.istp.FileChecks.empty_entry
   ~spacepy.pycdf.istp.VarBundle
   ~spacepy.pycdf.istp.VariableChecks.deltas
   ~spacepy.pycdf.istp.VariableChecks.empty_entry

Deprecations and removals
*************************
:mod:`~spacepy.pycdf` now warns if creating a new CDF file without
explicitly setting backward compatible or not backward compatible
(:meth:`~spacepy.pycdf.Library.set_backward`). The default is
still to make backward-compatible CDFs, but this will change in
0.3.0. Similarly it now warns if creating a time variable without
specifying a time type; the default is still to use EPOCH or
EPOCH16, but this will change to TIME_TT2000 in 0.3.0.

:func:`~spacepy.pybats.rim.fix_format` is now deprecated, as
:class:`~spacepy.pybats.rim.Iono` can now read these files directly.

Quaternion math functions have been moved to
:mod:`~spacepy.coordinates`; using the functions in
:mod:`~spacepy.toolbox` is deprecated.

:func:`~spacepy.toolbox.feq` is deprecated; numpy 1.7 added the equivalent
:func:`~numpy.isclose`.

The :class:`~spacepy.plot.spectrogram` class is now capitalized
(:class:`~spacepy.plot.Spectrogram`); the old, lower-case variant is
kept for compatibility but will be removed.

Dependency requirements
***********************
Not all dependencies are required for all functionality; see
:doc:`dependencies` for full details, including what functionality is
lost if a dependency is not installed.

numpy 1.10 is now required. (Many functions erroneously required it from 0.2.1, but this was not adequately documented.)

scipy 0.11 is now the minimum supported version of SciPy. (Again, this was erroneously required in 0.2.0 without appropriate documentation.)

Several dependencies without an established minimum version were tested.

As of 0.2.2, minimum supported versions of dependencies are:
  * CPython 2 2.7 or CPython 3 3.2
  * AstroPy 1.0
  * CDF 2.7
  * dateutil 1.4 (earlier may work)
  * ffnet 0.7 (earlier may work)
  * h5py 2.6 (earlier may work)
  * matplotlib 1.5
  * networkx 1.0 (earlier may work)
  * numpy 1.10
  * scipy 0.11

Major bugfixes
**************
Time conversions between time systems before 1961 now use the proper
number of leapseconds (0).

Many minor bugfixes.

Other changes
*************
Data sources for leapsecond files and :mod:`~spacepy.omni` Qin-Denton
files have been updated to provide current sources. If present,
entries in the :doc:`configuration file <configuration>` will still be
used instead. A (configurable) warning is issued for out-of-date leapsecond
files.

The representation of leap second intervals in time systems which
cannot directly represent them has been changed. Formerly times such
as 2008-12-31T23:59:60 were represented in e.g. UTC datetime as the
the beginning of the next day, e.g. 2009-01-01T00:00:00. They are
now represented by the last possible moment of the same day, e.g.
2008-12-31T23:59:59.999999. Fractional leapsecond counts are now rounded
to the integer instead of truncated; this rounding is applied to the total
TAI - UTC quantity not the individual increments of leap seconds. E.g
successive 0.2, 0.2, 0.2 leap seconds will result in 0, 0, and 1 new
leap seconds.

Similarly, leap seconds are now included in the fractional day
calculation of MJD, so MJD values around a leap second may be different
than in previous versions of SpacePy.

Most time systems are now converted to/from TAI rather than using
datetime. This may cause small differences with previous versions of
SpacePy, on order of a double precision. RDT and JD are particularly
affected for dates in the modern era. Time conversions around
leapseconds may also be different; in many cases they were undefined
in previous versions.

:meth:`~spacepy.time.Ticktock.now` and :meth:`~spacepy.time.Ticktock.today`
return times in UTC; in previous versions the value returned was local,
but was treated as UTC for all conversions (and thus inaccurate.)

See :mod:`~spacepy.time` for full discussion of leap seconds, time
resolution, and other conversion considerations.

0.2.1 (2019-10-02)
------------------

New features
************
The following module is new:

.. autosummary::
   ~spacepy.pycdf.istp

Deprecations and removals
*************************
None

Dependency requirements
***********************
No changes to minimum dependency versions.

As of 0.2.1, the minimum versions of dependencies are:
  * CPython 2 2.7 or CPython 3 3.2
  * CDF 2.7
  * matplotlib 1.5
  * numpy 1.4
  * scipy 0.10

Other dependencies have no established minimum. See
:doc:`dependencies` for full details.

Major bugfixes
**************
Fixed compliation of :mod:`~spacepy.irbempy` on several systems.

Other changes
*************
None of note.

0.2.0 (2019-06-22)
------------------

New features
************

Deprecations and removals
*************************
None

Dependency requirements
***********************
Support for Python 2.6 was removed; 2.7 is the only supported version
of Python 2.

As of 0.2.0, the minimum versions of dependencies are:
  * CPython 2 2.6 or CPython 3 3.2
  * CDF 2.7
  * matplotlib 1.5
  * numpy 1.4
  * scipy 0.10

Other dependencies have no established minimum. See
:doc:`dependencies` for full details.

Major bugfixes
**************
:meth:`~spacepy.toolbox.human_sort` was fixed for non-numeric inputs
(the normal case.) This had been broken since 0.1.6.

Many minor bugfixes as well.

Other changes
*************
Many updates to improve ease of installation, including Windows binary wheels.

0.1 Series
==========
See the CHANGELOG file in the source distribution for changes in the 0.1
release series.

0.1.6 (2016-09-08)
------------------

0.1.5 (2014-12-23)
------------------

0.1.4 (2013-05-21)
------------------

0.1.3 (2012-06-22)
------------------

0.1.2 (2012-05-25)
------------------

0.1.1 (2011-10-31)
------------------

0.1 (2011-08-24)
----------------
