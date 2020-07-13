=============
Release Notes
=============

This document presents the most notable user-level changes in each
release of SpacePy. The CHANGELOG file in the source distribution
contains more detail.

.. contents::
   :depth: 2
   :local:

0.2 Series
==========

0.2.2 (2020-xx-xx)
------------------

This will be the last release with full support for :doc:`Python 2 <py2k_eol>`.

New features
************
:mod:`~spacepy.irbempy` incorporates upstream IRBEMlib rev620. This
adds IGRF13 coefficients. :mod:`~spacepy.coordinates` and
:mod:`~spacepy.irbempy` now also support using all supported
coordinate systems as inputs to routines; if a routine does not
support an input system, it will be automatically converted.

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
0.3.0.

Quaternion math functions have been moved to
:mod:`~spacepy.coordinates`; using the functions in
:mod:`~spacepy.toolbox` is deprecated.

:func:`~spacepy.pybats.rim.fix_format` is now deprecated, as
:class:`~spacepy.pybats.rim.Iono` can now read these files directly.

Dependency requirements
***********************
numpy 1.10 is now required. (Many functions erroneously required it from 0.2.1, but this was not adequately documented.)

scipy 0.11 is now required. (Again, this was erroneously required in 0.2.0 without appropriate documentation.)

Several dependencies without an established minimum version were tested.

As of 0.2.2, the minimum versions of dependencies are:
  * CPython 2 2.7 or CPython 3 3.2
  * CDF 2.7
  * dateutil 1.4 (earlier may work)
  * ffnet 0.7 (earlier may work)
  * h5py 2.6 (earlier may work)
  * matplotlib 1.5
  * networkx 1.0 (earlier may work)
  * numpy 1.10
  * scipy 0.11

See :doc:`dependencies` for full details.

Major bugfixes
**************
None of note (but many minor ones).

Other changes
*************
Data sources for leapsecond files and :mod:`~spacepy.omni` Qin-Denton
files have been updated to provide current sources. If present,
entries in the :doc:`configuration file <configuration>` will still be
used instead.

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
None of note (but many minor ones).

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
