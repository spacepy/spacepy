=============
Release Notes
=============

This document presents user-visible changes in each release of SpacePy.

.. contents::
   :depth: 2
   :local:

0.8 Series
==========
0.8.0 (2025-xx-xx)
------------------

New features
************
`.pycdf.CDF.leapsecond_lastupdated` was added to allow reading and
setting the last update of the leapsecond table in a CDF.

Major bugfixes
**************
Several small issues in `~spacepy.ae9ap9` were fixed, relating to
plotting and reading files with non-MJD timestamps.

Other changes
*************
The IGRF provided with SpacePy (used in the default spacepy backend for
:doc:`coordinate <coordinates>` conversion) has been updated to the
newly-released `IGRF14 <https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/>`_.

The `IRBEM <https://github.com/PRBEM/IRBEM>`_ library bundled with SpacePy
has been updated to reflect the new version of IGRF along with recent updates
and bugfixes. This now reflects the upstream repository as of 2024-12-14 (commit
`8279aa9 <https://github.com/PRBEM/IRBEM/commit/8279aa9fd50c0a32d63efe762454b9ebf3b81ef6>`_).
This is used across `spacepy.irbempy` as well as in the (non-default)
irbem backend for coordinate conversion.


0.7 Series
==========
0.7.0 (2024-11-08)
------------------
Binary wheels are now provided for Linux on 64-bit ARM, intended for
Raspberry Pi. Installing :ref:`dependencies from system packages
<linux_dep_apt>` is recommended.

New features
************
`.help` now supports searching the documentation.

Dependency requirements
***********************
Numpy 2.0 is now fully supported. Several deprecations and errors
using 2.0 were fixed.

Numpy and f2py are no longer required to build SpacePy. This change
will not be noticable to most users. Binary wheels are no longer tied
to numpy or Python version. SpacePy will now build on numpy 1.26 and
later on Python 3.12. Numpy is still required to run SpacePy.

Sphinx 4.0 is now required to build the documentation; this is not
a concern for most users.

Support for Python 3.6 has been removed due to inability to
test. Python 3.7 is the oldest supported Python; as a result, Astropy
2.0 is the oldest supported Astropy (if using Astropy).

Deprecations and removals
*************************
`~spacepy.toolbox.timeout_check_call` is deprecated as redundant to using
`subprocess.check_call` with the ``timeout`` argument.

The ``spacepy.irbempy.irbempylib`` module has been removed. This was
the old internal interface to the IRBEM library and was not intended
for public use.

The ``data_assimilation``, ``radbelt``, and ``spacepy_EnKF`` modules
have been moved to the "sandbox". These modules were undertested and
minimally documented. They may be restored in the future with
appropriate testing and docs.

Major bugfixes
**************
`~spacepy.irbempy.find_Bmirror` now correctly returns one ``Bmirr`` per input
pitch angle instead of ignoring all but the first.

Other changes
*************
Operations on `~spacepy.datamodel.dmarray` which return a scalar value
will now return a numpy :std:term:`array scalar` rather than the base
Python type. This is consistent with the behavior of `~numpy.ndarray`.
`~spacepy.datamodel.dmarray` also supports assigning to its ``dtype``
and ``shape``. Together these changes should make ``dmarray`` a much
closer drop-in replacement for `~numpy.ndarray`; in particular, these
address known issues with conversion to masked arrays and passing to
`~matplotlib.pyplot.pcolormesh`.

Warnings issued by SpacePy are now associated with the line of the
calling code, not with the SpacePy code itself. The
``enable_deprecation_warning`` :doc:`configuration option <configuration>`
has been removed and SpacePy does not force display of deprecation
warnings, as warnings issued by SpacePy cannot be distinguished from
other warnings. It is recommended to occasionally run Python code with
the :option:`-Wd <-W>` option to enable deprecation warnings. See
:py:mod:`warnings` for more details.


0.6 Series
==========
0.6.0 (2024-04-25)
------------------
There are no changes to dependencies, minimum versions, or
installation process with this release.

New features
************
`.datamodel` supports output to Pandas
(`~.datamodel.ISTPContainer.toDataFrame`) and input from Pandas
(`~.datamodel.ISTPContainer.fromDataFrame`). This support requires
pandas 0.18.0, which is not installed by default with SpacePy.

`~.datamodel.dmarray` supports output to Astropy Quantities
(`~.datamodel.ISTPArray.toQuantity`) and input from them
(`~.datamodel.ISTPArray.fromQuantity`). This support requires Astropy
1.0, which is not installed by default with SpacePy.

`~.pybats.ImfInput` now has two new methods:
:meth:`~.pybats.ImfInput.calc_clock` which, calculates and stores the IMF
clock angle defined as the angle of the interplanetary magnetic field (IMF) in
the GSM Y-Z plane; and :meth:`~.pybats.ImfInput.calc_epsilon`, which calculates
the epsilon parameter, an approximation of power input into the magnetosphere.

0.5 Series
==========
0.5.0 (2024-03-11)
------------------
This release marks the end of all support for Python 2. SpacePy now
requires Python 3.6 or later. Minimum supported versions for other
dependencies were also increased; see :ref:`release_0_5_0_deps` for details.

The SpacePy team now delivers binary wheels for all supported Python
versions (3.6-3.12) for 64-bit Windows, Linux, and Mac. Most users
will not need to build SpacePy; see :doc:`install` for details.

Due to changes in f2py, SpacePy 0.5.0 will not build on Python 3.12
with numpy 1.26 or later. It will, however, run fine. The SpacePy team
provides binary wheels for Mac and Windows. Users who wish to build
from source on Python 3.12 should first install numpy (which will
likely be built from source)::

  pip install --no-build-isolation "numpy<1.26"
  pip install --no-build-isolation spacepy

.. _release_0_5_0_deps:

Dependency requirements
***********************
Not all dependencies are required for all functionality; see
:doc:`dependencies` for full details, including what functionality is
lost if a dependency is not installed.

The ``setup.py`` based install process is no longer supported; as such,
``pip`` and ``setuptools`` are now required. ``wheel`` is required if
building from source. The vast majority of modern Python distributions
already have these requirements.

Installing from ``pip`` normally installs all necessary dependencies.

The minimum supported version of all dependencies was updated in
SpacePy 0.5.0. Minimum versions are:

  * CPython 3.6
  * AstroPy 1.0
  * CDF 3.5
  * dateutil 2.1
  * h5py 2.10
  * matplotlib 3.1
  * numpy 1.15.1
  * scipy 1.0

The only supported compiler is the GNU gfortran compiler, aka
"gnu95". Support for the older g77 compiler, as well as the Portland
Group and Intel compilers, has been removed.

New features
************
`~.datamodel.readJSONheadedASCII` and `~.datamodel.readJSONMetadata`
now support reading from gzipped input files; filenames ending with
``.gz`` are assumed to be gzipped.

`.datamodel` now has explicit support for operations using
ISTP-compliant metadata; see `~.datamodel.dmarray` and
`~.datamodel.SpaceData` for details. In particular note the plotting
support with `~.datamodel.ISTPContainer.plot`.

`~.irbempy` now provides an interface to the SHIELDOSE-2 shielding
radiation dose model. This is implemented through the
`~.irbempy.Shieldose2` class.

`~.plot.spectrogram.simpleSpectrogram()` supports treating zeros as
valid data on log plots.

:mod:`.pycdf` :meth:`~.pycdf.Library.set_backward` now returns the prior
state of backward compatibility mode.

:func:`~.datamodel.toCDF` now supports writing backward-compatible CDFs.

`~.pybats.ImfInput` has been updated with a more flexible
quick-look plot capability. Variable naming has been unified to remove issues
with temperature and density naming.

`~.pybats.bats.Bats2d` can now calculate parallel and
perpendicular velocities with new `~.pybats.bats.Bats2d.calc_uperp`
and `~.pybats.bats.Bats2d.calc_upar` methods.

`~.pybats.bats.Bats2d` objects and the class that handles quad tree building
(`~.pybats.qotree.QTree`) now accepts a keyword argument to set the size
of each block: ``blocksize``. Default value is 8.

`~.toolbox.dictree` now supports returning the output instead of printing it.

Deprecations and removals
*************************
Since plot styles are no longer applied on import, importing
``spacepy.plot.apionly`` has no effect and is deprecated.

The new pip-based install method does not support a separate ``build``
step so `~spacepy_testing.add_build_to_path` is deprecated. Alternatives
include using an :ref:`editable install <install_--editable>` or installing
in a custom location using :ref:`\\\\\\-\\\\\\-prefix <install_--prefix>` and manually editing :envvar:`PYTHONPATH`.

:mod:`~spacepy.coordinates` no longer warns if a backend (irbem or
SpacePy) is not specified. The default backend is SpacePy, as it has
been since 0.4.0.

In `~spacepy.time`, the ``cls`` argument of
`~spacepy.time.Ticktock.update_items` has been removed (deprecated in 0.2.2).
The end-user normally does not call this function.

In `~spacepy.pybats.bats`, removed ``add_b_magsphere_legacy`` and
``add_b_magsphere_new``, which were deprecated in 0.2.0. Use
:meth:`~spacepy.pybats.bats.Bats2d.add_b_magsphere` (which is the same as
``add_b_magsphere_new``).

Major bugfixes
**************
:func:`~.datamodel.toCDF` handling of time types other than TT2000 has
been fixed.

Other changes
*************
`~.pybats.IdlFile` no longer sorts unstructured data from binary files; see
that documentation for details. Thanks Lutz Rastaetter.

:func:`~.datamodel.toCDF` now only accepts valid keyword arguments.

Plot styles are not automatically applied on import of :mod:`.plot`. Use
`.plot.style` directly to apply the desired style.

`~.pybats.bats.Bats2d` plot functions will no longer raise an exception when
trying to add a planet/inner boundary patch without an 'rbody' attribute
present. Rather, the patch will not be applied to the axes object.

0.4 Series
==========
0.4.1 (2022-09-15)
------------------
This minor release provides no changes in functionality, but fixes
installation problems on some systems. There is no need to upgrade
from a functioning 0.4.0 (and no harm in doing so).

Other changes
*************
Unicode characters were removed from the IRBEM sources, fixing
compilation problems for certain user locale settings.

The version of numpy used for building on Apple Silicon Mac was
updated.

Documentation on troubleshooting ``pip`` problems was improved.


0.4.0 (2022-09-07)
------------------
This release marks the end of support and/or fixes for bugs that cannot
be reproduced on Python 3. As with the previous release series, SpacePy
0.4.0 can still be built and installed "by hand" on Python 2, but no
Python 2 binaries are provided and this version will not install on Python 2
using ``pip``.

New features
************
The :mod:`~spacepy.LANLstar` module has been rewritten to use numpy to
evaluate the neural networks instead of relying on ``ffnet``. The
temporary removal of support for this module in SpacePy 0.3.0 has therefore
been lifted. The new implementation provides a slight performance increase
with no change in results or accuracy.

:class:`~spacepy.pycdf.istp.VarBundle` now supports output to and input from
:class:`~spacepy.datamodel.SpaceData` objects as well as
:class:`~spacepy.pycdf.CDF`.

Both :mod:`~spacepy.coordinates` backends now provide access to the TEME
coordinate system (as used by the SGP4 orbit propagator).

Deprecations and removals
*************************
The ``_nelems`` method of :class:`~spacepy.pycdf.Var` has been removed;
use the public interface :meth:`~spacepy.pycdf.Var.nelems`. (Deprecated
in 0.2.2).

:mod:`~spacepy.irbempy` ``get_sysaxes``, ``sph2car`` and ``car2sph``
were deprecated in SpacePy 0.2.2 and have been removed. In place
of the latter functions, :func:`~spacepy.coordinates.sph2car` and
:func:`~spacepy.coordinates.car2sph` should be used.

Major bugfixes
**************
The installer has been updated to address certain build issues,
particularly on Mac. The Mac :doc:`installation directions
<install_mac>` have been completely rewritten.

:mod:`~spacepy.pycdf` has been updated for Apple Silicon (ARM/M1);
Python 3.8 is required for this support.

:mod:`~spacepy.pycdf` contains a time conversion workaround for
versions of the NASA CDF library before 3.8.0.1. Non-integral epoch
values close to midnight would erroneously return the following day;
:meth:`~spacepy.pycdf.Library.epoch_to_datetime` now returns the
correct value on all CDF library versions.

The IRBEM backend for coordinate transformations has been updated to
correct the specification of transformations through the J2000 and TOD
systems, including correctly setting the GEI and TOD systems to be
equivalent. This may change results by a small amount. The IRBEM update
also traps a singularity at the South pole in the conversion to geodetic
(GDZ) coordinates.

Dependency requirements
***********************
:mod:`~spacepy.LANLstar` now uses a numpy-based implementation (based on
contributions from Aaron Hendry) so neither ``ffnet`` or ``networkx`` are
required to use it. These dependencies were removed in SpacePy 0.3.0, but
were still required for use of ``LANLstar``. Support for ``LANLstar`` is
reinstated in SpacePy 0.4.0.

Other changes
*************
:mod:`~spacepy.pycdf` no longer warns when defaulting to version 3 CDFs
and TIME_TT2000 time type if not specified; the warning was added in
0.2.2 and the default changed in 0.3.0. Use
:meth:`~spacepy.pycdf.Library.set_backward` to create version 2 CDFs and
explicitly specify a time type (e.g. with :meth:`~spacepy.pycdf.CDF.new`)
if TT2000 is not desired.

The IRBEM library bundled with SpacePy has been updated to reflect recent
updates and bugfixes, and reflects the upstream repository as of 2022-08-29
(commit dfb9d26).

0.3 Series
==========
0.3.0 (2022-04-27)
------------------
This release continues the phaseout of Python 2
support. No Python 2 binaries are provided, and 0.3.0 will not install
on Python 2 with ``pip``. Installation via ``setup.py`` from a source
distribution is still available.

This is the last release with Python 2 bugfix support. SpacePy 0.4.0
will make no attempt to maintain functionality for Python 2 and
SpacePy 0.5.0 will not function without Python 3.

Windows binaries are only provided as 64-bit wheels, installable with
``pip``, for Python 3.6 and later. Windows executable installers and
32-bit binaries are no longer provided.


New features
************
The :mod:`~spacepy.coordinates` module has been overhauled with a new,
Python-based backend. This provides comparable performance to the
existing :mod:`~spacepy.irbempy` backend with higher precision and
reduces the dependence on Fortran. By default, irbemlib will still be
built at installation time. The default backend remains IRBEM; in
0.4.0, this will switch to the new :mod:`~spacepy.ctrans` based
backend. The new :mod:`~spacepy.igrf` module is part of this support
but may be of interest on its own.

In accordance with a change from NASA, :mod:`~spacepy.pycdf` now
assumes strings in CDFs are UTF-8. It will no longer raise errors on
reading non-ASCII data from a CDF. See :ref:`pycdf_string_handling` in
the pycdf documentation for details.

:mod:`~spacepy.ae9ap9` now supports the new ephem model file format
(>=1.50.001) via :func:`~spacepy.ae9ap9.parseHeader`. The old file
format is deprecated.

Deprecations and removals
*************************
HTML documentation is no longer installed with
SpacePy. :func:`~spacepy.help` now opens the latest `online
documentation <https://spacepy.github.io/>`_. Offline documentation
are available separately (files named like ``spacepy-x.y.z-doc.zip``
and ``spacepy-x.y.z-doc.pdf``) and as part of the source distribution
(``spacepy-x.y.z.tar.gz`` or ``spacepy-x.y.z.zip``). These files can
be downloaded from SpacePy's `releases on GitHub
<https://github.com/spacepy/spacepy/releases>`_; the source can also
be found on `PyPI <https://pypi.org/project/spacepy/#files>`_.

``LANLstar`` requires `ffnet <http://ffnet.sourceforge.net/>`_, which
does not install properly with current `setuptools
<https://github.com/pypa/setuptools>`_ (version 58).  The SpacePy team
is working on replacing this dependency, but in the meantime
``LANLstar`` is unsupported and will require manually installing
``ffnet`` and `networkx <http://networkx.lanl.gov/>`_.

As mentioned above, :mod:`~spacepy.ae9ap9` support for the old ephem
model file format is deprecated.

Colourmaps have been removed from :class:`~spacepy.plot`. The same
colourmaps (``plasma`` and ``viridis``) have been available in
matplotlib since at least 1.5. (Deprecated in 0.2.3.)

The old name ``spectrogram`` for :class:`~spacepy.plot.spectrogram`
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

Dependency requirements
***********************
Due to the new backend, scipy is now required for
:mod:`~spacepy.coordinates` (even if using the old backend). 0.11
remains the minimum version.

Since ``LANLstar`` is not currently supported, ``ffnet`` and
``networkx`` are no longer treated as SpacePy dependencies.

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

On Windows, :mod:`~spacepy.pycdf` now looks in more locations for the
NASA CDF library. Newer versions of the library by default install to
a different location (``Program Files``). The DLL is also now placed
in the ``bin`` directory instead of ``lib``, so ``bin`` is searched
and the value of environment variable ``CDF_BIN`` in addition to
``lib`` and ``CDF_LIB``. The net effect should be to increase the
chance of successfully loading the library, with a small chance of
accidentally loading the wrong one.

The default data source for leapsecond files has been reverted from
NASA/MODIS to the USNO, as USNO data services are back online. If
present, entries in the :doc:`configuration file <configuration>` will
still be used instead of the default.

0.2 Series
==========

0.2.3 (2021-10-30)
------------------
This is the last release of the 0.2 series and the last with full
support for Python 2. Binary installers (including
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

The 0.2 series will be the last with full support for Python 2;
0.2.3 will likely be the last release. Binary installers
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

``spacepy.pybats.rim.fix_format()`` is now deprecated, as
:class:`~spacepy.pybats.rim.Iono` can now read these files directly.

Quaternion math functions have been moved to
:mod:`~spacepy.coordinates`; using the functions in
:mod:`~spacepy.toolbox` is deprecated.

``spacepy.toolbox.feq()`` is deprecated; numpy 1.7 added the equivalent
:func:`~numpy.isclose`.

The :class:`~spacepy.plot.spectrogram` class is now capitalized
(:class:`~spacepy.plot.spectrogram`); the old, lower-case variant is
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
