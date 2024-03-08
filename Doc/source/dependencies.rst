********************
SpacePy Dependencies
********************

SpacePy relies on several other pieces of software for complete
functionality.  :doc:`install` links to details on
installing the required software for each platform.

Unless otherwise noted, a dependency may be installed *after*
SpacePy, and the new functionality will be available the next
time SpacePy is imported.

Currently required versions are documented
here. :doc:`dep_versions` describes future support.

Hard Dependencies
=================
Without these packages installed, SpacePy will not
function. Installing via ``pip`` will normally install Python-based
dependencies automatically.

Python 3.6+
-----------

`Python <http://www.python.org/>`_ is the core language for SpacePy.

Required to install SpacePy.

NumPy 1.15.1+
-------------
`NumPy <http://numpy.scipy.org/>`_ provides the
high-performance array data structure used throughout SpacePy. Version
1.15.1 or later is required.

Required to install SpacePy. f2py is part of NumPy, but is sometimes
packaged separately; it is required (at installation time) if
:mod:`~spacepy.irbempy` is to be used.

Due to a numpy bug, numpy 1.15.0 is not supported. Use 1.15.1 or later.

On Python 3.9, numpy 1.18 or later is required.

dateutil
--------
`dateutil <https://dateutil.readthedocs.io/>`_ 2.1 or later is required.

C compiler
----------
If you are installing SpacePy from source, a working C compiler is
required. (Not necessary for binary wheels, provided on most
platforms.)

pip
---
`pip <https://pip.pypa.io/>`_ is now the expected installer for Python
packages. Almost any modern Python installation will include it.

setuptools
----------
`setuptools <https://setuptools.pypa.io>`_ is the installer backend. It
is commonly included and, if necessary, installable with ``pip``.

Soft Dependencies
=================
Without these packages, SpacePy will install, but certain features may
not be available. Usually an ImportError means a dependency is missing.

These are simply marked as dependencies in SpacePy metadata and thus
will be automatically installed when using dependency-resolving
methods such as pip.

.. _dependencies_scipy:

SciPy 1.0+
----------
`SciPy <http://www.scipy.org/>`_ provides several useful scientific
and numerical functions build on top of NumPy.  It is highly
recommended. The following modules may have limited functionality
without SciPy:

    * :mod:`~spacepy.coordinates`
    * :mod:`~spacepy.ctrans`
    * :mod:`~spacepy.empiricals`
    * :mod:`~spacepy.seapy`
    * :mod:`~spacepy.toolbox`


.. _dependencies_mpl:

matplotlib 3.1+
---------------
`matplotlib <http://matplotlib.sourceforge.net/>`_ is the preferred
plotting package for Python. It is highly recommended. Without it, you
will not be able to effectively visualize data, and the following
modules may have limited functionality or fail entirely:

    * :mod:`~spacepy.plot`
    * :mod:`~spacepy.poppy`
    * :mod:`~spacepy.pybats`
    * :mod:`~spacepy.radbelt`
    * :mod:`~spacepy.seapy`
    * :mod:`~spacepy.toolbox`

.. _dependencies_h5py:

h5py 2.10+
----------
`h5py <http://code.google.com/p/h5py/>`_ provides a Python interface to
HDF5 files. It is required for the HDF import/export capability of
:mod:`~spacepy.datamodel` and for use of the :mod:`~spacepy.omni` module.

.. _dependencies_cdf:

CDF 3.5+
--------
NASA's `CDF <http://cdf.gsfc.nasa.gov/>`_ library provides access to
Common Data Format files. It is required for :mod:`~spacepy.pycdf`,
and thus for the CDF import/export capability of
:mod:`~spacepy.datamodel`. The library is included with SpacePy binary
wheels, but a locally-installed library will be used instead if it is
found.

.. warning::
    If building SpacePy from source, the CDF library must be installed
    if pycdf support is needed. It can be installed after SpacePy.

.. _dependencies_fortran:

Fortran compiler
----------------
If installing from source, :mod:`~spacepy.irbempy` requires a Fortran
compiler. This is not required if using a binary wheel. The supported
compiler is the GNU compiler ``gfortran``.

If :mod:`~spacepy.irbempy` is to be used, the Fortran compiler (and
f2py) must be installed before SpacePy.

:mod:`~spacepy.coordinates` requires :mod:`~spacepy.irbempy` to use
the IRBEM-based backend, but the new CTrans-based backend can be used
without Fortran. See the :mod:`~spacepy.coordinates` documentation for
the ``use_irbem`` option.

.. _dependencies_astropy:

Astropy 1.0+
------------
:mod:`~spacepy.time` requires Astropy if conversion to/from
Astropy :class:`~astropy.time.Time` is desired.

:mod:`~spacepy.coordinates` requires Astropy if conversion to/from
Astropy :class:`~astropy.coordinates.SkyCoord` is desired.

Because Astropy is only required for these conversions, it is not
automatically installed with SpacePy.

Soft Dependency Summary
=======================

The following table summarizes, by SpacePy module, the functionality
that is *lost* if a soft dependency is not installed. If there is
nothing for a given dependency/module combination, the module is
unaffected by that dependency.

.. list-table:: SpacePy functionality lost without soft dependencies
   :header-rows: 1
   :stub-columns: 1

   * -
     - :ref:`CDF <dependencies_cdf>`
     - :ref:`Fortran compiler <dependencies_fortran>`
     - :ref:`h5py <dependencies_h5py>`
     - :ref:`matplotlib <dependencies_mpl>`
     - :ref:`SciPy <dependencies_scipy>`
     - :ref:`AstroPy <dependencies_astropy>`
   * - :mod:`~spacepy.coordinates`
     -
     - :class:`~spacepy.coordinates.Coords` IRBEM backend (except binaries)
     -
     -
     - :mod:`Entire module <spacepy.coordinates>`
     - * :meth:`~spacepy.coordinates.Coords.from_skycoord`
       * :meth:`~spacepy.coordinates.Coords.to_skycoord`
   * - :mod:`~spacepy.ctrans`
     -
     -
     -
     -
     - :mod:`Entire module <spacepy.ctrans>`
     -
   * - :mod:`~spacepy.datamodel`
     - * :meth:`~spacepy.datamodel.SpaceData.toCDF`
       * :func:`~spacepy.datamodel.fromCDF`
       * :func:`~spacepy.datamodel.toCDF`
       * (except binaries)
     -
     - * :meth:`~spacepy.datamodel.SpaceData.toHDF5`
       * :func:`~spacepy.datamodel.fromHDF5`
       * :func:`~spacepy.datamodel.toHDF5`
     -
     -
     -
   * - :mod:`~spacepy.empiricals`
     -
     -
     -
     -
     - * :func:`~spacepy.empiricals.vampolaPA`
       * :func:`~spacepy.empiricals.omniFromDirectionalFlux`
     -
   * - :mod:`~spacepy.irbempy`
     -
     - :mod:`Entire module <spacepy.irbempy>` (except binaries)
     -
     -
     -
     -
   * - :mod:`~spacepy.LANLstar`
     -
     -
     -
     -
     -
     -
   * - :mod:`~spacepy.omni`
     -
     -
     - :mod:`Entire module <spacepy.omni>`
     -
     -
     -
   * - :mod:`~spacepy.plot`
     -
     -
     -
     - :mod:`Entire module <spacepy.plot>`
     -
     -
   * - :mod:`~spacepy.poppy`
     -
     -
     -
     - * :meth:`~spacepy.poppy.PPro.assoc`
       * :meth:`~spacepy.poppy.PPro.plot`
       * :meth:`~spacepy.poppy.PPro.plot_mult`
       * :func:`~spacepy.poppy.plot_two_ppro`
     -
     -
   * - :mod:`~spacepy.pybats`
     -
     -
     -
     - * :meth:`~spacepy.pybats.bats.Bats2d.regrid`
       * :mod:`~spacepy.pybats.dgcpm`
       * :mod:`~spacepy.pybats.interact`
       * :mod:`~spacepy.pybats.kyoto`
       * :mod:`~spacepy.pybats.pwom`
       * :mod:`~spacepy.pybats.ram`
       * :mod:`~spacepy.pybats.rim`

       All plotting functions:

       * :func:`~spacepy.pybats.add_body`
       * :func:`~spacepy.pybats.add_planet`
       * :meth:`~spacepy.pybats.ImfInput.add_pram_bz`
       * :meth:`~spacepy.pybats.ImfInput.quicklook`
       * :meth:`~spacepy.pybats.bats.BatLog.add_dst_quicklook`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_b_magsphere`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_b_magsphere_legacy`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_body`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_comp_plot`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_contour`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_cont_shell`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_grid_plot`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_pcolor`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_planet`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_plot`
       * :meth:`~spacepy.pybats.bats.Bats2d.add_stream_scatter`
       * :meth:`~spacepy.pybats.bats.MagGridFile.add_ae_quicklook`
       * :meth:`~spacepy.pybats.bats.MagGridFile.add_contour`
       * :meth:`~spacepy.pybats.bats.MagGridFile.add_kp_quicklook`
       * :meth:`~spacepy.pybats.bats.MagGridFile.add_orbit_plot`
       * :meth:`~spacepy.pybats.quotree.QTree.plot_res`
       * :meth:`~spacepy.pybats.quotree.Branch.plotbox`
       * :meth:`~spacepy.pybats.quotree.Branch.plot_res`
       * :func:`~spacepy.pybats.trace2d.test_asymtote`
       * :func:`~spacepy.pybats.trace2d.test_dipole`
     -
     -
   * - :mod:`~spacepy.pycdf`
     - :mod:`Entire module <spacepy.pycdf>` (except binaries)
     -
     -
     -
     -
     -
   * - :mod:`~spacepy.radbelt`
     -
     -
     -
     - * :meth:`~spacepy.radbelt.RBmodel.plot`
       * :meth:`~spacepy.radbelt.RBmodel.plot_obs`
     -
     -
   * - :mod:`~spacepy.seapy`
     -
     -
     -
     - :mod:`Entire module <spacepy.seapy>`
     - * :func:`~spacepy.seapy.sea_signif`
     -
   * - :mod:`~spacepy.time`
     -
     -
     -
     -
     -
     - AstroPy support in :class:`~spacepy.time.Ticktock`
   * - :mod:`~spacepy.toolbox`
     -
     -
     -
     - * :func:`~spacepy.toolbox.tCommon`
       * :func:`~spacepy.toolbox.linspace` if using
         :class:`~datetime.datetime` inputs
       * :func:`~spacepy.toolbox.logspace` if using
         :class:`~datetime.datetime` inputs
     - * :func:`~spacepy.toolbox.dist_to_list`
       * :func:`~spacepy.toolbox.intsolve`
       * :func:`~spacepy.toolbox.poisson_fit`
     -
