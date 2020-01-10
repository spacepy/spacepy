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
here. :doc:`dep_versions`. describes future support.

Hard Dependencies
=================
Without these packages installed, SpacePy will not function.

Python 2.7+
-----------

`Python <http://www.python.org/>`_ is the core language for SpacePy.
Python 3 is strongly recommended; while SpacePy currently supports
Python 2.7, this support will be phased out over the course
of 2020. See :doc:`py2k_eol`.

Required to install SpacePy.

NumPy 1.4+
----------
`NumPy <http://numpy.scipy.org/>`_ provides the
high-performance array data structure used throughout SpacePy. Version
1.4 or later is required; 1.6 or later recommended.

Required to install SpacePy. f2py is part of NumPy, but is sometimes
packaged separately; it is required (at installation time) if
:mod:`~spacepy.irbempy` is to be used.

Due to a numpy bug, numpy 1.15.0 is not supported. Use 1.15.1 or later.

dateutil
--------
If you choose not to install :ref:`matplotlib <dependencies_mpl>`,
`dateutil <http://labix.org/python-dateutil>`_ is required.
(Installing matplotlib will fulfill this dependency.)

C compiler
----------
If you are installing SpacePy from source, a working C compiler
is required. (Not necessary for the Windows binary installer.)

Soft Dependencies
=================
Without these packages, SpacePy will install, but certain features may
not be available. Usually an ImportError means a dependency is missing.

These are simply marked as dependencies in SpacePy metadata and thus
will be automatically installed when using dependency-resolving
methods such as pip.

SciPy 0.10+
-----------
`SciPy <http://www.scipy.org/>`_ provides several useful scientific
and numerical functions build on top of NumPy.  It is highly
recommended. The following modules may have limited functionality
without SciPy:

    * :mod:`~spacepy.empiricals`
    * :mod:`~spacepy.seapy`
    * :mod:`~spacepy.toolbox`


.. _dependencies_mpl:

matplotlib 1.5.0+
-----------------
`matplotlib <http://matplotlib.sourceforge.net/>`_ is the preferred
plotting package for Python. It is highly recommended. Without it, you
will not be able to effectively visualize data, and the following
modules may have limited functionality:

    * :mod:`~spacepy.data_assimilation`
    * :mod:`~spacepy.plot`
    * :mod:`~spacepy.poppy`
    * :mod:`~spacepy.pybats`
    * :mod:`~spacepy.radbelt`
    * :mod:`~spacepy.time`
    * :mod:`~spacepy.toolbox`

ffnet
-----
`ffnet <http://ffnet.sourceforge.net/>`_ is a neural network package,
required for :mod:`~spacepy.LANLstar`.

networkx
--------
`networkx <http://networkx.lanl.gov/>`_ is a requirement for ffnet,
and thus :mod:`~spacepy.LANLstar`.

h5py
----
`h5py <http://code.google.com/p/h5py/>`_ provides a Python interface to
HDF5 files. It is required for the HDF import/export capability of
:mod:`~spacepy.datamodel` and for use of the :mod:`~spacepy.omni` module.

CDF
---
NASA's `CDF <http://cdf.gsfc.nasa.gov/>`_ library provides access to
Common Data Format files. It is required for :mod:`~spacepy.pycdf`,
and thus for the CDF import/export capability of
:mod:`~spacepy.datamodel`.

.. warning::
    Unlike the Python-based dependencies, the CDF library must be
    installed if pycdf support is needed; it will not be automatically
    installed.

Fortran compiler
----------------
If installing from source, :mod:`~spacepy.irbempy` requires a Fortran
compiler. (This is not required for the Windows binary installer).
Supported compilers are the GNU compiler ``gfortran``, the older GNU
compiler ``g77``, and the Portland Group PGI compiler.

If :mod:`~spacepy.irbempy` is to be used, the Fortran compiler (and
f2py) must be installed before SpacePy.
