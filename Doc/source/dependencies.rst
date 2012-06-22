********************
SpacePy Dependencies
********************

SpacePy relies on several other pieces of software for complete
functionality.  :doc:`install` links to details on
installing the required software for each platform.

Hard Dependencies
=================
Without these packages installed, SpacePy will not function.
The standard installer checks for these dependencies and will
not install without them.

Python 2.6+
-----------
`Python <http://www.python.org/>`_ is the core
language for SpacePy.  Python 2, version 2.6 or later is
required. Much of SpacePy works under Python 3, but it will not be
fully supported until all our dependencies work on Python 3. *Python 3
is not simply a "newer" Python 2; there are substantial differences in
the language*. See `Should I use Python 2 or Python 3?
<http://wiki.python.org/moin/Python2orPython3>`_.

NumPy 1.4+
----------
`NumPy <http://numpy.scipy.org/>`_ provides the
high-performance array data structure used throughout SpacePy. Version
1.4 or later is required; 1.6 or later recommended.

Soft Dependencies
=================
Without these packages, SpacePy will install, but certain features may
not be available. Usually an ImportError means a dependency is missing.

SciPy
-----
`SciPy <http://www.scipy.org/>`_ provides several useful scientific
and numerical functions build on top of NumPy.  It is highly
recommended. The following modules may have limited functionality
without SciPy:

    * :mod:`~spacepy.empiricals`
    * :mod:`~spacepy.seapy`
    * :mod:`~spacepy.toolbox`

matplotlib
----------
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
`h5py <http://code.google.com/p/h5py/>`_ provies a Python interface to
HDF5 files. It is required for the HDF import/export capability of
:mod:`~spacepy.datamodel`.

CDF
---
NASA's `CDF <http://cdf.gsfc.nasa.gov/>`_ library provides access to
Common Data Format files. It is required for :mod:`~spacepy.pycdf`,
and thus for the CDF import/export capability of
:mod:`~spacepy.datamodel`.
