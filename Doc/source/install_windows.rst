********************
Windows Installation
********************

The SpacePy team currently provides binary "wheels" via PyPI so it can
be installed on Windows without a compiler. Binaries are provided for
Python 3.6 through 3.10 in 64-bit variant only.
``pip install spacepy`` should find and install these binaries.

Our recommended (but not required) Python distribution is `Anaconda
<https://docs.anaconda.com/anaconda/>`_ running 64-bit
Python 3. Anaconda includes much of the scientific Python
stack. Another excellent distribution is `Canopy
<https://www.enthought.com/product/canopy/>`_.

If you prefer to install the dependencies some way other than pip, see
:ref:`win_dep_conda`.

.. contents::
   :local:

.. _windows_compiling:

Compiling
=========

If a binary wheel is not available for your version of Python, ``pip``
will try to compile SpacePy. The only supported compiler is
``mingw32``. Install it with::

  conda install m2w64-gcc-fortran libpython

This is also required if installing from a source distribution or git checkout.

:mod:`~spacepy.irbempy` requires Fortran to compile and the only
supported compiler is ``gnu95``; this is the default and provided
by ``m2w64-gcc-fortran``.

If you have difficulties, it may be useful to reference the `build
scripts
<https://github.com/spacepy/spacepy/tree/main/developer/scripts>`_
the SpacePy developers use.

.. _windows_CDF:

NASA CDF
========

SpacePy binary wheels ship with a copy of the NASA CDF library and
:mod:`~spacepy.pycdf` will use this copy if no other CDF libraries can
be found.

If you build SpacePy from source or wish to use a different version of
the library, you can download it from the `NASA CDF page
<https://cdf.gsfc.nasa.gov/html/sw_and_docs.html>`_. Binary
installers are available for Windows; be sure to pick the version
(32-bit or 64-bit) that matches your Python installation.

NASA CDF can be installed either before or after installing SpacePy.

.. _win_dep_conda:

Dependencies via conda
======================

Installation via ``pip`` will automatically install most Python
dependencies (but not the :ref:`NASA CDF library <windows_CDF>`).
They can also be installed from conda::

  conda install numpy python-dateutil scipy matplotlib h5py

Standalone dependencies
=======================

Most of the :doc:`dependencies` have Windows installers available via
their pages, but ``pip`` or ``conda`` are recommended instead.

Developers
==========

If you want to build the documentation yourself (rather than using the
documentation shipped with SpacePy), install sphinx and numpydoc. The
easiest way is via pip::

  pip install sphinx numpydoc

They are also available via conda::

  conda install sphinx numpydoc
