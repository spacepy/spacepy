********************
Windows Installation
********************

The SpacePy team currently provides binary "wheels" via PyPI so it can
be installed on Windows without a compiler. Binaries are provided for
Python 2.7 and 3.6 through 3.9 in 64-bit and 32-bit variants for each.
``pip install spacepy`` should find and install these binaries.

32-bit Windows binaries will be discontinued with SpacePy 0.3, so
installing 32-bit SpacePy will require :ref:`compiling <windows_compiling>`.

Our recommended (but not required) Python distribution is `Anaconda
<https://docs.anaconda.com/anaconda/>`_ running 64-bit
Python 3. Anaconda includes much of the scientific Python
stack. Another excellent distribution is `Canopy
<https://www.enthought.com/product/canopy/>`_.

You may need to install the dependencies some way other than pip; for
example, if you are running an earlier version of Python. The latest
version of many dependencies requires Python 3.6 and pip will not
install older versions to get around this. See :ref:`win_dep_conda`.

.. contents::
   :local:

.. _windows_fortran:

Fortran and ffnet
=================

`ffnet <http://ffnet.sourceforge.net/download.html>`_ is required for
:mod:`~spacepy.LANLstar`. It can be installed either before or after
SpacePy. Binary wheels are not provided, so a Fortran compiler is
required.

With Anaconda, the compiler and ffnet can be installed with::

  conda install m2w64-gcc-fortran libpython
  SET FC_VENDOR=gfortran
  pip install ffnet

The `FC_VENDOR` line is necessary because ffnet defaults to the
Microsoft compiler.

Some standalone binary installers (no ``pip`` support) are also
available on the `ffnet site
<http://ffnet.sourceforge.net/download.html>`_. These do not require a
compiler but support only a limited set of Python versions.

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
<https://github.com/spacepy/spacepy/tree/master/developer/scripts>`_
the SpacePy developers use.

.. _windows_CDF:

NASA CDF
========

:mod:`~spacepy.pycdf` requires the `NASA CDF library
<https://cdf.gsfc.nasa.gov/html/sw_and_docs.html>`_ . Binary
installers are available for Windows; be sure to pick the version
that matches your Python installation. The current 32-bit version
is `cdf37_1_0-setup-32.exe
<https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/windows/cdf37_1_0-setup-32.exe>`_;
for 64-bit, `cdf37_1_0-setup-64.exe
<https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/cdf37_1/windows/cdf37_1_0-setup-64.exe>`_.

This is a simple self-extracting installer that can be installed either before or after installing SpacePy.

Standalone installers
=====================

Standalone installers for Windows will be discontinued starting with
SpacePy 0.3. Binary wheels will still be provided and can be
downloaded from the same location and installed with ``pip``,
e.g. ``pip install spacepy-0.2.2-cp38-cp38-win_amd64.whl``.

Self-extracting and self-installing executables are also available for
download direct from `our github
<https://github.com/spacepy/spacepy/releases>`__. Be sure to choose the
file that matches your Python installation, both in Python version and
word size (64-bit vs. 32-bit.)
E.g. ``spacepy-0.2.0.win-amd64.py36.exe`` is the installer for Python
3.6 on 64-bit Windows and ``spacepy-0.2.0.win32.py27.exe`` is the
installer for Python 2.7 on 32-bit Windows.

Once downloaded, these can be installed without an internet connection.

If using these installers, the :doc:`dependencies` will not be
installed automatically.

.. _win_dep_conda:

Dependencies via conda
======================

Installation via ``pip`` will automatically install most Python
dependencies (but not :ref:`ffnet <windows_fortran>` or the :ref:`NASA
CDF library <windows_CDF>`). They can also be installed from conda::

  conda install numpy scipy matplotlib networkx h5py

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
