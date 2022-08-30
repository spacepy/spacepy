******************
Linux Installation
******************

Installation on Linux requires both a C and a Fortran compiler; a
recent GCC is recommended (the C compiler is likely included with your
distribution). On Debian and Ubuntu::
  
      sudo apt-get install gfortran

Once this is set up, ``pip install spacepy`` should Just Work. If
you're installing as a single user (not in a virtual environment) then
add the ``--user`` flag.

You will also need the :ref:`NASA CDF library <linux_CDF>` to use
:mod:`~spacepy.pycdf`.

Our recommended (but not required) Python distribution is `Anaconda
<https://docs.anaconda.com/anaconda/>`_ running 64-bit
Python 3. Anaconda includes much of the scientific Python
stack. Another excellent distribution is `Canopy
<https://www.enthought.com/product/canopy/>`_.

You may need to install the dependencies some way other than pip; for
example, if you are running an earlier version of Python. The latest
version of many dependencies requires Python 3.6 and pip will not
install older versions to get around this. See :ref:`linux_dep_conda`
and :ref:`linux_dep_apt`.

.. contents::
   :local:

.. _linux_dep_conda:

Dependencies via conda
======================

Installation via ``pip`` will automatically install most Python
dependencies (but not the :ref:`NASA CDF library <linux_CDF>`). They
can also be installed from conda::

  conda install numpy scipy matplotlib h5py

.. _linux_dep_apt:

Dependencies via system packages
================================

SpacePy usually works with the system Python on Linux. To install dependencies via the package manager on Debian or Ubuntu::

  sudo apt-get install python-dev python-h5py python-matplotlib python-numpy python-scipy

For Python 3, use::

  sudo apt-get install python3-dev python3-h5py python3-matplotlib python3-numpy python3-scipy

For other distributions, check :doc:`dependencies` and install by hand
or via your package manager. 

.. _linux_CDF:

CDF
===

It is recommended to install the ncurses library; on Ubuntu and Debian::

    sudo apt-get install ncurses-dev


Download the latest `CDF library <http://cdf.gsfc.nasa.gov/>`_. Choose
the file ending in ``-dist-all.tar.gz`` from the ``linux``
directory. Untar and cd into the resulting directory. Then build::

    make OS=linux ENV=gnu CURSES=yes FORTRAN=no UCOPTIONS=-O2 SHARED=yes all

Use ``CURSES=no`` if the curses library is not installed. (The
distribution-specific directions above will install curses.)

Install::

    sudo make install

This will install the library into the default location ``/usr/local/cdf``, where 
SpacePy can find it. If you choose to install elsewhere, see the CDF documentation, 
particularly the notes on the ``CDF_BASE`` and ``CDF_LIB`` environment variables. 
SpacePy uses these variables to find the library.

Compiling
=========

With the dependencies installed, SpacePy can be built from source.
This uses standard Python distutils.
You can always get the latest source code for SpacePy from our `github
repository <https://github.com/spacepy/spacepy>`_ and the latest
release from `PyPI <https://pypi.org/project/SpacePy/#files>`_

Build::

     python setup.py build

If this fails, specify a Fortran compiler::

    python setup.py build --fcompiler=gnu95

``python setup.py build --help-fcompiler`` will list options for
Fortran compilers. Currently available compilers are ``pg``,
``gnu95``, ``gnu``, ``intelem``, ``intel`` or ``none`` (to skip all
Fortran); ``gnu95`` (the GNU gfortran compiler) is recommended.

Install for one user::

    python setup.py install --user

If you're using conda, installation as user isn't recommended::

    python setup.py install

Or install for all users on the system::

    sudo python setup.py install

If you want to build the documentation yourself (rather than using the
documentation shipped with SpacePy), install sphinx and numpydoc. The
easiest way is via pip::

  pip install sphinx numpydoc

They are also available via conda::

  conda install sphinx numpydoc

Or the package manager:

  sudo apt-get install python-sphinx python-numpydoc

For Python 3:

  sudo apt-get install python3-sphinx python3-numpydoc

Raspberry Pi
============
SpacePy works on Raspberry Pi, using Raspberry Pi OS in 32-bit or
64-bit flavors. A few tips:

   * It is highly recommended to install all dependencies (numpy,
     etc.)  via the system package manager ``apt-get`` rather than
     pip, as prebuilt wheels are not generally available and compiling
     dependencies on the Pi can take a very long time::

      sudo apt-get install gfortran python3-numpy python3-scipy python3-h5py python3-matplotlib

   * Similarly, if installing SpacePy via pip, use the
     ``--no-build-isolation`` flag to use the system numpy.
