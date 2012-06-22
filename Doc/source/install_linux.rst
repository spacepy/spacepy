******************
Linux Installation
******************

Dependencies
============
To ensure that you have all dependencies for SpacePy satisfied, there are two
approaches. The first is to download a single software suite that installs a
comprehensive set of Python libraries; the second is to install the core 
dependencies yourself.

The Easy Way
------------
Either the `Enthought Python Distribution (EPD) <http://www.enthought.com/>`_
or the `Python(x,y) distribution <https://code.google.com/p/pythonxy-linux/>`_
will provide Python, numpy, scipy, matplotlib and a host of other useful
3rd-party libraries. In both cases you will still need to install :ref:`ffnet` 
by hand.

If you are installing by hand, follow the instructions for your linux
distribution below.

Debian and Ubuntu
-----------------
The following command will install most of the dependencies. It has
been checked for Ubuntu 11.10 and Debian 7.0 "wheezy"::

    sudo apt-get install python-dev python-numpy build-essential \
    python-scipy python-matplotlib python-networkx python-h5py \
    python-f2py gfortran ncurses-dev

You can also, of course, install the same packages via synaptic or
other package manager of your choice.

Since no packages are available for them, install :ref:`CDF`
and :ref:`ffnet` by hand.

Other distributions
-------------------
For other distributions, check :doc:`dependencies` and install by hand or via your package manager. Once you figure it out, please contact the SpacePy team so we can update this documentation.


.. _CDF:

CDF
---
Download the latest `CDF library <http://cdf.gsfc.nasa.gov/>`_. Choose
the file ending in ``-dist-all.tar.gz`` from the ``linux``
directory. Untar and cd into the resulting directory. Then build::

    make OS=linux ENV=gnu CURSES=yes FORTRAN=no UCOPTIONS=-O2 SHARED=yes all

Use ``CURSES=no`` if the curses library is not installed. (The
distribution-specific directions above will install curses.)

Install::

    sudo make install

This will install the library into the default location ``/usr/local/cdf``, where SpacePy can find it. If you choose to install elsewhere, see the CDF documentation, particularly the notes on the ``CDF_BASE`` and ``CDF_LIB`` environment variables. SpacePy uses these variables to find the library.


.. _ffnet:

ffnet
-----
Download the latest `ffnet module
<http://ffnet.sourceforge.net/install.html>`_. Compilation requires
f2py (from numpy), installed in the distribution-specific directions
above. Untar and cd into the resulting directory. Then build::

    python setup.py build

Either install just for one user::

    python setup.py install --user

Or install for all users on the system::

    sudo python setup.py install

Normally the correct Fortan compiler will be found; if compilation
fails, try specifying the older GNU compiler at the build step::

    python setup.py build --fcompiler=gnu

SpacePy
=======
With the dependencies installed, SpacePy is ready to build and install. This uses the same basic setup as ffnet (standard Python distutils).

Build::

     python setup.py build

If this fails, specify a fortran compiler::

    python setup.py build --fcompiler=gnu

(``python setup.py build --help-fcompiler`` will list options for
Fortan compilers.)

Install for one user::

    python setup.py install --user

Or install for all users on the system::

    sudo python setup.py install
