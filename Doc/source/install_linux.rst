******************
Linux Installation
******************

Dependencies
============
Installing SpacePy on linux using pip should automatically install the core
dependencies. Installing numpy prior to SpacePy should make the installation
process a little smoother.

The first is to download a single software suite that installs a
comprehensive set of Python libraries; the second is to install the core 
dependencies yourself.

The Easy Way
------------
Both `Anaconda <https://www/anaconda.com>`_ and 
`Enthought Canopy <https://www.enthought.com/>`_
will provide Python, numpy, scipy, matplotlib and a host of other useful
3rd-party libraries. In both cases you will still need to install 
:ref:`CDF <linux_CDF>` by hand.

Assuming you have a fortran compiler installed, you should be able to install
Spacepy using pip::

    pip install spacepy

If you're installing as a single user (not in a virtual environment)
then add the --user flag.

If you are installing by hand, follow the instructions for your linux
distribution below.

Debian and Ubuntu
-----------------
Installation on the Python distribution and scientific stack via the
system package manager is no longer recommended. Use of pip or conda
is recommended so that packages updates can be obtained in a timely
manner. However, your system may not come with some of the required
packages, like gfortran.

The following command should install gfortran on Ubuntu and Debian::

    sudo apt install gfortran

To make installing the NASA CDF library easier, and to enable the CDF
command line tools to be built::

    sudo apt install ncurses-dev

You can also, of course, install the same packages via synaptic or
other package manager of your choice.

Since no packages are available that provide it, install 
:ref:`CDF <linux_CDF>` by hand.

Other distributions
-------------------
For other distributions, check :doc:`dependencies` and install by hand
or via your package manager. 


.. _linux_CDF:

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

This will install the library into the default location ``/usr/local/cdf``, where 
SpacePy can find it. If you choose to install elsewhere, see the CDF documentation, 
particularly the notes on the ``CDF_BASE`` and ``CDF_LIB`` environment variables. 
SpacePy uses these variables to find the library.


.. _linux_ffnet:

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

Normally the correct Fortran compiler will be found; if compilation
fails, try specifying the older GNU compiler at the build step::

    python setup.py build --fcompiler=gnu

SpacePy
=======
With the dependencies installed, SpacePy is ready to build and install.
This uses the same basic setup as ffnet (standard Python distutils).
You can always get the latest source code for SpacePy from our
`github repository <https://github.com/spacepy/spacepy>`_.

Build::

     python setup.py build

If this fails, specify a Fortran compiler::

    python setup.py build --fcompiler=gnu95

(``python setup.py build --help-fcompiler`` will list options for
Fortran compilers.)

Install for one user::

    python setup.py install --user

If you're using conda, installation as user isn't recommended::

    python setup.py install

Or install for all users on the system::

    sudo python setup.py install
