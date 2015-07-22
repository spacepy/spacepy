******************
MacOS Installation
******************

Binary Installers
=================

OSX install has become a challenge.  With the Enthought transition to Canopy we
cannot figure out clean install directions for 3rd party packages and therefore
can no longer recommend using EPD for SpacePy.

Python
------
The Python that ships with OSX is not going to work.  We recommend Python 2.7
from http://www.python.org/download/

Numpy
-----
There is an OSX package at http://sourceforge.net/projects/numpy/files/NumPy/1.7.1/

Scipy
-----
There is an OSX package at http://sourceforge.net/projects/scipy/files/scipy/0.12.0/

Matplotlib
----------
There is an OSX package at http://matplotlib.org/downloads.html

.. _CDF:

CDF
---
If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_.   The default installation directory is recommended to
help SpacePy find the library.  There is an OSX install package.
    #. Get the package from ftp://cdaweb.gsfc.nasa.gov/pub/software/cdf/dist/cdf34_1/macosX/
    #. Unzip and install

.. _ffnet:

ffnet
-----
    easy_install --user ffnet

SpacePy
-------
With the dependencies installed, SpacePy is ready to build and install. This uses the same basic setup as ffnet (standard Python distutils).

Build::

     python setup.py build

If this fails, specify a FORTRAN compiler::

    python setup.py build --fcompiler=gnu

(``python setup.py build --help-fcompiler`` will list options for
FORTRAN compilers.)

Install for one user::

    python setup.py install --user

Or install for all users on the system::

    sudo python setup.py install


------------


Full MacPorts Install Instructions
====================================
These were current as of 13-May-2013, path especially go out of date quickly.
We recommend installing this way for the best results.

Dependencies
------------

    #. Install xcode from the app store
    #. Open xcode
    #. Xcode -> preferences -> Downloads -> install command line tools
    #. Install MacPorts (http://www.macports.org/install.php)
    #. sudo port -v selfupdate
    #. sudo xcode-select -switch /Applications/Xcode.app
    #. sudo port install gcc49 +gfortran
    #. sudo port select gcc mp-gcc49
    #. sudo port install python27
    #. sudo port select --set python python27
    #. sudo port install py27-ipython py27-scipy py27-numpy py27-matplotlib readline py27-h5py ipython-select
    #. sudo port select --set ipython ipython27
    #. Download spacepy source (or clone from git) (http://spacepy.lanl.gov/download.shtml)

INSTALL SPACEPY
---------------
    #. python setup.py install -–user  (in the spacepy directory after unzip)

POST INSTALL TWEAKING
---------------------
    #. Create .matplotlib/matplotlibrc
    #. Add:   backend      : MacOSX
    #. Add:   interactive : True






