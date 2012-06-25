******************
MacOS Installation
******************

The Easy Way
============

Download and install the `Enthought Python Distribution (EPD)
<http://www.enthought.com/>`_. The free version works fine. (If you are
considering purchasing a subscription, keep in mind that Enthought
supports numpy, scipy, and other scientific Python development)

The other dependencies (except ffnet) should be installable with ``easy_install``.
From a command prompt, run::

    easy_install h5py

Finally, install SpacePy using the installer tarball. 

.. _CDF:

CDF
---
If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_.   The default installation directory is recommended to
help SpacePy find the library.  There is an OSX install package.

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

Normally the correct Fortran compiler will be found; if compilation
fails, try specifying the older GNU compiler at the build step::

    python setup.py build --fcompiler=gnu


SpacePy
-------
With the dependencies installed, SpacePy is ready to build and install. This uses the same basic setup as ffnet (standard Python distutils).

Build::

     python setup.py build

If this fails, specify a Fortran compiler::

    python setup.py build --fcompiler=gnu

(``python setup.py build --help-fcompiler`` will list options for
Fortran compilers.)

Install for one user::

    python setup.py install --user

Or install for all users on the system::

    sudo python setup.py install


------------


Expert MacPorts Install Instructions
====================================
These were current as of 22-June-2012, path especially go out of date quickly.

REQUIRED OR STRONGLY RECOMMENDED
--------------------------------

    #. Install xcode from the app store
    #. Open xcode
    #. Xcode -> preferences -> Downloads -> install command line tools
    #. Install MacPorts (http://www.macports.org/install.php)
    #. Setup proxies if needed
    #. sudo port -v selfupdate
    #. sudo xcode-select -switch /Applications/Xcode.app
    #. sudo port install python27
    #. sudo port select --set python python27
    #. sudo port install py27-ipython py27-scipy py27-numpy py27-matplotlib readline py27-h5py ipython-select
    #. sudo port select --set ipython ipython27
    #. sudo port install py27-sphinx wget 
    #. sudo port select --set sphinx py27-sphinx
    #. wget ftp://cdaweb.gsfc.nasa.gov/pub/cdf/dist/cdf34_0/macosX/cdf34_0-setup_universal_binary.tar.gz
    #. tar -zxvf cdf34_0-setup_universal_binary.tar.gz 
    #. sudo /usr/sbin/installer -pkg CDF3400ub.pkg -target /

OPTIONAL (ALSO RECOMMENDED)
---------------------------
    #. sudo port install py27-spyder py27-xlrd py27-xlwt  py27-coverage py27-pyside (I have to run this exact command multiple times)
    #. install xquartz (http://xquartz.macosforge.org/landing/)
    #. wget http://xquartz.macosforge.org/downloads/SL/XQuartz-2.7.1.dmg
    #. hdiutil attach XQuartz-2.7.1.dmg 
    #. cd /Volumes/XQuartz-2.7.1/
    #. sudo /usr/sbin/installer -pkg XQuartz.pkg -target /
    #. cd ~
    #. hdiutil detach /Volumes/XQuartz-2.7.1
    #. sudo port install git-core
    #. Install numpydoc (http://pypi.python.org/pypi/numpydoc)
    #. wget http://pypi.python.org/packages/source/n/numpydoc/numpydoc-0.4.tar.gz#md5=e5bdd98f84f2bb220373819e20c27091
    #. tar -zxvf numpydoc-0.4.tar.gz
    #. cd numpydoc-0.4
    #. python setup.py install --user

INSTALL SPACEPY
---------------
    #. Download and install gfortran (http://gcc.gnu.org/wiki/GFortranBinaries)
    #. wget http://quatramaran.ens.fr/~coudert/gfortran/gfortran-4.6.2-x86_64-Lion.dmg
    #. hdiutil attach gfortran-4.6.2-x86_64-Lion.dmg
    #. cd /Volumes/gfortran-4.6.2-x86_64-Lion/
    #. sudo /usr/sbin/installer -pkg gfortran.pkg -target /
    #. Download source (or clone from git) (http://spacepy.lanl.gov/download.shtml)
    #. cd spacepy/
    #. python setup.py install –user
    #. python setup.py install --user --build-docs
    #. cd tests/
    #. python test_spacepy.py (there should be no errors or fails)

POST INSTALL TWEAKING
---------------------
    #. Create .matplotlib/matplotlibrc
    #. Add:   backend      : MacOSX
    #. Add:   interactive : True






