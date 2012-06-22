*******************
MacOS Installation
*******************

The Easy Way
===========

Download and install the `Enthought Python Distribution (EPD)
<http://www.enthought.com/>`_. The free version works fine. (If you are
considering purchasing a subscription, keep in mind that Enthought
supports numpy, scipy, and other scientific Python development, from
which we all benefit.) 

The other dependencies should be installable with ``easy_install``.
From a command prompt, run::

    easy_install ffnet h5py

Finally, install SpacePy using the installer EXE. Be sure to choose the
installer that matches your version of Python, either 2.6 or 2.7.



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




