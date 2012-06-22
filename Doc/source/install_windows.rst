********************
Windows Installation
********************

Although the SpacePy team makes every effort to keep code portable, and has
expended considerable effort on Windows specifically, Windows is not a
primary development platform for any of the developers. There may be some
limitations and/or bugs under Windows which are not present in OSX or Linux.
Bug reports and patches are always welcome.

The Easy Way
============

Download and install the `Enthought Python Distribution (EPD)
<http://www.enthought.com/>`_. The free version works fine. (If you are
considering purchasing a subscription, keep in mind that Enthought
supports numpy, scipy, and other scientific Python development, from
which we all benefit.) SpacePy ONLY supports 32-bit Python under
Windows; 32-bit Python should run fine on 64-bit Windows.

If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_. Again, the 32-bit installer is required, e.g.
CDF33_1-SETUP-32.EXE. The default installation directory is recommended to
help SpacePy find the library.

The other dependencies should be installable with ``easy_install``.
From a command prompt, run::

    easy_install ffnet h5py

Finally, install SpacePy using the installer EXE. Be sure to choose the
installer that matches your version of Python, either 2.6 or 2.7.


The Hard Way
============


This is a step-by-step guide to compiling and installing SpacePy from source.
The filenames listed for the dependencies are the latest at this writing.

Download and install the 32-bit `Python Windows installer
<http://python.org/download/>`_ .  SpacePy is developed against the
2.6 series, but either the 2.6 or 2.7 version should work.

Edit your Windows path: Right-click my computer, choose properties,
advanced, environment variables.  Edit the path variable: append
``;c:\python26;c:\python26\Scripts;c:\MinGW\bin`` (`more information
<http://docs.python.org/using/windows.html#finding-the-python-executable>`_).

Download the latest `mingwget <http://sourceforge.net/projects/mingw/files/Automated%20MinGW%20Installer/mingw-get/>`_. mingw32 is the only C and Fortran compiler supported by SpacePy on Windows. Unzip into C:\MinGW. At a command prompt (Start, Run, cmd) type::

      mingw-get install gcc g++ mingw32-make fortran

See also `the mingw docs <http://www.mingw.org/wiki/Getting_Started>`_

Download and install `numpy-1.6.1-win32-superpack-python2.6.exe
<http://sourceforge.net/projects/numpy/files/>`_.

Download and install `scipy-0.9.0-win32-superpack-python2.6.exe
<http://sourceforge.net/projects/scipy/files/>`_.

Download and install `matplotlib-1.0.1.win32-py2.6.exe
<http://matplotlib.sourceforge.net/>`_.

Download and install `setuptools-0.6c11.win32-py2.6.exe
<http://pypi.python.org/pypi/setuptools>`_.

From a command prompt, run::

    easy_install ffnet h5py sphinx numpydoc

If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_. Again, the 32-bit installer is required, e.g.
CDF33_1-SETUP-32.EXE. The default installation directory is recommended to
help SpacePy find the library.

Unzip the SpacePy source documentation. Open a command prompt in the
resulting directory and run::

    python setup.py install
