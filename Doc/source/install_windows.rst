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

Download the `Python(x,y) distribution 
<http://code.google.com/p/pythonxy/>`_. The standard plugins include
most of the dependencies required for SpacePy. Like SpacePy, this currently
only supports 32-bit code, but will run on 64-bit machines.

If you wish to use the lanlstar module, you will need `ffnet-0.7.1.win32-py2.7.exe <https://sourceforge.net/projects/ffnet/files/ffnet/ffnet-0.7.1/>`_.

If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_. Again, the 32-bit installer is required, e.g.
CDF34_1_0-SETUP-32.EXE. The default installation directory is recommended to
help SpacePy find the library.

Finally, install SpacePy using the installer EXE. Be sure to choose the
installer that matches your version of Python, either 2.6 or 2.7.


The Hard Way
============

Installing dependencies
-----------------------

This is a step-by-step guide to installing all of the spacepy dependencies
individiually. This is the easiest way to install all dependencies by hand;
if you are comfortable downloading, compiling, and installing Python modules
yourself, that will also work.

Download and install the 32-bit `Python Windows installer
<http://python.org/download/>`_ .  The 2.7 version is recommended since
most dependencies are now built against 2.7.

The following filenames are given for Python 2.7 and were the latest
as of this writing; download the appropriate file for your version of
Python and the latest version available of the required package. 

Download and install `numpy-1.7.1-win32-superpack-python2.7.exe
<http://sourceforge.net/projects/numpy/files/>`_.

Download and install `scipy-0.12.0-win32-superpack-python2.7.exe
<http://sourceforge.net/projects/scipy/files/>`_.

Download and install `matplotlib-1.2.1.win32-py2.7.exe
<http://matplotlib.org/downloads.html>`_.

Download and install `h5py-2.1.3.win32-py2.7.msi
<http://code.google.com/p/h5py/downloads/list>`_.

Download and install `networkx-1.7_py27.exe
<http://code.google.com/p/pythonxy/wiki/StandardPlugins>`_.

Download and install `ffnet-0.7.1.win32-py2.7.exe
<https://sourceforge.net/projects/ffnet/files/ffnet/ffnet-0.7.1/>`_.

If you wish to use CDF files, download and install the `NASA CDF library
<http://cdf.gsfc.nasa.gov/>`_. Again, the 32-bit installer is required, e.g.
CDF34_1_0-SETUP-32.EXE. The default installation directory is recommended to
help SpacePy find the library.

At this point, you can install SpacePy from the installer EXE, or continue
to build SpacePy from source.


Building and installing SpacePy
-------------------------------

Edit your Windows path: Right-click my computer, choose properties,
advanced, environment variables.  Edit the path variable: append
``;c:\python27;c:\python27\Scripts;c:\MinGW\bin`` (`more information
<http://docs.python.org/using/windows.html#finding-the-python-executable>`_).

Download the latest `mingwget <http://sourceforge.net/projects/mingw/files/Automated%20MinGW%20Installer/mingw-get/>`_. mingw32 is the only C and Fortran compiler supported by SpacePy on Windows. Unzip into C:\MinGW. At a command prompt (Start, Run, cmd) type::

      mingw-get install gcc g++ mingw32-make fortran

See also `the mingw docs <http://www.mingw.org/wiki/Getting_Started>`_

Create a file ``distutils.cfg`` in ``C:\Python27\Lib\distutils``
(change appropriately for Python 2.6). It should contain::

    [build]
    compiler=mingw32

Unzip the SpacePy source documentation. Open a command prompt in the
resulting directory and run::

    python setup.py install


Developers
==========

If you want to build the documentation yourself (rather than using the
documentation shipped with SpacePy), install sphinx and numpydoc. The
easiest way is via pip.

Download and run `distribute_setup.py
<https://pypi.python.org/pypi/distribute/#distribute-setup-py>`_.

Download and run `get-pip.py
<http://www.pip-installer.org/en/latest/installing.html#using-get-pip>`_.

Then::

    pip install sphinx numpydoc
