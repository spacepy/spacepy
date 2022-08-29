******************
MacOS Installation
******************

Unless otherwise noted, the commands in these instructions are run
from the MacOS terminal (command line).

Installation requires a working Python environment and compilers. The
two common ways to achieve this are via :ref:`conda <install_mac_conda>`
or via :ref:`MacPorts <install_mac_macports>`. As a weak recommendation
for choosing between them, conda may be better if your main focus is
running Python and you need conda's easy support for multiple Python
environments; MacPorts may be better if you want to use many of the other
open source tools provided in MacPorts.

Binary installers for SpacePy on Mac are in preparation for a future
release.

.. contents::
   :local:

.. _install_mac_conda:

Conda installation
==================
Our recommendation for `Anaconda
<https://docs.anaconda.com/anaconda/>`_ is to use Python 3 and 64-bit
binaries. Follow the directions to install conda and set up an
environment; a minimal setup can be had by downloading the latest
`miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_. Double-click the
miniconda pkg and run through the installation process, choosing
"install for me only". Then, at the terminal::

  source ~/opt/miniconda3/bin/activate

You will now be in an active conda environment. (Note: the
installation will modify your ``.zprofile`` file.)

Compiling under conda `requires the MacOS 10.9 SDK
<https://stackoverflow.com/questions/69236331/
conda-macos-big-sur-ld-unsupported-tapi-file-type-tapi-tbd-in-yaml-file/>`_
on Intel (x86_64) Macs and `11.0 on Apple Silicon
<https://conda-forge.org/blog/posts/2020-10-29-macos-arm64/>`_
(ARM/M1/M2). It can be downloaded `here
<https://github.com/phracker/MacOSX-SDKs/releases>`_ (choose
"MacOSX10.9.sdk.tar.xz" or "MacOSX11.0.sdk.tar.xz"). Uncompress it
into ``opt``, e.g.::

  sudo tar xf ~/Downloads/MacOSX10.9.sdk.tar.xz -C /opt

Install the Fortran compiler::

  conda install gfortran

If you do not have Xcode installed, you will be prompted with a
message like "The xcrun command requires the command line developer
tools." Accept the installation and allow it to finish before
continuing.

You may optionally install SpacePy dependencies via conda (otherwise
they will be installed via ``pip``)::

   conda install numpy scipy matplotlib h5py

Finally, install SpacePy::

  SDKROOT=/opt/MacOSX10.9.sdk pip install spacepy  # Intel
  SDKROOT=/opt/MacOSX11.0.sdk pip install spacepy  # ARM

If you're installing as a single user (not in a virtual environment) then
add the ``--user`` flag.

You will also need the :ref:`NASA CDF library <install_mac_cdf>` to use
:mod:`~spacepy.pycdf`.

.. _install_mac_macports:

MacPorts installation
=====================

You may install :ref:`the full Xcode suite <install_mac_xcode>` and
follow the `MacPorts guide <https://guide.macports.org/>`_; however,
these directions should suffice to install a working Python and
SpacePy.

Installing the Xcode command line tools is recommended before proceeding::

  xcode-select --install

Download the `MacPorts installer
<https://www.macports.org/install.php>`_ and double-click the pkg to
perform the installation. (Note this modifies your ``.zprofile``
environment file.)

Install Python and the needed compilers. You need to specify a
version; at this time, Python 3.9 and gcc 11 are reasonable choices::

  sudo port install gcc11  # Includes gfortran
  sudo port install python39
  sudo port install py39-pip
  # Installing the following is optional; pip will automatically install
  sudo port install py39-numpy py39-scipy py39-h5py py39-matplotlib

If you have not already installed the Xcode command line tools, you
will be prompted to do so. In that case, it is suggested to accept the
tools installation, and then quit the port command and restart once
the tools are installed.

To install via ``pip``, default versions of Python and gcc must be set::

  sudo port select --set python python39
  rehash #recalculate the pathing to not get system python
  sudo port select --set python3 python39
  sudo port select --set pip pip39
  sudo port select --set gcc mp-gcc11

Then you can install SpacePy::

  pip install spacepy

If you're installing as a single user (not in a virtual environment) then
add the ``--user`` flag.

You will also need the :ref:`NASA CDF library <install_mac_cdf>` to use
:mod:`~spacepy.pycdf`.

If you are installing from a source distribution, you can specify the
compiler at install time instead of using ``port select``::

  python3.9 setup.py install --fcompiler=gnu95 --f90exec=/opt/local/bin/gfortran-mp-11

.. _install_mac_cdf:

CDF
===

NASA provides `Mac binaries
<https://spdf.gsfc.nasa.gov/pub/software/cdf/dist/latest-release/macosx/>`_
of the CDF library. Download the file ending in ``binary_signed.pkg``
(e.g. ``CDF3_8_1-binary_signed.pkg``), double-click, and install per
the defaults.

.. _install_mac_xcode:

Xcode installation
==================
Installation of the full Xcode package is not required simply for
SpacePy; however, if you are interested in regular compiler use, it
may be useful. If you choose to install the full Xcode package,
perform these steps before installing conda or macports via the
directions above.

  * Create and log in to an Apple developer account at
    https://developer.apple.com/
  * Check the `Xcode release notes
    <https://developer.apple.com/documentation/xcode-release-notes/>`_
    to find the latest version of Xcode supported on your version of
    MacOS.
  * From the `more downloads
    <https://developer.apple.com/download/all/>`_ section of the Apple
    Developer site, search for and download that version of Xcode.
  * Double-click on the downloaded .xip file to open with the archive
    utility and extract the Xcode app.
  * Drag the resulting Xcode icon into Applications
  * From the `more downloads
    <https://developer.apple.com/download/all/>`_ section of the Apple
    Developer site, search for the Xcode command line tools for the
    same version of Xcode
  * Open the dmg file with the command line tools, open the resulting
    mounted disk image, and double-click the pkg file to install.

Proceed with the installation of conda or MacPorts and SpacePy
