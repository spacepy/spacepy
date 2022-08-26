******************
Installing SpacePy
******************

The simplest way from zero (no Python) to a working SpacePy setup is:

  1. Install the `Anaconda <https://docs.anaconda.com/anaconda/>`_ Python
     environment. Python 3 is strongly recommended (64-bit is recommended).
  2. ``pip install --upgrade spacepy``

If you already have a working Python setup, install SpacePy by:

  1. ``pip install --upgrade numpy``
  2. ``pip install --upgrade spacepy``

This will install a binary build of SpacePy if available (currently
only on Windows), otherwise it will attempt to compile. It will also
install most dependencies.

If you are familiar with installing Python packages, have particular
preferences for managing an installation, or if the above doesn't
work, refer to platform-specific instructions and the details
below.

For installing the NASA CDF library to support :mod:`~spacepy.pycdf`,
see the platform-specific instructions linked below.

The first time a user imports SpacePy, it automatically creates the
:doc:`configuration directory <configuration>`.

If you need further assistance, you can `open an issue
<https://github.com/spacepy/spacepy/issues>`_.

.. toctree::
    :maxdepth: 1

    dependencies
    install_linux
    install_mac
    install_windows

.. contents::
   :local:

SpacePy installs with the common Python distutils and pip.

The latest stable release is provided via `PyPI
<https://pypi.org/project/SpacePy/>`_ To install from PyPI, make sure
you have pip installed::

    pip install --upgrade spacepy

If you are installing for a single user, and are not working in a
virtual environment, add the ``--user`` flag when installing with pip.

Source releases are available from `PyPI
<https://pypi.org/project/SpacePy/#files>`__ and `our github
<https://github.com/spacepy/spacepy/releases>`__. Development versions
are on `github <https://github.com/spacepy/spacepy>`__. In addition to
downloading tarballs, the development version can be directly installed
with::

  pip install git+https://github.com/spacepy/spacepy

For source releases, after downloading and
unpacking, run (a virtual environment, such as a conda environment, is
recommended)::

    python setup.py install

or, to install for all users (not in a virtual environment)::

    sudo python setup.py install

or, to install for a single user (not in a virtual environment)::

    python setup.py install --user

If you do not have administrative privileges, or you will be
developing for SpacePy, we strongly recommend using virtual
environments.

To install in custom location, e.g.::

    python setup.py install --home=/n/packages/lib/python

Installs using ``setup.py`` do not require setuptools.

Troubleshooting
===============

irbempy
-------
The most common failures relate to compilation of the IRBEM
library. Unfortunately ``pip`` will hide these warnings, so they
manifest when running ``import spacepy.irbempy`` (or some other
component of SpacePy that uses irbempy).

The error ``ImportError: cannot import name 'irbempylib' from
partially initialized module 'spacepy.irbempy' (most likely due to a
circular import)`` means the IRBEM library did not compile at
all. This is most likely a compiler issue: either there is no Fortran
compiler, or, when using conda on :doc:`Mac <install_mac>`, the correct
SDK version has not been installed.

The error ``RuntimeError: module compiled against API version 0x10 but
this version of numpy is 0xe`` followed by ``ImportError:
numpy.core.multiarray failed to import`` means that the version of
numpy used at installation of SpacePy does not match that used at
runtime. Check that there is only one version of numpy installed. In
some cases ``pip`` will install another version of number to support
the build; this was fixed in SpacePy 0.4.0 but may recur.
