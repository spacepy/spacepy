******************
Installing SpacePy
******************

The simplest way from zero (no Python) to a working SpacePy setup is:

  1. Install the `Anaconda <https://docs.anaconda.com/anaconda/>`_ Python
     environment. 64-bit is recommended.
  2. ``pip install --upgrade spacepy``

If you already have a working Python setup, install SpacePy by:

  1. ``pip install --upgrade spacepy``

In most cases this will install a binary build of SpacePy which
includes runtime dependencies. Otherwise it will attempt to compile
from source.

If you are familiar with installing Python packages, have particular
preferences for managing an installation, or if the above doesn't
work, refer to platform-specific instructions and the details
below.

Depending on your Python environment, you may need to explicitly
specify Python 3 throughout these commands, e.g. ``pip3`` instead of
``pip``.

The first time a user imports SpacePy, it automatically creates the
:doc:`configuration directory <configuration>`.

If you need further assistance, you can `check our support discussions
<https://github.com/spacepy/spacepy/discussions/categories/support-and-q-a>`_
and `start a new discussion
<https://github.com/spacepy/spacepy/discussions/new?category=support-and-q-a>`_
if there are no discussions on your topic of interest.

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

    pip install .

or, to install for all users (not in a virtual environment)::

    sudo pip install .

or, to install for a single user (not in a virtual environment)::

    pip install . --user

If you do not have administrative privileges, or you will be
developing for SpacePy, we strongly recommend using virtual
environments.

To install in custom location, e.g.::

    pip install . --prefix /n/packages/lib/python

(See :ref:`\\\\\\-\\\\\\-prefix <install_--prefix>` documentation
and the related :option:`--target`, :option:`--root`).

The closest analogy to the old ``setup.py`` handling (no dependency
handling, no isolated build) is::

    pip install . --no-build-isolation --no-deps

This is recommended if dependencies are managed by the OS or manually.

If installing into the system Python version on Linux (and potentially
some other cases), you will need to pass the flag
:option:`--break-system-packages`. Despite the frightening name this is
usually safe, although in this case is recommended to use
:option:`--no-build-isolation` :option:`--no-deps` and manage the
dependencies using the system package manager or conda.

Documentation
=============
If you want to build the documentation yourself (rather than using the
online documentation), install sphinx and numpydoc. The easiest way is
via pip::

  pip install sphinx numpydoc

They are also available via conda::

  conda install sphinx numpydoc

Compiling
=========
With the dependencies installed, SpacePy can be built from source.
You can always get the latest source code for SpacePy from our `github
repository <https://github.com/spacepy/spacepy>`_ and the latest
release from `PyPI <https://pypi.org/project/SpacePy/#files>`__

Following the instructions above will compile before the installation,
if installing from source or a binary installer is not available. If
this fails, specify a Fortran compiler::

    pip install . --config-setting="--build-option=--fcompiler=gnu95"

The supported compiler is ``gnu95`` (the GNU gfortran compiler); ``none``
can be specified as a "compiler" to skip all Fortran. You can also specify
the full path to the Fortran 77 compiler with ``--f77exec`` and to the
Fortran 90 compiler with ``--f90exec``::

    pip install . --config-setting="--build-option=--fcompiler=gnu95" --config-setting="--build-option=--f77exec=/usr/bin/gfortran" --config-setting="--build-option=--f90exec=/usr/bin/gfortran"


Troubleshooting
===============

.. _install_pip_failures:

pip failures
------------
If ``pip`` completely fails to build, a common issue is a failure in
the isolated build environment that ``pip`` sets up. Usually this can
be addressed by installing numpy first and eschewing the separate
build environment::

  pip install numpy
  pip install spacepy --no-build-isolation

Another option is manually installing all dependencies (via ``pip``,
``conda``, or other means) and then installing from an unpacked source
release with::

  pip install --no-build-isolation --no-deps .

``pip`` suppresses detailed output from the build process. To
troubleshoot a failure to install, it is useful to write this detailed
output to a file using the ``--log`` option, e.g.::

  pip install spacepy --log=install_log.txt

Please include this log file if opening an issue related to installation.

``pip`` will also cache packages; unfortunately sometimes it will use
a cached package which is incompatible with the current
environment. In that case, try clearing the cache first, so all
locally-compiled packages are rebuilt::

  pip cache purge

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
SDK version has not been installed. This may also result from
:ref:`pip caching <install_pip_failures>`.

The error ``RuntimeError: module compiled against API version 0x10 but
this version of numpy is 0xe`` followed by ``ImportError:
numpy.core.multiarray failed to import`` means that the version of
numpy used at installation of SpacePy does not match that used at
runtime. Check that there is only one version of numpy installed. In
some cases ``pip`` will install another version of numpy to support
the build; try installing numpy separately first, and then using the
``--no-build-isolation`` flag to ``pip``.
