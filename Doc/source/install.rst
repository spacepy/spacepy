******************
Installing SpacePy
******************

The simplest way from zero to a working SpacePy setup is:

  1. Install the `Anaconda <https://docs.anaconda.com/anaconda/>`_ Python
     environment. Python 3 is strongly recommended (64-bit is recommended).
  2. ``pip install --upgrade spacepy``

If you are familiar with installing Python packages, have particular
preferences for managing an installation, or if the above doesn't
work, refer to platform-specific instructions and the details
below. (In particular, if you will use :mod:`~spacepy.LANLstar` on
Windows, see :doc:`install_windows`.)

For installing the NASA CDF library to support :mod:`~spacepy.pycdf`,
see the platform-specific instructions linked below.

If you need further assistance, you can `open an issue
<https://github.com/spacepy/spacepy/issues>`_.

.. toctree::
    :maxdepth: 1

    dependencies
    install_linux
    install_mac
    install_windows

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
are on `github <https://github.com/spacepy/spacepy>`__.


After downloading and
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
