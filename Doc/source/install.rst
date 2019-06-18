******************
Installing SpacePy
******************

SpacePy uses the standard Python distutils and setuptools to compile and install.
For detailed, platform-specific installation instructions, see:

.. toctree::
    :maxdepth: 1

    dependencies
    install_linux
    install_mac
    install_windows

For a list of dependencies, see :doc:`dependencies`.

Following are generic instructions.

The latest stable release is provided via the Python Package Index (PyPI).
To install from PyPI, make sure you have pip installed::

    pip install --upgrade spacepy

If you are installing for a single user, and are not working in a virtual environment,
add the --user flag when installing with pip.

To build from source (available at `our github <https://github.com/spacepy/spacepy>`_
just run (from a virtual environment, such as a conda environment)::

    python setup.py install

or, to install for all users (not in a virtual environment)::

    sudo python setup.py install

or, to install for a single user (not in a virtual environment)::

    python setup.py install --user

If you do not have administrative privileges, or you will be developing for SpacePy,
we strongly recommend using virtual environments.

To install in custom location, e.g.::

    python setup.py install --home=/n/packages/lib/python

It is also possible to select a specific compiler for installing the IRBEM-LIB library as part
of SpacePy. Currently the
following flags are supported: gnu95, gnu, pg. You can invoke these by using one of the
following commands below but not all of them are supported on all platforms:

* ``python setup.py install --fcompiler=pg``      #(will use pgi compiler suite)
* ``python setup.py install --fcompiler=gnu``    #(will use g77)
* ``python setup.py install --fcompiler=gnu95``   #(default option for using gfortran)
