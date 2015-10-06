******************
Installing SpacePy
******************

SpacePy uses the standard Python distutils to compile and install.
For detailed, platform-specific installation instructions, see:

.. toctree::
    :maxdepth: 1

    dependencies
    install_linux
    install_mac
    install_windows

For a list of dependencies, see :doc:`dependencies`.

Following are generic instructions.

Option 1) to install it in a standard location (depending on your system)::

    python setup.py install

or::

    sudo python setup.py install

or::

    python setup.py install --user

If you do not have administrative privileges, or you will be developing for SpacePy,
the latter is recommended.

Option 2) to install in custom location, e.g.::

    python setup.py install --home=/n/packages/lib/python

It is also possible to select a specific compiler for installing the IRBEM-LIB library as part
of SpacePy. Currently the
following flags are supported: gnu95, gnu, pg. You can invoke these by using one of the
following commands below but not all of them are supported on all platforms:

* ``python setup.py install --fcompiler=pg``      #(will use pgi compiler suite)
* ``python setup.py install --fcompiler=gnu``    #(will use g77)
* ``python setup.py install --fcompiler=gnu95``   #(default option for using gfortran)
