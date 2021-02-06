==========
Unit tests
==========


.. contents::
   :local:

The spacepy_testing module
==========================

.. module:: spacepy_testing

The ``spacepy_testing`` module contains utilities for assistance with
testing SpacePy. It is not installed as part of SpacePy and is thus
only importable from the unit test scripts themselves (which are in the
same directory). All unit test scripts import this module, and should
do so before importing any spacepy modules to ensure pathing is correct.

On import, :func:`~spacepy_testing.add_build_to_path` is run so that
the ``build`` directory is added to the Python search path. This means
the tests run against the latest build, not the installed version. Remove
the build directory to run against the installed version instead. The build
directory does not completely separate out Python versions, so removing
the build directory (and rebuilding) is recommended when switching Python
versions.

Build the package before running the tests::

  python setup.py build

.. contents::
   :local:

Classes
-------

|

.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    assertWarns
    assertDoesntWarn

Functions
---------

|

.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    add_build_to_path

Data
----

|

.. autosummary::
   :toctree: autosummary

   datadir
   testsdir

--------------------------

:Release: |version|
:Doc generation date: |today|
