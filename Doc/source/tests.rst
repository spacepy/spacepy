==========
Unit tests
==========


.. contents::
   :local:

Unit tests are in the ``tests`` directory. Individual scripts can be
run from this directory, or ``test_all.py`` will run all tests in all
scripts.

spacepy must be installed before running the tests. To avoid affecting an
installed version while testing, use a separate :mod:`venv` or conda
environment, or install in a custom location using
:ref:`\\\\\\-\\\\\\-prefix <install_--prefix>` and manually edit
:envvar:`PYTHONPATH`.

Using an :ref:`editable install <install_--editable>` may be useful to save
time with repeated installs while editing.

The spacepy_testing module
==========================

.. module:: spacepy_testing

The ``spacepy_testing`` module contains utilities for assistance with
testing SpacePy. It is not installed as part of SpacePy and is thus
only importable from the unit test scripts themselves (which are in the
same directory). All unit test scripts import this module.

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
