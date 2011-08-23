#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
pyCDF: Python interface to NASA's Common Data Format library

This package provides a Python interface to the Common Data Format (CDF)
library used for many NASA missions, available at http://cdf.gsfc.nasa.gov/.
It is targeted at Python 2.6+ and should work without change on either
Python 2 or Python 3.

The interface is intended to be quite 'pythonic' rather than reproducing the
C interface. To open or close a CDF and access its variables, see the :py:class:`pycdf.CDF`
class documentation. Accessing data within the variables is via the :py:class:`pycdf.Var`
class. The :py:class:`pycdf.lib` object (of type :py:class:`pycdf.lib`) provides access to some routines
that affect the functionality of the library in general.

The base CDF library must be properly installed in order to use this package.
pycdf will search for the library in this order:

    1. A directory named by the environment variable CDF_LIB (which should be set if using the definitions file provided with the CDF library).
    2. A subdirectory C{lib} in a directory named by the environment variable CDF_BASE.
    3. The standard system library search path.

If pycdf has trouble finding the library, try setting CDF_LIB before importing
the module, e.g. if the library is in CDF/lib in the user's home directory::
  import os
  os.putenv("CDF_LIB", "~/CDF/lib")
  import pycdf

.. currentmodule:: spacepy.pycdf

The :py:mod:`pycdf.const` module contains constants useful for accessing the underlying library.

.. NOTE... there is an error with this reference

Authors: Jon Niehof
Institution: Los Alamos National Laboratory
Contact: jniehof@lanl.gov


Copyright ©2010 Los Alamos National Security, LLC.

"""

from ._pycdf import *
