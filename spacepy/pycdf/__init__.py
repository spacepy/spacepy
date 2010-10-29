#!/usr/bin/env python

"""pyCDF: Python interface to NASA's Common Data Format library

This package provides a Python interface to the Common Data Format (CDF)
library used for many NASA missions, available at U{http://cdf.gsfc.nasa.gov/}.
It is targeted at Python 2.6+ and should work without change on either
Python 2 or Python 3.

The interface is intended to be quite 'pythonic' rather than reproducing the
C interface. To open or close a CDF and access its variables, see the L{CDF}
class documentation. Accessing data within the variables is via the L{Var}
class.

The base CDF library must be properly installed in order to use this package.
pycdf will search for the library in this order:
  1. A directory named by the environment variable CDF_LIB (which should be
     set if using the definitions file provided with the CDF library).
  2. A subdirectory C{lib} in a directory named by the environment variable
     CDF_BASE.
  3. The standard system library search path.
If pycdf has trouble finding the library, try setting CDF_LIB before importing
the module, e.g. if the library is in CDF/lib in the user's home directory::
  import os
  os.putenv("CDF_LIB", "~/CDF/lib")
  import pycdf

The L{const} module contains constants useful for accessing the underlying library.
"""

__version__ = "0.6"
__author__ = "Jonathan Niehof <jniehof@lanl.gov>"

from ._pycdf import *
