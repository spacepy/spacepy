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
L{Library} provides details on how the library is found.

The L{const} module contains constants useful for accessing the underlying library.
"""

__version__ = "0.2"
__author__ = "Jonathan Niehof <jniehof@lanl.gov>"

from ._pycdf import *
