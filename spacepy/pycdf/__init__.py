#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This package provides a Python interface to the Common Data Format (CDF)
library used for many NASA missions, available at http://cdf.gsfc.nasa.gov/.
It is targeted at Python 2.6+ and should work without change on either
Python 2 or Python 3.

The interface is intended to be 'pythonic' rather than reproducing the
C interface. To open or close a CDF and access its variables, see the :class:`CDF`
class. Accessing data within the variables is via the :class:`Var`
class. The :data:`lib` object provides access to some routines
that affect the functionality of the library in general. The
:mod:`~spacepy.pycdf.const` module contains constants useful for accessing
the underlying library.


The CDF C library must be properly installed in order to use this package.
The CDF distribution provides scripts meant to be called in a user's
login scripts, ``definitions.B`` for bash and ``definitions.C`` for C-shell
derivatives. (See the installation instructions which come with the CDF library.)
These will set environment variables specifying the location
of the library; pycdf will respect these variables if they are set. Otherwise
it will search the standard system library path and the default installation
locations for the CDF library.

If pycdf has trouble finding the library, try setting ``CDF_LIB`` before importing
the module, e.g. if the library is in ``CDF/lib`` in the user's home directory:
    
>>> import os
>>> os.putenv("CDF_LIB", "~/CDF/lib")
>>> from spacepy import pycdf
    
If this works, make the environment setting permanent. Note that on OSX,
using plists to set the environment may not carry over to Python terminal
sessions; use ``.cshrc`` or ``.bashrc`` instead.

.. currentmodule:: spacepy.pycdf


Authors: Jon Niehof

Institution: Los Alamos National Laboratory

Contact: jniehof@lanl.gov


Copyright 2010-2012 Los Alamos National Security, LLC.

"""

from ._pycdf import *

__contact__ = 'Jon Niehof, jniehof@lanl.gov'
