#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Toolbox of various functions and generic utilities.

Authors: Steve Morley, Jon Niehof, Brian Larsen, Josef Koller, Dan Welling
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov, jniehof@lanl.gov, balarsen@lanl.gov, jkoller@lanl.gov, dwelling@lanl.gov
Los Alamos National Laboratory

Copyright 2010-2012 Los Alamos National Security, LLC.
"""
#from __future__ import absolute_import

from .toolbox import *
try:
    from _toolbox import hypot
except ImportError:
    #probably py3k, use slow python versions for now
    import numpy
    def hypot(*args):
        """
        Compute sqrt(vals[0] **2 + vals[1] **2 ...), i.e.
        n-dimensional hypoteneuse
        """
        return numpy.sqrt(numpy.inner(args, args))


__contact__ = 'Brian Larsen, balarsen@lanl.gov'
