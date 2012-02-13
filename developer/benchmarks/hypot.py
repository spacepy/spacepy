
from __future__ import division

import matplotlib.pyplot as plt 
import timeit
import os
import ctypes
import math

import numpy as np
import spacepy.toolbox as tb

# this is a git access module useful in here.
import dulwich


lib = ctypes.CDLL('hypot.so')
lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long]
lib.hypot.restype = ctypes.c_double


from benchmarker import Benchmarker
#==============================================================================
# Compare them all
#==============================================================================
lib = ctypes.CDLL('hypot.so')
lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long]
lib.hypot.restype = ctypes.c_double

def ctypes_t(n):
    data = np.arange(1, n , dtype=ctypes.c_double)
    ans = lib.hypot(data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

def python_t(n):
    data = np.arange(1, n , dtype=ctypes.c_double)
    ans = math.sqrt(sum([v**2 for v in data]))

def numpy_t(n):
    data = np.arange(1, n , dtype=ctypes.c_double)
    ans = np.sqrt(np.inner(data, data))

def extension_t(n):
    data = np.arange(1, n , dtype=ctypes.c_double)
    ans = tb.hypot(data)

loop = 1000
for bm in Benchmarker(width=25, cycle=10, extra=1):
    bm.run(ctypes_t, loop)
    bm.run(python_t, loop)
    bm.run(numpy_t, loop)
    bm.run(extension_t, loop)




