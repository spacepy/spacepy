
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
from benchmarker import Benchmarker

import matplotlib.pyplot as plt

lib = ctypes.CDLL('hypot.so')
lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long]
lib.hypot.restype = ctypes.c_double


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

ans = {}
ans['ctypes_t'] = []
ans['python_t'] = []
ans['numpy_t'] = []
ans['extension_t'] = []

for loop in tb.logspace(3, 1e7, 30):
    print "loop", loop
    for bm in Benchmarker(width=25, cycle=5, extra=1):
        bm.run(ctypes_t, loop)
        bm.run(python_t, loop)
        bm.run(numpy_t, loop)
        bm.run(extension_t, loop)
    for result in bm.results:
        ans[result.label].append((loop, result.real))

for key in ans:
    plt.loglog(zip(*ans[key])[0], zip(*ans[key])[1],'o-',  label=key[:-2], lw=2 )

rep = dulwich.repo.Repo(os.path.abspath('../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.xlabel('hypot points')
plt.ylabel('Real execution time')
plt.legend(loc='upper left', shadow=True, fancybox=True)
plt.savefig('hypot_bench.png')

