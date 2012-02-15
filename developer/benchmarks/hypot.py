
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


#==============================================================================
# Compare them all
#==============================================================================
max_size = 1e6
data = np.arange(1, max_size , dtype=ctypes.c_double)

lib = ctypes.CDLL('hypot.so')
lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long]
lib.hypot_c.restype = ctypes.c_double

def ctypes_t(n, data):
#    ans = lib.hypot_c(data[0:n].ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
#                      len(data[0:n]))
    ans = lib.hypot_c(data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                      len(data))

def python_t(n, data):
#    ans = math.sqrt(sum([v**2 for v in data[0:n]]))
    ans = math.sqrt(sum([v**2 for v in data]))

def numpy_t(n, data):
#    ans = np.sqrt(np.inner(data[0:n], data[0:n]))
    ans = np.sqrt(np.inner(data, data))

def extension_t(n, data):
#    ans = tb.hypot(data[0:n])
    ans = tb.hypot(data)

ans = {}
ans['ctypes_t'] = []
ans['python_t'] = []
ans['numpy_t'] = []
ans['extension_t'] = []

for loop in tb.logspace(3, max_size, 20):
    print "loop", loop
    for bm in Benchmarker(width=25, cycle=5, extra=1):
        data = np.arange(1, loop, dtype=ctypes.c_double)
        bm.run(ctypes_t, loop, data)
        bm.run(python_t, loop, data)
        bm.run(numpy_t, loop, data)
        bm.run(extension_t, loop, data)
    for result in bm.results:
        ans[result.label].append((loop, result.real))

#==============================================================================
# plot up the times
#==============================================================================
plt.figure()
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

#==============================================================================
# plot up the ratios
#==============================================================================
plt.figure()
for key in ans:
    plt.loglog(np.asarray(zip(*ans[key])[0]), np.asarray(zip(*ans['python_t'])[1])/np.asarray(zip(*ans[key])[1]),'o-',  label=key[:-2], lw=2 )

rep = dulwich.repo.Repo(os.path.abspath('../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.xlabel('hypot points')
plt.ylabel('Speedup factor')
plt.legend(loc='upper left', shadow=True, fancybox=True)
plt.savefig('hypot_bench_ratio.png')


