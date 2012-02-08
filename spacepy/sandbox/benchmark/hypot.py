
from __future__ import division

import matplotlib.pyplot as plt 
import timeit
import os
import ctypes

import numpy as np
import spacepy.toolbox as tb

# this is a git access module useful in here.
import dulwich


lib = ctypes.CDLL('hypot.so')
lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long]
lib.hypot.restype = ctypes.c_double



ans = {}
for val in xrange(1, 6*3):
    ctypes_t = timeit.Timer("lib.hypot(data.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))", 
                            setup="import ctypes; lib = ctypes.CDLL('hypot.so'); " +
                            "lib.hypot_c.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_long];" +
                            "lib.hypot.restype = ctypes.c_double; " + 
                            "import numpy as np;" + 
                            "data = np.arange(1, " + str(10**(val//3)) + ", dtype=ctypes.c_double)")
    numpy_t = timeit.Timer("np.sqrt(np.sum(data*data))", 
                           setup="import numpy as np; import ctypes;" + 
                           "data = np.arange(1, " + str(10**(val//3)) + ", dtype=ctypes.c_double)")
    bnh = [numpy_t.timeit(1000)/ctypes_t.timeit(1000) for v in range(100)]    
    ans[10**(val//3)] = [np.mean(bnh), np.std(bnh)]
    tb.progressbar(val, 1, len(xrange(1, 6*3)), 'Progress')
    
x = np.array(ans.keys())
x.sort()
y = np.asarray([ans[v][0] for v in x])
yerr = np.asarray([ans[v][1] for v in x])
plt.errorbar(x, y, yerr, color='b', fmt='o-')
plt.xscale('log')
plt.yscale('log')
plt.axhline(0)
plt.xlabel('hypot points')
plt.ylabel('numpy/ctypes')
# get the title
rep = dulwich.repo.Repo(os.path.abspath('../../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.savefig('hypot_bench.png')


## CONCLUSION:





