

import matplotlib.pyplot as plt 
import timeit
import os

import numpy as np
import spacepy.toolbox as tb

# this is a git access module useful in here.
import dulwich

ans = {}
for val in xrange(1, 9):
    tb_t = timeit.Timer("tb.linspace(1.5, 10.5, " + str(10**(val/2)) + ")", setup="import spacepy.toolbox as tb")
    np_t = timeit.Timer("np.linspace(1.5, 10.5, " + str(10**(val/2)) + ")", setup="import numpy as np")
    bnh = [np_t.timeit(10000)/tb_t.timeit(10000) for v in range(50)]    
    ans[10**(val/2)] = [np.mean(bnh), np.std(bnh)]
    tb.progressbar(val, 1, len(xrange(1, 9)), 'Progress')
    
x = np.array(ans.keys())
x.sort()
y = np.asarray([ans[v][0] for v in x])
yerr = np.asarray([ans[v][1] for v in x])
plt.errorbar(x, y, yerr, color='b', fmt='o-')
plt.xscale('log')
plt.xlabel('linspace points')
plt.ylabel('tb.linspace/np.linspace')
# get the title
rep = dulwich.repo.Repo(os.path.abspath('../../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.savefig('linspace_bench.png')


## CONCLUSION:
#   the toolbox version of linspace is faster than the numpy version for all tested range sizes
#   when the size is small it is ~16x an when it is large ~2x





