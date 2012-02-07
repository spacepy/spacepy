

import matplotlib.pyplot as plt 
import timeit

import numpy as np
import spacepy.toolbox as tb

# this is a git access module useful in here.
import dulwich

ans = {}
for val in range(1, 9):
    tb_t = timeit.Timer("tb.linspace(1.5, 10.5, " + str(10**(val/2)) + ")", setup="import spacepy.toolbox as tb")
    np_t = timeit.Timer("np.linspace(1.5, 10.5, " + str(10**(val/2)) + ")", setup="import numpy as np")
    val = [np_t.timeit(10000)/tb_t.timeit(10000) for v in range(50)]    
    ans[10**(val/2)] = [np.mean(val), np.std(val)]
    
x = np.array(ans.keys())
x.sort()
y = np.asarray([ans[v][0] for v in x])
yerr = np.asarray([ans[v][1] for v in x])
plt.errorbar(x, y, yerr, color='b')
plt.xscale('log')
plt.xtitle('linspace points')
plt.ytitle('tb.linspace/np.linspace')
# get the title

