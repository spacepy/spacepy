
from __future__ import division

import os
import ctypes
import datetime

from matplotlib.dates import date2num, num2date
import numpy as np
import spacepy.toolbox as tb
import spacepy.time as spt
import matplotlib.pyplot as plt

# this is a git access module useful in here.
import dulwich
from benchmarker import Benchmarker



#==============================================================================
# Compare them all
#==============================================================================
max_size = 1e6
date_data = np.arange(1, max_size , dtype=ctypes.c_double)
date_data = np.asarray([datetime.datetime(2000, 1, 1) + datetime.timedelta(minutes=v) for v in date_data])

def date2num_mpl(data):
    ans = date2num(data)

def date2num_spt(data):
    ans = spt.date2num(data)

ans = {}
ans['date2num_mpl'] = []
ans['date2num_spt'] = []

for loop in tb.logspace(1, max_size, 20):
    print "loop", loop
    for bm in Benchmarker(width=25, cycle=5, extra=1):
        data = date_data[:loop]
        bm.run(date2num_mpl, data)
        bm.run(date2num_spt, data)
    for result in bm.results:
        ans[result.label].append((loop, result.real))

#==============================================================================
# plot up the times
#==============================================================================
plt.figure()
for key in ans:
    plt.loglog(zip(*ans[key])[0], zip(*ans[key])[1],'o-',  label=key[:-2], lw=2 )

rep = dulwich.repo.Repo(os.path.abspath('../../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.xlabel('array size')
plt.ylabel('Real execution time')
plt.legend(loc='upper left', shadow=True, fancybox=True)
plt.savefig('date2num_bench' + title + '.png')


#==============================================================================
# Compare them all
#==============================================================================
max_size = 1e6
date_data = np.arange(1, max_size , dtype=ctypes.c_double) + 730120.0

def num2date_mpl(data):
    ans = num2date(data)

def num2date_spt(data):
    ans = spt.num2date(data)

ans = {}
ans['num2date_mpl'] = []
ans['num2date_spt'] = []

for loop in tb.logspace(1, max_size, 20):
    print "loop", loop
    for bm in Benchmarker(width=25, cycle=5, extra=1):
        data = date_data[:loop]
        bm.run(num2date_mpl, data)
        bm.run(num2date_spt, data)
    for result in bm.results:
        ans[result.label].append((loop, result.real))

#==============================================================================
# plot up the times
#==============================================================================
plt.figure()
for key in ans:
    plt.loglog(zip(*ans[key])[0], zip(*ans[key])[1],'o-',  label=key[:-2], lw=2 )

rep = dulwich.repo.Repo(os.path.abspath('../../../'))
refs = rep.get_refs()
title = refs['HEAD']
plt.title(title)
plt.xlabel('array size')
plt.ylabel('Real execution time')
plt.legend(loc='upper left', shadow=True, fancybox=True)
plt.savefig('num2date_bench' + title + '.png')
