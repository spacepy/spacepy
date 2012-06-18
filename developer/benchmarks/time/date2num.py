
from __future__ import division

import os
import datetime

from matplotlib.dates import date2num as mpl_stock_date2num
from mpl_dates import date2num as mpl_date2num
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
date_data = np.arange(1, max_size , dtype=np.float)
date_data = np.asarray([datetime.datetime(2000, 1, 1) + datetime.timedelta(minutes=v) for v in date_data])

def date2num_mpl(data):
    ans = mpl_stock_date2num(data)

def date2num_spt(data):
    tt = spt.Ticktock(data)
    ans = tt.JD - 1721424.5

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
