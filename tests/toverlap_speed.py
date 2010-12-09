#!/usr/bin/env python
"""Tests for the speed of tOverlap"""

import datetime
import random
import timeit
from spacepy import toolbox

n_iter = 500

#Very simple call to make sure everything loaded/init before timing
dt_a = [datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 2)]
dt_b = [datetime.datetime(1999, 12, 31), datetime.datetime(2000, 1, 1)]
toolbox.tOverlap(dt_a, dt_b, presort=False)
toolbox.tOverlap(dt_a, dt_b, presort=True)

dt1 = datetime.datetime(2000, 11, 12)
dt_a = [dt1 + datetime.timedelta(hours=offset) for offset in range(10000)]
dt_b = [dt1 + datetime.timedelta(hours=offset)
        for offset in range(-5000, 5000)]
timing = timeit.timeit('toolbox.tOverlap(dt_a, dt_b)',
                       'from __main__ import toolbox, dt_a, dt_b',
                       number=n_iter)
print('Sorted arrays took ' + str(timing) + ' seconds.')

timing = timeit.timeit('toolbox.tOverlap(dt_a, dt_b, presort=True)',
                       'from __main__ import toolbox, dt_a, dt_b',
                       number=n_iter)
print('Sorted arrays, sorted algorithm took ' + str(timing) + ' seconds.')

random.seed(0)
random.shuffle(dt_a)
random.shuffle(dt_b)
timing = timeit.timeit('toolbox.tOverlap(dt_a, dt_b)',
                       'from __main__ import toolbox, dt_a, dt_b',
                       number=n_iter)
print('Scrambled arrays took ' + str(timing) + ' seconds.')



