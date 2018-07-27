#!/usr/bin/env python

"""Test assigning to zEntries and type guessing...hoping to speed this up"""

import itertools
import os
import os.path
import time

from spacepy import pycdf

if os.path.exists('foo.cdf'):
    os.remove('foo.cdf')

with pycdf.CDF('foo.cdf', '') as f:
    for i in range(100):
        f.new('DOUBLE_{0}'.format(i), type=pycdf.const.CDF_DOUBLE)
        f.new('INT_{0}'.format(i), type=pycdf.const.CDF_INT4)
        f.new('STR_{0}'.format(i), type=pycdf.const.CDF_CHAR,
                      n_elements=10)
    for i in range(100):
        for aname in ['{0}_{1}'.format(*p) for p in itertools.product(
                ('VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'),
                range(10))]:
            f['DOUBLE_{0}'.format(i)].attrs.new(
                aname, type=pycdf.const.CDF_DOUBLE)
            f['INT_{0}'.format(i)].attrs.new(
                aname, type=pycdf.const.CDF_INT4)
            f['STR_{0}'.format(i)].attrs.new(
                aname, type=pycdf.const.CDF_CHAR)
        for aname in ('DEPEND_0', 'DEPEND_1', 'DEPEND_2'):
            for varname in ('DOUBLE', 'INT', 'STR'):
                f['{0}_{1}'.format(varname, i)].attrs.new(
                    aname, type=pycdf.const.CDF_CHAR)
    t = time.time()
    f['STR_0'].attrs['stuff'] = 1
    print(time.time() - t)
