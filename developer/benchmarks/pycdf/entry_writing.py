#!/usr/bin/env python

"""Write a bunch of zEntries to test for speed

The actual test case was for CHAR entries on FLOAT vars with types
specified, so that's why that's used here. That has 1492 variables with
18 attrs each.

32s initially. Down to 11.3 at this point.
"""

import os

import spacepy.pycdf


with spacepy.pycdf.CDF('foo.cdf', create=True) as f:
    for i in range(1500):
        f.new('VAR_{:03d}'.format(i), type=spacepy.pycdf.const.CDF_FLOAT)
    for i in range(1500):
        var = f['VAR_{:03d}'.format(i)]
        for j in range(20):
            var.attrs.new('ATTR_{:02d}'.format(j), 'foo',
                          type=spacepy.pycdf.const.CDF_CHAR)

os.remove('foo.cdf')
