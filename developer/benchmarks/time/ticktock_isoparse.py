#!/usr/bin/env python
"""Speed test for ISO parsing in ticktock"""

import datetime
import time

import dateutil.parser
import numpy
import spacepy.time
import spacepy.toolbox

t1 = time.time()
dt = spacepy.toolbox.linspace(datetime.datetime(2012, 1, 1),
                              datetime.datetime(2012, 1, 2), 140000)
#dts = numpy.asarray([v.isoformat() for v in dt])
dts = numpy.asarray([v.strftime('%Y-%m-%dT%H:%M:%S.%f') for v in dt])
print(time.time() - t1)

t1 = time.time()
t = spacepy.time.Ticktock(dts, 'ISO').UTC
print(time.time() - t1)

t1 = time.time()
t = numpy.asarray([dateutil.parser.parse(v) for v in dts])
print(time.time() - t1)

t1 = time.time()
t = spacepy.time.Ticktock(dts, 'ISO', isoformat='%Y-%m-%dT%H:%M:%S.%f').UTC
print(time.time() - t1)

