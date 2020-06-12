#!/usr/bin/env python

import matplotlib.colors
import matplotlib.pyplot
import numpy
import spacepy.datamanager
import spacepy.toolbox

def j(e, t, a):
    return e ** -2 * (1 / (90 * numpy.sqrt(2 * numpy.pi))) \
        * numpy.exp(-0.5 * (
            (a - 90 + 90 * numpy.sin(t / 573.))
            / 90.) ** 2)

def pa(d, t):
    return (d * 22.5 + t / 6 * (2 * (d % 2) - 1)) % 180
