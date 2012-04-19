#!/usr/bin/env python
"""Compare various means of vectorizing epoch conversions"""

import datetime
import timeit

import numpy
from spacepy import pycdf

n_epochs = 10000
n_iter = 5

startepoch = pycdf.lib.datetime_to_epoch(datetime.datetime.now())
epochs = numpy.arange(startepoch, startepoch + n_epochs)
command = 'vepoch_to_datetime(epochs)'
setup = 'from __main__ import epochs, vepoch_to_datetime'

print('Tests of EPOCH')
#case 1: vectorize
vepoch_to_datetime = numpy.vectorize(pycdf.lib.epoch_to_datetime)
timing = timeit.timeit(command, setup, number=n_iter)
print('vectorize: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))

#case 2: ufunc
vepoch_to_datetime = numpy.frompyfunc(pycdf.lib.epoch_to_datetime, 1, 1)
timing = timeit.timeit(command, setup, number=n_iter)
print('ufunc: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))

#case 3: ndindex
def vepoch_to_datetime(inarray):
    outarray = numpy.empty(inarray.shape, dtype='O')
    for idx in numpy.ndindex(inarray.shape):
        outarray[idx] == pycdf.lib.epoch_to_datetime(inarray[idx])
    return outarray
timing = timeit.timeit(command, setup, number=n_iter)
print('ndindex: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))

#case 4: ndenumerate
def vepoch_to_datetime(inarray):
    outarray = numpy.empty(inarray.shape, dtype='O')
    for idx, x in numpy.ndenumerate(inarray):
        outarray[idx] == pycdf.lib.epoch_to_datetime(x)
    return outarray
timing = timeit.timeit(command, setup, number=n_iter)
print('ndenumerate: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))


##EPOCH16
print('')
print('Tests of EPOCH16')
dt = datetime.datetime.now()
epochs = numpy.empty((n_epochs, 2), dtype='float64')
for i in range(n_epochs):
    epochs[i] = numpy.array(
        pycdf.lib.datetime_to_epoch16(dt +
                                      datetime.timedelta(seconds=i)))
setup = 'from __main__ import epochs, vepoch16_to_datetime'

#case 1: vectorize
vepoch16_to_datetime = numpy.vectorize(pycdf.lib.epoch16_to_datetime)
command = 'vepoch16_to_datetime(epochs[:,0], epochs[:,1])'
timing = timeit.timeit(command, setup, number=n_iter)
print('vectorize: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))

#case 2: ufunc
vepoch16_to_datetime = numpy.frompyfunc(pycdf.lib.epoch16_to_datetime, 2, 1)
command = 'vepoch16_to_datetime(epochs[:,0], epochs[:,1])'
timing = timeit.timeit(command, setup, number=n_iter)
print('ufunc: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))

#case 3: ndindex
def vepoch16_to_datetime(inarray):
    outarray = numpy.empty(inarray.shape[0:-1], dtype='O')
    for idx in numpy.ndindex(outarray.shape):
        outarray[idx] = pycdf.lib.epoch16_to_datetime(*inarray[idx])
    return outarray
command = 'vepoch16_to_datetime(epochs)'
timing = timeit.timeit(command, setup, number=n_iter)
print('ndindex: {0} loops of {1} took {2} seconds.'.format(
    n_iter, n_epochs, timing))
