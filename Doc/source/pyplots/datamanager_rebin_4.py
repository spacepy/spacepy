#!/usr/bin/env python

from datamanager_rebin_common import *

times = numpy.arange(0., 7200, 300)
alpha = pa(numpy.arange(8).reshape(1, -1), times.reshape(-1, 1))
energies = numpy.logspace(0, 3, 10)
flux = j(numpy.reshape(energies, (1, 1, -1)),
                     numpy.reshape(times, (-1, 1, 1)),
                     numpy.expand_dims(alpha, -1))
spacepy.plot.simpleSpectrogram(times, numpy.arange(8),
     flux[..., 0], cb=False, ylog=False)
matplotlib.pyplot.xlim(0, 7200)
matplotlib.pyplot.ylim(-0.5, 7.5)
matplotlib.pyplot.ylabel('Detector')
matplotlib.pyplot.xlabel('Time (sec)')
matplotlib.pyplot.title('Flux at 1 MeV')
matplotlib.pyplot.show()
