#!/usr/bin/env python

from datamanager_rebin_common import *

times = numpy.arange(0., 7200, 5)
alpha = numpy.arange(0, 181., 2)
# Add a dimension so the flux is a 2D array
flux = j(1., numpy.expand_dims(times, 1),
             numpy.expand_dims(alpha, 0))
spacepy.plot.simpleSpectrogram(times, alpha, flux, cb=False,
                               ylog=False)
matplotlib.pyplot.xlim(0, 7200)
matplotlib.pyplot.ylim(0, 180)
matplotlib.pyplot.ylabel('Pitch angle (deg)')
matplotlib.pyplot.xlabel('Time (sec)')
matplotlib.pyplot.title('Flux at 1 MeV')
matplotlib.pyplot.show()
