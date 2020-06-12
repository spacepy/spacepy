#!/usr/bin/env python

from datamanager_rebin_common import *

times = numpy.arange(0., 7200, 5)
energies = numpy.logspace(0, 3, 50)
flux = j(numpy.expand_dims(energies, 0),
         numpy.expand_dims(times, 1), 90.)
matplotlib.pyplot.pcolormesh(
    spacepy.toolbox.bin_center_to_edges(times),
    spacepy.toolbox.bin_center_to_edges(energies),
    flux.transpose(), norm=matplotlib.colors.LogNorm())
matplotlib.pyplot.yscale('log')
matplotlib.pyplot.ylabel('Energy (MeV)')
matplotlib.pyplot.xlabel('Time (sec)')
matplotlib.pyplot.title('Flux at 90 degrees')
matplotlib.pyplot.show()
