#!/usr/bin/env python

from datamanager_rebin_common import *

times = numpy.arange(0., 7200, 300)
alpha = pa(numpy.arange(8).reshape(1, -1), times.reshape(-1, 1))
energies = numpy.logspace(0, 3, 10)
flux = j(numpy.reshape(energies, (1, 1, -1)),
                     numpy.reshape(times, (-1, 1, 1)),
                     numpy.expand_dims(alpha, -1))
matplotlib.pyplot.pcolormesh(
    spacepy.toolbox.bin_center_to_edges(times),
    spacepy.toolbox.bin_center_to_edges(energies),
    flux[:, 0, :].transpose(), norm=matplotlib.colors.LogNorm())
matplotlib.pyplot.yscale('log')
matplotlib.pyplot.ylabel('Energy (MeV)')
matplotlib.pyplot.xlabel('Time (sec)')
matplotlib.pyplot.title('Flux in detector 0')
matplotlib.pyplot.show()
