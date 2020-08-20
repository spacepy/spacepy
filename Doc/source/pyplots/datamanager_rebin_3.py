#!/usr/bin/env python

from datamanager_rebin_common import *

times = numpy.arange(0., 7200, 5)
lines = matplotlib.pyplot.plot(
    times, pa(numpy.arange(8).reshape(1, -1), times.reshape(-1, 1)),
    marker='o', ms=1, linestyle='')
matplotlib.pyplot.legend(lines, ['Detector {}'.format(i) for i in range(8)],
                         loc='best')
matplotlib.pyplot.xlim(0, 7200)
matplotlib.pyplot.ylim(0, 180)
matplotlib.pyplot.xlabel('Time (sec)')
matplotlib.pyplot.ylabel('Pitch angle (deg)')
matplotlib.pyplot.title('Measured pitch angle by detector')
matplotlib.pyplot.show()
