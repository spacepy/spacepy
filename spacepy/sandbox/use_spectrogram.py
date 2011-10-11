
import datetime

import numpy as np

import spacepy.toolbox as tb
import spacepy.pycdf as pycdf
from spacepy.datamodel import SpaceData, dmarray

cdf = pycdf.CDF('poes_n15_short.cdf')

plotTime = np.asarray([datetime.datetime(2000, 11, 1), datetime.datetime(2000, 12, 1)])

rng = tb.tOverlapHalf(plotTime, cdf['EPOCH'], presort=True)
ind1 = rng[0]
ind2 = rng[-1]

data = SpaceData()
data['Epoch'] = dmarray(cdf['EPOCH'][ind1:ind2])
data['lValue'] = dmarray(cdf['lValue'][ind1:ind2])
data['mepOmni'] = dmarray(cdf['mepOmni'][ind1:ind2])[:,0]
#data['Epoch'] = cdf['EPOCH'][ind1:ind2]
#data['lValue'] = cdf['lValue'][ind1:ind2]


cdf.close()

avgt = datetime.timedelta(minutes=300)              #time scale to average over in minutes
Tbins = np.asarray([plotTime[0] + stp * avgt for
                    stp in xrange(np.long(np.ceil(plotTime.ptp().total_seconds()/avgt.total_seconds())))])

minL = 1                #min L value to plot (>1)
maxL = 10               #max L value to plot (<20)
step = 0.25             #L step size to use (.1, .2, .25, .5, 1, 2, ...), #smaller step sizes take considerably
Lbins = np.arange(minL, maxL, step)

data['Epoch'].attrs['bins'] = Tbins
data['lValue'].attrs['bins'] = Lbins

kwargs = {}

kwargs['bins'] = [Tbins, Lbins]
kwargs['ylim'] = [Lbins[0], Lbins[-1]]
kwargs['variables'] = ['Epoch', 'lValue', 'mepOmni']
#kwargs['zlim'] = [0, np.max(data['mepOmni'])]

from spacepy.plot import spectrogram

a = spectrogram(data, **kwargs)
a.plot()


#from pylab import *
#figure()
#bb = np.ma.masked_less_equal(a['spectrogram']['spectrogram'], 0)
#pcolormesh(a['spectrogram']['xedges'], a['spectrogram']['yedges'], np.ma.log10(bb))
