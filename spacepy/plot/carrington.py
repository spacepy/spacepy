#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Module for plotting data by Carrington or Bartels rotation

Authors: Steve Morley
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov
Los Alamos National Laboratory

Copyright 2011-2015 Los Alamos National Security, LLC.

.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    solarRotationPlot
"""
import datetime as dt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
from matplotlib.colors import LogNorm

import spacepy.datamodel as dm
import spacepy.toolbox as tb
import spacepy.empiricals as emp
import spacepy.time as spt

__contact__ = 'Steve Morley, smorley@lanl.gov'

def solarRotationPlot(ticks, data, targ_ax=None, rtype='bartels', nbins=27):
    '''Plots a 1-D time series as a Carrington or Bartels plot

    '''
    SRlength = {'bartels': 27.0, 'carrington':27.2753}
    if not targ_ax:
        #setup axes
        fig, targ_ax = plt.subplots()

    #data are 1-D, get range of SR numbers for Y-axis
    min_sr = emp.getSolarRotation(ticks[0],  rtype=rtype)
    max_sr = emp.getSolarRotation(ticks[-1], rtype=rtype)
    n_sr = 1+max_sr-min_sr
    sr_values = np.linspace(min_sr, max_sr, n_sr)

    #get bin size in time
    start_time = emp.getSolarRotation(min_sr,  rtype=rtype, reverse=True)[0] #convert back from Solar Rotation to date, so we get the start time of the rotation
    end_time = emp.getSolarRotation(max_sr+1,  rtype=rtype, reverse=True)[0] #add 1 to last solar rotation number then convert back to date
    binned_times = spt.tickrange(start_time, end_time, SRlength[rtype.lower()]/nbins)

    #now bin data on new time grid -- TODO: use windowMean from toolbox??
    digital = np.digitize(ticks.RDT, binned_times.RDT)
    bin_means = np.zeros(len(binned_times)-1)
    bin_means.fill(np.nan)
    for i in range(1, len(binned_times)):
        try:
            bin_means[i] = data[digital == i].mean()
        except (ZeroDivisionError, IndexError): #ZDE happens when no data for that bin, IE is when bins extend beyond data
            pass

    plot_view = np.array(bin_means).reshape((n_sr, nbins))
    im = targ_ax.imshow(plot_view, interpolation='nearest', vmin=data.min(), vmax=data.max(), extent=[1,nbins, max_sr,min_sr], aspect=100)
    plt.colorbar(im, ax=targ_ax)
    targ_ax.set_xlabel('Day of Rotation')
    targ_ax.set_ylabel('{0} Rotation number'.format(rtype.title()))
    #targ_ax.pcolormesh(plot_view)
    
    return
