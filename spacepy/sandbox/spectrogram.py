#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
code to create and plot generic spectrograms for space science - these are different
than the standard spectrogram in other fields in that there is no fft used


Authors: Brian Larsen and Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov, smorley@lanl.gov
Los Alamos National Laboratory

Copyright 2011 Los Alamos National Security, LLC.
"""

import datetime

import numpy as np
import matplotlib as mpl
from matplotlib.dates import date2num, num2date

import spacepy.datamodel as dm
import spacepy.toolbox as tb

class spectrogram(dm.SpaceData):
    """
    This class generates and then contains the data binned into he spectrogram

    It is meant to be used on arbitrary data series.  The first series "x" is
    plotted on the abscissa and second series "y" is plotted on the ordinate and
    the third series "z" is plotted in color.

    The series are not passed in independently but instead inside a
    spacepy.datamodel.SpaceData container.  Helper routines are provided to
    facilitate the creation of the SpaceData container if the data are not in the format.
    """

    ## NOTE this will need to set the sphinx var autoclass_content to "both"

    def __init__(self, data, **kwargs):
        """
        Parameters
        ==========
        data : spacepy.datamodel.SpaceData
            The data for the spectrogram, the variables to be used default to
            "Epoch" for x, "Energy" for y, and "Flux" for z.  Other names are
            specified using the 'variables' keyword.  All keywords override .attrs
            contents.

        Other Parameters
        ================
        variables : list
            keyword containing the names of the variables to use for the spectrogram
            the list is a list of the SpaceData keys in x, y, z, order
        bins : list
            if the name "bins" is not specified in the .attrs of the dmarray variable
            this specifies the bins for each vairable in a
            [[xbins], [ybins]] format
        xlim : list
            if the name "lim" is not specified in the .attrs of the dmarray variable
            this specifies the limit for the x variable [xlow, xhigh]
        ylim : list
            if the name "lim" is not specified in the .attrs of the dmarray variable
            this specifies the limit for the y variable [ylow, yhigh]
        zlim : list
            if the name "lim" is not specified in the .attrs of the dmarray variable
            this specifies the limit for the z variable [zlow, zhigh]
        # plotting
        #zlog : bool
        #    if the name "log" is not specified in the .attrs of the dmarray variable
        #    this specifies if the colorbar should be log (default=True)

        TODO
        ====
        Allow for the input of a list of SpaceData objecys for different sats
        """
        ## setup a default dictionary to step through to set values from kwargs
        self.settings = {}
        self.settings['variables'] = ['Epoch', 'Energy', 'Flux']
        self.settings['bins'] = None  # this is the linspace over the range with sqrt() of the len bins
        self.settings['xlim'] = None
        self.settings['ylim'] = None
        self.settings['zlim'] = None

        # if the key exists in kwargs replace setting with it, otherwise its an error
        for key in kwargs:
            if key not in self.settings:
                raise(KeyError('Invalid keyword specified ' + str(key)))
            self.settings[key] = kwargs[key]

        # check to see if the variables are in the spacedata
        for var in self.settings['variables']:
            if not var in data:  # TODO could check other capitalization
                raise(ValueError(str(var) + ' not found in the input data' ))

        # set default limits
        if self.settings['xlim'] == None:
            self.settings['xlim'] = (np.min(data[self.settings['variables'][0]]),
                                np.max(data[self.settings['variables'][0]]))
        if self.settings['ylim'] == None:
            self.settings['ylim'] = (np.min(data[self.settings['variables'][1]]),
                                np.max(data[self.settings['variables'][1]]))
        if self.settings['zlim'] == None:
            self.settings['zlim'] = (np.min(data[self.settings['variables'][2]]),
                                np.max(data[self.settings['variables'][2]]))

        # set default bins)()
        if self.settings['bins'] == None:
            # use the toolbox version of linspace so it works on dates
            forcedate = [False] * 2
            if isinstance(data[self.settings['variables'][0]][0], datetime.datetime):
                forcedate[0] = True
            if isinstance(data[self.settings['variables'][1]][0], datetime.datetime):
                forcedate[1] = True
            self.settings['bins'] = [np.asarray(tb.linspace(self.settings['xlim'][0],
                                             self.settings['xlim'][1],
                                             np.sqrt(len(data[self.settings['variables'][0]])), forcedate=forcedate[0])),
                np.asarray(tb.linspace(self.settings['ylim'][0],
                                             self.settings['ylim'][1],
                                             np.sqrt(len(data[self.settings['variables'][1]])), forcedate=forcedate[1]))]

        for key in data: # TODO can this be done easier?
            self[key] = data[key]
            try:
                self[key].attrs = data[key].attrs
            except:
                pass
        # do the spectrogram
        self._doSpec()


    def _doSpec(self):
        """
        Method operates on the inout data to bin up the spectrogram and adds it
        to the spectrogram class data
        """
        # this is here for in the future when we take a list a SpaceData objects
        sz = (self.settings['bins'][1].shape[0]-1, self.settings['bins'][0].shape[0]-1)
        overall_sum = np.zeros(sz, dtype=np.double)
        overall_count = np.zeros(sz, dtype=np.long)

        # the valid range for the histograms
        _range = [self.settings['xlim'], self.settings['ylim']]

        # if x/y is a time need to convert it to numbers (checking first element)
        var_time = [False]*2
        for ivar, var in enumerate(self.settings['variables']):
            if isinstance(self[var][0], datetime.datetime):
                self.settings['bins'][ivar] = mpl.dates.date2num(self.settings['bins'][ivar])
                try: #weird issue with arrays and date2num
                    _range[ivar] = mpl.dates.date2num(_range[ivar].tolist())
                except AttributeError:
                    try: # ugg if this is not an array this breaks
                        _range[ivar] = [mpl.dates.date2num(val.tolist()) for val in _range[ivar]]
                    except AttributeError:
                        _range[ivar] = [mpl.dates.date2num(val) for val in _range[ivar]]
                var_time[ivar] = True

        plt_data = np.vstack((self[self.settings['variables'][0]], self[self.settings['variables'][1]]))

        for ival, val in enumerate(var_time):
            if val:
                plt_data[ival] = date2num(plt_data[ival])

        if plt_data.dtype.name == 'object':  # why was this an abject
            plt_data = plt_data.astype(float)

        # go through and get rid of bad counts
        zdat = np.ma.masked_outside(self[self.settings['variables'][2]], self.settings['zlim'][0], self.settings['zlim'][1])
        zind = ~zdat.mask
        # get the number in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0, zind],
                                plt_data[1, zind],
                                bins = self.settings['bins'],
                                range = _range,
                                )
        # this is here for in the future when we take a list a SpaceData objects
        np.add(overall_count, H.transpose(), overall_count)

        # get the sum in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0,zind],
                                plt_data[1, zind],
                                bins = self.settings['bins'],
                                range = _range,
                                weights = zdat.data[zind]
                                )
        np.add(overall_sum, H.transpose(), overall_sum)

        overall_count = np.ma.masked_array(overall_count, overall_count == 0)
        data = np.ma.divide(overall_sum, overall_count)

        ## for plotting
        #ind0 = data.data <= 0
        #data[ind0] = np.ma.masked # add to the mask
        #data = np.ma.log10(data)
        self['spectrogram'] = dm.SpaceData()
        self['spectrogram'].attrs = self.settings
        self['spectrogram']['spectrogram'] = dm.dmarray(data,
            attrs={'name':str(self.settings['variables'][0]) + ' ' + str(self.settings['variables'][1]),
                   'xedges':'xedges',
                   'yedges':'yedges',})
        self['spectrogram']['xedges'] = dm.dmarray(xedges,
            attrs={'name':str(self.settings['variables'][0]),
                   'lim':[self.settings['xlim']],})
        self['spectrogram']['yedges'] = dm.dmarray(yedges,
            attrs={'name':str(self.settings['variables'][1]),
                   'lim':[self.settings['ylim']],})



