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
        # is variables in kwargs?
        if not 'variables' in kwargs:
            self._variables = ['Epoch', 'Energy', 'Flux']
        else:
            self._variables = kwargs['variables']
            # check to see if the variables are in the spacedata
            for var in kwargs['variables']:
                if not var in data:  # TODO could check other capitalization
                    raise(ValueError(str(var) + ' not found in the input data' ))
        for key in data: # TODO can this be done easier?
            self[key] = data[key]
            try:
                self[key].attrs = data[key].attrs
            except:
                pass
        if not 'bins' in kwargs:
            self.bins = []
            try:
                for var in self._variables:
                    self.bins.append(self[var].attrs['bins'])
            except (AttributeError, KeyError):
                raise(ValueError(str(var) + ' has no bins specified and no "bins" in .attrs'))
        else:
            self.bins = kwargs['bins']
        # check for xlim in kwargs
        if not 'xlim' in kwargs:
            try:
                self.xlim = self[self._variables[0]].attrs['lim']
            except (AttributeError, KeyError):
                self.xlim = (np.min(self[self._variables[0]]),
                             np.max(self[self._variables[0]]))
        else:
            self.xlim = kwargs['xlim']
        # check for ylim in kwargs
        if not 'ylim' in kwargs:
            try:
                self.ylim = self[self._variables[1]].attrs['lim']
            except (AttributeError, KeyError):
                self.ylim = (np.min(self[self._variables[1]]),
                             np.max(self[self._variables[1]]))
        else:
            self.ylim = kwargs['ylim']
        if not 'zlim' in kwargs:
            try:
                self.zlim = self[self._variables[2]].attrs['lim']
            except (AttributeError, KeyError):
                self.zlim = (np.min(self[self._variables[2]]),
                             np.max(self[self._variables[2]]))
        else:
            self.zlim = kwargs['zlim']
        # do the spectrogram
        self._doSpec()


    def _doSpec(self):
        """
        Method operates on the inout data to bin up the spectrogram and adds it
        to the spectrogram class data
        """
        # this is here for in the future when we take a list a SpaceData objects
        overall_sum = np.zeros((self.bins[1].shape[0]-1, self.bins[0].shape[0]-1), dtype=np.double)
        overall_count = np.zeros((self.bins[1].shape[0]-1, self.bins[0].shape[0]-1), dtype=np.long)

        # the valid range for the histograms
        _range = [self.xlim, self.ylim]

        # if x/y is a time need to convert it to numbers (checking first element)
        var_time = [False]*2
        for ivar, var in enumerate(self._variables):
            if isinstance(self[var][0], datetime.datetime):
                self.bins[ivar] = mpl.dates.date2num(self.bins[ivar])
                try: #weird issue with arrays and date2num
                    _range[ivar] = mpl.dates.date2num(_range[ivar].tolist())
                except AttributeError:
                    try: # ugg if this is not an array this breaks
                        _range[ivar] = [mpl.dates.date2num(val.tolist()) for val in _range[ivar]]
                    except AttributeError:
                        _range[ivar] = [mpl.dates.date2num(val) for val in _range[ivar]]
                var_time[ivar] = True

        plt_data = np.vstack((self[self._variables[0]], self[self._variables[1]]))

        for ival, val in enumerate(var_time):
            if val:
                plt_data[ival] = date2num(plt_data[ival])

        if plt_data.dtype.name == 'object':  # why was this an abject
            plt_data = plt_data.astype(float)

        # go through and get rid of bad counts
        zdat = np.ma.masked_outside(self[self._variables[2]], self.zlim[0], self.zlim[1])
        zind = ~zdat.mask
        # get the number in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0, zind],
                                plt_data[1, zind],
                                bins = self.bins,
                                range = _range,
                                )
        # this is here for in the future when we take a list a SpaceData objects
        np.add(overall_count, H.transpose(), overall_count)

        # get the sum in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0,zind],
                                plt_data[1, zind],
                                bins = self.bins,
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
        self['spectrogram']['spectrogram'] = dm.dmarray(data,
            attrs={'name':str(self._variables[0]) + ' ' + str(self._variables[1]),
                   'xedges':'xedges',
                   'yedges':'yedges',})
        self['spectrogram']['xedges'] = dm.dmarray(xedges,
            attrs={'name':str(self._variables[0]),
                   'lim':[self.xlim],})
        self['spectrogram']['yedges'] = dm.dmarray(yedges,
            attrs={'name':str(self._variables[1]),
                   'lim':[self.ylim],})



