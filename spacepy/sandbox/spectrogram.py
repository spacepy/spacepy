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
        self.specSettings = {}
        self.specSettings['variables'] = ['Epoch', 'Energy', 'Flux']
        self.specSettings['bins'] = None  # this is the linspace over the range with sqrt() of the len bins
        self.specSettings['xlim'] = None
        self.specSettings['ylim'] = None
        self.specSettings['zlim'] = None

        # if the key exists in kwargs replace setting with it, otherwise its an error
        for key in kwargs:
            if key not in self.specSettings:
                raise(KeyError('Invalid keyword specified ' + str(key)))
            self.specSettings[key] = kwargs[key]

        # check to see if the variables are in the spacedata
        for var in self.specSettings['variables']:
            if not var in data:  # TODO could check other capitalization
                raise(ValueError(str(var) + ' not found in the input data' ))

        # set default limits
        if self.specSettings['xlim'] == None:
            self.specSettings['xlim'] = (np.min(data[self.specSettings['variables'][0]]),
                                np.max(data[self.specSettings['variables'][0]]))
        if self.specSettings['ylim'] == None:
            self.specSettings['ylim'] = (np.min(data[self.specSettings['variables'][1]]),
                                np.max(data[self.specSettings['variables'][1]]))
        if self.specSettings['zlim'] == None:
            self.specSettings['zlim'] = (np.min(data[self.specSettings['variables'][2]]),
                                np.max(data[self.specSettings['variables'][2]]))

        # set default bins)()
        if self.specSettings['bins'] == None:
            # use the toolbox version of linspace so it works on dates
            forcedate = [False] * 2
            if isinstance(data[self.specSettings['variables'][0]][0], datetime.datetime):
                forcedate[0] = True
            if isinstance(data[self.specSettings['variables'][1]][0], datetime.datetime):
                forcedate[1] = True
            self.specSettings['bins'] = [np.asarray(tb.linspace(self.specSettings['xlim'][0],
                                             self.specSettings['xlim'][1],
                                             np.sqrt(len(data[self.specSettings['variables'][0]])), forcedate=forcedate[0])),
                np.asarray(tb.linspace(self.specSettings['ylim'][0],
                                             self.specSettings['ylim'][1],
                                             np.sqrt(len(data[self.specSettings['variables'][1]])), forcedate=forcedate[1]))]

        for key in data: # TODO can this be done easier?
            self[key] = data[key]
            try:
                self[key].attrs = data[key].attrs
            except:
                pass
        # do the spectrogram
        self._computeSpec()

    def _computeSpec(self):
        """
        Method operates on the inout data to bin up the spectrogram and adds it
        to the spectrogram class data
        """
        # this is here for in the future when we take a list a SpaceData objects
        sz = (self.specSettings['bins'][1].shape[0]-1, self.specSettings['bins'][0].shape[0]-1)
        overall_sum = np.zeros(sz, dtype=np.double)
        overall_count = np.zeros(sz, dtype=np.long)

        # the valid range for the histograms
        _range = [self.specSettings['xlim'], self.specSettings['ylim']]

        # if x/y is a time need to convert it to numbers (checking first element)
        var_time = [False]*2
        for ivar, var in enumerate(self.specSettings['variables']):
            if isinstance(self[var][0], datetime.datetime):
                self.specSettings['bins'][ivar] = mpl.dates.date2num(self.specSettings['bins'][ivar])
                try: #weird issue with arrays and date2num
                    _range[ivar] = mpl.dates.date2num(_range[ivar].tolist())
                except AttributeError:
                    try: # ugg if this is not an array this breaks
                        _range[ivar] = [mpl.dates.date2num(val.tolist()) for val in _range[ivar]]
                    except AttributeError:
                        _range[ivar] = [mpl.dates.date2num(val) for val in _range[ivar]]
                var_time[ivar] = True

        plt_data = np.vstack((self[self.specSettings['variables'][0]], self[self.specSettings['variables'][1]]))

        for ival, val in enumerate(var_time):
            if val:
                plt_data[ival] = date2num(plt_data[ival])

        if plt_data.dtype.name == 'object':  # why was this an abject
            plt_data = plt_data.astype(float)

        # go through and get rid of bad counts
        zdat = np.ma.masked_outside(self[self.specSettings['variables'][2]], self.specSettings['zlim'][0], self.specSettings['zlim'][1])
        zind = ~zdat.mask
        try:
            if len(zind) == len(zdat):
                pass
        except TypeError: # no len to a scalar
            zind = np.asarray([zind]*len(zdat))

        # get the number in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0, zind],
                                plt_data[1, zind],
                                bins = self.specSettings['bins'],
                                range = _range,
                                )
        # this is here for in the future when we take a list a SpaceData objects
        np.add(overall_count, H.transpose(), overall_count)

        # get the sum in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0,zind],
                                plt_data[1, zind],
                                bins = self.specSettings['bins'],
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
        self['spectrogram'].attrs = self.specSettings
        self['spectrogram']['spectrogram'] = dm.dmarray(data,
            attrs={'name':str(self.specSettings['variables'][0]) + ' ' + str(self.specSettings['variables'][1]),
                   'xedges':'xedges',
                   'yedges':'yedges',})
        self['spectrogram']['xedges'] = dm.dmarray(xedges,
            attrs={'name':str(self.specSettings['variables'][0]),
                   'lim':[self.specSettings['xlim']],})
        self['spectrogram']['yedges'] = dm.dmarray(yedges,
            attrs={'name':str(self.specSettings['variables'][1]),
                   'lim':[self.specSettings['ylim']],})

    def plot(self, fignum=None, **kwargs):
        """
        Plot the spectrogram

        TODO: **kwargs is just passed straight through o pcolormesh for now
        """
        # go through the passed in kwargs to plot and look at defaults
        import matplotlib.pyplot as plt

        plotSettings_keys = ('title', 'x_label', 'y_label', 'DateFormatter', 
                             'zlim', 'nocolorbar', 'colorbar_label')
        for key in kwargs:
            if key not in plotSettings_keys:
                raise(ValueError('Invalid keyword argument to plot(), "' + key + '"'))

        self.plotSettings = {}
        if 'title' in kwargs:
            self.plotSettings['title'] = kwargs['title']
        else:
            self.plotSettings['title'] = ''
        if 'x_label' in kwargs:
            self.plotSettings['x_label'] = kwargs['x_label']
        else:
            self.plotSettings['x_label'] = ''
        if 'y_label' in kwargs:
            self.plotSettings['y_label'] = kwargs['y_label']
        else:
            self.plotSettings['y_label'] = ''

        if 'DateFormatter' in kwargs:
            self.plotSettings['DateFormatter'] = kwargs['DateFormatter']
        else:
            self.plotSettings['DateFormatter'] = mpl.dates.DateFormatter("%d %b %Y")
        if 'zlim' in kwargs:
            self.plotSettings['zlim'] = kwargs['zlim']
        else:
            self.plotSettings['zlim'] = self.specSettings['zlim']
        if 'nocolorbar' in kwargs:
            self.plotSettings['nocolorbar'] =  kwargs['nocolorbar']
        else:
            self.plotSettings['nocolorbar'] = False
        if 'colorbar_label' in kwargs:
            self.plotSettings['colorbar_label'] = kwargs['colorbar_label']
        else:
            self.plotSettings['colorbar_label'] = ''


        if fignum == None:
            fig = plt.figure()
        else:
            fig = plt.figure(fignum)
        ax = fig.add_subplot(111)
        bb = np.ma.masked_outside(self['spectrogram']['spectrogram'], *self.plotSettings['zlim'])
        pcm = ax.pcolormesh(self['spectrogram']['xedges'], self['spectrogram']['yedges'], np.ma.log10(bb))
        time_ticks = self._set_ticks_to_time(ax)
        ax.set_title(self.plotSettings['title'])
        ax.set_xlabel(self.plotSettings['x_label'])
        ax.set_ylabel(self.plotSettings['y_label'])
        if not self.plotSettings['nocolorbar']:
            if 'cmap' in kwargs:
                self.plotSettings['cmap'] = kwargs['cmap']
            else:
                self.plotSettings['cmap'] = mpl.cm.rainbow
            #cb = mpl.colorbar.ColorbarBase(ax, cmap=self.plotSettings['cmap'],
                                   #norm=norm,
            #                       orientation='vertical')
            cb = mpl.pyplot.colorbar(pcm)
            #cb = plt.colorbar(np.ma.log10(bb))
            cb.set_label(self.plotSettings['colorbar_label'])
        return fig

    def _set_ticks_to_time(self, axis):
        """
        given the axis change the ticks to times
        """
        import matplotlib.dates
        ticks = axis.get_xticks()
        axis.set_xticklabels(mpl.dates.num2date(ticks))
        timeFmt = mpl.dates.DateFormatter("%d %b %Y")
        axis.xaxis.set_major_formatter(timeFmt)
        axis.get_figure().autofmt_xdate()




