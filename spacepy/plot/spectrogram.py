#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create and plot generic 'spectrograms' for space science.
This is not a signal processing routine and does not apply
Fourier transforms (or similar) to the data. The functionality
provided here is the binning (and averaging) of multi-dimensional
to provide a 2D output map of some quantity as a function of two
parameters. An example would be particle data from a satellite mission:
electron flux, at a given energy, can be binned as a function of
both time and McIlwain L, then plotted as a 2D color-map,
colloquially known as a spectrogram.

In many other settings 'spectrogram' refers to a transform of data
from the time domain to the frequency domain, and the subsequent plotting
of some quantity (e.g., power spectral density) as a function of time and
frequency. To approximate this functionality for, e.g., time-series magnetic field
data you would first calculate a the power spectral density and then use
:class:`Spectrogram` to rebin the data for visualization.

Authors: Brian Larsen and Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov, smorley@lanl.gov
Los Alamos National Laboratory

Copyright 2011 Los Alamos National Security, LLC.

.. rubric:: Class
.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    Spectrogram

.. rubric:: Function
.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    simpleSpectrogram
"""

import bisect
import copy
import datetime

import numpy as np
import matplotlib
import matplotlib.axes
import matplotlib.collections as mcoll
from matplotlib.dates import date2num, num2date
import matplotlib.colors
from matplotlib.colors import LogNorm
import matplotlib.transforms as mtransforms
import matplotlib.pyplot as plt

import spacepy.datamodel as dm
import spacepy.toolbox as tb
from  . import utils as spu

__contact__ = 'Brian Larsen, balarsen@lanl.gov'

__all__ = ['Spectrogram', 'simpleSpectrogram']

class Spectrogram(dm.SpaceData):
    """
    This class rebins data to produce a 2D data map that can be plotted as a spectrogram

    It is meant to be used on arbitrary data series.  The first series "x" is
    plotted on the abscissa and second series "y" is plotted on the ordinate and
    the third series "z" is plotted in color.

    The series are not passed in independently but instead inside a
    :class:`~spacepy.datamodel.SpaceData` container.

    Parameters
    ==========
    data : :class:`~spacepy.datamodel.SpaceData`
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
        this specifies the bins for each variable in a
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
    extended_out : bool (optional)
        if this is True add more information to the output data model (default True)

    Notes
    =====
    Helper routines are planned to
    facilitate the creation of the SpaceData container if the data are not in the format.

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> import numpy as np
    >>> import spacepy.plot as splot
    >>> sd = dm.SpaceData()
    >>> sd['radius'] = dm.dmarray(2*np.sin(np.linspace(0,30,500))+4, attrs={'units':'km'})
    >>> sd['day_of_year'] = dm.dmarray(np.linspace(74,77,500))
    >>> sd['1D_dataset'] = dm.dmarray(np.random.normal(10,3,500)*sd['radius'])
    >>> spec = splot.Spectrogram(sd, variables=['day_of_year', 'radius', '1D_dataset'])
    >>> ax = spec.plot()

    .. autosummary::

        ~Spectrogram.plot

    .. automethod:: Spectrogram.plot

    """

#    TODO
#    ====
#    Allow for the input of a list of SpaceData objects for different sats
#    Make "subclasses" that allow for data to be passed in directly avoiding the data model


    def __init__(self, data, **kwargs):
        """
        """
        super(Spectrogram, self).__init__()
        ## setup a default dictionary to step through to set values from kwargs
        self.specSettings = dm.SpaceData()
        self.specSettings['variables'] = ['Epoch', 'Energy', 'Flux']
        self.specSettings['bins'] = None  # this is the linspace over the range with sqrt() of the len bins
        self.specSettings['xlim'] = None
        self.specSettings['ylim'] = None
        self.specSettings['zlim'] = None
        self.specSettings['extended_out'] = True
        self.specSettings['axisDates'] = False

        # if the key exists in kwargs replace setting with it, otherwise its an error
        for key in kwargs:
            if key not in self.specSettings:
                raise(KeyError('Invalid keyword specified ' + str(key)))
            self.specSettings[key] = kwargs[key]

        # check to see if the variables are in the spacedata
        for var in self.specSettings['variables']:
            if not var in data:  # TODO could check other capitalization
                raise(KeyError('"{0}" not found in the input data'.format(var) ))

        # if the variables are empty error and quit
        if len(data[self.specSettings['variables'][0]]) == 0:
            raise(ValueError('No {0} data passed in'.format(self.specSettings['variables'][0])))
        if len(data[self.specSettings['variables'][1]]) == 0:
            raise(ValueError('No {0} data passed in'.format(self.specSettings['variables'][1])))
        if len(data[self.specSettings['variables'][2]]) == 0:
            raise(ValueError('No {0} data passed in'.format(self.specSettings['variables'][2])))

        # set limits, keywords override those in the data
        #if (self.specSettings['xlim'] is None) and (self.specSettings['bins'] is None):
        if self.specSettings['xlim'] is None:
            try:
                if 'lim' in data[self.specSettings['variables'][0]].attrs:
                    self.specSettings['xlim'] = (data[self.specSettings['variables'][0]].attrs['lim'][0],
                                                data[self.specSettings['variables'][0]].attrs['lim'][1])
                else:
                    dum1 = np.min(data[self.specSettings['variables'][0]]).tolist() #TODO: tolist here is a workaround for a bug in
                    dum2 = np.max(data[self.specSettings['variables'][0]]).tolist() #      datamodel's min method
                    self.specSettings['xlim'] = (dum1, dum2)
            except AttributeError: # was a numpy array not dmarray
                self.specSettings['xlim'] = (np.min(data[self.specSettings['variables'][0]]),
                                            np.max(data[self.specSettings['variables'][0]]))
        if self.specSettings['ylim'] is None:
            try:
                if 'lim' in data[self.specSettings['variables'][1]].attrs:
                    self.specSettings['ylim'] = (data[self.specSettings['variables'][1]].attrs['lim'][0],
                                                data[self.specSettings['variables'][1]].attrs['lim'][1])
                else:
                    self.specSettings['ylim'] = (np.min(data[self.specSettings['variables'][1]]),
                                                np.max(data[self.specSettings['variables'][1]]))
            except AttributeError: # was a numpy array not dmarray
                self.specSettings['ylim'] = (np.min(data[self.specSettings['variables'][1]]),
                                            np.max(data[self.specSettings['variables'][1]]))
        if self.specSettings['zlim'] is None:
            try:
                if 'lim' in data[self.specSettings['variables'][2]].attrs:
                    self.specSettings['zlim'] = (data[self.specSettings['variables'][2]].attrs['lim'][0],
                                                data[self.specSettings['variables'][2]].attrs['lim'][1])
                else:
                    self.specSettings['zlim'] = (np.min(data[self.specSettings['variables'][2]]),
                                                np.max(data[self.specSettings['variables'][2]]))
            except AttributeError: # was a numpy array not dmarray
                self.specSettings['zlim'] = (np.min(data[self.specSettings['variables'][2]]),
                                            np.max(data[self.specSettings['variables'][2]]))

        # are the axes dates?
        forcedate = [False] * 2
        if isinstance(data[self.specSettings['variables'][0]][0], datetime.datetime):
            forcedate[0] = True
        if isinstance(data[self.specSettings['variables'][1]][0], datetime.datetime):
            forcedate[1] = True
        self.specSettings['axisDates'] = forcedate

        # set default bins
        if self.specSettings['bins'] is None:
            # since it is not set by keyword was it set in the datamodel?
            attr_bins = ['bins' in data[var].attrs for var in self.specSettings['variables']]
            if dm.dmarray(attr_bins[0:2]).all():
                self.specSettings['bins'] = [dm.dmarray(data[self.specSettings['variables'][0]].attrs['bins']),
                                             dm.dmarray(data[self.specSettings['variables'][1]].attrs['bins']),]
                # TODO this is not a hard extension to doing one with bins and one default
            else:
                # use the toolbox version of linspace so it works on dates
                self.specSettings['bins'] = [dm.dmarray(tb.linspace(self.specSettings['xlim'][0],
                                                 self.specSettings['xlim'][1],
                                                 int(np.sqrt(len(data[self.specSettings['variables'][0]]))))),
                    dm.dmarray(tb.linspace(self.specSettings['ylim'][0],
                                                 self.specSettings['ylim'][1],
                                                 int(np.sqrt(len(data[self.specSettings['variables'][1]])))))]

        # copy all the used keys
        for key in self.specSettings['variables']:
            self[key] = data[key]
            try:
                self[key].attrs = data[key].attrs
            except:
                pass
        # do the spectrogram
        self._computeSpec()

    def _computeSpec(self):
        """
        Method operates on the input data to bin up the spectrogram and adds it
        to the Spectrogram class data
        """
        # this is here for in the future when we take a list a SpaceData objects
        sz = (self.specSettings['bins'][1].shape[0]-1, self.specSettings['bins'][0].shape[0]-1)
        overall_sum = dm.dmarray(np.zeros(sz, dtype=np.double))
        overall_count = dm.dmarray(np.zeros(sz, dtype=np.int_))

        # the valid range for the histograms
        _range = [self.specSettings['xlim'], self.specSettings['ylim']]

        # if x/y is a time need to convert it to numbers (checking first element)
        var_time = [False]*2
        for ivar, var in enumerate(self.specSettings['variables']):
            if isinstance(self[var][0], datetime.datetime):
                try:
                    self.specSettings['bins'][ivar] = matplotlib.dates.date2num(self.specSettings['bins'][ivar])
                except AttributeError: # it is already changed to date2num
                    pass
                try: #weird issue with arrays and date2num
                    _range[ivar] = matplotlib.dates.date2num(_range[ivar].tolist())
                except AttributeError:
                    try: # ugg if this is not an array this breaks
                        _range[ivar] = [matplotlib.dates.date2num(val.tolist()) for val in _range[ivar]]
                    except AttributeError:
                        _range[ivar] = [matplotlib.dates.date2num(val) for val in _range[ivar]]
                var_time[ivar] = True

        # ok not as a dmarray since it is local now
        plt_data = np.vstack((self[self.specSettings['variables'][0]], self[self.specSettings['variables'][1]]))

        for ival, val in enumerate(var_time):
            if val:
                plt_data[ival] = date2num(plt_data[ival])

        if plt_data.dtype.name == 'object':  # why was this an abject
            plt_data = plt_data.astype(float)

        # go through and get rid of "bad" counts
        zdat = np.ma.masked_outside(self[self.specSettings['variables'][2]],
                                    self.specSettings['zlim'][0],
                                    self.specSettings['zlim'][1])
        zind = ~zdat.mask
        # ma has the annoying feature of if all the masks are the same just giving one value
        try:
            if len(zind) == len(zdat):
                pass
        except TypeError: # no len to a scalar
            # ok not as a dmarray since it is local now
            zind = np.asarray([zind]*len(zdat))

        # get the number in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0, zind],
                                plt_data[1, zind],
                                bins = self.specSettings['bins'],
                                range = _range,
                                )
        # this is here for in the future when we take a list a SpaceData objects
        np.add(overall_count,
               np.require(H.transpose(), dtype=overall_count.dtype),
               overall_count)

        # get the sum in each bin
        H, xedges, yedges = np.histogram2d(plt_data[0, zind],
                                plt_data[1, zind],
                                bins = self.specSettings['bins'],
                                range = _range,
                                weights = zdat.data[zind]
                                )
        np.add(overall_sum, H.transpose(), overall_sum)

        overall_count = np.ma.masked_array(overall_count, overall_count == 0)
        # Explicitly ensure this array owns its mask
        overall_count.unshare_mask()
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
        if self.specSettings['extended_out']:
            self['spectrogram']['count'] = overall_count
            self['spectrogram']['sum'] = overall_sum

    def add_data(self, data):
        """
        Add another SpaceData with same keys, etc. to Spectrogram instance

        Examples
        --------
        >>> import spacepy.datamodel as dm
        >>> import numpy as np
        >>> import spacepy.plot as splot
        >>> sd = dm.SpaceData()
        >>> sd['radius'] = dm.dmarray(2*np.sin(np.linspace(0,30,500))+4, attrs={'units':'km'})
        >>> sd['day_of_year'] = dm.dmarray(np.linspace(74,77,500))
        >>> sd['1D_dataset'] = dm.dmarray(np.random.normal(10,3,500)*sd['radius'])
        >>> sd2 = dm.dmcopy(sd)
        >>> sd2['radius'] = dm.dmarray(2*np.cos(np.linspace(0,30,500))+4, attrs={'units':'km'})
        >>> sd2['1D_dataset'] = dm.dmarray(np.random.normal(10,3,500)*sd2['radius'])
        >>> spec = splot.Spectrogram(sd, variables=['day_of_year', 'radius', '1D_dataset'])
        >>> spec.add_data(sd2)
        >>> ax = spec.plot()
        """
        if not self.specSettings['extended_out']:
            raise(NotImplementedError('Cannot add data to a Spectrogram unless "extended_out" was True on initial creation'))
        b = Spectrogram(data, **self.specSettings)
        # if they are both masked keep them that way
        mask = self['spectrogram']['count'].mask & b['spectrogram']['count'].mask
        # turn off the mask by setting sum and count to zero where it masked for self an be sure b mask doesn't it self
        self['spectrogram']['count'][self['spectrogram']['count'].mask] = 0
        b['spectrogram']['count'][b['spectrogram']['count'].mask] = 0
        # put the mask back they are both bad
        self['spectrogram']['count'].mask = mask
        b['spectrogram']['count'].mask = mask
        self['spectrogram']['count'] += b['spectrogram']['count']
        self['spectrogram']['sum'] += b['spectrogram']['sum']
        self['spectrogram']['spectrogram'][...] = np.ma.divide(self['spectrogram']['sum'], self['spectrogram']['count'])

    def plot(self, target=None, loc=111, figsize=None, **kwargs):
        """
        Plot the spectrogram

        Other Parameters
        ================
        title : str
            plot title (default '')
        xlabel : str
            x axis label (default '')
        ylabel : str
            y axis label (default '')
        colorbar_label : str
            colorbar label (default '')
        DateFormatter : matplotlib.dates.DateFormatter
            The formatting to use on the dates on the x-axis (default matplotlib.dates.DateFormatter("%d %b %Y"))
        zlog : bool
            plot the z variable on a log scale (default True)
        cmap : matplotlib Colormap
            colormap instance to use
        colorbar : bool
            plot the colorbar (default True)
        axis : matplotlib axis object
            axis to plot the spectrogram to
        zlim : np.array
            array like 2 element that overrides (interior) the spectrogram zlim (default Spectrogram.specSettings['zlim'])
        figsize : tuple (optional)
            tuple of size to pass to figure(), None does the default
        """
        # go through the passed in kwargs to plot and look at defaults
        import matplotlib.pyplot as plt
        plotSettings_keys = ('title', 'xlabel', 'ylabel', 'DateFormatter',
                             'zlim', 'colorbar', 'colorbar_label', 'zlog',
                             'xlim', 'ylim', 'figsize', 'cmap')
        for key in kwargs:
            if key not in plotSettings_keys:
                raise(KeyError('Invalid keyword argument to plot(), "' + key + '"'))

        self.plotSettings = dm.SpaceData()
        for key in ['title', 'xlabel', 'ylabel']:
            if key in kwargs:
                self.plotSettings[key] = kwargs[key]
            else:
                self.plotSettings[key] = ''
        if 'zlog' in kwargs:
            self.plotSettings['zlog'] = kwargs['zlog']
        else:
            self.plotSettings['zlog'] = True

        if 'DateFormatter' in kwargs:
            self.plotSettings['DateFormatter'] = kwargs['DateFormatter']
        else:
            self.plotSettings['DateFormatter'] = matplotlib.dates.DateFormatter("%d %b %Y")
        if 'zlim' in kwargs:
            self.plotSettings['zlim'] = kwargs['zlim']
        elif 'zlim' in self.specSettings:
            self.plotSettings['zlim'] = self.specSettings['zlim']
        else:
            self.plotSettings['zlim'] = self.specSettings['zlim']
        if 'colorbar' in kwargs:
            self.plotSettings['colorbar'] =  kwargs['colorbar']
        else:
            self.plotSettings['colorbar'] = True
        if 'colorbar_label' in kwargs:
            self.plotSettings['colorbar_label'] = kwargs['colorbar_label']
        else:
            self.plotSettings['colorbar_label'] = ''
        if 'xlim' in kwargs:
            self.plotSettings['xlim'] = kwargs['xlim']
        elif 'xlim' in self.specSettings:
            self.plotSettings['xlim'] = self.specSettings['xlim']
        else:
            self.plotSettings['xlim'] = None
        if 'ylim' in kwargs:
            self.plotSettings['ylim'] = kwargs['ylim']
        elif 'ylim' in self.specSettings:
            self.plotSettings['ylim'] = self.specSettings['ylim']
        else:
            self.plotSettings['ylim'] = None

        fig, ax = spu.set_target(target, loc=loc, figsize=figsize)

        bb = np.ma.masked_outside(self['spectrogram']['spectrogram'], *self.plotSettings['zlim'])
        if 'cmap' in kwargs:
            self.plotSettings['cmap'] = kwargs['cmap']
        else:
            self.plotSettings['cmap'] = matplotlib.cm.rainbow

        if self.plotSettings['zlog']:
            pcm = ax.pcolormesh(self['spectrogram']['xedges'], self['spectrogram']['yedges'], np.asarray(bb),
                                norm=LogNorm(vmin=self.plotSettings['zlim'][0], vmax=self.plotSettings['zlim'][1]),
                                cmap=self.plotSettings['cmap'])
        else:
            pcm = ax.pcolormesh(self['spectrogram']['xedges'], self['spectrogram']['yedges'], np.asarray(bb),
                                cmap=self.plotSettings['cmap'], vmin=self.plotSettings['zlim'][0],
                                vmax=self.plotSettings['zlim'][1])

        if self.specSettings['axisDates'][0]:
            time_ticks = self._set_ticks_to_time(ax, 'x')
        elif self.specSettings['axisDates'][1]:
            time_ticks = self._set_ticks_to_time(ax, 'y')
        ax.set_title(self.plotSettings['title'])
        ax.set_xlabel(self.plotSettings['xlabel'])
        ax.set_ylabel(self.plotSettings['ylabel'])
        if self.plotSettings['ylim'] != None:
            ax.set_ylim(self.plotSettings['ylim'])
        if self.plotSettings['xlim'] != None:
            ax.set_xlim(self.plotSettings['xlim'])

        if self.plotSettings['colorbar']:
            cb =plt.colorbar(pcm, ax=ax)
            cb.set_label(self.plotSettings['colorbar_label'])
        return ax

    def vslice(self, value):
        """
        slice a spectrogram at a given position along the x axis, maintains
        variable names from spectrogram

        Parameters
        ==========
        value : float or datetime.datetime
            the value to slice the spectrogram at

        Returns
        =======
        out : datamodel.SpaceData
            spacedata containing the slice
        """
        # using bisect find the index of the spectrogram to use
        if isinstance(value, datetime.datetime):
            value = date2num(value)
        ind = bisect.bisect_right(self['spectrogram']['xedges'], value)
        ans = dm.SpaceData()
        ans[self['spectrogram'].attrs['variables'][1]] = tb.bin_edges_to_center(self['spectrogram']['yedges'])
        ans['yedges'] = self['spectrogram']['yedges'].copy()
        ans['xedges'] = self['spectrogram']['xedges'][ind:ind+2].copy()
        ans[self['spectrogram'].attrs['variables'][2]] = self['spectrogram']['spectrogram'][:,ind:ind+1]
        return ans

    def hslice(self, value):
        """
        slice a spectrogram at a given position along the y axis, maintains
        variable names from spectrogram

        Parameters
        ==========
        value : float or datetime.datetime
            the value to slice the spectrogram at

        Returns
        =======
        out : datamodel.SpaceData
            spacedata containing the slice
        """
        # using bisect find the index of the spectrogram to use
        if isinstance(value, datetime.datetime):
            value = date2num(value)
        ind = bisect.bisect_right(self['spectrogram']['yedges'], value)
        ans = dm.SpaceData()
        ans[self['spectrogram'].attrs['variables'][0]] = tb.bin_edges_to_center(self['spectrogram']['xedges'])
        ans['yedges'] = self['spectrogram']['yedges'][ind:ind+2].copy()
        ans['xedges'] = self['spectrogram']['xedges'].copy()
        ans[self['spectrogram'].attrs['variables'][2]] = self['spectrogram']['spectrogram'][ind, :]
        return ans

    def _set_ticks_to_time(self, axis, xy):
        """
        given the axis change the ticks to times
        """
        timeFmt = self.plotSettings['DateFormatter']
        if xy == 'x':
            ticks = axis.get_xticks()
            axis.set_xticklabels(matplotlib.dates.num2date(ticks))
            axis.xaxis.set_major_formatter(timeFmt)
            axis.get_figure().autofmt_xdate()
        elif xy == 'y':
            ticks = axis.get_yticks()
            axis.set_yticklabels(matplotlib.dates.num2date(ticks))
            axis.yaxis.set_major_formatter(timeFmt)
            axis.get_figure().autofmt_ydate()

    def __str__(self):
        return "<Spectrogram object>"

    __repr__ = __str__




def simpleSpectrogram(*args, **kwargs):
    """
    Plot a spectrogram given Z or X,Y,Z. This is a wrapper around pcolormesh()
    that can handle Y being a 2d array of time dependent bins. Like in the
    Van Allen Probes HOPE and MagEIS data files.

    Parameters
    ==========
    *args : 1 or 3 arraylike

        Call Signatures::

            simpleSpectrogram(Z, **kwargs)
            simpleSpectrogram(X, Y, Z, **kwargs)

    Other Parameters
    ================
    zlog : bool
        Plot the color with a log colorbar (default: True)
    ylog : bool
        Plot the Y axis with a log scale (default: True)
    alpha : scalar (0-1)
        The alpha blending value (default: None)
    cmap : string
        The name of the colormap to use (default: system default)
    vmin : float
        Minimum color value (default: Z.min(), if log non-zero min)
    vmax : float
        Maximum color value (default: Z.max())
    ax : matplotlib.axes
        Axes to plot the spectrogram on (default: None - new axes)
    cb : bool
        Plot a colorbar (default: True)
    cbtitle : string
        Label to go on the colorbar (default: None)
    zero_valid : bool
        Treat zero as a valid value on zlog plots and use same color as
        other under-minimum values. No effect with linear colorbar.
        (default: False, draw as fill)

        .. versionadded:: 0.5.0

    Returns
    =======
    ax : matplotlib.axes._subplots.AxesSubplot
        Matplotlib axes object that the plot is on
    """
    if len(args) not in [1,3]:
        raise(TypeError("simpleSpectrogram, takes Z or X, Y, Z"))

    if len(args) == 1: # just Z in, makeup an X and Y
        Z = np.ma.masked_invalid(np.asarray(args[0])) # make sure not a dmarray or VarCopy
        X = np.arange(Z.shape[0]+1)
        Y = np.arange(Z.shape[1]+1)
    else:
        # we really want X and Y to be one larger than Z
        X = np.asarray(args[0]) # make sure not a dmarray or VarCopy
        Y = np.asarray(args[1]) # make sure not a dmarray or VarCopy
        Z = np.ma.masked_invalid(np.asarray(args[2]))
        if X.shape[0] == Z.shape[0]: # same length, expand X
            X = tb.bin_center_to_edges(X) # hopefully evenly spaced
        if len(Y.shape) == 1: # 1d, just use as axis
            if Y.shape[0] == Z.shape[1]: # same length, expand Y
                Y = tb.bin_center_to_edges(Y) # hopefully evenly spaced
        elif len(Y.shape) == 2:
            Y_orig = Y
            # 2d this is time dependent and thus need to overplot several
            Y = tb.unique_columns(Y, axis=1)
            if Y.shape[1] == Z.shape[1]: # hopefully evenly spaced
                Y_uniq = Y
                Y_tmp = np.empty((Y.shape[0], Y.shape[1]+1), dtype=Y.dtype)
                for ii in range(Y.shape[0]):
                    Y_tmp[ii] = tb.bin_center_to_edges(Y[ii])
                Y = Y_tmp

    # deal with all the default keywords
    zlog    = kwargs.pop('zlog', True)
    zero_valid = kwargs.pop('zero_valid', False)
    ylog    = kwargs.pop('ylog', True)
    alpha   = kwargs.pop('alpha', None)
    cmap    = kwargs.pop('cmap', None)
    vmin    = kwargs.pop('vmin', np.min(Z) if not zlog else np.min(Z[np.nonzero(Z)]))
    vmax    = kwargs.pop('vmax', np.max(Z))
    ax      = kwargs.pop('ax', None)
    cb      = kwargs.pop('cb', True)
    cbtitle = kwargs.pop('cbtitle', None)

    # the first case is that X, Y are 1d and Z is 2d, just make the plot
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        fig = ax.get_figure()
    if zlog and zero_valid:
        Z[Z == 0] = vmin / 2.
    Z = Z.filled(0)
    Norm = matplotlib.colors.LogNorm if zlog else matplotlib.colors.Normalize
    if len(Y.shape) > 1:
        for i, yy in enumerate(Y_uniq):
            # which indices in X have these values
            ind = (yy == Y_orig).all(axis=1)
            Y_tmp = np.zeros((Y_orig.shape[0], Y.shape[1]), dtype=Y.dtype)
            Y_tmp[ind] = Y[i]
            Z_tmp = np.zeros_like(Z)
            Z_tmp[ind] = Z[ind]
            pc = ax.pcolormesh(X, Y_tmp[ind][0], Z_tmp.T,
                               norm=Norm(vmin=vmin, vmax=vmax), cmap=cmap,
                               alpha=alpha)
    else:
        pc = ax.pcolormesh(X, Y, Z.T, norm=Norm(vmin=vmin, vmax=vmax),
                           cmap=cmap,
                           alpha=alpha)
    if ylog: ax.set_yscale('log')

    if cb: # add a colorbar
        cb_ = fig.colorbar(pc)
        if cbtitle is not None:
            cb_.set_label(cbtitle)
    return ax
