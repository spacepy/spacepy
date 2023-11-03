#!/usr/bin/python
# -*- coding: utf-8 -*-

"""SeaPy -- Superposed Epoch in Python.

This module contains superposed epoch class types and a variety of functions
for using on superposed epoch objects.
Each instance must be initialized with (assuming import seapy as se):

>>> obj = se.Sea(data, times, epochs)

To perform a superposed epoch analysis

>>> obj.sea()

To plot

>>> obj.plot()

If multiple SeaPy objects exist, these can be combined into a single object

>>> objdict = seadict([obj1, obj2],['obj1name','obj2name'])

and then used to create a multipanel plot

>>> multisea(objdict)


For two-dimensional superposed epoch analyses, initialize an Sea2d() instance

>>> obj = se.Sea2d(data, times, epochs, y=[4., 12.])

All object methods are the same as for the 1D object. Also, the multisea()
function should accept both 1D and 2D objects, even mixed together. Currently,
the plot() method is recommended for 2D SEA.


--++-- By Steve Morley --++--

smorley@lanl.gov
Los Alamos National Laboratory

Copyright 2010 Los Alamos National Security, LLC.
"""
import warnings

import numpy as np
import numpy.ma as ma
import numbers
import datetime as dt
import spacepy.toolbox as tb
from spacepy import help
import spacepy.time as spt
import spacepy.datamodel as dm
import spacepy.plot as spplt
import matplotlib.pyplot as plt
from matplotlib.dates import date2num, num2date

__contact__ = 'Steve Morley, smorley@lanl.gov'

class SeaBase(object):
    """SeaPy Superposed epoch analysis base class

    Do not use directly -- subclass it!

    """
    def __init__(self, data, times, epochs, **kwargs):
        self._kwargs = kwargs
        self._times = times
        self._epochs = epochs
        self.data = np.asarray(data, dtype=float)
        if isinstance(times, spt.Ticktock):
            self.times=times.UTC
        else:
            self.times = times
        td = np.diff(times)
        noncontig = (not np.isclose(td, td[0]).all()) if np.issubdtype(td.dtype, np.inexact)\
              else (len(np.unique(td)) > 1)
        nonmon = ((np.array([t.total_seconds() for t in td])
                if isinstance(td[0], dt.timedelta) else td) < 0).any()
        if nonmon or noncontig:
            warnings.warn('Input time not {}; results are unlikely to be valid.'.format(
                          ('monotonic or contiguous' if nonmon else 'contiguous')
                          if noncontig else 'monotonic'))
        if len(epochs) > len(times) / 2:
            warnings.warn('Too many epochs; results are unlikely to be valid.')
        if isinstance(epochs, spt.Ticktock):
            self.epochs = epochs.UTC
        else:
            self.epochs = epochs
        self.verbose = kwargs['verbose']
        if type(kwargs['delta']) == dt.timedelta:
            t_delt = kwargs['delta'].days + kwargs['delta'].seconds/86400
            self.delta = t_delt
        else:
            self.delta = kwargs['delta']
        if type(kwargs['window']) == dt.timedelta:
            self.window = kwargs['window'].days + kwargs['window'].seconds/86400
        else:
            self.window = kwargs['window']
        self.window = np.ceil(float(self.window)/float(self.delta))
        if kwargs['window'] != self.window:
            warnings.warn(
                'Window size changed to {0} (points) to fit resolution ({1})'.format(
                self.window, self.delta))
        self.bound_type = None

    def __str__(self):
        """Define String Representation of Sea object"""

        strhead = 'Superposed Epoch Object:'
        strtail = 'Data array - {0}; #Epochs - {1}'.format(self.data.shape, len(self.epochs))
        return ' '.join((strhead, strtail))

    __repr__ = __str__

    def __len__(self):
        """Calling len(obj) will return the number of epochs."""

        return len(self.epochs)

    def restoreepochs(self):
        """Replaces epoch times stored in obj.badepochs in the epochs attribute

        """
        try:
            dum = self.badepochs
        except AttributeError:
            raise AttributeError('No bad epochs to restore')

        self.epochs = np.union1d(self.badepochs,self.epochs)
        if self.verbose: print('Bad epochs restored to epochs attribute')

        return None

    def _timeepoch(self, delt):
        #check type of time input and throw error message
        el1,ep1 = self.times[0], self.epochs[0]
        el1num = isinstance(el1, numbers.Number) or (type(el1)==np.float64)
        ep1num = isinstance(ep1, numbers.Number) or (type(ep1)==np.float64)
        if isinstance(el1, dt.datetime) and (type(el1) == type(self.epochs[0])):
            #both time and epochs are datetime objects
            #convert to serial time
            if self.verbose:
                print('converting to serial time')
            dum = date2num(self.epochs)
            t_epoch = np.array(dum)
            dum = date2num(self.times)
            time = np.array(dum)
            ser_flag = False
        elif el1num and ep1num:
            #time is serial, do nothing
            if self.verbose:
                print('time is serial')
            t_epoch = np.array(self.epochs, dtype=float)
            time = np.array(self.times, dtype=float)
            ser_flag = True
        else:
            raise ValueError('Time and Epochs must be consistently typed (numeric/datetime)')

        lose0 = np.where(t_epoch > time[-1]-(self.window*self.delta))
        lose1 = np.where(t_epoch < time[0]+(self.window*self.delta))
        if len(lose0)>0 and len(lose1)>0:
            linds = np.union1d(lose0[0],lose1[0])
            if len(linds)>0:
                if self.verbose:
                    print('sea(): {0} out-of-range epochs moved to badepochs attribute'.format(len(linds)))
            if ser_flag:
                self.badepochs = t_epoch[linds]
            else:
                self.badepochs = num2date(t_epoch[linds])
            ##TODO: replace following two lines with calls to tOverlapHalf
            keep0 = np.where(t_epoch <= time[-1]-(self.window*self.delta))
            keep1 = np.where(t_epoch >= time[0]+(self.window*self.delta))
            kinds = np.intersect1d(keep0[0],keep1[0])

            self.epochs = t_epoch[kinds]
            t_epoch = t_epoch[kinds]
            if len(self.epochs)==0:
                raise RuntimeError('No valid epochs for data supplied')
        return time, t_epoch

    def random(self, n=None):
        '''Return a new Sea object, of same dimensionality, with a set of random epochs
        '''
        if not n:
            n = len(self._epochs)

        sttime = self._epochs[0]
        entime = self._epochs[-1]
        try:
            dum = sttime + np.pi
            cflag = False
        except TypeError:
            sttime = date2num(sttime)
            entime = date2num(entime)
            cflag = True
        finally:
            repochs = (entime - sttime) * np.random.random_sample(n) + sttime
        if cflag:
            repochs = num2date(repochs)
        newclass = type(self)
        new = newclass(self.data, self._times, sorted(repochs), **self._kwargs)
        return new



class Sea(SeaBase):
    """SeaPy Superposed epoch analysis object

    Initialize object with data, times, epochs, window (half-width) and delta (optional).
    'times' and epochs should be in some useful format
    Includes method to perform superposed epoch analysis of input data series

    Parameters
    ==========
    data : array_like
        list or array of data
    times : array_like
        list of datetime objects (or list of serial times). Must be contiguous
        (constant cadence) and monotonically increasing.

        .. versionchanged:: 0.5.0
           Issues a warning for non-contiguous/non-monotonic times.

    epochs : array_like
        list of datetime objects (or serial times) for zero epochs in SEA.
        For a suitable SEA, this should be substantially shorter than
        ``times``.

        .. versionchanged:: 0.5.0
           Issues a warning for too many epochs; arbitrarily defined as
           more than half the number of times.

    window : datetime.timedelta
        size of the half-window for the SEA (can also be given as serial time)
    delta : datetime.timedelta
        resolution of the input data series, which must be uniform (can also be
        given as serial time)

    Notes
    =====
    Output can be nicely plotted with :py:meth:`plot`, or for multiple objects
    use the :py:func:`multisea` function

    .. currentmodule:: spacepy.seapy
    .. autosummary::
        ~Sea.sea
        ~Sea.plot
    .. automethod:: sea
    .. automethod:: plot
    """
    def __init__(self, data, times, epochs, window=3., delta=1., verbose=True):
        super(Sea, self).__init__(data, times, epochs, \
              window=window, delta=delta, verbose=verbose)

    def sea(self, **kwargs):
        """Method called to perform superposed epoch analysis on data in object.

        Uses object attributes obj.data, obj.times, obj.epochs, obj.delta,
        obj.window, all of which must be available on instantiation.

	Other Parameters
        ================
        storedata : boolean
	    saves matrix of epoch windows as obj.datacube (default = False)
	quartiles : list
	    calculates the quartiles as the upper and lower bounds (and is default);
        ci : float
	    will find the bootstrapped confidence intervals of ci_quan at the ci percent level (default=95)
        mad : float
	    will use +/- the median absolute deviation for the bounds;
        ci_quan : string
	    can be set to 'median' (default) or 'mean'

        Notes
        =====
        A basic plot can be raised with :meth:`plot`
        """
        #check this hasn't already been done
        #TODO: find out why doing two .sea() calls back-to-back fails 2nd time
        if hasattr(self, 'semedian') or hasattr(self, 'semean'):
            return None

        #check defaults
        defaults = {'storedata': True, 'quartiles': True, 'ci': False,
                    'mad': False, 'ci_quan': 'median'}
        for default in defaults:
            if default not in kwargs:
                kwargs[default] = defaults[default]

        #ensure all input is np array
        delt = float(self.delta)
        if isinstance(self.data, np.ndarray):
            y = self.data
        else:
            y = np.asarray(self.data, dtype=float)

        if kwargs['ci']:
            kwargs['quartiles'], kwargs['mad'] = False, False
        if kwargs['mad']:
            kwargs['quartiles'], kwargs['ci'] = False, False

        time, t_epoch = self._timeepoch(delt)

        #build SEA matrix and perform analysis
        wind = int(self.window)
        m = int(2*wind + 1)
        n = len(t_epoch)
        y_sea = np.zeros((n,m), dtype=float)
        blankslice = np.zeros([m], dtype=float)
        for i in range(n):
            dif = np.abs(time-t_epoch[i])
            j = np.where(dif == np.min(dif))
            stpt = j[0][0]-wind
            enpt = j[0][0]+wind+1
            sea_slice = blankslice.copy()
            if stpt < 0: #fix for bad epochs not correctly moved to badepochs attr #TODO: make badepochs robust or do all checking here
                sea_slice[0:abs(stpt)] = np.NaN
                sea_slice[abs(stpt):] = y[0:enpt]
            elif enpt >= len(y):
                tmpslice = y[stpt:]
                sea_slice[:len(tmpslice)] = tmpslice
                sea_slice[len(tmpslice):] = np.NaN
            else:
                sea_slice = y[stpt:enpt]

            y_sea[i,0:] = sea_slice

        #find SEA mean, median and percentiles - exclude NaNs (or badval)
        try:
            badval = kwargs['badval']
        except KeyError:
            badval = np.nan
            y_sea_m = ma.masked_where(np.isnan(y_sea), y_sea)
        else:
            y_sea_m = ma.masked_values(y_sea, badval)
        self.semean = [np.mean(y_sea_m[:,i].compressed()) for i in range(m)]
        self.semedian = [np.median(y_sea_m[:,i].compressed()) for i in range(m)]
        self.semean, self.semedian = np.array(self.semean), np.array(self.semedian)
        self.bound_low = np.zeros((m,1))
        self.bound_high = np.zeros((m,1))

        if kwargs['quartiles']:
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                qul = np.percentile(dum, (25,75))
                self.bound_low[i], self.bound_high[i] = qul[0], qul[1]
                self.bound_type = 'quartiles'
        elif kwargs['ci']: #bootstrapped confidence intervals (95%)
            funcdict = {'mean': np.mean,
                        'median': np.median}
            try:
                if isinstance(kwargs['ci'], bool):
                    raise ValueError #fall through to default case
                else:
                    ci_level = float(kwargs['ci'])
            except ValueError:
                ci_level = 95
            from spacepy.poppy import boots_ci
            if hasattr(kwargs['ci_quan'], "__call__"): #ci_quan is a function
                ci_func = kwargs['ci_quan']
            else:
                ci_func = funcdict[kwargs['ci_quan']]
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                self.bound_low[i], self.bound_high[i] = \
                     boots_ci(dum, 800, ci_level, ci_func)
                self.bound_type = 'ci'
        elif kwargs['mad']: #median absolute deviation
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                spread_mad = tb.medAbsDev(dum)
                self.bound_low[i] = self.semedian[i]-spread_mad
                self.bound_high[i] = self.semedian[i]+spread_mad
                self.bound_type = 'mad'

        self.x = np.linspace(-1.*self.window*self.delta, self.window*self.delta, \
         len(self.semedian))
        if kwargs['storedata']:
            self.datacube = y_sea_m
            if self.verbose:
                print('sea(): datacube added as new attribute')

        if self.verbose:
            print('Superposed epoch analysis complete')


    #def sea_norm(self, epoch2, nbins=100, storedata=False, quartiles=True, ci=False, mad=False,
    #ci_quan='median'):
        #"""Method for normalized superposed epoch analysis (creates new object)

        #Examples
        #Assuming an already instantiated Sea object, xampl

        #>>> xampl_norm = xampl.sea_norm(epoch2)

        #Note that epoch2 must be a paired set with xampl.epochs

        #Inputs:
        #Uses object attributes obj.data, obj.times, obj.epochs, obj.delta.
        #all of which must be available on instantiation.
        #Also requires second set of epochs.

        #Other Parameters
        #nbins (default = 100) - number of bins for normalized SEA
        #storedata (default = False) - saves matrix of epoch windows as obj.datacube
        #quartiles calculates the quartiles as the upper and lower bounds (and is default);
        #ci will find the bootstrapped confidence intervals (and requires ci_quan to be set);
        #mad will use +/- the median absolute deviation for the bounds;
        #ci_quan can be set to 'median' or 'mean'

        #Output is a new SeaPy object with all Sea methods.
        #"""
        #import numpy as np
        #import numpy.ma as ma

        ##ensure all input is np array
        #delt = float(self.delta)
        #y = np.array(self.data, dtype=float)

        #if ci:
            #quartiles, mad = False, False
        #if mad:
            #quartiles, ci = False, False

        #time,t_epoch = self._timeepoch(delt)

        ##build SEA matrix and perform analysis
        #wind = self.window
        #m = int(2*wind + 1)
        #n = len(t_epoch)
        #y_sea = np.zeros((n,m), dtype=float)

        #for i in range(n):
            #t_a, t_b = abs(time-t_epoch1[i]), abs(time-t_epoch2[i])
            #j = np.where(t_a == np.min(t_a))
            #k = np.where(t_b == np.min(t_b))
            #norm_ax = np.linspace(t_epoch1[i], t_epoch2[i], nbins+1)
            #norm_event = np.interp(norm_ax, time[j[1]:k[1]], y[j[1]:k[1]]) #piecewise linear
            #norm_event[0], norm_event[-1] = y[j[0]], y[k[0]] #Fix endpoints
            #y_sea[i,0:] = norm_event

        ##for i in range(n):
            ##dif = np.abs(time-t_epoch[i])
            ##j = np.where(dif == np.min(dif))
            ##sea_slice = y[j[0][0]-wind:j[0][0]+wind+1]
            ##y_sea[i,0:] = sea_slice

        ##Instantiate new SeaPy object
        #outobj = Sea(self.data, self.times, self.epochs, self.window, self.delta)
        #print 'New Sea object created'

        ##find SEA mean, median and percentiles - exclude NaNs
        #y_sea_m = ma.masked_where(np.isnan(y_sea), y_sea)
        #outobj.semean = [np.mean(dum.compressed()) for dum in y_sea_m]
        #outobj.semedian = [np.median(dum.compressed()) for dum in y_sea_m]
        #outobj.bound_low = np.zeros((m,1))
        #outobj.bound_high = np.zeros((m,1))

        #if quartiles:
            #for i in range(m):
                #dum = np.sort(y_sea_m[:,i].compressed())
                #qul = np.percentile(dum, (25,75))
                #outobj.bound_low[i], outobj.bound_high[i] = qul[0], qul[1]
        #elif ci: #bootstrapped confidence intervals (95%)
            #from spacepy.poppy import boots_ci
            #if ci_quan == 'mean':
                #ci_func = lambda x: np.mean(x)
            #else:
                #ci_func = lambda x: np.median(x)
            #for i in range(m):
                #dum = np.sort(y_sea_m[:,i].compressed())
                #outobj.bound_low[i], outobj.bound_high[i] = \
                    #boots_ci(dum, 800, 95, ci_func)
        #elif mad: #median absolute deviation
            #from spacepy.toolbox import medabsdev
            #for i in range(m):
                #dum = np.sort(y_sea_m[:,i].compressed())
                #spread_mad = medabsdev(data)
                #outobj.bound_low[i] = outobj.semedian[i]-spread_mad
                #outobj.bound_high[i] = outobj.semedian[i]+spread_mad

        #outobj.x = np.linspace(0, nbins, len(outobj.semedian))
        #if storedata:
            #outobj.datacube = y_sea_m
            #print 'sea_norm(): datacube added to new object as attribute'

        #return 'Superposed epoch analysis complete'

    def plot(self, xquan = 'Time Since Epoch', yquan='', xunits='',
                yunits='', epochline=True, usrlimy=[], show=True, target=None,
                loc=111, figsize=None, dpi=None, transparent=True, color='#7F7FFF'):
        """Method called to create basic plot of superposed epoch analysis.

        Parameters
	==========
        Uses object attributes created by the obj.sea() method.

        Other Parameters
        ================
        xquan : str
	    (default = 'Time since epoch' ) - x-axis label.
	yquan : str
	    default None - yaxus label
	xunits : str
	    (default = None) - x-axis units.
        yunits : str
	    (default = None) - y-axis units.
        epochline : boolean
	    (default = True) - put vertical line at zero epoch.
        usrlimy : list
	    (default = []) - override automatic y-limits on plot.
        transparent : boolean
	    (default True): make patch for low/high bounds transparent
        color : str
            Color to use for the patch if not transparent.
            (default #7F7FFF, a medium blue)

        Notes
        =====
        If both quan and units are supplied, axis label will read
        'Quantity Entered By User [Units]'
        """
        try:
            assert hasattr(self, 'semedian') or hasattr(self, 'semean')
        except AssertionError:
            raise ValueError('No superposed epoch results to plot')

        if len(xunits)<1:
            xlstr = '{0}'.format(xquan)
        else:
            xlstr = '{0} [{1}]'.format(xquan, xunits)
        if len(yquan)>=1 and len(yunits)>=1:
            ylstr = '{0} [{1}]'.format(yquan, yunits)
        elif len(yquan)>=1 and len(yunits)<1:
            ylstr = '{0}'.format(yquan)
        else:
            ylstr = ''

        fig, ax0 = spplt.utils.set_target(target, loc=loc, figsize=figsize)

        if transparent:
            ax0.fill_between(self.x, self.bound_low.ravel(), self.bound_high.ravel(),
                             edgecolor='none', facecolor=color, interpolate=True, alpha=0.25)
        else:
            ax0.fill_between(self.x, self.bound_low.ravel(), self.bound_high.ravel(),
                             edgecolor='none', facecolor=color, interpolate=True)
        ax0.plot(self.x, self.semedian, 'k-', lw=2.0)
        ax0.plot(self.x, self.semean, 'r--', lw=1.25)
        plt.xlabel(xlstr)
        plt.ylabel(ylstr)

        if usrlimy:
            ax0.set_ylim(usrlimy)

        if epochline:
            ax0.axvline(0, 0, 1, color='k', ls=':')

        if show:
            plt.show()
            return None
        else:
            return ax0


class Sea2d(SeaBase):
    """SeaPy 2D Superposed epoch analysis object

    Initialize object with data (n element vector), times(y*n array), 
    epochs, window (half-width), delta (optional), and 
    y (two-element vector with max and min of y;optional)
    'times' and epochs should be in some useful format
    Includes method to perform superposed epoch analysis of input data series

    Parameters
    ==========
    data : array_like
        2-D array of data (0th dimension is quantity y, 1st dimension is time)
    times : array_like
        list of datetime objects (or list of serial times). Must be contiguous
        (constant cadence) and monotonically increasing.

        .. versionchanged: 0.5.0
           Issues a warning for non-contiguous/non-monotonic times.

    epochs : array_like
        list of datetime objects (or serial times) for zero epochs in SEA.
        For a suitable SEA, this should be substantially shorter than
        ``times``.

        .. versionchanged: 0.5.0
           Issues a warning for too many epochs; arbitrarily defined as
           more than half the number of times.

    window : datetime.timedelta
        size of the half-window for the SEA (can also be given as serial time)
    delta : datetime.timedelta
        resolution of the input data series, which must be uniform (can also be
        given as serial time)

    Notes
    =====
    Output can be nicely plotted with :meth:`plot`, or for multiple
    objects use the :func:`multisea` function


    .. currentmodule:: spacepy.seapy
    .. autosummary::
        ~Sea2d.sea
        ~Sea2d.plot
    .. automethod:: sea
    .. automethod:: plot
    """
    def __init__(self, data, times, epochs, window=3., delta=1., verbose=False, y=[]):
        super(Sea2d, self).__init__(data, times, epochs, window=window, \
              delta=delta, verbose=False, y=[])

        if y:
            self.y = np.linspace(y[0], y[1], data.shape[0]+1)
        else:
            self.y = np.linspace(0, data.shape[0]-1, data.shape[0]+1)

    def sea(self, storedata=False, quartiles=True, ci=False, mad=False,
        ci_quan='median', nmask=1, **kwargs):
        """Perform 2D superposed epoch analysis on data in object

        Uses object attributes obj.data, obj.times, obj.epochs, obj.delta,
        obj.window, all of which must be available on instantiation.

        Other Parameters
        ================
        storedata : boolean
            saves matrix of epoch windows as obj.datacube (default = False)
        quartiles : list
            calculates the inter-quartile range to show the spread
            (and is default);
        ci : float
            will find the bootstrapped confidence interval
            (and requires ci_quan to be set)
        mad : float
            will use the median absolute deviation for the spread;
        ci_quan : string
            can be set to 'median' or 'mean'

        Notes
        =====
        A basic plot can be raised with :meth:`plot`
        """
        #ensure all input is np array or correct form
        delt = float(self.delta)
        y = np.asarray(self.data, dtype=float)
        time,t_epoch = self._timeepoch(delt)
        if not nmask:
            nmask = 0 #set mask to exclude none

        #build SEA matrix and perform analysis
        wind = int(self.window)
        l = y.shape[0]
        m = 2*int(wind) + 1
        n = len(t_epoch)
        y_sea = np.zeros((l,m,n), dtype=float)
        for i in range(n):
            dif = np.abs(time-t_epoch[i])
            j = np.where(dif == np.min(dif))
            sea_slice = y[:,j[0][0]-wind:j[0][0]+wind+1]
            y_sea[:,:,i] = sea_slice

        #find SEA mean, median and percentiles - exclude NaNs
        #y_sea_m = ma.masked_where(y_sea < 0., y_sea)
        try:
            badval = kwargs['badval']
        except KeyError:
            badval = np.nan
            y_sea_m = ma.masked_where(np.isnan(y_sea), y_sea)
        else:
            y_sea_m = ma.masked_values(y_sea, badval)
        #now get SEA quantities
        self.semean, self.semedian, self.countmask = np.empty((l,m)), np.empty((l,m)), np.empty((l,m))
        yj=0
        for ti in range(int(m)):
            for yj in range(l):
                self.semean[yj,ti] = np.mean(y_sea_m[yj,ti,:].compressed()) #np.mean(y_sea_m, axis=2)
                self.semedian[yj,ti] = np.median(y_sea_m[yj,ti,:].compressed())
                self.countmask[yj,ti] = len(y_sea_m[yj,ti,:].compressed())
        self.semean = ma.masked_where(self.countmask<nmask, self.semean)#(np.isnan(self.semean), self.semean)
        self.semedian = ma.masked_where(self.countmask<nmask, self.semedian)#(np.isnan(self.semedian), self.semedian)
        self.perc25 = np.zeros((l,m))
        self.perc50 = np.zeros((l,m))
        self.perc75 = np.zeros((l,m))

        #maybe use median absolute deviation?
        #p75 = np.around(75*n/100)
        #p25 = np.around(25*n/100)
        #sea_sort = y_sea_m.copy()
        #sea_sort.sort(axis=0)
        ##self.perc25 = sea_sort[p25,0:]
        ##self.perc75 = sea_sort[p75,0:]

        self.x = np.linspace(-1.*self.window*delt,self.window*delt, \
        self.semedian.shape[1]+1)

        if storedata:
            self.datacube = y_sea_m
            print('sea(): datacube added as new attribute')

        if self.verbose:
            print('Superposed epoch analysis complete')

    def plot(self, xquan = 'Time Since Epoch', yquan='', xunits='',
                yunits='', zunits='', epochline=True, usrlimy=[],
                show=True, zlog=True, figsize=None, dpi=300):
        """Method called to create basic plot of 2D superposed epoch analysis.

        Uses object attributes created by :meth:`sea`.

        Other Parameters
        ================
        x(y)quan : str
            x(y)-axis label.  (default = 'Time since epoch' (None))
        x(y/z)units : str
            x(y/z)-axis units. (default = None (None))
        epochline : boolean
            put vertical line at zero epoch. (default = True)
        usrlimy : list
            override automatic y-limits on plot. (default = [])
        show : boolean
            shows plot; set to false to output plot object to variable
            (default = True)
        figsize : tuple
            (width, height) in inches
        dpi : int
            figure resolution in dots per inch (default=300)

        Notes
        =====
        If both quan and units are supplied, axis label will read
        'Quantity Entered By User [Units]'
        """
        try:
            dum = self.semedian
        except AttributeError:
            return 'Error: No superposed epoch results to plot'

        from matplotlib.colors import LogNorm

        if len(xunits)<1:
            xlstr = '%s' % xquan
        else:
            xlstr = '%s [%s]' % (xquan, xunits)
        if len(yquan)>=1 and len(yunits)>=1:
            ylstr = '%s [%s]' % (yquan, yunits)
        elif len(yquan)>=1 and len(yunits)<1:
            ylstr = '%s' % yquan
        else:
            ylstr = ''

        usry = self.y
        if usrlimy:
            st = (self.y>=usrlimy[0]).nonzero()
            st = st[0][0]
            en = (self.y<=usrlimy[1]).nonzero()
            en = en[0][-1]
            pld = self.semedian[st:en+1,:]
            usry = self.y[st:en+1]

        fig = plt.figure(figsize=figsize)
        ax0 = fig.add_subplot(111)
        if zlog:
            cax = ax0.pcolor(self.x, usry, pld, norm=LogNorm(vmin=np.nanmin(pld), vmax=np.nanmax(pld)))
        else:
            cax = ax0.pcolor(self.x, usry, pld)

        plt.xlabel(xlstr)
        plt.ylabel(ylstr)

        if usrlimy:
            ax0.set_ylim(usrlimy)

        if epochline:
            ax0.axvline(0, 0, 1, color='k', ls=':', lw=1.5)

        ax0.set_xlim([self.x[0],self.x[-1]])
        hc1 = plt.colorbar(cax)
        if len(zunits)>=1:
            hc1.set_label(zunits)

        if show:
            plt.show()
            return None
        else:
            return ax0


def seadict(objlist, namelist):
    """Function to create dictionary of SeaPy.Sea objects.

    Parameters
    ==========
        - objlist: List of Sea objects.
        - namelist: List of variable labels for input objects.

    Other Parameters
    ================
    namelist = List containing names for y-axes.

    """
    try:
        assert type(objlist) == \
        list, 'seadict(): Inputs must be in lists'
    except AssertionError as args:
        raise ValueError('{0} -- {1}'.format(args.__class__.__name__, args))

    nobj, nname = len(objlist), len(namelist)
    if nobj!=nname:
        raise ValueError('seadict(): Lengths of object list and names list must be equal')
    else:
        pass

    outdict = {}
    for i  in range(nobj):
        outdict[namelist[i]] = objlist[i]

    return outdict


def multisea(dictobj, n_cols=1, epochline=True, usrlimx=[], usrlimy=[],
    xunits='', show=True, zunits='', zlog=True, figsize=None):
    """Function to create multipanel plot of superposed epoch analyses.

    Parameters
    ==========
    Dictionary of Sea objects (from superposedepoch.seadict()).

    Other Parameters
    ================
        - epochline (default = True) - put vertical line at zero epoch.
        - usrlimy (default = []) - override automatic y-limits on plot (same for all plots).
        - show (default = True) - shows plot; set to false to output plot object to variable
        - x/zunits - Units for labeling x and z axes, if required
        - figsize - tuple of (width, height) in inches
        - dpi (default=300) - figure resolution in dots per inch
        - n_cols - Number of columns: not yet implemented.

    Returns
    =======
    Plot of input object median and bounds (ci, mad, quartiles - see sea()).
    If keyword 'show' is False, output is a plot object.

    """
    import matplotlib.ticker as tik

    keys = list(dictobj.keys())
    pld, ply, plx, ylab = [], [], [], []
    qup, qlo = [], []
    keylist = sorted(keys)
    for key in keylist:
        tmp = dictobj[key].semedian
        pld.append(tmp) #pld (PLot Data)
        plx.append(dictobj[key].x) #plx (PLot X)
        qup.append(dictobj[key].bound_high)
        qlo.append(dictobj[key].bound_low)
        ylab.append(key) #ylab (Y-axis LABel)
        if tmp.ndim==2:
            ply.append(dictobj[key].y) #ply (Plot Y)

    if n_cols > 1:
        print('multisea(): Multiple column output not yet implemented')
        n_fig = len(dictobj) #force single column output for now
        fignum = list(range(1,n_fig+1))
        fig = plt.figure()
    else:
        n_fig = len(dictobj)
        fignum = list(range(1,n_fig+1))
        fig = plt.figure(figsize=figsize)

    for i in range(n_fig):
        #loop over each subplot
        #create subplot panel number
        dum = str(n_fig)+'1'+str(fignum[i])
        ax = fig.add_subplot(dum)
        if pld[i].ndim==2: #if semean/median 2D, do colour-plot
            if zlog:
                cax = ax0.pcolor(plx[i], ply[i], pld[i], norm=LogNorm(vmin=np.nanmin(pld), vmax=np.nanmax(pld)))
            else:
                cax = ax.pcolor(plx[i],ply[i],pld[i])
            plt.ylabel(ylab[i])
            hc1 = plt.colorbar(cax)
            if len(zunits)>=1:
                hc1.set_label(zunits)
        else: #else do line plot
            ax_b = dictobj[keylist[i]].plot(xquan='', yquan=ylab[i], epochline=epochline, show=False, target=ax, loc=dum)

        if i==n_fig-1:
            #Add x-label to bottom panel
            if xunits:
                xlab = 'Time since epoch [' + xunits + ']'
            else:
                xlab = 'Time since epoch'
            plt.xlabel(xlab)

        if usrlimx:
            plt.xlim(usrlimx)

        if usrlimy:
            plt.ylim(usrlimy)

        majorticky=tik.MaxNLocator(7)
        ax.yaxis.set_major_locator(majorticky)

    if show:
        plt.show()
        return None
    else:
        return fig


def readepochs(fname, iso=False, isofmt="%Y-%m-%dT%H:%M:%S"):
    """Read epochs from text file assuming YYYY MM DD hh mm ss format

    Parameters
    ==========
    Filename (include path)

    Other Parameters
    ================
    iso (default = False), read in ISO date format
    isofmt (default is YYYY-mm-ddTHH:MM:SS, code is %Y-%m-%dT%H:%M:%S)

    Returns
    =======
    epochs (type=list)

    """
    if iso==False:
        from numpy import loadtxt, zeros

        intime = loadtxt(fname)
        dum = zeros([intime.shape[0],6])
        for i in range(intime.shape[1]):
            dum[:,i] = intime[:,i]

        epochs=[]
        for i in range(dum.shape[0]):
            dtobj = dt.datetime(int(dum[i,0]), int(dum[i,1]), \
            int(dum[i,2]), int(dum[i,3]), int(dum[i,4]), int(dum[i,5]))
            epochs.append(dtobj)
    else:
        fh = open(fname,'r')
        intime = fh.readlines()
        epochs=[]
        for d in intime:
            d = d.strip()
            dtobj = dt.datetime.strptime(d, isofmt)
            epochs.append(dtobj)
    return epochs


def sea_signif(obj1, obj2, test='KS', show=True, xquan = 'Time Since Epoch',
        yquan='', xunits='', yunits='', epochline=True, usrlimy=[]):
    """Test for similarity between distributions at each lag in two 1-D SEAs

    Parameters
    ==========
    obj1 : Sea
        First instance for comparison
    obj2 : Sea
        Second instance for comparison


    Other Parameters
    ================
    test
        (default = 'KS') Test to apply at each lag:
        KS is 2-smaple Kolmogorov-Smirnov; U is Mann-Whitney U-test
    show
        (default = True)
    xquan
        (default = 'Time since epoch' (None)) - x-axis label.
    yquan
        (default = 'Time since epoch' (None)) - y-axis label.
    xunits
        (default = None (None)) - x-axis units.
    yunits
        (default = None (None)) - y-axis units.
    epochline
        (default = True) - put vertical line at zero epoch.
    usrlimy
        (default = []) - override automatic y-limits on plot.

    Examples
    ========
    >>> obj1 = seapy.Sea(data1, times1, epochs1)
    >>> obj2 = seapy.Sea(data2, times2, epochs2)
    >>> obj1.sea(storedata=True)
    >>> obj2.sea(storedata=True)
    >>> seapy.sea_signif(obj1, obj2)

    """
    try:
        assert isinstance(obj1, Sea)
        assert isinstance(obj2, Sea)
    except:
        raise TypeError("Inputs must both be seapy.Sea() instances of same length")

    keylist, S, prob = ['KS','U'], [], []
    try:
        assert test.upper() in keylist
    except:
        raise NotImplementedError("Test "+self.dtype+" not implemented, only "+str(keylist))

    from scipy.stats import ks_2samp, mannwhitneyu

    if test.upper() == 'KS':
        try:
            for x in range(obj1.datacube.shape[1]):
                tD, tprobKS2 = ks_2samp(obj1.datacube[:,x].compressed(),
                    obj2.datacube[:,x].compressed())
                S.append(tD)
                prob.append(tprobKS2)
        except:
            raise AttributeError('''KS significance testing requires datacube attribute\n
                Run sea() with kwarg storedata''')

    if test.upper() == 'U':
        for x,y in zip(obj1.datacube.transpose(), obj2.datacube.transpose()):
            tU, tprobU = mannwhitneyu(x, y)
            S.append(tU)
            prob.append(tprobU*2)

    if len(xunits)<1:
        xlstr = '{0}'.format(xquan)
    else:
        xlstr = '{0} [{1}]'.format(xquan, xunits)
    if len(yquan)>=1 and len(yunits)>=1:
        ylstr = '{0} [{1}]'.format(yquan, yunits)
    elif len(yquan)>=1 and len(yunits)<1:
        ylstr = '{0}'.format(yquan)
    else:
        ylstr = ''

    #set up plot panels
    fig = plt.figure()
    ax0 = fig.add_subplot(211)
    pos = ax0.get_position()
    ax0.set_position([pos.bounds[0], 0.4, pos.bounds[2], 0.55])
    ax1 = fig.add_subplot(212, position=[pos.bounds[0], 0.1, pos.bounds[2], 0.25])
    #overlay superposed epochs
    ax0.plot(obj1.x, obj1.semedian, lw=1.5)
    ax0.plot(obj1.x, obj1.bound_high, color='royalblue', ls='--')
    ax0.plot(obj1.x, obj1.bound_low, color='royalblue', ls='--')
    ax0.plot(obj2.x, obj2.semedian, color='crimson', lw=1.5)
    ax0.plot(obj2.x, obj2.bound_high, color='crimson', ls='--')
    ax0.plot(obj2.x, obj2.bound_low, color='crimson', ls='--')
    if epochline:
        ax0.axvline(0, 0, 1, color='k', ls=':')
    ax0.set_ylabel(ylstr)
    #plot prob in lower panel
    ax1.plot(obj1.x, prob, color='crimson', ls='-', lw=1.5, drawstyle='steps-mid')
    ax1.plot([obj1.x[0], obj1.x[-1]], [0.05, 0.05], color='royalblue', ls=':')
    if epochline:
        ax1.axvline(0, 0, 1, color='k', ls=':')
    ax1.set_xlabel(xlstr)
    ax1.set_ylabel('Prob. of H0')

    if show:
        plt.show()
        return None
    else:
        return fig

        #now make an output object to store the test stats and prob.
    #return None
