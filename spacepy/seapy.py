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

smorley@lanl.gov/morley_steve@hotmail.com,
Los Alamos National Laboratory, ISR-1,
PO Box 1663, Los Alamos, NM 87545
"""
from __future__ import division
import numpy as np
import spacepy.toolbox as tb
from spacepy import help

__author__ = 'Steve Morley (smorley@lanl.gov/morley_steve@hotmail.com)'

class Sea(object):
    """SeaPy Superposed epoch analysis object

    Initialize object with data, times, epochs, window (half-width) and delta (optional).
    'times' and epochs should be in some useful format
    Includes method to perform superposed epoch analysis of input data series

    Output can be nicely plotted with plot method, or for multiple objects
    use the seamulti function
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    def __init__(self, data, times, epochs, window=3., delta=1.):
        import datetime as dt
        from matplotlib.dates import date2num
        self.data = np.array(data, dtype=float)
        self.times = times
        self.epochs = epochs
        if type(delta)==dt.timedelta:
            #t_delt = date2num(times[1])-date2num(times[0])
            t_delt = delta.days + delta.seconds/86400
            self.delta = t_delt
        else:
            self.delta = delta
        if type(window)==dt.timedelta:
            self.window = window.days + window.seconds/86400
        else:
            self.window = window
        self.window = np.ceil(float(self.window)/float(self.delta))
        if window!=self.window:
            print 'Window size changed to %s (points) to fit resolution' % self.window, '(%s)' % self.delta
        self.bound_type = None
    
    def __str__(self):
        """Define String Representation of Sea object"""
        
        return """Superposed Epoch Object:
        Data array - %s ; Epochs - %d ; Window - %d
        """ % (self.data.shape, len(self.epochs), self.window)
    
    __repr__ = __str__
    
    def __len__(self):
        """Calling len(obj) will return the number of epochs."""
        
        return len(self.epochs)
    
    def restoreepochs(self):
        """Replaces epoch times stored in obj.badepochs in the epochs attribute
        
        Quite why you'd want this feature I don't know.
        """
        
        try:
                dum = self.badepochs
        except AttributeError:
                return 'No bad epochs to restore'
        
        from numpy import union1d
        self.epochs = union1d(self.badepochs,self.epochs)
        return 'Bad epochs restored to obj.epochs attribute'
    
    def _timeepoch(self,delt):
        import datetime as dt
        from matplotlib.dates import date2num,num2date
        #check type of time input and throw error message
        el1,ep1 = self.times[0], self.epochs[0]
        el1num = ((type(el1) == int) or (type(el1) == float) or \
        (type(el1) == long) or (type(el1)==np.float64))
        ep1num = ((type(ep1) == int) or (type(ep1) == float) or \
        (type(ep1) == long) or (type(ep1)==np.float64))
        if (type(el1) == dt.datetime) and (type(el1) == type(self.epochs[0])):
            #both time and epochs are datetime objects
            #convert to serial time
            print '''converting to serial time'''
            dum = date2num(self.epochs)
            t_epoch = np.array(dum)
            dum = date2num(self.times)
            time = np.array(dum)
            ser_flag=False
        elif el1num and ep1num:
            #time is serial, do nothing
            print '''time is serial'''
            t_epoch = np.array(self.epochs, dtype=float)
            time = np.array(self.times, dtype=float)
            ser_flag=True
        else:
            return 'Error: Time and Epochs must be consistently typed (numeric/datetime)'
        
        lose0 = np.where(t_epoch > time[-1]-(self.window*self.delta))
        lose1 = np.where(t_epoch < time[0]+(self.window*self.delta))
        if len(lose0)>0 and len(lose1)>0:
            linds = np.union1d(lose0[0],lose1[0])
            if len(linds)>0:
                print 'sea(): %s out-of-range epochs moved to badepochs attribute'\
                % len(linds)
            if ser_flag:
                self.badepochs = t_epoch[linds]
            else:
                self.badepochs = num2date(t_epoch[linds])
            keep0 = np.where(t_epoch <= time[-1]-(self.window*self.delta))
            keep1 = np.where(t_epoch >= time[0]+(self.window*self.delta))
            kinds = np.intersect1d(keep0[0],keep1[0])
            #print 'first %s' % time[0], 'lose %s ' % t_epoch[linds], 'keep %s ' % t_epoch[kinds]
            if ser_flag: #if serial transform flagged, output to datetime obj.
                t_epoch = date2num(t_epoch[kinds])
                time = date2num(self.times)
            else:
                self.epochs = t_epoch[kinds]
            t_epoch = t_epoch[kinds]
            
        return time,t_epoch
    
    def sea(self, storedata=False, quartiles=True, ci=False, mad=False, 
        ci_quan='median'):
        """Method called to perform superposed epoch analysis on data in object.
        
        Inputs:
        =======
        Uses object attributes obj.data, obj.times, obj.epochs, obj.delta, obj.window,
        all of which must be available on instatiation.
        
        Optional keyword(s):
        ====================
            - storedata (default = False) - saves matrix of epoch windows as obj.datacube
            - quartiles calculates the quartiles as the upper and lower bounds (and is default);
            - ci will find the bootstrapped confidence intervals (and requires ci_quan to be set);
            - mad will use +/- the median absolute deviation for the bounds;
            - ci_quan can be set to 'median' or 'mean'
        
        A basic plot can be raised with the obj.plot() method
        """
        import numpy.ma as ma
        
        #ensure all input is np array
        delt = float(self.delta)
        y = np.array(self.data, dtype=float)
        
        if ci:
            quartiles, mad = False, False
        if mad:
            quartiles, ci = False, False
        
        time,t_epoch = self._timeepoch(delt)
            
        #build SEA matrix and perform analysis
        wind = self.window
        m = int(2*wind + 1)
        n = len(t_epoch)
        y_sea = np.zeros((n,m), dtype=float)
        for i in range(n):
            dif = np.abs(time-t_epoch[i])
            j = np.where(dif == np.min(dif))
            sea_slice = y[j[0][0]-wind:j[0][0]+wind+1]
            y_sea[i,0:] = sea_slice
        
        #find SEA mean, median and percentiles - exclude NaNs
        y_sea_m = ma.masked_where(np.isnan(y_sea), y_sea)
        self.semean = [np.mean(y_sea_m[:,i].compressed()) for i in range(m)]
        self.semedian = [np.median(y_sea_m[:,i].compressed()) for i in range(m)]
        self.semean, self.semedian = np.array(self.semean), np.array(self.semedian)
        self.bound_low = np.zeros((m,1))
        self.bound_high = np.zeros((m,1))
        
        if quartiles:
            from matplotlib.mlab import prctile
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                qul = prctile(dum,p=(25,75))
                self.bound_low[i], self.bound_high[i] = qul[0], qul[1]
        elif ci: #bootstrapped confidence intervals (95%)
            from spacepy.poppy import boots_ci
            if ci_quan == 'mean':
                ci_func = lambda x: np.mean(x)
            else:
                ci_func = lambda x: np.median(x)
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                self.bound_low[i], self.bound_high[i] = boots_ci(dum, 800, 95, ci_func)
        elif mad: #median absolute deviation
            from spacepy.toolbox import medAbsDev
            for i in range(m):
                dum = np.sort(y_sea_m[:,i].compressed())
                spread_mad = medAbsDev(data)
                self.bound_low[i] = self.semedian[i]-spread_mad
                self.bound_high[i] = self.semedian[i]+spread_mad
        
        self.x = np.linspace(-1.*self.window*self.delta, self.window*self.delta, \
         len(self.semedian))
        if storedata:
            self.datacube = y_sea_m
            print 'sea(): datacube added as new attribute'
        
        return 'Superposed epoch analysis complete'
        
    #def sea_norm(self, epoch2, nbins=100, storedata=False, quartiles=True, ci=False, mad=False, 
    #ci_quan='median'):
        #"""Method for normalized superposed epoch analysis (creates new object)
        
        #Example:
        #Assuming an already instantiated Sea object, xampl
        
        #>>> xampl_norm = xampl.sea_norm(epoch2)
        
        #Note that epoch2 must be a paired set with xampl.epochs
        
        #Inputs:
        #Uses object attributes obj.data, obj.times, obj.epochs, obj.delta.
        #all of which must be available on instatiation.
        #Also requires second set of epochs.
        
        #Optional keyword(s):
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
            #from matplotlib.mlab import prctile
            #for i in range(m):
                #dum = np.sort(y_sea_m[:,i].compressed())
                #qul = prctile(dum,p=(25,75))
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
                yunits='', epochline=False, usrlimy=[], show=True, 
                figsize=None, dpi=None):
        """Method called to create basic plot of superposed epoch analysis.
        
        Inputs:
        =======
        Uses object attributes created by the obj.sea() method.
        
        Optional keyword(s):
        ====================
            - x(y)quan (default = 'Time since epoch' (None)) - x(y)-axis label.
            - x(y)units (default = None (None)) - x(y)-axis units.
            - epochline (default = False) - put vertical line at zero epoch.
            - usrlimy (default = []) - override automatic y-limits on plot.
        
        If both ?quan and ?units are supplied, axis label will read
        'Quantity Entered By User [Units]'
        """
        try:
            dum = self.semedian
        except AttributeError:
            return 'Error: No superposed epoch results to plot'
        
        from spacepy.toolbox import makePoly
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        
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
        
        if show==True and dpi==None:
            dpi=80
        elif show==False and dpi==None:
            dpi=300
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax0 = fig.add_subplot(111)
        
        #ax0.plot(self.x, self.semean, 'b-', lw=1.5)
        #plt.hold(True)
        #ax0.plot(self.x, self.semedian, 'k-', lw=2.5)
        #ax0.plot(self.x, self.bound_low, 'k-.', lw=2)
        #ax0.plot(self.x, self.bound_high, 'k-.', lw=2)
        
        poly = makePoly(self.x, self.bound_low.ravel(), self.bound_high.ravel(), alpha=0.25)
        ax0.add_patch(poly)
        plt.hold(True)
        ax0.plot(self.x, self.semedian, 'k-', lw=2.0)
        ax0.plot(self.x, self.semean, 'r--', lw=1.25)
        
        plt.xlabel(xlstr)
        plt.ylabel(ylstr)
        
        if usrlimy:
            ax0.set_ylim(usrlimy)
        
        if epochline:
            yr = ax0.get_ylim()
            if yr[0] < 0:
                yrlo = yr[0]+yr[0]
            else:
                yrlo = yr[0]-yr[0]
            if yr[1] < 0:
                yrhi = yr[1]-yr[1]
            else:
                yrhi = yr[1]+yr[1]
            ax0.plot([0,0], [yrlo,yrhi], 'k:', lw=1)
            plt.ylim(yr)
        
        if show:
            plt.show()
            return None
        else:
            return fig


class Sea2d(Sea):
    """SeaPy 2D Superposed epoch analysis object
        
    Initialize object with data, times, epochs, window (half-width),
    delta (optional), and y (two-element vector with max and min of y;optional)
    'times' and epochs should be in some useful format
    Includes method to perform superposed epoch analysis of input data series
            
    Output can be nicely plotted with plot method, or for multiple
    objects use the seamulti function
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    def __init__(self, data, times, epochs, window=3., delta=1., y=[]):
        import datetime as dt
        self.data = np.array(data, dtype=float)
        self.times = times
        self.epochs = epochs
        if type(delta)==dt.timedelta:
            #t_delt = date2num(times[1])-date2num(times[0])
            t_delt = delta.days + delta.seconds/86400
            self.delta = t_delt
        else:
            self.delta = delta
        if type(window)==dt.timedelta:
            self.window = window.days + window.seconds/86400
        else:
            self.window = window
        self.window = np.ceil(float(self.window)/float(self.delta))
        if window!=self.window:
            print 'Window size changed to %s (points) to fit resolution' % self.window, '(%s)' % self.delta

        if y:
            self.y = np.linspace(y[0], y[1], data.shape[0]+1)
        else:
            self.y = np.linspace(0, data.shape[0]-1, data.shape[0]+1)

    def sea(self, storedata=False, quartiles=True, ci=False, mad=False, 
        ci_quan='median', nmask=1):
        """Method called to perform 2D superposed epoch analysis on data in object.
        
        Inputs:
        =======
        Uses object attributes obj.data, obj.times, obj.epochs, obj.delta, obj.window,
        all of which must be available on instatiation.
        
        Optional keyword(s):
        ====================
            - storedata (default = False) - saves matrix of epoch windows as obj.datacube
            - quartiles calculates the interquartile range to show the spread (and is default);
            - ci will find the bootstrapped confidence interval (and requires ci_quan to be set);
            - mad will use the median absolute deviation for the spread;
            - ci_quan can be set to 'median' or 'mean'
        
        A basic plot can be raised with the obj.plot() method
        """
        
        import numpy.ma as ma
        import datetime as dt
        from matplotlib.dates import date2num
        
        #ensure all input is np array or correct form
        delt = float(self.delta)
        y = np.array(self.data, dtype=float)
        time,t_epoch = self._timeepoch(delt)
        if not nmask:
            nmask = 0 #set mask to exclude none
        
        #build SEA matrix and perform analysis
        wind = self.window
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
        y_sea_m = ma.masked_where(y_sea < 0., y_sea)	
        #now get SEA quantities
        self.semean, self.semedian, self.countmask = np.empty((l,m)), np.empty((l,m)), np.empty((l,m))
        yj=0
        for ti in xrange(int(m)):
            for yj in xrange(l):
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
            print 'sea(): datacube added as new attribute'
        
        return 'Superposed epoch analysis complete'
            
    def plot(self, xquan = 'Time Since Epoch', yquan='', xunits='',
                yunits='', zunits='', epochline=False, usrlimy=[], 
                show=True, zlog=True, figsize=None, dpi=300):
        """Method called to create basic plot of 2D superposed epoch analysis.
        
        Inputs:
        =======
        Uses object attributes created by the obj.sea() method.
        
        Optional keyword(s):
        ====================
            - x(y)quan (default = 'Time since epoch' (None)) - x(y)-axis label.
            - x(y/z)units (default = None (None)) - x(y/z)-axis units.
            - epochline (default = False) - put vertical line at zero epoch.
            - usrlimy (default = []) - override automatic y-limits on plot.
            - show (default = True) - shows plot; set to false to output plot object to variable
            - figsize - tuple of (width, height) in inches
            - dpi (default=300) - figure resolution in dots per inch
            
        If both ?quan and ?units are supplied, axis label will read
        'Quantity Entered By User [Units]'
        """
        try:
            dum = self.semedian
        except AttributeError:
            return 'Error: No superposed epoch results to plot'
        
        import matplotlib as mpl
        import matplotlib.pyplot as plt
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
            #ax0.set_ylim(usrlimy)
        
        if show==True and dpi==None:
            dpi=80
        elif show==False and dpi==None:
            dpi=300
            
        fig = plt.figure(figsize=figsize, dpi=dpi)
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
            yr = ax0.get_ylim()
            if yr[0] < 0:
                yrlo = yr[0]+yr[0]
            else:
                yrlo = yr[0]-yr[0]
            if yr[1] < 0:
                yrhi = yr[1]-yr[1]
            else:
                yrhi = yr[1]+yr[1]
            ax0.plot([0,0], [yrlo,yrhi], 'k:', lw=2)
            plt.ylim(yr)
            
        ax0.set_xlim([self.x[0],self.x[-1]])
        hc1 = plt.colorbar(cax)
        if len(zunits)>=1:
            hc1.set_label(zunits)
            
        if show:
            plt.show()	
            return None
        else:
            return fig


def seadict(objlist, namelist):
    """Function to create dictionary of SeaPy.Sea objects.
    
    Inputs:
    =======
        - objlist: List of Sea objects.
        - namelist: List of variable labels for input objects.
        
    Optional keyword(s):
    ====================
    namelist = List containing names for y-axes.
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    try:
        assert type(objlist) == \
        list, 'seadict(): Inputs must be in lists'
    except AssertionError, args:
        return '%s -- %s' % (args.__class__.__name__, args)
    
    nobj, nname = len(objlist), len(namelist)
    if nobj!=nname:
        return 'seadict(): Lengths of object list and names list must be equal'
    else:
        pass
    
    outdict = {}
    for i  in xrange(nobj):
        outdict[namelist[i]] = objlist[i]
            
    return outdict
    

def multisea(dictobj, n_cols=1, epochline=False, usrlimx=[], usrlimy=[], 
    xunits='', show=True, zunits='', zlog=True, figsize=None, dpi=300):
    """Function to create multipanel plot of superposed epoch analyses.
    
    Inputs:
    =======
    Dictionary of Sea objects (from superposedepoch.seadict()).
    
    Optional keyword(s):
    ====================
        - epochline (default = False) - put vertical line at zero epoch.
        - usrlimy (default = []) - override automatic y-limits on plot (same for all plots).
        - show (default = True) - shows plot; set to false to output plot object to variable
        - x/zunits - Units for labelling x and z axes, if required
        - figsize - tuple of (width, height) in inches
        - dpi (default=300) - figure resolution in dots per inch
        - n_cols - Number of columns: not yet implemented.
    
    Output:
    =======
    Plot of input object median and bounds (ci, mad, quartiles - see sea()).
    If keyword 'show' is False, output is a plot object.
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    import matplotlib.ticker as tik
    import matplotlib.pyplot as plt
    
    keys = dictobj.keys()
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

    if show==True and dpi==None:
        dpi=80
    elif show==False and dpi==None:
        dpi=300
    
    if n_cols > 1:
        print 'multisea(): Multiple column output not yet implemented'
        n_fig = len(dictobj) #force single column output for now
        fignum = range(1,n_fig+1)
        fig = plt.figure()
    else:
        n_fig = len(dictobj)
        fignum = range(1,n_fig+1)
        fig = plt.figure(figsize=figsize, dpi=dpi)
    
    for i in xrange(n_fig):
        #loop over each subplot
        #create subplot panel number
        dum = str(n_fig)+'1'+str(fignum[i])
        #print dum
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
            ax.plot(plx[i],pld[i], 'k-', lw=2)
            plt.hold(True)
            ax.plot(plx[i],qup[i], 'k-', lw=1.5)
            ax.plot(plx[i],qlo[i], 'k-', lw=1.5)
            plt.ylabel(ylab[i])
        
        if i==n_fig-1:
            #Add x-label to bottom panel
            if xunits:
                xlab = 'Time since epoch [' + xunits + ']'
            else:
                xlab = 'Time since epoch'
            plt.xlabel(xlab)
        
        if epochline:
            yr = plt.ylim()
            if yr[0] < 0:
                yrlo = yr[0]+yr[0]
            else:
                yrlo = yr[0]-yr[0]
            if yr[1] < 0:
                yrhi = yr[1]-yr[1]
            else:
                yrhi = yr[1]+yr[1]
            ax.plot([0,0], [yrlo,yrhi], 'k:', lw=1)
            plt.ylim(yr)
        
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
    
    Input:
    ======
    Filename (include path)
    
    Optional inputs:
    ================
    iso (default = False), read in ISO date format
    isofmt (default is YYYY-mm-ddTHH:MM:SS, code is %Y-%m-%dT%H:%M:%S)
    
    Output:
    =======
    epochs (type=list)
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    import datetime as dt
    
    if iso==False:
        from numpy import loadtxt, zeros

        intime = loadtxt(fname)
        dum = zeros([intime.shape[0],6])
        for i in xrange(intime.shape[1]):
            dum[:,i] = intime[:,i]
        
        epochs=[]
        for i in xrange(dum.shape[0]):
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
    
    Inputs:
    =======
    Two seapy.Sea() instances for comparison
    
    Optional keyword(s):
    ====================
        - show (default = True)
        - x(y)quan (default = 'Time since epoch' (None)) - x(y)-axis label.
        - x(y)units (default = None (None)) - x(y)-axis units.
        - epochline (default = True) - put vertical line at zero epoch.
        - usrlimy (default = []) - override automatic y-limits on plot.
    
    Example:
    ========
    >>> obj1 = seapy.Sea(data1, times1, epochs1)
    >>> obj2 = seapy.Sea(data2, times2, epochs2)
    >>> obj1.sea(storedata=True)
    >>> obj2.sea(storedata=True)
    >>> D, probH0 = seapy.sea_signif(obj1, obj2)
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    """
    try:
        assert isinstance(obj1, Sea)
        assert isinstance(obj2, Sea)
    except:
        raise Exception("Inputs must both be seapy.Sea() instances of same length")
    
    keylist, S, prob = ['KS','U'], [], []
    try:
        assert test.upper() in keylist
    except:
        raise TypeError("Test "+self.dtype+" not implemented, only "+str(keylist))
    
    from scipy.stats import ks_2samp, mannwhitneyu
    
    if test.upper() == 'KS':
        try:
            for x in xrange(obj1.datacube.shape[0]):
                tD, tprobKS2 = ks_2samp(obj1.datacube[:,x].compressed(), 
                    obj2.datacube[:,x].compressed())
                S.append(tD)
                prob.append(tprobKS2)        
        except:
            raise Exception('''KS significance testing requires datacube attribute\n
                Run sea() with kwarg storedata''')
    
    if test.upper() == 'U':
        for x,y in zip(obj1.semedian, obj2.semedian):
            tU, tprobU = mannwhitneyu(x, y)
            S.append(tU)
            prob.append(tprobU*2)
    
    def add_epochline(ax):
        yr = ax.get_ylim()
        if yr[0] < 0:
            yrlo = yr[0]+yr[0]
        else:
            yrlo = yr[0]-yr[0]
        if yr[1] < 0:
            yrhi = yr[1]-yr[1]
        else:
            yrhi = yr[1]+yr[1]
        
        return (yr, yrlo, yrhi)
    
    if show:
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
            
        import matplotlib.pyplot as plt
        #set up plot panels
        fig = plt.figure()
        ax0 = fig.add_subplot(211)
        pos = ax0.get_position()
        ax0.set_position([pos.bounds[0], 0.4, pos.bounds[2], 0.55])
        ax1 = fig.add_subplot(212, position=[pos.bounds[0], 0.1, pos.bounds[2], 0.25])
        #overlay superposed epochs
        ax0.plot(obj1.x, obj1.semedian, 'b-', lw=1.5)
        plt.hold(True)
        ax0.plot(obj1.x, obj1.bound_high, 'b--')
        ax0.plot(obj1.x, obj1.bound_low, 'b--')
        ax0.plot(obj2.x, obj2.semedian, 'r-', lw=1.5)
        ax0.plot(obj2.x, obj2.bound_high, 'r--')
        ax0.plot(obj2.x, obj2.bound_low, 'r--')
        if epochline:
            yr, yrlo, yrhi = add_epochline(ax0)
            ax0.plot([0,0], [yrlo,yrhi], 'k:', lw=1)
            ax0.set_ylim(yr)
        ax0.set_ylabel(ylstr)
        #plot prob in lower panel
        ax1.plot(obj.x, prob, 'r-', lw=1.5, drawstyle='steps-mid')
        ax1.plot([obj.x[0], obj.x[-1]], [0.05, 0.05], 'b-')
        if epochline:
            yr, yrlo, yrhi = add_epochline(ax1)
            ax1.plot([0,0], [yrlo,yrhi], 'k:', lw=1)
            ax1.set_ylim(yr)
        ax1.set_xlabel(xlstr)
        ax1.set_ylabel('Prob. of H0')
        
        plt.show()
        
    return S, prob
    