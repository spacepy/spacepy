#!/usr/bin/python
# -*- coding: utf-8 -*-

"""PoPPy -- Point Processes in Python.

This module contains point process class types and a variety of functions for 
association analysis. The routines given here grew from work presented by 
Morley and Freeman (Geophysical Research Letters, 34, L08104, doi:10.1029/ 
2006GL028891, 2007), which were originally written in IDL. This module is
intended for application to discrete time series of events to assess 
statistical association between the series and to calculate confidence limits.
Any mis-application or mis-interpretation by the user is the user's own fault.

>>> import datetime as dt
>>> import spacepy.time as spt

Since association analysis is rather computationally expensive, this example
shows timing.

>>> t0 = dt.datetime.now()
>>> onsets = spt.Ticktock(onset_epochs, 'CDF')
>>> ticksR1 = spt.Ticktock(tr_list, 'CDF')

Each instance must be initialized
        
>>> lags = [dt.timedelta(minutes=n) for n in xrange(-400,401,2)]
>>> halfwindow = dt.timedelta(minutes=10)
>>> pp1 = poppy.PPro(onsets.UTC, ticksR1.UTC, lags, halfwindow)

To perform association analysis

>>> pp1.assoc()
Starting association analysis                     
calculating association for series of length [3494, 1323] at 401 lags
>>> t1 = dt.datetime.now()
>>> print("Elapsed:  " + str(t1-t0))
Elapsed:  0:35:46.927138

Note that for calculating associations between long series at a large number of
lags is SLOW!!

To plot

>>> pp1.plot(dpi=80)
Error: No confidence intervals to plot - skipping

To add 95% confidence limits (using 4000 bootstrap samples)

>>> pp1.aa_ci(95, n_boots=4000)

The plot method will then add the 95% confidence intervals as a semi-
transparent patch.


--++-- By Steve Morley --++--

smorley@lanl.gov,
Los Alamos National Laboratory, ISR-1,
PO Box 1663, MS D466, Los Alamos, NM 87545
"""

import numpy as np
from matplotlib.mlab import prctile

from spacepy import help
__author__ = 'Steve Morley, Los Alamos National Lab (smorley@lanl.gov)'


#Try to pull in the C version. Assumption is that if you import this module,
#you want to do some association analysis, so the overhead in the import
#is OK.
try:
    import ctypes
    libpoppy = ctypes.CDLL('./libpoppy.so')
    dptr = ctypes.POINTER(ctypes.c_double)
    ulptr = ctypes.POINTER(ctypes.c_ulong)
    libpoppy.boots.restype = None
    libpoppy.boots.argtypes = [dptr, dptr, ctypes.c_ulong,
                               ctypes.c_ulong, ctypes.c_ulong,
                               ctypes.c_int]
    libpoppy.assoc.restype = None
    libpoppy.assoc.argtypes = [dptr, dptr, dptr, dptr,
                               ctypes.POINTER(ctypes.c_long),
                               ctypes.c_long, ctypes.c_long,
                               ctypes.c_long]
    libpoppy.aa_ci.restype = None
    libpoppy.aa_ci.argtypes=[ulptr, ulptr,
                             ctypes.c_ulong, ctypes.c_ulong, ctypes.c_ulong,
                             ctypes.c_ulong, ctypes.c_int]
    if hasattr(libpoppy, 'aa_ci_threaded'):
        libpoppy.aa_ci_threaded.restype = None
        libpoppy.aa_ci_threaded.argtypes=[ulptr, ulptr,
                                          ctypes.c_ulong, ctypes.c_ulong,
                                          ctypes.c_ulong,
                                          ctypes.c_ulong, ctypes.c_int,
                                          ctypes.c_int]
except:
    import bisect
    libpoppy = None


class PPro(object):
    """PoPPy point process object

    Initialize object with series1 and series2. These should be timeseries of
    events, given as lists, arrays, or lists of datetime objects.
    Includes method to perform association analysis of input series

    Output can be nicely plotted with plot method
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov
    Jonathan Niehof, Los Alamos National Lab, jniehof@lanl.gov
    """
    
    def __init__(self, process1, process2, lags=None, winhalf=None):
        self.process1 = process1
        self.process2 = process2
        self.lags = lags
        self.winhalf = winhalf
    
    def __str__(self):
        """String Representation of PoPPy object"""
        
        try:
            pk = max(self.assoc_total)
        except AttributeError:
            pk = 'N/A'
            
        try:
            asy = self.asympt_assoc
        except:
            asy = 'N/A'

        # TODO this is broken not enough values for the format statement
        return """Point Process Object:
        Points in process #1 - %d ; Points in process #2 - %d
        Peak association number - %s ; Asymptotic association - %s
        """ % (len(self.process1), len(self.process2))
    
    __repr__ = __str__
    
    def __len__(self):
        """Calling len(obj) will return the number of points in process 1"""
        
        return len(self.process1)
        
    def swap(self):
        """Swaps process 1 and process 2"""
        self.process1, self.process2 = self.process2, self.process1
        return None
        
    def assoc(self, u=None, h=None):
        """Perform association analysis on input series

        @param u: the time lags to use
        @type u: list
        @param h: association window half-width, same type as L{process1}
        """
        
        #check for existence of lags and winhalf
        try:
            if u:
                self.lags = u
            assert self.lags
            if h:
                self.winhalf = h
            assert self.winhalf
            print('calculating association for series of length %s at %d lags' \
                % ([len(self.process1), len(self.process2)], len(self.lags)))
        except:
            return 'assoc error: attributes lags and winhalf must be populated'

        import matplotlib as mpl
        from matplotlib import mlab
        if libpoppy == None:
            dtype = 'int64'
        else:
            dtype = 'int' + str(ctypes.sizeof(ctypes.c_long) * 8)
        
        ##Method 1 - use tb.tOverlap
        #create list for association number
        self.n_assoc = np.zeros((len(self.lags), len(self.process1)),
                                order='C', dtype=dtype)
        p2 = sorted(self.process2)
        starts = [t - self.winhalf for t in self.process1]
        stops = [t + self.winhalf for t in self.process1]

        if libpoppy == None:
            nss_list = list(range(len(self.process1)))
            for ilag in xrange(len(self.lags)):
                self.n_assoc[ilag, :] = [
                    bisect.bisect_right(p2, stops[nss] + self.lags[ilag]) - 
                    bisect.bisect_left(p2, starts[nss] + self.lags[ilag])
                    for nss in nss_list]
        else:
            p2 = (ctypes.c_double * len(p2))(*p2)
            starts = (ctypes.c_double * len(starts))(*starts)
            stops = (ctypes.c_double * len(stops))(*stops)
            lags = (ctypes.c_double * len(self.lags))(*self.lags)
            n_assoc = self.n_assoc.ctypes.data_as(ctypes.POINTER(ctypes.c_long))
            libpoppy.assoc(p2, starts, stops, lags, n_assoc,
                           len(p2), len(self.process1), len(self.lags))
        self.assoc_total = np.sum(self.n_assoc, axis=1)
        pul = mlab.prctile_rank(self.lags, p=(20,80))
        valsL = self.assoc_total[pul==0]
        valsR = self.assoc_total[pul==2]
        self.asympt_assoc = np.mean([np.mean(valsL), np.mean(valsR)])
        
        return None
    
    def plot(self, figsize=None, dpi=80, asympt=True, show=True, norm=True):
        """Create basic plot of association analysis.
        
        Uses object attributes created by the L{assoc} method and,
        optionally, L{aa_ci}.

        @param figsize: passed through to matplotlib.pyplot.figure
        @param dpi: passed through to matplotlib.pyplot.figure
        @param asympt: True to overplot the line of asymptotic association
                       number
        @type asympt: bool
        @param show: Show the plot? (if false, will create without showing)
        @type show: bool
        @param norm: Normalize plot to the asymptotic association number
        @type norm: bool
        """
        try:
            dum = self.n_assoc
        except AttributeError:
            return 'Error: No association analysis results to plot'
        
        import datetime as dt
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        from spacepy.toolbox import makePoly
        
        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax0 = fig.add_subplot(111)
        
        #fix such that self.lags (a timedelta) are converted to a time
        if type(self.lags[0]) == dt.timedelta:
            self.x = [i.seconds/60 + i.days*1440. for i in self.lags]
        else:
            self.x = self.lags

        ci = None
        if norm:
            if hasattr(self, 'ci'):
                ci = [[j / self.asympt_assoc for j in self.ci[i]]
                    for i in [0,1]]
            asympt_assoc = 1.0
            assoc_total = [assoc / self.asympt_assoc
                           for assoc in self.assoc_total]
        else:
            try:
                ci = self.ci
            except AttributeError:
                pass
            asympt_assoc = self.asympt_assoc
            assoc_total = self.assoc_total
        
        if asympt:
            ax0.plot([self.x[0], self.x[-1]], [asympt_assoc]*2, 'r--', lw=1.5)
        if ci != None:
            polyci = makePoly(self.x, ci[0], ci[1])
            ax0.add_patch(polyci)
        else:
            print('Error: No confidence intervals to plot - skipping')
        
        ax0.plot(self.x, assoc_total, 'b-', lw=1.5)
        
        if show:
            plt.show()
            return None
        else:
            return fig
    
    def aa_ci(self, inter, n_boots=1000, seed=None, threads=None):
        """Get bootstrap confidence intervals for association number
        
        Requires input of desired confidence interval, e.g.,
        >>> obj.aa_ci(95)
        
        Upper and lower confidence limits are added to the ci attribute
        """
        ci_low = np.zeros([len(self.lags)])
        ci_high = np.zeros([len(self.lags)])

        if libpoppy == None:
            if seed != None:
                np.random.seed(seed)
            for i in range(len(self.lags)):
                ci_low[i], ci_high[i] = boots_ci(self.n_assoc[i, :],
                                                 n_boots, inter, np.add.reduce,
                                                 seed=None)
        else:
            perc_low = (100.-inter)/2. #set confidence interval
            perc_high = inter + perc_low
            n_assoc_p = self.n_assoc.ctypes.data_as(
                ctypes.POINTER(ctypes.c_ulong))
            dtype = 'int' + str(ctypes.sizeof(ctypes.c_long) * 8)
            assoc_totals = np.zeros([len(self.lags), n_boots], dtype=dtype)
            assoc_totals_p = assoc_totals.ctypes.data_as(
                ctypes.POINTER(ctypes.c_ulong))
            if seed == None:
                seed = 0
                clock_seed = ctypes.c_int(1)
            else:
                clock_seed = ctypes.c_int(0)
            if threads != None and hasattr(libpoppy, 'aa_ci_threaded'):
                libpoppy.aa_ci_threaded(n_assoc_p, assoc_totals_p,
                                        len(self.lags), len(self.process1),
                                        n_boots, seed, clock_seed,
                                        threads)
            else:
                libpoppy.aa_ci(n_assoc_p, assoc_totals_p, len(self.lags),
                               len(self.process1), n_boots, seed, clock_seed)
            for i in range(len(self.lags)):
                ci_low[i], ci_high[i] = prctile(
                    assoc_totals[i, :], p=(perc_low,perc_high))
        self.ci = [ci_low, ci_high]
        return None

    def dump(self, f):
        """Dump this object out to a file

        @param f: open file (or file-like object) to dump to.
        """
        f.write('process1:' + str(len(self.process1)) + '\n')
        for val in self.process1:
            f.write(float(val).hex())
            f.write('\n')
        f.write('process2:' + str(len(self.process2)) + '\n')
        for val in self.process2:
            f.write(float(val).hex())
            f.write('\n')
        if self.lags != None:
            f.write('lags:' + str(len(self.lags)) + '\n')
            for val in self.lags:
                f.write(float(val).hex())
                f.write('\n')
        if self.winhalf != None:
            f.write('winhalf:' + float(self.winhalf).hex() + '\n')
        if hasattr(self, 'n_assoc'):
            f.write('asympt_assoc:' + float(self.asympt_assoc).hex() + '\n')
            assoc_shape = self.n_assoc.shape
            f.write('n_assoc:' + str(assoc_shape[0]) + ',' +
                    str(assoc_shape[1]) + '\n')
            for val in self.n_assoc.flat:
                f.write(float(val).hex() + '\n')
        if hasattr(self, 'ci'):
            f.write('ci:' + str(len(self.ci[0])) + '\n')
            for val in self.ci[0]:
                f.write(float(val).hex() + '\n')
            for val in self.ci[1]:
                f.write(float(val).hex() + '\n')

    @classmethod
    def load(cls, f):
        """Read a L{PPro} object from a file

        @param f: open file (or file-like object) to read from.
        @return: object stored in L{f}
        @rtype: L{PPro}
        """
        line = f.readline()
        (name, length) = line.split(':')
        if name != 'process1':
            raise ValueError('Unable to parse line ' + line)
        length = int(length)
        p1 = [0.0] * length
        for i in range(length):
            line = f.readline()
            p1[i] = float.fromhex(line)
            
        line = f.readline()
        (name, length) = line.split(':')
        if name != 'process2':
            raise ValueError('Unable to parse line ' + line)
        length = int(length)
        p2 = [0.0] * length
        for i in range(length):
            line = f.readline()
            p2[i] = float.fromhex(line)

        line = f.readline()
        (name, length) = line.split(':')
        if name == 'lags':
            length = int(length)
            lags = [0.0] * length
            for i in range(length):
                line = f.readline()
                lags[i] = float.fromhex(line)
            line = f.readline()
            (name, length) = line.split(':')
        else:
            lags = None
        if name == 'winhalf':
            winhalf = float.fromhex(length)
            line = f.readline()
            (name, length) = line.split(':')
        else:
            winhalf = None
        pop = PPro(p1, p2, lags, winhalf)

        if name == 'asympt_assoc':
            pop.asympt_assoc = float.fromhex(length)
            line = f.readline()
            (name, length) = line.split(':')

        if name == 'n_assoc':
            size1, size2 = length.split(',')
            size1 = int(size1)
            size2 = int(size2)
            pop.n_assoc = np.zeros([size1, size2])
            x = pop.n_assoc.flat
            for i in range(size1 * size2):
                line = f.readline()
                x[i] = float.fromhex(line)
            pop.assoc_total = np.sum(pop.n_assoc, axis=1)
            line = f.readline()
            (name, length) = line.split(':')

        if name == 'ci':
            length = int(length)
            pop.ci = [np.zeros([length]),
                      np.zeros([length])]
            for i in range(2):
                for j in range(length):
                    line = f.readline()
                    pop.ci[i][j] = float.fromhex(line)
        return pop

    
#Functions outside class
def boots_ci(data, n, inter, func, seed=None):
    """Construct bootstrap confidence interval
    
    The bootstrap is a statistical tool that uses multiple samples derived from
    the original data (called surrogates) to estimate a parameter of the 
    population from which the sample was drawn. This assumes that the sample is
    randomly drawn and hence is representative of the underlying distribution.
    The benefit of the bootstrap is that it is non-parametric and can be 
    applied in situations where there is reasonable doubt about the 
    characteristics of the underlying distribution. This routine uses the boot-
    strap for its most common application - the estimation of confidence 
    intervals.
    
    Example:
    ========
    >>> data, n = numpy.random.lognormal(mean=5.1, sigma=0.3, size=3000), 4000.
    >>> myfunc = lambda x: numpy.median(x)
    >>> ci_low, ci_high = poppy.boots_ci(data, n, 95, myfunc)
    >>> ci_low, numpy.median(data), ci_high
    (163.96354196633686, 165.2393331896551, 166.60491435416566) iter. 1
    ... repeat
    (162.50379144492726, 164.15218265100233, 165.42840588032755) iter. 2
    
    For comparison:
    ===============
    >>> data = numpy.random.lognormal(mean=5.1, sigma=0.3, size=90000)
    >>> numpy.median(data)
    163.83888237895815
    
    Note that the true value of the desired quantity may lie outside the
    95% confidence interval one time in 20 realizations. This occurred
    for the first iteration here.
    
    For the lognormal distribution, the median is found exactly by taking
    the exponential of the "mean" parameter. Thus here, the theoretical
    median is 164.022 (6 s.f.) and this is well captured by the above
    bootstrap confidence interval.
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov and
    Jonathan Niehof, Los Alamos National Lab, jniehof@lanl.gov

    @param data: data to bootstrap
    @type data: sequence
    @param n: number of surrogate series to select, i.e. number of bootstrap
              iterations.
    @type n: int
    @param inter: desired percentage confidence interval
    @type inter: numerical
    @param func: Function to apply to each surrogate series
    @type func: callable
    @param seed: Optional seed for the random number generator. If not
                 specified, numpy generator will not be reseeded;
                 C generator will be seeded from the clock.
    @type seed: int
    @return: L{inter} percent confidence interval on value derived from
             L{func} applied to the population sampled by L{data}.
    @rtype: sequence of float
    """
    perc_low = (100.-inter)/2. #set confidence interval
    perc_high = inter + perc_low
    
    n_els = len(data)
    if n_els <= 2:
        return np.nan, np.nan
    surr_quan = np.empty([n])

    if libpoppy == None:
        from numpy.random import randint
        if seed != None:
            np.random.seed(seed)
        ran_el = randint(n_els, size=[n, n_els])
        for i in range(int(n)): #compute n bootstrapped series
            surr_ser = [data[rec] for rec in ran_el[i, :]] #resample w/ replace
            surr_quan[i] = func(surr_ser) #get desired quantity from surrogates
    else:
        n = int(n)
        data = (ctypes.c_double * n_els)(*data)
        surr_ser = (ctypes.c_double * (n * n_els))()
        if seed == None:
            seed = 0
            clock_seed = ctypes.c_int(1)
        else:
            clock_seed= ctypes.c_int(0)
        libpoppy.boots(surr_ser, data,
                       ctypes.c_ulong(n), ctypes.c_ulong(n_els),
                       ctypes.c_ulong(seed), clock_seed)
        surr_quan = [func(surr_ser[i * n_els:(i  + 1) * n_els]) for i in xrange(n)]
    pul = prctile(surr_quan, p=(perc_low,perc_high)) #get confidence interval
    return pul[0], pul[1]
