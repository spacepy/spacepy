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

from spacepy import help
__author__ = 'Steve Morley, Los Alamos National Lab (smorley@lanl.gov)'
    
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
        import numpy as np
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
        temp2 = self.process1
        temp1 = self.process2
        self.process1 = temp1
        self.process2 = temp2
        
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

        import numpy as np
        import matplotlib as mpl
        from matplotlib import mlab
        import bisect
        
        ##Method 1 - use tb.tOverlap
        #create list for association number
        self.n_assoc = np.zeros((len(self.lags), len(self.process1)))
        p2 = sorted(self.process2)
        starts = [t - self.winhalf for t in self.process1]
        stops = [t + self.winhalf for t in self.process1]
        nss_list = list(range(len(self.process1)))
        for ilag in xrange(len(self.lags)):
            self.n_assoc[ilag, :] = [
                bisect.bisect_right(p2, stops[nss] + self.lags[ilag]) - 
                bisect.bisect_left(p2, starts[nss] + self.lags[ilag])
                for nss in nss_list]
        self.assoc_total = np.sum(self.n_assoc, axis=1)
        pul = mlab.prctile_rank(self.lags, p=(20,80))
        valsL = self.assoc_total[pul==0]
        valsR = self.assoc_total[pul==2]
        self.asympt_assoc = np.mean([np.mean(valsL), np.mean(valsR)])
        
        return None
    
    def plot(self, figsize=None, dpi=80, asympt=True, show=True):
        """Method called to create basic plot of association analysis.
        
        Inputs:
        =======
        Uses object attributes created by the obj.assoc() method.
        
        Optional keyword(s):
        ====================
        asympt (default = True) - add line of asymptotic assoc. number
        
        To Do:
        ======
        Add normalization keyword for plotting as n(u,h)/n_asympt
        """
        try:
            dum = self.n_assoc
        except AttributeError:
            return 'Error: No association analysis results to plot'
        
        import numpy as np
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
        
        if asympt:
            ax0.plot([self.x[0], self.x[-1]], [self.asympt_assoc]*2, 'r--', lw=1.5)
        try:
            dum = self.ci
            polyci = makePoly(self.x, self.ci[0], self.ci[1])
            ax0.add_patch(polyci)
        except AttributeError:
            print('Error: No confidence intervals to plot - skipping')
        
        ax0.plot(self.x, self.assoc_total, 'b-', lw=1.5)
        
        if show:
            plt.show()
            return None
        else:
            return fig
    
    def aa_ci(self, inter, n_boots=1000):
        """Get bootstrap confidence intervals for association number
        
        Requires input of desired confidence interval, e.g.,
        >>> obj.aa_ci(95)
        
        Upper and lower confidence limits are added to the ci attribute
        """
        import numpy as np
        
        aa_fun = lambda x: np.add.reduce(x)
        ci_low, ci_high = np.zeros([len(self.lags)]), np.zeros([len(self.lags)])
        for i in range(len(self.lags)):
            ci_low[i], ci_high[i] = boots_ci(self.n_assoc[i, :],
                                             n_boots, inter, aa_fun)
        self.ci = [ci_low, ci_high]
        return None

    
#Functions outside class
def boots_ci(data, n, inter, func):
    """Construct bootstrap confidence interval - caution: slow!
    
    Input:
    ======
        - n is number of surrogates;
        - data is data array (1D);
        - inter is desired confidence interval (e.g. 95%);
        - func is a user-defined function (lambda)
    
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
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov and
    Jonathan Niehof, Los Alamos National Lab, jniehof@lanl.gov
    """
    
    import numpy as np
    from numpy.random import randint
    from matplotlib.mlab import prctile

    perc_low = (100.-inter)/2. #set confidence interval
    perc_high = inter + perc_low
    
    n_els = len(data)
    if n_els <= 2:
        return np.nan, np.nan
    surr_quan = np.empty([n])
    ran_el = randint(0, n_els, size=[n, n_els])
    for i in range(int(n)): #compute n bootstrapped series
        surr_ser = [data[rec] for rec in ran_el[i, :]] #resample w/ replacemen
        surr_quan[i] = func(surr_ser) #get desired quantity from surrogates
    pul = prctile(surr_quan, p=(perc_low,perc_high)) #get confidence interval
    return pul[0], pul[1]
