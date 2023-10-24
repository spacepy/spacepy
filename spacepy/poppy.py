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

>>> lags = [dt.timedelta(minutes=n) for n in range(-400,401,2)]
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


Authors: Steve Morley and Jon Niehof
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov, jniehof@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.
"""

import bisect
import ctypes
import sys
import warnings
import datetime as dt

import numpy as np

from spacepy import help
from spacepy import lib
import spacepy.toolbox as tb
import spacepy.datamodel as dm

__contact__ = 'Steve Morley, smorley@lanl.gov'

class PPro(object):
    """PoPPy point process object

    Initialize object with series1 and series2. These should be timeseries of
    events, given as lists, arrays, or lists of datetime objects.
    Includes method to perform association analysis of input series

    Output can be nicely plotted with :py:meth:`plot`.

    .. currentmodule:: spacepy.poppy
    
    .. autosummary::
        ~PPro.aa_ci
        ~PPro.assoc
        ~PPro.assoc_mult
        ci
        conf_above
        ~PPro.plot
        ~PPro.plot_mult
        ~PPro.swap

    .. automethod:: aa_ci
    .. automethod:: assoc
    .. automethod:: assoc_mult

    .. attribute:: ci
    
       Contains the upper and lower confidence limits for the association
       number as a function of lag. The first element is the array of lower
       limits; the second, the array of upper limits. Not available until
       after calling :meth:`aa_ci`.

    .. attribute:: conf_above
    
       Contains the confidence that the association number, as a function
       of lag, is above the asymptotic association number. (The confidence
       of being below is 100 - ``conf_above``.)  Not available until
       after calling :meth:`aa_ci`.

    .. automethod:: plot
    .. automethod:: plot_mult
    .. automethod:: swap
    """
    #NB: P2 is the "master" timescale, P1 gets shifted by lags
    #Add lag to p1 to reach p2's timescale, subtract lag from p2 to reach p1's

    ci = None
    """Upper and lower confidence limits for the association number"""

    conf_above = None
    """Confidence that the association number is above the asymptotic"""

    def __init__(self, process1, process2, lags=None, winhalf=None, verbose=False):
        self.process1 = process1
        self.process2 = process2
        self.lags = lags
        self.winhalf = winhalf
        self.verbose = verbose

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
        """ % (len(self.process1), len(self.process2), pk, asy)

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

        Parameters
        ==========
        u : list, optional
            the time lags to use
        h :
            association window half-width, same type as process1
        """

        #check for existence of lags and winhalf
        if u is not None:
            self.lags = u
        if h is not None:
            self.winhalf = h
        if self.lags is None or self.winhalf is None:
            if self.verbose:
                print('assoc error: attributes lags and winhalf must be populated')
            return
        if self.verbose:
            print('calculating association for series of length %s at %d lags'
                  % ([len(self.process1), len(self.process2)], len(self.lags)))

        import matplotlib as mpl
        import matplotlib.dates as mpd
        dtype = 'int' + str(ctypes.sizeof(ctypes.c_long) * 8
                            if lib.have_libspacepy else 64)

        ##Method 1 - use tb.tOverlap
        #create list for association number
        self.n_assoc = np.empty((len(self.lags), len(self.process1)),
                                order='C', dtype=dtype)
        p2 = sorted(self.process2)
        p1 = sorted(self.process1)

        if lib.have_libspacepy == False:
            winhalf = self.winhalf
            lags = self.lags
            starts = [t - self.winhalf for t in p1]
            stops = [t + self.winhalf for t in p1]
            nss_list = list(range(len(self.process1)))
            for ilag in range(len(self.lags)):
                last_idx = [bisect.bisect_right(p2, stops[nss] + self.lags[ilag])
                            for nss in nss_list]
                first_idx = [bisect.bisect_left(p2, starts[nss] + self.lags[ilag])
                             for nss in nss_list]
                self.n_assoc[ilag, :] = [last_idx[i] - first_idx[i] for i in nss_list]
        else:
            def fracday(dd):
                '''turn a timedelta into a fractional day'''
                #NOTE: If 2.6 compatibility is ever dropped, can use dd.total_seconds()/86400
                return dd.days + dd.seconds/86400.0 + dd.microseconds/(86400.0 * 10**6)
            #test for datetime input and try to convert to numeric
            if isinstance(p1[0], dt.datetime):
                p1 = [mpd.date2num(pp) for pp in p1]
            if isinstance(p2[0], dt.datetime):
                p2 = [mpd.date2num(pp) for pp in p2]
            #TODO: convert lags and winhalf from timedelta to fractional days
            if isinstance(self.lags[0], dt.timedelta):
                lags = [fracday(dd) for dd in self.lags]
            else:
                lags = self.lags
            if isinstance(self.winhalf, dt.timedelta):
                winhalf = fracday(self.winhalf)
            else:
                winhalf = self.winhalf
            p2 = (ctypes.c_double * len(p2))(*p2)
            p1 = (ctypes.c_double * len(p1))(*p1)
            #prctile_rank dies on some versions if this is a ctypes array
            lags = np.require(lags, dtype=ctypes.c_double, requirements='C')
            n_assoc = self.n_assoc.ctypes.data_as(lib.lptr)
            lib.assoc(p2, p1, lags.ctypes.data_as(lib.dptr), n_assoc,
                      winhalf, len(p2), len(p1), len(self.lags))
        left20perc = int(np.round((len(lags)*0.2)))
        right20perc = int(np.round((len(lags)*0.8)))
        self.assoc_total = np.sum(self.n_assoc, axis=1)
        valsL = self.assoc_total[:left20perc]
        valsR = self.assoc_total[right20perc:]
        self.asympt_assoc = np.mean([np.mean(valsL), np.mean(valsR)])

        self.expected = np.empty([len(self.lags)], dtype='float64')
        for i in range(len(self.lags)):
            start_time = max([p1[0] + lags[i], p2[0]]) - winhalf
            stop_time = min([p1[-1] + lags[i], p2[-1]]) + winhalf
            if start_time != stop_time:
                n1 = bisect.bisect_right(p1, stop_time - lags[i]) - \
                     bisect.bisect_left(p1, start_time - lags[i])
                n2 = bisect.bisect_right(p2, stop_time) - \
                     bisect.bisect_left(p2, start_time)
                self.expected[i] = 2.0 * n1 * n2 * winhalf / \
                                       (stop_time - start_time)
            else:
                self.expected[i] = 0.0
        return None

    def assoc_mult(self, windows, inter=95, n_boots=1000, seed=None):
        """Association analysis w/confidence interval on multiple windows

        Using the time sequence and lags stored in this object, perform
        full association analysis, including bootstrapping of confidence
        intervals, for every listed window half-size

        Parameters
        ==========
        windows : sequence
            window half-size for each analysis
        inter : float, optional
            desired confidence interval, default 95
        n_boots : int, optional
            number of bootstrap iterations, default 1000
        seed : int, optional
            Random number generator seed. It is STRONGLY
            recommended not to specify (i.e. leave None) to permit
            multithreading.

        Warnings
        ========
        This function is likely to take a LOT of time.

        Returns
        =======
        out : three numpy array
            Three numpy arrays, (windows x lags), containing (in order)
            low values of confidence interval, high values of ci,
            percentage confidence above the asymptotic association number
        """
        n_lags = len(self.lags)
        ci_low = np.empty([len(windows), n_lags])
        ci_high = np.empty([len(windows), n_lags])
        percentiles = np.empty([len(windows), n_lags])
        for i in range(len(windows)):
            window = windows[i]
            self.assoc(h=window)
            self.aa_ci(inter, n_boots, seed)
            ci_low[i, :] = self.ci[0]
            ci_high[i, :] = self.ci[1]
            percentiles[i, :] = self.conf_above
        return (ci_low, ci_high, percentiles)

    def plot_mult(self, windows, data, min=None, max=None, cbar_label=None,
                  figsize=None, dpi=80,
                  xlabel='Lag', ylabel='Window Size'):
        """Plots a 2D function of window size and lag

        Parameters
        ==========
        windows : list
            list of window sizes (y axis)
        data : list
            list of data, dimensioned (windows x lags)
        min : float, optional
            clip L{data} to this minimum value
        max : float, optional
            clip L{data} to this maximum value
        """
        import matplotlib.pyplot as plt

        x = np.array(tb.bin_center_to_edges(self.lags))
        y = np.array(tb.bin_center_to_edges(windows))

        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax0 = fig.add_subplot(111)
        cax = ax0.pcolormesh(x, y, data, vmin=min, vmax=max,
                             shading='flat')
        ax0.set_xlim((x[0], x[-1]))
        ax0.set_ylim((y[0], y[-1]))
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if cbar_label is None:
            cbar_label = '{} confident above asymptotic association'.format(
                r'\%' if plt.rcParams['text.usetex'] else r'%')
        plt.colorbar(cax, fraction=0.05).set_label(cbar_label)
        return fig

    def plot(self, figsize=None, dpi=80, asympt=True, show=True, norm=True,
             xlabel='Time lag', xscale=None, ylabel=None, title=None, transparent=True):
        """Create basic plot of association analysis.

        Uses object attributes created by :py:meth:`assoc` and,
        optionally, :py:meth:`aa_ci`.

        Parameters
        ==========
        figsize : , optional
            passed through to matplotlib.pyplot.figure
        dpi : int, optional
            passed through to matplotlib.pyplot.figure
        asympt : boolean, optional
            True to overplot the line of asymptotic association number
        show : boolean, optional
            Show the plot? (if false, will create without showing)
        norm : boolean, optional
            Normalize plot to the asymptotic association number
        title : string, optional
            label/title for the plot
        xlabel : string, optional
            label to put on the X axis of the resulting plot
        xscale : float, optional
            scale x-axis by this factor (e.g. 60.0 to convert seconds to minutes)
        ylabel : string, optional
            label to put on the Y axis of the resulting plot
        transparent : boolean, optional
            make c.i. patch transparent (default)
        """
        try:
            dum = self.n_assoc
        except AttributeError:
            return 'Error: No association analysis results to plot'

        import datetime as dt
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=figsize, dpi=dpi)
        ax0 = fig.add_subplot(111)

        #fix such that self.lags (a timedelta) are converted to a time
        if type(self.lags[0]) == dt.timedelta:
            x = [i.seconds/60 + i.days*1440. for i in self.lags]
        else:
            x = self.lags
        if xscale != None:
            x = [i / xscale for i in x]
        ax0.set_xlim((min(x), max(x)))

        ci = self.ci
        if norm:
            if self.ci is not None:
                ci = [[j / self.asympt_assoc for j in self.ci[i]]
                    for i in [0,1]]
            asympt_assoc = 1.0
            assoc_total = [assoc / self.asympt_assoc
                           for assoc in self.assoc_total]
        else:
            asympt_assoc = self.asympt_assoc
            assoc_total = self.assoc_total

        if ci is not None:
            if transparent:
                ax0.fill_between(x, ci[0], ci[1],
                                 edgecolor='none', facecolor='blue', alpha=0.5)
            else:
                ax0.fill_between(x, ci[0], ci[1],
                                 edgecolor='none', facecolor='#ABABFF')
        else:
            print('Error: No confidence intervals to plot - skipping')

        ax0.plot(x, assoc_total, 'b-', lw=1.0)
        if asympt:
            ax0.plot([x[0], x[-1]], [asympt_assoc]*2, 'r--', lw=1.0)
        if ylabel is None:
            if norm:
                plt.ylabel(
                    'Normalized Association Number n(u, h={0}) / n({1}, h={0})'.format(
                    self.winhalf,
                    r'$\mathrm{u\rightarrow\infty}$'
                    if plt.rcParams['text.usetex'] else 'u->Inf'))
            else:
                plt.ylabel('Association Number n(u, h={0})'.format(
                    self.winhalf))
        else:
            plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        if title is not None:
            plt.title(title)

        if show:
            plt.show()
            return None
        return fig

    def aa_ci(self, inter, n_boots=1000, seed=None):
        """Get bootstrap confidence intervals for association number

        Requires input of desired confidence interval, e.g.:
            >>> obj.aa_ci(95)

        Upper and lower confidence limits are added to :attr:`~PPro.ci`.

        After calling, :attr:`~PPro.conf_above` will contain the confidence
        (in percent) that
        the association number at that lag is *above* the asymptotic
        association number. (The confidence of being below is 100 - conf_above)
        For minor variations in conf_above to be meaningful, a *large* number
        of bootstraps is required. (Rougly, 1000 to be meaningful to the
        nearest percent; 10000 to be meaningful to a tenth of a percent.) A
        conf_above of 100 usually indicates an insufficient sample size to
        resolve, *not* perfect certainty.

        Note also that a 95% chance of being above indicates an exclusion
        from the *90%* confidence interval!

        Parameters
        ==========
        inter : float
            percentage confidence interval to calculate
        n_boots : int, optional
            number of bootstrap iterations to run
        seed : int, optional
            seed for the random number generator. If not specified,
            Python code will use numpy's RNG and its current seed;
            C code will seed from the clock.

        Warnings
        ========
        If ``seed`` is specified, results may not be reproducible between
        systems with different sizes for C long type. Note that 64-bit
        Windows uses a 32-bit long and so results will be the same between
        64 and 32-bit Windows, but not between 64-bit Windows and other 64-bit
        operating systems. If ``seed`` is not specified, results are
        not reproducible anyhow.
        """
        lags = self.lags

        ci_low = np.empty([len(lags)])
        ci_high = np.empty([len(lags)])
        conf_above = np.empty([len(lags)])
        long_size = ctypes.sizeof(ctypes.c_long) * 8

        if seed is not None:
            np.random.seed(seed)
            # numpy random seeds must always be 32-bit
            seed_size = long_size if lib.have_libspacepy else 32
            minseed = -2 ** (seed_size - 1)
            maxseed = 2 ** (seed_size - 1) - 1
            #randint used to be system-size signed integer only.
            #so used that and cast to the required unsigned later
            #cast does not lose entropy: negative numbers map to high positives.
            #For reproducibility, keep doing that even though dtype
            #kwarg now available.
            lag_seeds = np.random.randint(minseed, maxseed, [len(lags)])
            if not lib.have_libspacepy:
                lag_seeds = np.require(lag_seeds, np.int32)
            newtype = np.dtype('u' + str(lag_seeds.dtype))
            lag_seeds = np.require(lag_seeds, dtype=newtype)
        if not lib.have_libspacepy:
            for i in range(len(lags)):
                if seed != None:
                    np.random.seed(lag_seeds[i])
                ci_low[i], ci_high[i], conf_above[i]= boots_ci(
                    self.n_assoc[i, :],
                    n_boots, inter, np.add.reduce,
                    seed=None, target=self.asympt_assoc)
        else:
            perc_low = (100.-inter)/2. #set confidence interval
            perc_high = inter + perc_low
            dtype = 'int' + str(long_size)
            assoc_totals = np.empty([len(lags), n_boots],
                                    dtype=dtype, order='C')
            if seed is None:
                clock_seed = ctypes.c_int(1)
                lag_seeds = np.empty([len(lags)], dtype=dtype)
            else:
                clock_seed = ctypes.c_int(0)
            def thread_targ(start, size):
                lib.aa_ci(
                    self.n_assoc[start:start+size].ctypes.data_as(lib.ulptr),
                    assoc_totals[start:start+size].ctypes.data_as(lib.ulptr),
                    size, len(self.process1), n_boots,
                    lag_seeds[start:start+size].ctypes.data_as(lib.ulptr),
                    clock_seed)
            tb.thread_job(len(lags), 0, thread_targ)
            for i in range(len(lags)):
                assoc_totals[i, :].sort()
                ci_low[i], ci_high[i] = np.percentile(
                    assoc_totals[i, :], (perc_low,perc_high))
                conf_above[i] = 100.0 - value_percentile(assoc_totals[i, :],
                                                         self.asympt_assoc)
        self.ci = [ci_low, ci_high]
        self.conf_above = conf_above


#Functions outside class

def plot_two_ppro(pprodata, pproref, ratio=None, norm=False,
                  title=None, xscale=None, figsize=None, dpi=80,
                  ylim=[None, None], log=False, xticks=None, yticks=None):
    """Overplots two PPro objects

    Parameters
    ==========
    pprodata : PPro
        first point process to plot (in blue)
    pproref : PPro
        second process to plot (in red)
    ratio : float
        multiply L{pprodata} by this ratio before plotting,
        useful for comparing processes of different magnitude
    norm : boolean
        normalize everything to L{pproref}, i.e. the association
        number for L{pproref} will always plot as 1.
    title : string
        title to put on the plot
    xscale : float
        scale x-axis by this factor (e.g. 60.0 to convert
        seconds to minutes)
    figsize :
        passed through to matplotlib.pyplot.figure
    dpi : int
        passed through to matplotlib.pyplot.figure
    ylim : list
        [minimum, maximum] values of y for the axis
    log : bollean
        True for a log plot
    xticks : sequence or float
        if provided, a list of tickmarks for the X axis
    yticks : sequance or float
        if provided, a list of tickmarks for the Y axis
    """
    import matplotlib.pyplot as plt
    if ratio is None:
        ratio = float(pproref.asympt_assoc) / pprodata.asympt_assoc
    lags = pproref.lags
    nlags = len(lags)
    assert lags[:] == pprodata.lags[:]
    if xscale is not None:
        lags = [float(i) / xscale for i in lags]
    fig = plt.figure(figsize=figsize, dpi=dpi)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    ax0 = fig.add_subplot(111)
    ax0.set_xlim((lags[0], lags[-1]))
    if norm:
        scaleddata = [ratio *
                      float(pprodata.assoc_total[i]) / pproref.assoc_total[i]
                      for i in range(nlags)]
        scaledhi = [ratio *
                    float(pprodata.ci[1][i]) / pproref.assoc_total[i]
                    for i in range(nlags)]
        scaledlo = [ratio *
                    float(pprodata.ci[0][i]) / pproref.assoc_total[i]
                    for i in range(nlags)]
        scaledref = [1.0 for i in range(nlags)]
        refhi = [float(pproref.ci[1][i])  / pproref.assoc_total[i]
                 for i in range(nlags)]
        reflo = [float(pproref.ci[0][i])  / pproref.assoc_total[i]
                 for i in range(nlags)]
    else:
        scaleddata = [ratio * float(pprodata.assoc_total[i])
                      for i in range(nlags)]
        scaledhi = [ratio * float(pprodata.ci[1][i])
                    for i in range(nlags)]
        scaledlo = [ratio * float(pprodata.ci[0][i])
                    for i in range(nlags)]
        scaledref = pproref.assoc_total
        refhi = pproref.ci[1]
        reflo = pproref.ci[0]
    ax0.fill_between(lags, reflo, refhi,
                     edgecolor='none', facecolor='#FFABAB', interpolate=True)
    ax0.fill_between(lags, scaledlo, scaledhi,
                     edgecolor='none', facecolor='#ABABFF', interpolate=True)
    bottom = np.fromiter((max([scaledlo[i], reflo[i]]) for i in range(nlags)),
                         np.float64, count=-1)
    top = np.fromiter((min([scaledhi[i], refhi[i]]) for i in range(nlags)),
                      np.float64, count=-1)
    ax0.fill_between(lags, bottom, top, where=(bottom <= top),
                     edgecolor='none', facecolor='#AB81D5',
                     interpolate=True)
    ax0.plot(lags, scaleddata, lw=1.0)
    ax0.plot(lags, scaledref, 'r--', lw=1.0)
    ax0.set_ylim(bottom = 0 if ylim[0] is None else ylim[0])
    if ylim[1] is not None:
        ax0.set_ylim(top=ylim[1])
    if log:
        ax0.set_yscale('log', nonposy='clip')
    if xticks:
        ax0.set_xticks(xticks)
    if yticks:
        ax0.set_yticks(yticks)
    if title:
        plt.title(title)


def boots_ci(data, n, inter, func, seed=None, target=None, sample_size=None, usepy=False, nretvals=1):
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

    Examples
    ========
    >>> data, n = numpy.random.lognormal(mean=5.1, sigma=0.3, size=3000), 4000.
    >>> myfunc = lambda x: numpy.median(x)
    >>> ci_low, ci_high = poppy.boots_ci(data, n, 95, myfunc)
    >>> ci_low, numpy.median(data), ci_high
    (163.96354196633686, 165.2393331896551, 166.60491435416566) iter. 1
    ... repeat
    (162.50379144492726, 164.15218265100233, 165.42840588032755) iter. 2

    For comparison

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

    Parameters
    ==========
    data : array like
        data to bootstrap
    n : int
        number of surrogate series to select, i.e. number of bootstrap iterations.
    inter : numerical
        desired percentage confidence interval
    func : callable
        Function to apply to each surrogate series
    sample_size : int
        number of samples in the surrogate series, default
        length of L{data}. This will change the statistical
        properties of the bootstrap and should only be used
        for good reason!
    seed : int
        Optional seed for the random number generator. If not
        specified, numpy generator will not be reseeded;
        C generator will be seeded from the clock.
    target : same as data
        a 'target' value. If specified, will also calculate
        percentage confidence of being at or above this value.
    nretvals : int
        number of return values from input function

    Returns
    =======
    out : sequence of float
        inter percent confidence interval on value derived from
        func applied to the population sampled by data.
        If target is specified, also the percentage confidence of
        being above that value.
    """
    perc_low = (100.-inter)/2. #set confidence interval
    perc_high = inter + perc_low

    n_els = len(data)
    if n_els <= 2:
        return np.nan, np.nan if target is None else np.nan, np.nan, np.nan
    if sample_size is None:
        sample_size = n_els
    if not lib.have_libspacepy or usepy:
        surr_quan = np.empty([n, nretvals] if nretvals > 1 else [n])
        if seed is not None:
            np.random.seed(seed)
        ran_el = np.random.randint(n_els, size=[n, sample_size])
        for i in range(int(n)): #compute n bootstrapped series
            surr_ser = np.array([data[rec] for rec in ran_el[i, :]]) #resample w/ replace
            surr_quan[i] = func(surr_ser) #get desired quantity from surrogates
        if len(surr_quan.shape) == 1:
            surr_quan.sort()
        else:
            surr_quan = surr_quan[surr_quan.argsort(axis=0)[:,0]]
    else:
        n = int(n)
        data = (ctypes.c_double * n_els)(*data)
        surr_ser = (ctypes.c_double * (n * sample_size))()
        if seed is None:
            seed = 0
            clock_seed = ctypes.c_int(1)
        else:
            clock_seed= ctypes.c_int(0)
        lib.boots(surr_ser, data,
                  ctypes.c_ulong(n), ctypes.c_ulong(n_els),
                  ctypes.c_ulong(sample_size),
                  ctypes.c_ulong(seed), clock_seed)
        gen = [func(surr_ser[i * sample_size:(i  + 1) * sample_size])
               for i in range(n)]
        #Can't sort if more than one return value
        surr_quan = sorted(gen) if nretvals == 1 else np.stack(gen)
    #get confidence interval
    if nretvals>1:
        pul = np.empty((2, nretvals))
        for nn in range(nretvals):
            pul[:,nn] = np.percentile(surr_quan[:,nn], (perc_low,perc_high))
    else:
        pul = np.percentile(surr_quan, (perc_low,perc_high))
    if target is None:
        return pul[0], pul[1]
    vp = value_percentile(surr_quan, target)
    return pul[0], pul[1], 100.0 - vp


def value_percentile(sequence, target):
    """Find the percentile of a particular value in a sequence

    Parameters
    ==========
    sequence : sequence
        a sequence of values, sorted in ascending order
    target : same type as sequence
         a target value

    Returns
    =======
    out : float
        the percentile of target in sequence
    """
    min = bisect.bisect_left(sequence, target) #first value >= x
    max = bisect.bisect_right(sequence, target, lo=min) #first value > x
    #min:max, in Python syntax, will include all values that are target
    count = len(sequence)
    if max == 0: #everything is bigger
        return 0.0
    if min == count: #everything is smaller
        return 100.0
    #Index 0 is 0th percentile (by definition); count-1 is 100th
    #so index n is (n * 100.0) / (count - 1)
    #find a 'virtual' index
    if min == max: #Not equal to any, between min-1 and min
        #interpolate between two points based on value.
        idx = min - 1 + float(target - sequence[min - 1]) / \
              (sequence[min] - sequence[min - 1])
    else: #Equal to one or more values (min...max-1), average them
        idx = (min + max - 1) / 2
    return idx * 100.0 / (count - 1)


def applyRefractory(process1, period):
    '''Apply a refractory period to an input discrete event time sequence

    All events in the refractory period are removed from the point process.

    Parameters
    ==========
    process1 : iterable
        an iterable of datetimes, or a spacepy.time.Ticktock
    period : datetime.timedelta
        length of refractory period

    Returns
    =======
    keep : iterable
        returns pruned set of datetimes with same type as input
        NOTE: array subclasses will be lost
    '''
    import spacepy.time as spt
    if isinstance(process1, spt.Ticktock):
        #SpacePy Ticktock
        p1 = dm.dmcopy(process1.UTC).tolist()
        tickt = True
    elif isinstance(process1[0], dt.datetime):
        #iterable containing datetimes
        p1 = dm.dmcopy(process1)
        try:
            p1 = p1.tolist()
            wasArr = True
        except:
            wasArr = False
        tickt = False
    else:
        raise NotImplementedError('Input process must be a list/array of datetimes, or a spacepy.time.Ticktock')

    try:
        assert period.seconds
    except AssertionError:
        raise AttributeError('period must be a datetime.timedelta')

    done = len(p1)<2
    keep, discard = [], []
    while not done:
        t1 = p1[0]
        t2 = t1 + period
        inds = tb.tOverlapHalf([t1, t2], p1[1:])
        for idx in inds:
            discard.append(p1.pop(idx + 1))
        keep.append(p1.pop(0)) # put test element into keep array
        done = len(p1) < 2

    if tickt:
        return spt.Ticktock(keep)
    return np.array(keep) if wasArr else keep
