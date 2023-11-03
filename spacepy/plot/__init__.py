#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
plot: SpacePy plotting routines

This package aims to make getting publication ready plots easier.
It provides classes and functions for different types of plot 
(e.g. Spectrogram, levelPlot), for helping make plots more cleanly 
(e.g. set_target, dual_half_circle), and for making plots convey 
information more cleanly, with less effort 
(e.g. applySmartTimeTicks, style).

This plot module now provides style sheets. For most standard plotting 
we recommend the *default* style sheet (aka *spacepy*). To apply the
default plot style, use the following::

    import spacepy.plot as splot
    splot.style()

.. versionchanged:: 0.5.0
   Plot styles are no longer applied automatically.

Different plot types may not work well with this style, so we have provided
alternatives. For polar plots, spectrograms, or anything with larger blocks 
of color, it may be better to use one of the alternatives::

    import spacepy.plot as splot
    splot.style('altgrid') # inverts background from default so it's white
    splot.style('polar') # designed for filled polar plots
    splot.revert_style() # put the style back to matplotlib defaults

Authors: Brian Larsen and Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov


Copyright 2011-2016 Los Alamos National Security, LLC.

.. autosummary::
   :template: clean_function.rst
   :toctree: autosummary

   add_logo
   annotate_xaxis
   applySmartTimeTicks
   available
   collapse_vertical
   dual_half_circle
   levelPlot
   plot
   revert_style
   set_target
   shared_ylabel
   solarRotationPlot
   Spectrogram
   style
   timestamp
   add_arrows

Most of the functionality in the plot module is made available directly 
through the *plot* namespace. However, the plot module does contain
several submodules listed below

.. autosummary::
    :toctree: autosummary
    :template: clean_module.rst

    carrington
    spectrogram
    utils

"""
from collections.abc import Mapping
import os
import warnings
import numpy as np
import matplotlib as mpl
from matplotlib.patches import Wedge
from matplotlib import __version__ as mplv
import matplotlib.pyplot as plt
from spacepy import __path__ as basepath
from spacepy.datamodel import dmcopy
import spacepy.datamanager as dman
from .. import config
from .spectrogram import *
from .utils import *
from .carrington import *

def plot(*args, **kwargs):
    '''Convenience wrapper for matplotlib's plot function

    As with matplotlib's plot function, *args* is a variable length
    argument, allowing for multiple *x*, *y* pairs, each with optional 
    format string. For full details, see matplotlib.pyplot.plot

    Other Parameters
    ----------------
    smartTimeTicks : boolean
        If True then use applySmartTimeTicks to set x-axis labeling
    figsize : array-like, 2 elements
        Set figure size directly on call to plot, (width, height)
    **kwargs : other keywords
        Other keywords to pass to matplotlib.pyplot.plot

    '''
    if 'smartTimeTicks' in kwargs:
        sTT = kwargs['smartTimeTicks']
        del kwargs['smartTimeTicks']
    else:
        sTT = False
    if 'figsize' not in kwargs:
        kwargs['figsize'] = (10,6)
        fig = plt.figure(figsize=kwargs['figsize'])
        del kwargs['figsize']
    pobj = plt.plot(*args, **kwargs)
    if sTT is True:
        ax = plt.gca()
        applySmartTimeTicks(ax, args[0], dolimit=True, dolabel=False)
    return pobj

def available(returnvals=False):
    '''List the available plot styles provided by spacepy.plot

    Note that some of the available styles have multiple aliases.
    To apply an available style, use `spacepy.plot.style`.
    '''
    spacepystyle = os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy.mplstyle')
    spacepyaltstyle = os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_altgrid.mplstyle')
    polarstyle = os.path.join('{0}'.format(basepath[0]), 'data', 'spacepy_polar.mplstyle')
    lookdict = {'default': spacepystyle,
                'spacepy': spacepystyle,
                'spacepy_altgrid': spacepyaltstyle,
                'altgrid': spacepyaltstyle,
                'spacepy_polar': polarstyle,
                'polar': polarstyle
               }
    if returnvals:
        return lookdict
    else:
        return list(lookdict.keys())

def style(look=None, cmap='plasma'):
    '''
    Apply SpacePy's matplotlib style settings from a known style sheet.

    Parameters
    ----------
    look : str, optional
        Name of style. For a list of available style names, see
        `spacepy.plot.available`. If not specified, will use default
        ``"spacepy"`` style.
    '''
    lookdict = available(returnvals=True)
    usestyle = lookdict.get(look, lookdict['default'])
    plt.style.use(usestyle)
    mpl.rcParams['image.cmap'] = cmap

#save current rcParams before applying spacepy style
oldParams = {key: dmcopy(val) for key, val in mpl.rcParams.items()}

def revert_style():
    '''Revert plot style settings to those in use prior to importing spacepy.plot
    '''
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for key in oldParams:
            try:
                mpl.rcParams[key] = oldParams[key]
            except ValueError:
                pass

def dual_half_circle(center=(0,0), radius=1.0,
                     sun_direction='right', ax=None, colors=('w','k'),
                     **kwargs):
    """
    Plot two half circles to a plot with the specified face colors and
    rotation. This is normal to use to denote the sun direction in
    magnetospheric science plots.

    Other Parameters
    ----------------
    center : array-like, 2 elements
        Center in data coordinates of the circles, default (0,0)
    radius : float
        Radius of the circles, defualt 1.0
    sun_direction : string or float
        The rotation direction of the first (white) circle. Options are ['down',
        'down right', 'right', 'up left', 'up right', 'up', 'down left', 'left']
        or an angle in degrees counter-clockwise from up. Default right.
    ax : matplotlib.axes
        Axis to plot the circles on.
    colors : array-like, 2 elements
        The two colors for the circle fill. The First number is the light and
        second is the dark.
    **kwargs : other keywords
        Other keywords to pass to matplotlib.patches.Wedge

    Returns
    -------
    out : tuple
        Tuple of the two wedge objects

    Examples
    --------
    >>> import spacepy.plot
    >>> spacepy.plot.dual_half_circle()

    """
    sun_dict = {'left':90,
                'right':-90,
                'up':0,
                'down':180,
                'up right':-45,
                'up left':45,
                'down right':45+180,
                'down left':-45+180,
                }
    try:
        angle = sun_dict[sun_direction]
    except KeyError:
        try:
            angle = float(sun_direction)
        except ValueError:
            raise(ValueError("Sun_direction was not understood, must be a float or {0}".format(sun_dict.keys())))
                
    theta1, theta2 = angle, angle + 180
    if ax is None:
        ax = plt.gca()
        
    w1 = Wedge(center, radius, theta1, theta2, fc=colors[0],transform=ax.transData._b, **kwargs)
    w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], transform=ax.transData._b,**kwargs)
    for wedge in [w1, w2]:
        ax.add_artist(wedge)
    return (w1, w2)

def levelPlot(data, var=None, time=None, levels=(3, 5), target=None, colors=None, **kwargs):
    """
    Draw a step-plot with up to 5 levels following a color cycle (e.g. Kp index "stoplight")

    Parameters
    ----------
    data : array-like, or dict-like
        Data for plotting. If dict-like, the key providing an array-like 
        to plot must be given to var keyword argument.

    Other Parameters
    ----------------
    var    : string
        Name of key in dict-like input that contains data
    time   : array-like or string
        Name of key in dict-like that contains time, or arraylike of datetimes
    levels : array-like, up to 5 levels
        Breaks between levels in data that should be shown as distinct colors
    target : figure or axes
        Target axes or figure window
    colors : array-like
        Colors to use for the color sequence (if insufficient colors, will use as a cycle)
    **kwargs : other keywords
        Other keywords to pass to spacepy.toolbox.binHisto

    Returns
    -------
    binned : tuple
        Tuple of the binned data and bins

    Examples
    --------
    >>> import spacepy.plot as splot
    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> tt = spt.tickrange('2012/09/28','2012/10/2', 3/24.)
    >>> omni = om.get_omni(tt)
    >>> splot.levelPlot(omni, var='Kp', time='UTC', colors=['seagreen', 'orange', 'crimson'])
    """
    #assume dict-like/key-access, before moving to array-like
    if var is not None:
        try:
            usearr = data[var]
        except KeyError:
            raise KeyError('Key "{1}" not present in data'.format(var))
    else:
        #var is None, so make sure we don't have a dict-like
        if not isinstance(data, Mapping):
            usearr = np.asarray(data)
        else:
            raise TypeError('Data appears to be dict-like without a key being given')
    tflag = False
    if time is not None:
        from scipy.stats import mode
        try:
            times = data[time]
        except (KeyError, ValueError, IndexError):
            times = time
        try:
            times = matplotlib.dates.date2num(times)
            tflag = True
        except AttributeError:
            #the x-data are a non-datetime
            times = np.asarray(time)
        #now add the end-point
        stepsize, dum = mode(np.diff(times), axis=None)
        times = np.hstack([times, times[-1]+stepsize])
    else:
        times = np.asarray(range(0, len(usearr)+1))
    if not colors:
        if len(levels)<=3:
            #traffic light colours that are distinct to protanopes and deuteranopes
            colors = ['lime', 'yellow', 'crimson', 'saddlebrown']
        else:
            colors = matplotlib.rcParams['axes.color_cycle']
    else:
        try:
            assert len(colors) > len(levels)
        except AssertionError:
            #cycle the given colors, if not enough are given
            colors = list(colors)*int(1+len(levels)/len(colors))
    if 'alpha' not in kwargs:
        kwargs['alpha']=0.75
    if 'legend' not in kwargs:
        legend = False
    else:
        legend = kwargs['legend']
        del kwargs['legend']
    fig, ax = set_target(target)
    subset = np.asarray(dmcopy(usearr))

    def fill_between_steps(ax, x, y1, **kwargs):
        y2 = np.zeros_like(y1)
        stepsxx = x.repeat(2)[1:-1]
        stepsyy = y1.repeat(2)
        y2 = np.zeros_like(stepsyy)
        ax.fill_between(stepsxx, stepsyy, y2, **kwargs)
        if mpl.__version__<'1.5.0':
            #pre-v1.5.0, need to manually add an artist for the legend
            p = plt.Rectangle((0, 0), 0, 0, **kwargs)
            ax.add_patch(p)
    
    #below threshold 1
    idx = 0
    inds = usearr>levels[0]
    subset[inds] = np.nan
    kwargs['label'] = u'≤{0}'.format(levels[idx])
    fill_between_steps(ax, times, subset, color=colors[0], zorder=30, **kwargs)
    #for each of the "between" thresholds
    for idx in range(1,len(levels)):
        subset = np.asarray(dmcopy(usearr))
        inds = np.bitwise_or(usearr<=levels[idx-1], usearr>levels[idx])
        subset[inds] = np.nan
        kwargs['label'] = u'>{0},≤{1}'.format(levels[idx-1], levels[idx])
        fill_between_steps(ax, times, subset, color=colors[idx], zorder=30-(idx*2), **kwargs)
    #last
    idx += 1
    try:
        inds = usearr<=levels[idx-1]
        subset = np.asarray(dmcopy(usearr))
        subset[inds] = np.nan
        kwargs['label'] = '>{0}'.format(levels[-1])
        fill_between_steps(ax, times, subset, color=colors[idx], zorder=30-(idx*2), **kwargs)
    except:
        pass

    #if required, set x axis to times
    if tflag:
        try:
            applySmartTimeTicks(ax, data[time])
        except (IndexError, KeyError):
            #using data array to index, so should just use time
            applySmartTimeTicks(ax, time)
        ax.grid(False, which='minor') #minor grid usually looks bad on these...

    if legend:
        ncols = len(levels)+1
        if ncols > 3: ncols = ncols//2
        ax.legend(loc='upper left', ncol=ncols)

    return ax
