"""
spacepy.plot.utils

various utilitiy routines for plotting and plot related activities

Authors: Brian Larsen

Institution: Los Alamos National Laboratory

Contact: balarsen@lanl.gov

Copyright 2011 Los Alamos National Security, LLC.
"""

__contact__ = 'Brian Larsen: balarsen@lanl.gov'

from pylab import gca, gcf, show, close

def print_clicks(fig=None, ax=None):
    """
    given a figure print out the values of the clicks
           
    Other Parameters
    ================
    fig : matplotlib.figure.Figure (optional)
        a figure instance to use for the clicking if not specified grabbed from gcf()
    ax : matplotlib.axes.AxesSubplot (optional)
        if an axis is specified then use that axis, otherwise grabbed from gca()
        
    """
    if fig == None:
        fig = gcf()
    if ax == None:
        ax = gca()
    def onclick(event):
        if event.button == 1: # only use left clicks
            print '{}, {}'.format(event.xdata, event.ydata)
            return [event.xdata, event.ydata]
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    show()
