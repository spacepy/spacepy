"""
spacepy.plot.utils

various utility routines for plotting and plot related activities

Authors: Brian Larsen

Institution: Los Alamos National Laboratory

Contact: balarsen@lanl.gov

Copyright 2011 Los Alamos National Security, LLC.
"""

__contact__ = 'Brian Larsen: balarsen@lanl.gov'

import bisect
import datetime

import matplotlib.pyplot as plt
import matplotlib.dates
import numpy
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


#TODO: make this accept a pre-existing plot
#TODO: do something with the Y values, maybe snap-to-data
class EventClicker(object):
    """
    .. codeauthor:: Jon Niehof <jniehof@lanl.gov>

    Present a time series for 

    Returns
    =======

    out : array
        2-D array of x values clicked on. First dimension is sized n_phases.
    """
    _colors = ['k', 'r', 'g']
    _styles = ['solid', 'dashed', 'dotted']

    def __init__(self, x, y, n_phases=1, interval=None):
        """
        Create a clicker 
        
        Parameters
        ==========
        x : sequence
            a sequence (e.g. list, 1-D numpy array) of the X values to plot
        y : sequence
            Y values to plot

        Other Parameters
        ================
        n_phases : int (optional, default 1)
            number of phases to an event, i.e. number of subevents to mark.
            E.g. for a storm where one wants the onset and the minimum, set
            n_phases to 2 and double click on the onset, then minimum, and
            then the next double-click will be onset of the next storm.

        interval : same as elements of x (optional)
            Size of the X window to show. This should be in units that can
            be added to/subtracted from individual elements of x (e.g.
            timedelta is x is a series of datetime.) Defaults to showing
            1/20th of the total
        """
        self._x_is_datetime = isinstance(x[0], datetime.datetime)
        self.x = x
        self.y = y
        self.n_phases = n_phases
        if interval == None:
            if self._x_is_datetime:
                self.interval = (self.x[-1] - self.x[0]) / 20
            else:
                self.interval = (self.x[-1] - self.x[0]) / 20.0
            self._autointerval = True
        else:
            self.interval = interval
            self._autointerval = False
            self._intervalcount = 0
            self._intervaltotal = None

    def analyze(self):
        """Iterate over the elements"""
        self._lastclick_x = None
        self._lastclick_y = None
        self._lastclick_button = None
        self._curr_phase = 0
        self._events = numpy.array([[self.x[0]] * self.n_phases])
        
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.plot(self.x, self.y)
        self._relim(self.x[0])
        self._cid1 = self.fig.canvas.mpl_connect('button_press_event', self._onclick)
        self._cid2 = self.fig.canvas.mpl_connect('close_event', self._onclose)
        plt.show()

    def get_eventlist(self):
        """Get back the list of events"""
        return self._events[0:-1]

    def _add_event_phase(self, xval, yval):
        """Add a phase of the event"""
        if self._x_is_datetime:
            xval = matplotlib.dates.num2date(xval).replace(tzinfo=None)
        self.ax.axvline(
            xval,
            color=self._colors[self._curr_phase % len(self._colors)],
            ls=self._styles[self._curr_phase / len(self._colors) % len(self._styles)])
        self._events[-1, self._curr_phase] = xval
        self._curr_phase += 1
        if self._curr_phase >= self.n_phases:
            self._curr_phase = 0
            if self._autointerval:
                if self._events.shape[0] > 2:
                    self._intervalcount += 1
                    self._intervaltotal += (self._events[-1, 0] - self._events[-2, 0])
                    self.interval = self._intervaltotal / self._intervalcount
                elif self._events.shape[0] == 2:
                    self._intervalcount = 1
                    self._intervaltotal = self._events[1, 0] - self._events[0, 0]
                    self.interval = self._intervaltotal
                
            self._events.resize((self._events.shape[0] + 1,
                                 self.n_phases
                                 ))
            self._relim(xval)
        else:
            self.fig.canvas.draw()

    def _delete_event_phase(self):
        """Delete the most recent phase of the event"""
        if self._curr_phase == 0:
            if self._events.shape[0] > 1:
                del self.ax.lines[-1]
                self._events.resize((self._events.shape[0] - 1,
                                     self.n_phases
                                     ))
                self._curr_phase = self.n_phases - 1
        else:
            del self.ax.lines[-1]
            self._curr_phase -= 1
        self.fig.canvas.draw()

    def _onclick(self, event):
        """Handle a click"""
        # a doubleclick gives us two IDENTICAL click events, same X and Y
        if event.xdata == self._lastclick_x and \
               event.ydata == self._lastclick_y and \
               event.button == self._lastclick_button:
            if event.button == 1:
                self._add_event_phase(event.xdata, event.ydata)
            else:
                self._delete_event_phase()
            self._lastclick_x = None
            self._lastclick_y = None
            self._lastclick_button = None
        else:
            self._lastclick_x = event.xdata
            self._lastclick_y = event.ydata
            self._lastclick_button = event.button

    def _onclose(self, event):
        """Handle the window closing"""
        self.fig.canvas.mpl_disconnect(self._cid1)
        self.fig.canvas.mpl_disconnect(self._cid2)
        
    def _relim(self, left_x):
        """Reset the limits based on a particular X value"""
        if self._x_is_datetime:
            xmin = left_x - self.interval/10
            xmax = left_x + self.interval + self.interval/10
        else:
            xmin = left_x - 0.1 * self.interval
            xmax = left_x + 1.1 * self.interval
        idx_l = bisect.bisect_left(self.x, xmin)
        idx_r = bisect.bisect_right(self.x, xmax)
        if idx_l >= len(self.y):
            idx_l = len(self.y) - 1
        ymin = min(self.y[idx_l:idx_r])
        ymax = max(self.y[idx_l:idx_r])
        self.ax.set_xlim(xmin, xmax)
        self.ax.set_ylim(ymin, ymax)
        self.ax.autoscale_view()
        self.fig.canvas.draw()
