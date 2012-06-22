"""
Utility routines for plotting and related activities

Authors: Jonathan Niehof

Institution: Los Alamos National Laboratory

Contact: jniehof@lanl.gov

Copyright 2012 Los Alamos National Security, LLC.

.. currentmodule:: spacepy.plot.utils

Classes
-------

.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    EventClicker

Functions
---------

.. autosummary::
    :toctree: autosummary

    applySmartTimeTicks
    printfig
    smartTimeTicks
    timestamp
"""

__contact__ = 'Jonathan Niehof: jniehof@lanl.gov'

import bisect
import datetime

import matplotlib.dates
import matplotlib.pyplot as plt
import numpy

class EventClicker(object):
    """
    Presents a provided figure (normally a time series) and provides
    an interface to mark events shown in the plot. The user interface
    is explained in :meth:`analyze` and results are returned
    by :meth:`get_events`

    Other Parameters
    ================
    ax : maplotlib.axes.AxesSubplot
        The subplot to display and grab data from. If not provided, the
        current subplot is grabbed from gca() (Lookup of the current
        subplot is done when :meth:`analyze` is called.)
        
    n_phases : int (optional, default 1)
        number of phases to an event, i.e. number of subevents to mark.
        E.g. for a storm where one wants the onset and the minimum, set
        n_phases to 2 and double click on the onset, then minimum, and
        then the next double-click will be onset of the next storm.

    interval : (optional)
        Size of the X window to show. This should be in units that can
        be added to/subtracted from individual elements of x (e.g.
        timedelta if x is a series of datetime.) Defaults to showing
        the entire plot.

    auto_interval : boolean (optional)
        Automatically adjust interval based on the average distance
        between selected events. Default is True if interval is not
        specified; False if interval is specified.

    auto_scale : boolean (optional, default True):
        Automatically adjust the Y axis to match the data as the X
        axis is panned.

    ymin : (optional, default None)
        If auto_scale is True, the bottom of the autoscaled Y axis will
        never be above ymin (i.e. ymin will always be shown on the plot).
        This prevents the autoscaling from blowing up very small features
        in mostly flat portions of the plot. The user can still manually
        zoom in past this point. The autoscaler will always zoom out to
        show the data.

    ymax : (optional, default None)
        Similar to ymin, but the top of the Y axis will never be below ymax.

    line : matplotlib.lines.Line2D (optional)
        Specify the matplotlib line object to use for autoscaling the
        Y axis. If this is not specified, the first line object on the
        provided subplot will be used. This should usually be correct.

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> import numpy
    >>> import matplotlib.pyplot as plt
    >>> x = numpy.arange(630) / 100.0 * numpy.pi
    >>> y = numpy.sin(x)
    >>> clicker = spacepy.plot.utils.EventClicker(
    ... n_phases=2, #Two picks per event
    ... interval=numpy.pi * 2) #Display one cycle at a time
    >>> plt.plot(x, y)
    >>> clicker.analyze() #Double-click on max and min of each cycle; close
    >>> e = clicker.get_events()
    >>> peaks = e[:, 0, 0] #x value of event starts
    >>> peaks -= 2 * numpy.pi * numpy.floor(peaks / (2 * numpy.pi)) #mod 2pi
    >>> max(numpy.abs(peaks - numpy.pi / 2)) < 0.2 #Peaks should be near pi/2
    True
    >>> troughs = e[:, 1, 0] #x value of event ends
    >>> troughs -= 2 * numpy.pi * numpy.floor(troughs / (2 * numpy.pi))
    >>> max(numpy.abs(peaks - 3 * numpy.pi / 2)) < 0.2 #troughs near 3pi/2
    True
    >>> d = clicker.get_events_data() #snap-to-data of events
    >>> peakvals = d[:, 0, 1] #y value, snapped near peaks
    >>> max(peakvals) <= 1.0 #should peak at 1
    True
    >>> min(peakvals) > 0.9 #should click near 1
    True
    >>> troughvals = d[:, 1, 1] #y value, snapped near peaks
    >>> max(troughvals) <= -0.9 #should click near -1
    True
    >>> min(troughvals) <= -1.0 #should bottom-out at -1
    True

    .. autosummary::
         ~EventClicker.analyze
         ~EventClicker.get_events
         ~EventClicker.get_events_data

    .. codeauthor:: Jon Niehof <jniehof@lanl.gov>
    .. automethod:: analyze
    .. automethod:: get_events
    .. automethod:: get_events_data
    """
    _colors = ['k', 'r', 'g']
    _styles = ['solid', 'dashed', 'dotted']

    def __init__(self, ax=None, n_phases=1, interval=None, auto_interval=None,
                 auto_scale=True, ymin=None, ymax=None, line=None):
        """Initialize EventClicker

        Other Parameters
        ================
        ax : maplotlib.axes.AxesSubplot
            The subplot to display and grab data from. If not provided, the
            current subplot is grabbed from gca() (Lookup of the current
            subplot is done when :meth:`analyze` is called.)
        
        n_phases : int (optional, default 1)
            number of phases to an event, i.e. number of subevents to mark.
            E.g. for a storm where one wants the onset and the minimum, set
            n_phases to 2 and double click on the onset, then minimum, and
            then the next double-click will be onset of the next storm.

        interval : (optional)
            Size of the X window to show. This should be in units that can
            be added to/subtracted from individual elements of x (e.g.
            timedelta if x is a series of datetime.) Defaults to showing
            the entire plot.

        auto_interval : boolean (optional)
            Automatically adjust interval based on the average distance
            between selected events. Default is True if interval is not
            specified; False if interval is specified.

        auto_scale : boolean (optional, default True):
            Automatically adjust the Y axis to match the data as the X
            axis is panned.

        ymin : (optional, default None)
            If auto_scale is True, the bottom of the autoscaled Y axis will
            never be above ymin (i.e. ymin will always be shown on the plot).
            This prevents the autoscaling from blowing up very small features
            in mostly flat portions of the plot. The user can still manually
            zoom in past this point. The autoscaler will always zoom out to
            show the data.

        ymax : (optional, default None)
            Similar to ymin, but the top of the Y axis will never be below ymax.

        line : matplotlib.lines.Line2D (optional)
            Specify the matplotlib line object to use for autoscaling the
            Y axis. If this is not specified, the first line object on the
            provided subplot will be used. This should usually be correct.
        """
        self.n_phases = n_phases
        self.interval = interval
        self._autointerval = auto_interval
        self._autoscale = auto_scale
        self._intervalcount = 0
        self._intervaltotal = None
        self._events = None
        self._data_events = None #snap-to-data version of events
        self._ymax = ymax
        self._ymin = ymin
        self._line = line
        self.ax = None

    def analyze(self):
        """
        Displays the figure provided and allows the user to select events.

        All matplot lib controls for zooming, panning, etc. the figure
        remain active.

        Double left click
            Mark this point as an event phase. One-phase events are the
            simplest: they occur at a particular time. Two-phase events
            have two times associated with them; an example is any event
            with a distinct start and stop time. In that case, the first
            double-click would mark the beginning, the second one, the end;
            the next double-click would mark the beginning of the next event.
            Each phase of an event is annotated with a vertical line on the
            plot; the color and line style is the same for all events, but
            different for each phase.

            After marking the final phase of an event, the X axis will scroll
            and zoom to place that phase near the left of the screeen and
            include one full interval of data (as defined in the constructor).
            The Y axis will be scaled to cover the data in that X range.
            
        Double right click or delete button
            Remove the last marked event phase. If an entire event (i.e., the
            first phase of an event) is removed, the X axis will be scrolled
            left to the previous event and the Y axis will be scaled to cover
            the data in the new range.
            
        Space bar
            Scroll the X axis by one interval. Y axis will be scaled to cover
            the data.

        When finished, close the figure window (if necessary) and call
        :meth:`get_events` to get the list of events.
        """
        self._lastclick_x = None
        self._lastclick_y = None
        self._lastclick_button = None
        self._curr_phase = 0
        
        if self.ax is None:
            self.ax = plt.gca()
        self.fig = self.ax.get_figure()
        lines = self.ax.get_lines()
        if self._line is None:
            if len(lines) > 0:
                self._line = lines[0]
        else:
            if not self._line in lines:
                self._line = None
                
        if self._line is None:
            self._xdata = None
            self._ydata = None
            self._autoscale = False
            self._x_is_datetime = False
        else:
            self._xdata = self._line.get_xdata()
            self._ydata = self._line.get_ydata()
            self._x_is_datetime = isinstance(self._xdata[0],
                                             datetime.datetime)
            if self._x_is_datetime:
                self._xydata = numpy.column_stack(
                    (matplotlib.dates.date2num(self._xdata), self._ydata))
            else:
                self._xydata = numpy.column_stack((self._xdata, self._ydata))
            if self._ymin is None: #Make the clipping comparison always fail
                self._ymin = numpy.nanmax(self._ydata)
            if self._ymax is None:
                self._ymax = numpy.nanmin(self._ydata)

        if self._autointerval is None:
            self._autointerval = self.interval is None
        if self.interval is None:
            (left, right) = self.ax.get_xaxis().get_view_interval()
            if self._x_is_datetime:
                right = matplotlib.dates.num2date(right)
                left = matplotlib.dates.num2date(left)
            self.interval = right - left

        if not self._xdata is None:
            self._relim(self._xdata[0])
        else:
            self._relim(self.ax.get_xaxis().get_view_interval()[0])
        self._cids = []
        self._cids.append(self.fig.canvas.mpl_connect('button_press_event', self._onclick))
        self._cids.append(self.fig.canvas.mpl_connect('close_event', self._onclose))
        self._cids.append(self.fig.canvas.mpl_connect('key_press_event', self._onkeypress))
        plt.show()

    def get_events(self):
        """Get back the list of events.

        Call after :meth:`analyze`.

        Returns
        =======
        out : array
            3-D array of (x, y) values clicked on.
            Shape is (n_events, n_phases, 2), i.e. indexed by event
            number, then phase of the event, then (x, y).
        """
        if self._events is None:
            return None
        else:
            return self._events[0:-1].copy()

    def get_events_data(self):
        """Get a list of events, "snapped" to the data.

        For each point selected as a phase of an event, selects the point
        from the original data which is closest to the clicked point. Distance
        from point to data is calculated based on the screen distance, not
        in data coordinates.

        Note that this snaps to data points, not to the closest point on the
        line between points.

        Call after :meth:`analyze`.

        Returns
        =======        
        out : array
            3-D array of (x, y) values in the data which are closest to each
            point clicked on. Shape is (n_events, n_phases, 2), i.e. indexed
            by event number, then phase of the event, then (x, y).
        """
        if self._data_events is None:
            return None
        else:
            return self._data_events[0:-1].copy()

    def _add_event_phase(self, xval, yval):
        """Add a phase of the event"""
        self.ax.axvline(
            xval,
            color=self._colors[self._curr_phase % len(self._colors)],
            ls=self._styles[self._curr_phase / len(self._colors) % len(self._styles)])
        if not self._xydata is None:
            point_disp = self.ax.transData.transform(
                numpy.array([[xval, yval]])
                )[0]
            data_disp = self.ax.transData.transform(self._xydata)
            idx = numpy.argmin(numpy.sum(
                (data_disp - point_disp) ** 2, axis=1
                ))
            if self._data_events is None:
                self._data_events = numpy.array(
                    [[[self._xdata[0], self._ydata[0]]] * self.n_phases])
            self._data_events[-1, self._curr_phase] = \
                                  [self._xdata[idx], self._ydata[idx]]
        if self._x_is_datetime:
            xval = matplotlib.dates.num2date(xval)
        if self._events is None:
            self._events = numpy.array([[[xval, yval]] * self.n_phases])
        self._events[-1, self._curr_phase] = [xval, yval]
        self._curr_phase += 1
        if self._curr_phase >= self.n_phases:
            self._curr_phase = 0
            if self._autointerval:
                if self._events.shape[0] > 2:
                    self._intervalcount += 1
                    self._intervaltotal += (self._events[-1, 0, 0] - self._events[-2, 0, 0])
                    self.interval = self._intervaltotal / self._intervalcount
                elif self._events.shape[0] == 2:
                    self._intervalcount = 1
                    self._intervaltotal = self._events[1, 0, 0] - self._events[0, 0, 0]
                    self.interval = self._intervaltotal
                
            self._events.resize((self._events.shape[0] + 1,
                                 self.n_phases, 2
                                 ))
            self._data_events.resize((self._data_events.shape[0] + 1,
                                      self.n_phases, 2
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
                                     self.n_phases, 2
                                     ))
                if not self._data_events is None:
                    self._data_events.resize((self._data_events.shape[0] - 1,
                                              self.n_phases, 2
                                              ))
                self._curr_phase = self.n_phases - 1
        else:
            del self.ax.lines[-1]
            self._curr_phase -= 1
        if self._curr_phase == 0 and self._events.shape[0] > 1:
            self._relim(self._events[-2, -1, 0])
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
        for cid in self._cids:
            self.fig.canvas.mpl_disconnect(cid)

    def _onkeypress(self, event):
        """Handle a keypress"""
        if event.key == ' ':
            rightside = self.ax.xaxis.get_view_interval()[1]
            if self._x_is_datetime:
                rightside = matplotlib.dates.num2date(rightside)
            self._relim(rightside)
        if event.key == 'delete':
            self._delete_event_phase()
        
    def _relim(self, left_x):
        """Reset the limits based on a particular X value"""
        if self._x_is_datetime:
            xmin = left_x - self.interval/10
            xmax = left_x + self.interval + self.interval/10
        else:
            xmin = left_x - 0.1 * self.interval
            xmax = left_x + 1.1 * self.interval
        if self._autoscale:
            idx_l = bisect.bisect_left(self._xdata, xmin)
            idx_r = bisect.bisect_right(self._xdata, xmax)
            if idx_l >= len(self._ydata):
                idx_l = len(self._ydata) - 1
            ymin = numpy.nanmin(self._ydata[idx_l:idx_r])
            ymax = numpy.nanmax(self._ydata[idx_l:idx_r])
            if ymin > self._ymin:
                ymin = self._ymin
            if ymax < self._ymax:
                ymax = self._ymax
            ydiff = (ymax - ymin) / 10
            ymin -= ydiff
            ymax += ydiff
            self.ax.set_xlim(xmin, xmax)
            self.ax.set_ylim(ymin, ymax)
            self.ax.autoscale_view()
        else:
            self.ax.set_xlim(xmin, xmax)
            self.ax.autoscale_view(scalex=True, scaley=False)
        self.fig.canvas.draw()

def applySmartTimeTicks(ax, time, dolimit = True):
    """
    Given an axis 'ax' and a list/array of datetime objects, 'time',
    use the smartTimeTicks function to build smart time ticks and
    then immediately apply them to the given axis.  The first and
    last elements of the time list will be used as bounds for the
    x-axis range.

    The range of the 'time' input value will be used to set the limits
    of the x-axis as well.  Set kwarg 'dolimit' to False to override
    this behavior.

    Parameters
    ==========
    ax : matplotlib.pyplot.Axes
        A matplotlib Axis object.
    time : list
        list of datetime objects
    dolimit : boolean (optional)
        The range of the 'time' input value will be used to set the limits
        of the x-axis as well. Setting this overrides this behavior.

    See Also
    ========
    smartTimeTicks
    """
    Mtick, mtick, fmt = smartTimeTicks(time)
    ax.xaxis.set_major_locator(Mtick)
    ax.xaxis.set_minor_locator(mtick)
    ax.xaxis.set_major_formatter(fmt)
    if dolimit:
        ax.set_xlim([time[0], time[-1]])

def smartTimeTicks(time):
    """
    Returns major ticks, minor ticks and format for time-based plots

    smartTimeTicks takes a list of datetime objects and uses the range
    to calculate the best tick spacing and format.  Returned to the user
    is a tuple containing the major tick locator, minor tick locator, and
    a format string -- all necessary to apply the ticks to an axis.

    It is suggested that, unless the user explicitly needs this info,
    to use the convenience function applySmartTimeTicks to place the
    ticks directly on a given axis.

    Parameters
    ==========
    time : list
        list of datetime objects

    Returns
    =======
    out : tuple
        tuple of Mtick - major ticks, mtick - minor ticks, fmt - format

    See Also
    ========
    applySmartTimeTicks
    """
    from matplotlib.dates import (MinuteLocator, HourLocator,
                                  DayLocator, DateFormatter)
    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if nHours < 1:
        Mtick = MinuteLocator(byminute = [0,15,30,45])
        mtick = MinuteLocator(byminute = list(range(60)), interval = 5)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick = MinuteLocator(byminute = [0,30])
        mtick = MinuteLocator(byminute = list(range(60)), interval = 10)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick = HourLocator(byhour = list(range(24)), interval = 2)
        mtick = MinuteLocator(byminute = [0,15,30,45])
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 24:
        Mtick = HourLocator(byhour = [0,3,6,9,12,15,18,21])
        mtick = HourLocator(byhour = list(range(24)))
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick = HourLocator(byhour = [0,6,12,18])
        mtick = HourLocator(byhour = list(range(24)))
        fmt = DateFormatter('%H:%M UT')
    else:
        Mtick = DayLocator(bymonthday = list(range(1,32)))
        mtick = HourLocator(byhour = [0,6,12,18])
        fmt = DateFormatter('%d %b')

    return (Mtick, mtick, fmt)

def printfig(fignum, saveonly=False, pngonly=False, clean=False, filename=None):
    """save current figure to file and call lpr (print).

    This routine will create a total of 3 files (png, ps and c.png) in the
    current working directory with a sequence number attached. Also, a time
    stamp and the location of the file will be imprinted on the figure. The
    file ending with c.png is clean and no directory or time stamp are
    attached (good for PowerPoint presentations).

    Parameters
    ----------
    fignum : integer
        matplotlib figure number
    saveonly : boolean (optional)
        True (don't print and save only to file)  False (print and save)
    pngolny : boolean (optional)
        True (only save png files and print png directly) False (print ps file, and generate png, ps; can be slow)
    clean : boolean (optional)
        True (print and save only clean files without directory info) False (print and save directory location as well)
    filename : string (optional)
        None (If specified then the filename is set and code does not use the sequence number)

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> import matplotlib.pyplot as plt
    >>> p = plt.plot([1,2,3],[2,3,2])
    >>> spacepy.plot.utils.printfig(1, pngonly=True, saveonly=True)
    """
    import matplotlib.pyplot as plt

    try:
        nfigs = len(fignum)
    except:
        nfigs = 1
        fignum = [fignum]

    for ifig in fignum:
        # active this figure
        plt.figure(ifig)

        if filename == None:
            # create a filename for the figure
            cwd = os.getcwd()
            num = len(glob.glob('*.png'))
            fln = cwd+'/figure_'+str(num)
        else:
            fln = filename
        # truncate fln if too long
        if len(fln) > 60:
            flnstamp = '[...]'+fln[-60:]
        else:
            flnstamp = fln

        # save a clean figure without timestamps
        if clean == True:
            plt.savefig(fln+'_clean.png')
            plt.savefig(fln+'_clena.ps')

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")
        # add the filename to the figure for reference
        plt.figtext(0.01, 0.01, timestamp+flnstamp+'.png', rotation='vertical', va='bottom', size=8)

        # now save the figure to this filename
        if pngonly == False:
            plt.savefig(fln+'.ps')

        plt.savefig(fln+'.png')

        # send it to the printer
        if saveonly != True:
            if pngonly == False:
                os.popen('lpr '+fln+'.ps')
            else:
                os.popen('lpr '+fln+'.png')
    return

def timestamp(position=[1.003, 0.01], size='xx-small', draw=True, **kwargs):
    """
    print a timestamp on the current plot, vertical lower right

    Parameters
    ==========
    position : list
        position for the timestamp
    size : string (optional)
        text size
    draw : Boolean (optional)
        call draw to make sure it appears
    kwargs : keywords
        other keywords to axis.annotate

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> from pylab import plot, arange
    >>> plot(arange(11))
    [<matplotlib.lines.Line2D object at 0x49072b0>]
    >>> spacepy.plot.utils.timestamp()
    """
    from matplotlib.pyplot import gca, draw
    now = datetime.datetime.now()
    strnow = now.strftime("%d%b%Y %H:%M")
    ax=gca()
    ax.annotate(strnow, position, xycoords='axes fraction', rotation='vertical', size=size, va='bottom',  **kwargs)
    if draw:
        draw()
