"""
Utility routines for plotting and related activities

Authors: Jonathan Niehof, Steven Morley, Daniel Welling

Institution: Los Alamos National Laboratory

Contact: jniehof@lanl.gov

Copyright 2012-2014 Los Alamos National Security, LLC.

.. currentmodule:: spacepy.plot.utils

Classes
-------

.. autosummary::
    :template: clean_class.rst
    :toctree:

    EventClicker

Functions
---------

.. autosummary::
    :toctree:

    add_logo
    annotate_xaxis
    applySmartTimeTicks
    collapse_vertical
    printfig
    set_target
    shared_ylabel
    show_used
    smartTimeTicks
    timestamp
    add_arrows
"""

__contact__ = 'Jonathan Niehof: jniehof@lanl.gov'

import bisect
import datetime
import itertools

import matplotlib
import matplotlib.axis
import matplotlib.dates
import matplotlib.image
import matplotlib.patches
try:
    import matplotlib.pyplot as plt
except RuntimeError:
    pass
import numpy

__all__ = ['add_logo', 'annotate_xaxis', 'applySmartTimeTicks', 'collapse_vertical', 'filter_boxes',
           'smartTimeTicks', 'get_biggest_clear', 'get_clear', 'get_used_boxes', 'EventClicker',
           'set_target', 'shared_ylabel', 'show_used', 'timestamp', 'add_arrows']

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
    >>> max(numpy.abs(troughs - 3 * numpy.pi / 2)) < 0.2 #troughs near 3pi/2
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

    >>> import spacepy.plot.utils
    >>> import spacepy.time
    >>> import datetime
    >>> import matplotlib.pyplot as plt
    >>> import numpy
    >>> t = spacepy.time.tickrange('2019-01-01', #get a range of days
    ...                            '2019-12-31',
    ...                            deltadays=datetime.timedelta(days=1))
    >>> y = numpy.linspace(0, 100, 1001)
    >>> seconds = t.TAI - t.TAI[0]
    >>> seconds = numpy.asarray(seconds) #normal ndarray so reshape (in meshgrid) works
    >>> tt, yy = numpy.meshgrid(seconds, y) #use TAI to get seconds
    >>> z = 1 + (numpy.exp(-(yy - 20)**2 / 625) #something like a spectrogram
    ...          * numpy.sin(1e-7 * numpy.pi**2 * tt)**2) #pi*1e7 seconds per year
    >>> plt.pcolormesh(t.UTC, y, z)
    >>> clicker = spacepy.plot.utils.EventClicker(n_phases=1)
    >>> clicker.analyze() #double-click on center of peak; close
    >>> events = clicker.get_events() #returns an array of the things clicked
    >>> len(events) == 10 #10 if you click on the centers, including the last one
    True
    >>> clicker.get_events_data() is None #should be nothing
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
            self._xydata = None
            self._autoscale = False
            self._x_is_datetime = isinstance(self.ax.xaxis.converter,
                                             matplotlib.dates.DateConverter)
            if self._x_is_datetime:
                try:
                    self._tz = self.ax.xaxis.get_major_formatter()._tz
                except AttributeError:
                    self._tz = None
        else:
            self._xdata = self._line.get_xdata()
            self._ydata = self._line.get_ydata()
            self._x_is_datetime = isinstance(self._xdata[0],
                                             datetime.datetime)
            if self._x_is_datetime:
                self._xydata = numpy.column_stack(
                    (matplotlib.dates.date2num(self._xdata), self._ydata))
                self._tz = self._xdata[0].tzinfo
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
                #For naive datetimes, smash explicitly back to naive;
                #otherwise it's just replacing with the same.
                right = matplotlib.dates.num2date(right, tz=self._tz).replace(
                    tzinfo=self._tz)
                left = matplotlib.dates.num2date(left, tz=self._tz).replace(
                    tzinfo=self._tz)
            self.interval = right - left

        if not self._xdata is None:
            self._relim(self._xdata[0])
        elif self._x_is_datetime:
            #Handle the case of no xdata but we are using a datetime, such as
            #spectrum data (that's not a line, hence no xdata):
            self._relim(matplotlib.dates.num2date(self.ax.get_xaxis()\
                .get_view_interval()[0]).replace(tzinfo=self._tz))
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
            ls=self._styles[self._curr_phase // len(self._colors) % len(self._styles)])
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
            xval = matplotlib.dates.num2date(xval, tz=self._tz).replace(
                tzinfo=self._tz)
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
            if self._data_events is not None:
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
                rightside = matplotlib.dates.num2date(rightside, tz=self._tz)\
                    .replace(tzinfo=self._tz)
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


def annotate_xaxis(txt, ax=None):
    """
    Write text in-line and to the right of the x-axis tick labels

    Annotates the x axis of an :class:`~matplotlib.axes.Axes` object with text
    placed in-line with the tick labels and immediately to the right of the
    last label. This is formatted to match the existing tick marks.

    Parameters
    ==========
    txt : str
        The annotation text.

    Other Parameters
    ================
    ax : matplotlib.axes.Axes
        The axes to annotate; if not specified, the
        :func:`~matplotlib.pyplot.gca` function will be used.

    Returns
    =======
    out : matplotlib.text.Text
        The :class:`~matplotlib.text.Text` object for the annotation.

    Notes
    =====
    The annotation is placed *immediately* to the right of the last tick label.
    Generally the first character of ``txt`` should be a space to allow some
    room.

    Calls :func:`~matplotlib.pyplot.draw` to ensure tick marker locations are
    up to date.

    Examples
    ========
    .. plot::
        :include-source:

        >>> import spacepy.plot.utils
        >>> import matplotlib.pyplot as plt
        >>> import datetime
        >>> times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i)
        ...     for i in range(0, 48, 3)]
        >>> plt.plot(times, range(16))
        [<matplotlib.lines.Line2D object at 0x0000000>]
        >>> spacepy.plot.utils.annotate_xaxis(' UT') #mark that times are UT
        <matplotlib.text.Text object at 0x0000000>
    """
    if ax is None:
        ax = plt.gca()
    #For some reason the last one is sometimes null, so search for non-null
    t = next((t for t in ax.get_xticklabels()[::-1] if t.get_text()), None)
    if not t:
        plt.draw() #force a redraw, try again
        t = next((t for t in ax.get_xticklabels()[::-1] if t.get_text()), None)
        if not t:
            return
    transform = t.get_transform()
    pos = transform.inverted().transform(t.get_window_extent())
    left = pos[1, 0] #line up the left of annotation with right of existing
    bottom = pos[0, 1] #for some reason bottom matches better than top
    props = dict((p, getattr(t, 'get_' + p)()) for p in
         ['color', 'family', 'size', 'style', 'variant', 'weight'])
    return ax.text(left, bottom, txt, transform=transform,
                   ha='left', va='bottom', **props)


def applySmartTimeTicks(ax, time, dolimit=True, dolabel=False):
    """
    Given an axis *ax* and a list/array of datetime objects, *time*,
    use the smartTimeTicks function to build smart time ticks and
    then immediately apply them to the given axis.  The first and
    last elements of the time list will be used as bounds for the
    x-axis range.

    The range of the *time* input value will be used to set the limits
    of the x-axis as well.  Set kwarg 'dolimit' to False to override
    this behavior.

    Parameters
    ==========
    ax : matplotlib.pyplot.Axes
        A matplotlib Axis object.
    time : list
        list of datetime objects
    dolimit : boolean (optional)
        The range of the *time* input value will be used to set the limits
        of the x-axis as well. Setting this overrides this behavior.
    dolabel : boolean (optional)
        Sets autolabeling of the time axis with "Time from" time[0]

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
    if dolabel:
        ax.set_xlabel('Time from {0}'.format(time[0].isoformat()))
    return True

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
                                  DayLocator, MonthLocator, YearLocator,
                                  DateFormatter)
    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if deltaT.total_seconds()<600:
        Mtick = MinuteLocator(byminute=list(range(0,60,2)) )
        mtick = MinuteLocator(byminute=list(range(60)), interval=1)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < .5:
        Mtick = MinuteLocator(byminute=list(range(0,60,5)) )
        mtick = MinuteLocator(byminute=list(range(60)), interval=5)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 1:
        Mtick = MinuteLocator(byminute = [0,15,30,45])
        mtick = MinuteLocator(byminute = list(range(60)), interval = 5)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 2:
        Mtick = MinuteLocator(byminute=[0,15,30,45])
        mtick = MinuteLocator(byminute=list(range(60)), interval=5)
        fmt =  DateFormatter('%H:%M UT')
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
    elif deltaT.days < 8:
        Mtick = DayLocator(bymonthday=list(range(32)))
        mtick = HourLocator(byhour=list(range(0,24,2)))
        fmt =  DateFormatter('%d %b')
    elif deltaT.days < 15:
        Mtick = DayLocator(bymonthday=list(range(2,32,2)))
        mtick = HourLocator(byhour=[0,6,12,18])
        fmt =  DateFormatter('%d %b')
    elif deltaT.days < 32:
        Mtick = DayLocator(bymonthday=list(range(5,35,5)))
        mtick = HourLocator(byhour=[0,6,12,18])
        fmt =  DateFormatter('%d %b')
    elif deltaT.days < 60:
        Mtick = MonthLocator()
        mtick = DayLocator(bymonthday=list(range(5,35,5)))
        fmt =  DateFormatter('%d %b')
    elif deltaT.days < 731:
        Mtick = MonthLocator()
        mtick = DayLocator(bymonthday=15)
        fmt =  DateFormatter('%b %Y')
    else:
        Mtick = YearLocator()
        mtick = MonthLocator(bymonth=7)
        fmt =  DateFormatter('%Y')
    return(Mtick, mtick, fmt)


def set_target(target, figsize=None, loc=None, polar=False):
    '''
    Given a *target* on which to plot a figure, determine if that *target*
    is **None** or a matplotlib figure or axes object.  Based on the type
    of *target*, a figure and/or axes will be either located or generated.
    Both the figure and axes objects are returned to the caller for further
    manipulation.  This is used in nearly all *add_plot*-type methods.

    Parameters
    ==========
    target : object
        The object on which plotting will happen.

    Other Parameters
    ================
    figsize : tuple
        A two-item tuple/list giving the dimensions of the figure, in inches.
        Defaults to Matplotlib defaults.
    loc : integer
        The subplot triple that specifies the location of the axes object.
        Defaults to matplotlib default (111).
    polar : bool
        Set the axes object to polar coodinates.  Defaults to **False**.

    Returns
    =======
    fig : object
      A matplotlib figure object on which to plot.

    ax : object
      A matplotlib subplot object on which to plot.

    Examples
    ========
    >>> import matplotlib.pyplot as plt
    >>> from spacepy.pybats import set_target
    >>> fig = plt.figure()
    >>> fig, ax = set_target(target=fig, loc=211)

    '''
    # Is target an axes?  Make no new items.
    if isinstance(target, plt.Axes):
        ax  = target
        fig = ax.figure
    else: # Make a new axis
        # Make a new figure if target isn't one
        fig = target if isinstance(target, plt.Figure) \
              else plt.figure(figsize=figsize)
        if loc is None:
            ax = fig.add_subplot(polar=polar)
            if ax is None: # matplotlib <3.1, no default subplot position
                ax = fig.add_subplot(111, polar=polar)
        else:
            ax = fig.add_subplot(loc, polar=polar)
    return fig, ax


def collapse_vertical(combine, others=(), leave_axis=False):
    """
    Collapse the vertical spacing between two or more subplots.

    Useful for a multi-panel plot where most subplots should have
    space between them but several adjacent ones should not (i.e.,
    appear as a single plot.) This function will remove all the
    vertical space between the subplots listed in ``combine`` and
    redistribute the space between all of the subplots in both
    ``combine`` and ``others`` in proportion to their current size,
    so that the relative size of the subplots does not change.

    Parameters
    ==========
    combine : sequence
        The :class:`~matplotlib.axes.Axes` objects (i.e. subplots)
        which should be placed together with no vertical space.

    Other Parameters
    ================
    others : sequence
        The :class:`~matplotlib.axes.Axes` objects (i.e. subplots)
        which will keep their vertical spacing, but will be expanded
        with the space taken away from between the elements of ``combine``.
    leave_axis : bool
        If set to true, will leave the axis lines and tick marks between
        the collapsed subplots. By default, the axis line ("spine") is
        removed so the two subplots appear as one.

    Notes
    =====
    This function can be fairly fragile and should only be used for fairly
    simple layouts, e.g., a one-column multi-row plot stack.

    This may require some clean-up of the y axis labels, as they are likely
    to overlap.

    Examples
    ========
    .. plot::
        :include-source:

        >>> import spacepy.plot.utils
        >>> import matplotlib.pyplot as plt
        >>> fig = plt.figure()
        >>> #Make three stacked subplots
        >>> ax0 = fig.add_subplot(311)
        >>> ax1 = fig.add_subplot(312)
        >>> ax2 = fig.add_subplot(313)
        >>> ax0.plot([1, 2, 3], [1, 2, 1]) #just make some lines
        [<matplotlib.lines.Line2D object at 0x0000000>]
        >>> ax1.plot([1, 2, 3], [1, 2, 1])
        [<matplotlib.lines.Line2D object at 0x0000000>]
        >>> ax2.plot([1, 2, 3], [1, 2, 1])
        [<matplotlib.lines.Line2D object at 0x0000000>]
        >>> #Collapse space between top two plots, leave bottom one alone
        >>> spacepy.plot.utils.collapse_vertical([ax0, ax1], [ax2])
    """
    combine = tuple(combine)
    others = tuple(others)
    #bounding box for ALL subplots/axes
    boxes = dict(((ax, ax.get_position()) for ax in combine + others))
    #vertical sizes
    sizes = dict(((ax, boxes[ax].ymax - boxes[ax].ymin) for ax in boxes))
    #Fraction of vertical for each?
    vtotal = float(sum(sizes.values()))
    for s in sizes:
        sizes[s] /= vtotal
    #get the ones to combine in top-to-bottom order
    c_sort = sorted(combine, key=(lambda x: boxes[x].ymax), reverse=True)
    #and figure out how much space we're going to take away
    v_additional = float(0)
    for i in range(len(c_sort) - 1):
        v_additional += (boxes[c_sort[i]].ymin - boxes[c_sort[i + 1]].ymax)
    #get EVERYTHING in top-to-bottom order
    all_sort = sorted(combine + others, key=(lambda x: boxes[x].ymax),
                      reverse=True)
    shift = 0.0 #how far UP to move each plot
    for i, ax in enumerate(all_sort):
        bb = boxes[ax]
        pos = [bb.xmin, bb.ymin, bb.xmax - bb.xmin, bb.ymax - bb.ymin]
        pos[3] += (sizes[ax] * v_additional) #expand vertically
        shift -= sizes[ax] * v_additional #everything shifted down by expansion
        pos[1] += shift  #slide the bottom by total shift
        ax.set_position(pos)
        if ax in combine: #This subplot is participating in combination
            #Combining with one below?
            if i < len(all_sort) - 1 and all_sort[i + 1] in combine:
                plt.setp(ax.get_xticklabels(), visible=False) #no labels
                if not leave_axis:
                    #no bottom ticks
                    ax.tick_params(axis='x', which='both', bottom=False)
                    ax.spines['bottom'].set_visible(False) #no axis line
                #ALSO slide everything up by difference between these two
                shift += (bb.ymin - boxes[all_sort[i + 1]].ymax)
            #combining with one above?
            if i > 0 and all_sort[i - 1] in combine:
                if not leave_axis:
                    #no top ticks
                    ax.tick_params(axis='x', which='both', top=False)
                    ax.spines['top'].set_visible(False) #no axis line


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


def shared_ylabel(axes, txt, *args, **kwargs):
    """
    Create a ylabel that spans several subplots

    Useful for a multi-panel plot where several subplots have the
    same units/quantities on the y axis.

    Parameters
    ==========
    axes : list
        The :class:`~matplotlib.axes.Axes` objects (i.e. subplots)
        which should share a single label
    txt : str
        The label to place in the middle of all the `axes` objects.

    Other Parameters
    ================
    Additional arguments and keywords are passed through to
    :meth:`~matplotlib.axes.Axes.set_ylabel`

    Returns
    =======
    out : matplotlib.text.Text
        The :class:`~matplotlib.text.Text` object for the label.

    Notes
    =====
    This function can be fairly fragile and should only be used for fairly
    simple layouts, e.g., a one-column multi-row plot stack.

    The label is associated with the bottommost subplot in ``axes``.

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> #Make three stacked subplots
    >>> ax0 = fig.add_subplot(311)
    >>> ax1 = fig.add_subplot(312)
    >>> ax2 = fig.add_subplot(313)
    >>> ax0.plot([1, 2, 3], [1, 2, 1]) #just make some lines
    [<matplotlib.lines.Line2D object at 0x0000000>]
    >>> ax1.plot([1, 2, 3], [1, 2, 1])
    [<matplotlib.lines.Line2D object at 0x0000000>]
    >>> ax2.plot([1, 2, 3], [1, 2, 1])
    [<matplotlib.lines.Line2D object at 0x0000000>]
    >>> #Create a green label across all three axes
    >>> spacepy.plot.utils.shared_ylabel([ax0, ax1, ax2],
    ... 'this is a very long label that spans all three axes', color='g')
    """
    fig = axes[0].get_figure() #better all be the same!
    #these are in Figure coordinate space
    #transform to display coords for sorting
    boxes = dict(((ax, fig.transFigure.transform(ax.get_position()))
                  for ax in axes))
    #top-to-bottom by upper edge
    top = sorted(axes, key=(lambda x: boxes[x][1, 1]), reverse=True)[0]
    #bottom-to-top by lower edge
    bottom = sorted(axes, key=(lambda x: boxes[x][0, 1]))[0]
    #get the TOP of the TOP subplot in axes coordinates of BOTTOM subplot
    top_in_bottom = bottom.transAxes.inverted().transform( #into bottom coords
        boxes[top]) #into display coords from fig
    bottom_in_bottom = bottom.transAxes.inverted().transform( #into bottom
        boxes[bottom]) #into display coords from fig
    #The mean of bottom-of-bottom and top-of-top, in bottom coords
    middle = (top_in_bottom[1, 1] + bottom_in_bottom[0, 1]) / 2
    bottom.set_ylabel(txt, *args, **kwargs)
    lbl = bottom.get_yaxis().get_label()
    lbl.set_verticalalignment('center')
    lbl.set_y(middle)
    return lbl


def timestamp(position=(1.003, 0.01), size='xx-small', draw=True, strnow=None,
              rotation='vertical', ax=None, **kwargs):
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
    if strnow is None:
            now = datetime.datetime.now()
            strnow = now.strftime("%d%b%Y %H:%M")
    if ax is None:
        ax=plt.gca()
    ann=ax.annotate(strnow, position,
                    xycoords='axes fraction', rotation=rotation,
                    size=size, va='bottom',  **kwargs)
    if draw:
        plt.draw()
    return ann

def _used_boxes_helper(obj, renderer=None):
    """Recursively-called helper function for get_used_boxes. Internal."""
    boxes = []
    if hasattr(obj, 'get_renderer_cache'): #I know how to render myself
        renderer = obj.get_renderer_cache()
    #Axis objects are weird, go for the tick/axis labels directly
    if isinstance(obj, matplotlib.axis.Axis):
        boxes = [tl.get_window_extent() for tl in obj.get_ticklabels()
                 if tl.get_text()]
        if obj.get_label().get_text():
            boxes.append(obj.get_label().get_window_extent())
    #Base size on children, *unless* there are none
    elif hasattr(obj, 'get_children') and obj.get_children():
        for child in obj.get_children():
            res = _used_boxes_helper(child, renderer)
            if res is None: #Child can't find its size, just use own bounds
                boxes = []
                break
            boxes.extend(res)
    if boxes: #found details from children
        return boxes
    #Nothing from children, try own bounds
    try:
        return [obj.get_window_extent()]
    except (TypeError, RuntimeError): #need a renderer
        if not renderer is None:
            return [obj.get_window_extent(renderer)]
        else: #I can't figure out my size!
            return None

def get_used_boxes(fig=None):
    """Go through all elements of a figure and find the "boxes" they occupy,
    in figure coordinates. Mostly helper for add_logo
    """
    plt.draw() #invoke the renderer to figure everything out
    if fig is None:
        fig = plt.gcf()
    #Get rid of double-nesting, and don't include top-level z-order 1
    #(background rectangle) OR anything completely degenerate (point only)
    boxes = [box for child in fig.get_children()
             if (child.get_zorder() == 0 or child.get_zorder() > 1)
             for box in _used_boxes_helper(child)
             if (box.xmin != box.xmax or box.ymin != box.ymax)
             ]
    #Transform to figure
    boxes = [fig.transFigure.inverted().transform(b) for b in boxes]
    return [b for b in boxes if numpy.isfinite(b).all()]


def filter_boxes(boxes):
    """From a list of boxes, exclude those that are completely contained by another"""
    #Filter exact overlap (any box before this one have same bounds?)
    boxes= [b for i, b in enumerate(boxes) if i==0 or not max(
        [(b[0][0] == other[0][0] and b[1][0] == other[1][0] and
          b[0][1] == other[0][1] and b[1][1] == other[1][1])
         for other in boxes[0:i]])]
    #and filter "completely enclosed"
    return [b for b in boxes
            if not max( #Is this contained in ANY other box? If so, drop it.
                [(b[0][0] >= other[0][0] and b[0][1] >= other[0][1] and
                    b[1][0] <= other[1][0] and b[1][1] <= other[1][1])
                for other in boxes if not other is b] #don't compare to self
                )]


def get_clear(boxes, pos='br'):
    """Take a list of boxes which *obstruct* the plot, i.e., don't overplot

    Return a list of boxes which are "clear".

    Mostly a helper for add_logo

    pos is where to look for the clear area:
    br: bottom right
    bl: bottom left
    tl: top left
    tr: top right
    """
    pos = pos.lower()
    assert(pos in ('br', 'bl', 'tl', 'tr'))
    clear = []
    if pos[1] == 'l': #sort obstructing boxes on left edge
        sboxes = sorted(boxes, key=lambda b: b[0][0])
    else: #sort on right edge, descending (work in from right edge)
        sboxes = sorted(boxes, key=lambda b: b[1][0], reverse=True)
    if pos[0] == 't':  #There's a clear space across the top of everything
        top = max([b[1][1] for b in sboxes])
        if top < 1.0: # there is space at the top
            clear.append(numpy.array([[0.0, top], [1.0, 1.0]]))
    else: #clear space across bottom of everything
        bottom = min([b[0][1] for b in sboxes])
        if bottom > 0.0: # there is space at the bottom
            clear.append(numpy.array([[0.0, 0.0], [1.0, bottom]]))
    #default corners
    left = 0.0
    right = 1.0
    bottom = 0.0
    top = 1.0
    #Work in from left or right edge, and avoid all boxes that we've
    #reached so far
    for i, box in enumerate(sboxes):
        if pos[0] == 't':
            #bottom of clear zone is top of every box from here to edge
            bottom = 0.0 if i == 0 else max([b[1][1] for b in sboxes[0:i]])
        else:
            #top of clear zone is bottom of every box from here to edge
            top = 1.0 if i == 0 else min([b[0][1] for b in sboxes[0:i]])
        if pos[1] == 'l':
            right = box[0][0] #right edge of clear zone is the left of this box
        else:
            left = box[1][0] #left of clear zone is right of this obstructing box
        clearbox = numpy.array([[left, bottom], [right, top]])
        clearbox = numpy.clip(clearbox, 0, 1)
        clear.append(clearbox)
    return filter_boxes(clear) #and remove overlaps


def get_biggest_clear(boxes, fig_aspect=1.0, img_aspect=1.0):
    """Given a list of boxes with clear space, figure aspect ratio (width/height),
    and image aspect ratio (width/height), return the largest clear space
    that maintains the aspect ratio of the image

    Mostly a helper for add_logo
    """
    def effective_width(box):
        """Returns "effective" width of the box"""
        width = box[1][0] - box[0][0]
        height = box[1][1] - box[0][1]
        #If figure is wide, each unit of height is smaller than unit of width
        real_height = height / fig_aspect #in width units
        #Box aspect ratio, corrected for figure. Is it "taller" than image?
        if width / real_height <= img_aspect:
            #yes, so the width is the limiter
            return width
        else:
            #no, take the height, correct for figure aspect, and find the
            #width the image would have at this height and its aspect ratio.
            return real_height * img_aspect
    return sorted(boxes, key=effective_width, reverse=True)[0]


def add_logo(img, fig=None, pos='br', margin=0.05):
    """
    Add an image (logo) to one corner of a plot.

    The provided image will be placed in a corner of the plot and sized
    to maintain its aspect ratio and be as large as possible without
    overlapping any existing elements of the figure. Thus this should
    be the last call in constructing a figure.

    Parameters
    ==========
    img : str or numpy.ndarray
        The image to place on the figure. If a string, assumed to be a
        filename to be read with :func:`~matplotlib.image.imread`; if
        a numpy array, assumed to be the image itself
        (in a simliar format).

    Other Parameters
    ================
    fig : matplotlib.figure.Figure
        The figure on which to place the logo; if not specified, the
        :func:`~matplotlib.pyplot.gcf` function will be used.

    pos : str
        The position to place the logo.
        br: bottom right;
        bl: bottom left;
        tl: top left;
        tr: top right

    margin : float
        Margin to include on each side of figure, as a fraction of the larger
        dimension of the figure (width or height). Default is 0.05 (5%).

    Returns
    =======
    (axes, axesimg) : tuple of Axes and AxesImage
        The :class:`~matplotlib.axes.Axes` object created to hold the iamge,
        and the :class:`~matplotlib.image.AxesImage` object for the image
        itself.

    Notes
    =====
    Calls :func:`~matplotlib.pyplot.draw` to ensure locations are
    up to date.

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot([1, 2, 3], [1, 2, 1])
    [<matplotlib.lines.Line2D at 0x00000000>]
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot([1, 2, 3], [2, 1, 2])
    [<matplotlib.lines.Line2D at 0x00000000>]
    >>> spacepy.plot.utils.add_logo('logo.png', fig)
    (<matplotlib.axes.Axes at 0x00000000>,
     <matplotlib.image.AxesImage at 0x00000000>)
    """
    #Consider an alpha keyword (for watermarking)
    #Do something about margin/padding
    pos = pos.lower()
    assert(pos in ('br', 'bl', 'tl', 'tr'))
    if not hasattr(img, 'size'):
        img = matplotlib.image.imread(img)
    if fig is None:
        fig = plt.gcf()
    width = float(img.shape[1])
    height = float(img.shape[0])
    #Margin can change effective aspect ratio of figure
    margin_px = round(margin * (width if width > height else height))
    img_aspect = (width + 2 * margin_px) / (height + 2 * margin_px)
    fig_aspect = float(fig.get_figwidth()) / fig.get_figheight()
    clear_boxes = get_clear(filter_boxes(get_used_boxes(fig)), pos)
    clear_box = get_biggest_clear(clear_boxes, fig_aspect, img_aspect)
    box_width = clear_box[1][0] - clear_box[0][0]
    box_height = clear_box[1][1] - clear_box[0][1]
    if box_width / box_height * fig_aspect > img_aspect: #img uses full height
        box_width = box_height / fig_aspect * img_aspect
    else: #img uses full width, correct the height
        box_height = box_width / img_aspect * fig_aspect
    #default corners
    left = 0.0
    bottom = 0.0
    if pos[0] == 't':
        bottom = 1.0 - box_height
    if pos[1] == 'r':
        left = 1.0 - box_width
    ax = fig.add_axes([left, bottom, box_width, box_height])
    ax.axis('off')
    axesimg = ax.imshow(img)
    #Resize to include the margin
    for getter, setter in ((ax.get_xlim, ax.set_xlim),
                           (ax.get_ylim, ax.set_ylim)):
        limits = getter()
        orientation = 1 if limits[1] > limits[0] else -1 #upside-down?
        setter((limits[0] - margin_px * orientation,
                limits[1] + margin_px * orientation))
    return (ax, axesimg)


def show_used(fig=None):
    """
    Show the areas of a figure which are used/occupied by plot elements.

    This function will overplot each element of a plot with a rectangle
    showing the full bounds of that element, to see for example
    the margins and such used by a text label.

    Other Parameters
    ================
    fig : matplotlib.figure.Figure
        The figure to mark up; if not specified, the
        :func:`~matplotlib.pyplot.gcf` function will be used.

    Notes
    =====
    Calls :func:`~matplotlib.pyplot.draw` to ensure locations are
    up to date.

    Returns
    =======
    boxes : list of Rectangle
        The :class:`~matplotlib.patches.Rectangle` objects used for the overplot.

    Examples
    ========
    >>> import spacepy.plot.utils
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax0 = fig.add_subplot(211)
    >>> ax0.plot([1, 2, 3], [1, 2, 1])
    [<matplotlib.lines.Line2D at 0x00000000>]
    >>> ax1 = fig.add_subplot(212)
    >>> ax1.plot([1, 2, 3], [2, 1, 2])
    [<matplotlib.lines.Line2D at 0x00000000>]
    >>> spacepy.plot.utils.show_used(fig)
    [<matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>,
     <matplotlib.patches.Rectangle at 0x0000000>]
    """
    colors = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
    rects = []
    if fig is None:
        fig = plt.gcf()
    boxes = get_used_boxes(fig)
#    boxes = [b for b in get_used_boxes(fig)
#             if b[0, 0] != b[1, 0] and b[0, 1] != b[1, 1]]
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis('off')
    for b, c in zip(boxes, colors):
        rects.append(ax.add_patch(matplotlib.patches.Rectangle(
            b[0], b[1, 0] - b[0, 0], b[1, 1] - b[0, 1],
                     alpha=0.3, figure=fig, axes=ax, ec='none', fc=next(colors),
                     fill=True)))
    return rects

def add_arrows(lines, n=3, size=12, style='->', dorestrict=False,
               positions=False):
    '''
    Add directional arrows along a plotted line.  Useful for plotting flow
    lines, magnetic field lines, or other directional traces.

    *lines* can be either :class:`~matplotlib.lines.Line2D`, a list or tuple
    of lines, or a :class:`~matplotlib.collections.LineCollection` object.

    For each line, arrows will be added using 
    :meth:`~matplotlib.axes.Axes.annotate`.  Arrows will be spread evenly
    over the line using the number of points in the line as the metric for
    spacing.  For example, if a line has 120 points and 3 arrows are requested,
    an arrow will be added at point number 30, 60, and 90.
    Arrow color and alpha is obtained from the parent line.

    Parameters
    ==========
    lines : :class:`~matplotlib.lines.Line2D`, a list/tuple, or :class:`~matplotlib.collections.LineCollection`
        A single line or group of lines on which to place the arrows.
        Arrows inherent color and transparency (alpha) from the line on which
        they are placed.


    Other Parameters
    ================
    n : integer
        Number of arrows to add to each line; defaults to 3.

    size : integer
        The size of the arrows in points.  Defaults to 12.

    style : string
        Set the style of the arrow via :class:`~matplotlib.patches.ArrowStyle`,
        e.g. '->' (default)

    dorestrict : boolean
        If True (default), only points along the line within the current
        limits of the axes will be considered when distributing arrows.

    positions : Nx2 array
        N must be the number of lines provided via the argument *lines*.
        If provided, only one arrow will be placed per line.
        *positions* sets the explicit location of each arrow for each line as
        X-Y space in Axes coordinates.

    Returns
    =======
    None

    Notes
    =====
    The algorithm works by dividing the line in to *n*+1 segments and
    placing an arrow between each segment, endpoints excluded.  Arrows span
    the shortest distance possible, i.e., two adjacent points along a line.
    For lines that are long spatially but sparse in points, the arrows will
    have long tails that may extend beyond axes bounds.  For explicit positions,
    the arrow is placed at the point on the curve closest to that position
    and the exact position is not always attainable.  A maximum number of arrows
    equal to one-half of the number of points in a line per line will be
    created, so not all lines will receive *n* arrows.

    Examples
    ========
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> x = np.arange(1, 10, .01)
    >>> y = np.sin(x)
    >>> line = plt.plot(x,y)
    >>> add_arrows(line, n=15, style='-|>')

    '''
    from numpy import floor, argmin
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    # Convert our line, lines, or LineCollection into
    # a series of arrays with the x/y data.
    # Also, grab the axes, colors, and line alpha.
    def _parse_line_list(lines): #handles lists of lines
        data = [x.get_xydata() for x in lines]
        cols = [x.get_color()  for x in lines]
        alph = [x.get_alpha()  for x in lines]
        ax   = lines[0].axes
        return data, cols, alph, ax

    if hasattr(lines, 'axes'):  # Matplotlib-like objects:
        try: # Collection-like: use "get paths" to extract points into list
            ax   = lines.axes
            data = [x.vertices for x in lines.get_paths()]
            cols = lines.get_colors()
            alph = len(data) * [lines.get_alpha()]
        except AttributeError: # Line2D-like: put in list and parse
            data, cols, alph, ax = _parse_line_list( [lines] )
    elif isinstance(lines, (tuple, list)): # List/tuple of Line-like objects:
        try:
            data, cols, alph, ax = _parse_line_list( lines )
        except AttributeError:
            raise ValueError('Non-Line2d-like item found in list of lines.')
    else:
        raise ValueError('Unknown input type for lines.')

    # Check to make sure that we have as many colors as lines.
    # With LineCollections, this isn't always the case!
    if len(data)>1 and len(cols)==1:
        cols=list(cols)*len(data)

    # Grab plot limits (sometimes needed):
    xlim, ylim = ax.get_xlim(), ax.get_ylim()

    if not isinstance(positions, bool):
        # Explicitly set positions of arrows, one per line.
        for l, c, a, p in zip(data, cols, alph, positions):
            # Get x-y points of line:
            x, y = l[:,0], l[:,1]
            # Get point on line closest to desired position.
            # Ensure we have at least 1 point after to set direction.
            i = min(argmin(abs(x - p[0])), x.size - 2)
            j = min(argmin(abs(y - p[1])), y.size - 2)

            # Annotate, matching color/alpha:
            ax.annotate('',xytext=(x[i], y[j]), xy=(x[i+1], y[j+1]), alpha=a,
                        arrowprops=dict(arrowstyle=style,color=c),size=size)
        # Nothing else to do at this point.
        return

    # Get positions for arrows and add them:
    for l, c, a in zip(data, cols, alph):
        # Get x-y points of line:
        x, y = l[:,0], l[:,1]
        # Restrict to axes limits as necessary:
        if dorestrict:
            loc = (x>=xlim[0])&(x<=xlim[1])&(y>=ylim[0])&(y<=ylim[1])
            x, y = x[loc], y[loc]

        # Get size of lines.  Skip those with too few points:
        npts = x.size
        if npts <= 3: continue

        # If number of arrows is greater than the number of
        # points along the line/2, reduce the number of arrows.
        # Guarantee at least one arrow per line.
        n_now = max(1, min(n, int(floor(npts/2))) )

        # Place evenly-spaced arrows:
        for iArr in range(0,n_now):
            # Get location as fraction of no. of points along line:
            i  = int( floor((iArr+1) * npts/(n_now+1)) )
            # Place an arrow:
            ax.annotate('',xytext=(x[i], y[i]), xy=(x[i+1], y[i+1]), alpha=a,
                        arrowprops=dict(arrowstyle=style,color=c),size=size)
