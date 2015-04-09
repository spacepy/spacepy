#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. module:: spacepy.datamanager

The datamanager classes and functions are useful for locating the correct
data file for a particular day and manipulating data and subsets in a generic
way.

Authors: Jon Niehof

Institution: University of New Hampshire

Contact: Jonathan.Niehof@unh.edu

Copyright 2015


About datamanager
-----------------


Examples
--------

"""

import datetime
import os.path
import re

import numpy


class DataManager(object):
    """
    THIS CLASS IS NOT YET COMPLETE, doesn't do much useful.

    Will have to do something that allows the config file to specify regex
    and other things, and then just the directory to be changed (since
    regex, etc. 

    Parameters
    ==========
    directories : list
        A list of directories that might contain the data
    file_fmt : string
        Regular expression that matches the files desired. Will also recognize
        strftime parameters %w %d %m %y %Y %H %M %s %j %U %W, all zero-pad.
        https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
        Can have subdirectory reference, but separator should be unix-style,
        with no leading slash.

    period : string
        Size of file; can be a number followed by one of d, m, y, H, M, s.
        Anything else is assumed to be "irregular" and files treated
        as if there are neither gaps nor overlaps in the sequence.
        If not specified, will be assumed to match one count of the
        smallest unit in the format string.

    Examples
    =======
    """
    def __init__(self, directories, file_fmt, descend=False, period=None):
        """
        """
        #Duplicate from class docstring. Consider autoclass_content = "both"
        #Convert to system path format
        directories = [os.path.expandvars(os.path.expanduser(
            os.path.normpath(d))) for d in directories]
        file_fmt = os.path.normpath(file_fmt)
        #The period matching might go into RePath
        if period is None:
            period = next(
                ('1{0}'.format(i) for i in 'sMHdmyY' if i in file_fmt),
                None)
        elif not re.match(r'^\d+[yYmdHMs]$', period):
            period = None #irregular, fun!
        self.directories = directories
        self.period = period
        self.file_fmt = RePath(file_fmt)

    def get_filename(self, dt):
        """
        Returns the filename corresponding to a particular point in time
        """
        if self.period:
            this_re = datetime.datetime.strftime(dt, self.file_fmt)
        else:
            raise NotImplementedError

    def _files_matching(self, file_re):
        """
        Return all the files matching a particular regular expression
        """
        #Use os.walk. If descend is False, only continue for matching
        #the re to this point. If True, compare branch to entire re but
        #walk everything
        for d in self.directories:
            for (dirpath, dirnames, filenames) in \
                os.walk(d, topdown=True, followlinks=True):
                #dirpath is FULL DIRECTORY to this point
                relpath = dirpath[len(d):]
                if not self.descend:
                    #Prune dirnames based on whether they match the re
                    pass

class RePath(object):
    """
    A path based on a regular expression + time format

    Parameters
    ==========
    expression : string
        Regular expression that matches the files desired. Will also recognize
        strftime parameters %w %d %m %y %Y %H %M %s %j %U %W, all zero-pad.
        https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior
        Can have subdirectory reference, but separator should be unix-style,
        with no leading slash.
        Matching is normally done against entire string, anchors should NOT
        be included.
    """
    fmt_to_re = { 'd': r'[0-3]\d', 'm': r'[0-2]\d', 'y': r'\d\d', 'Y': r'\d{4}',
                  'H': r'[0-5]\d', 'M': r'[0-5]\d', 's': r'[0-6]\d',
                  'j': r'[0-3]\d\d', 'U': r'[0-5]\d', 'W': r'[0-5]\d'
                  }
    def __init__(self, expression):
        self.file_fmt = expression
        self.file_re = re.sub(r'(?!(?:%%)+)%([wdmyYHMsjUW])',
                             lambda x: self.fmt_to_re(x.group(1)),
                              expression)
        self.file_fmt_split = self.path_split(self.file_fmt)
        self.file_re_split = self.path_split(self.file_re)

    def match(self, string, dt=None):
        """
        Matches a string against the entire path, optionally anchored
        at a particular date/time

        Other Parameters
        ================
        dt : datetime.datetime
            The time to specifically match; otherwise matches all files.
        """
        return re.match('^' + (datetime.datetime.strftime(dt, self.file_fmt)
                               if dt else self.file_re) + '$',
                        string)

    def match_end(self, string, dt=None):
        """
        Matches a string against the end of the path, optionally anchored
        at a particular date/time.
        Note this matches the last several elements of the path, NOT just
        the last several characters, i.e., foo/bar will not match oo/bar

        Other Parameters
        ================
        dt : datetime.datetime
            The time to specifically match; otherwise matches all files.

        Returns
        =======
        out : re.MatchObject
        """
        return re.match(
            '^'+
            self.path_slice(datetime.datetime.strftime(dt, self.file_fmt)
                            if dt else self.file_re,
                            -len(self.path_split(string)), 999)
            + '$', string)

    def match_start(self, string, dt=None):
        """
        Matches a string against the start of the path, optionally anchored
        at a particular date/time.
        Note this matches the first several elements of the path, NOT just
        the last several characters, i.e., foo/ba will not match foo/bar

        Other Parameters
        ================
        dt : datetime.datetime
            The time to specifically match; otherwise matches all files.
        """
        return re.match(
            '^'+
            self.path_slice(datetime.datetime.strftime(dt, self.file_fmt)
                            if dt else self.file_re,
                            0, len(self.path_split(string)))
            + '$', string)
    
    @staticmethod
    def path_split(path):
        """
        Break a path apart into a list for each path element.

        Parameters
        ==========
        path : str
            Path to split

        Returns
        =======
        out : list of str
            One path element (directory or file) per item
        """
        res = []
        while path:
            path, tail = os.path.split(path)
            res.insert(0, tail)
        return res

    @staticmethod
    def path_slice(path, start, stop=None, step=None):
        """
        Slice a path by elements, as in getitem or []

        Parameters
        ==========
        path : str
            Path to slice
        start : int
            First path element to return, or the only element if
            stop and step are not specified.

        Other Parameters
        ================
        stop : int
            First path element to NOT return, i.e., one past the last
        step : int
            Increment between each.

        Returns
        =======
        out : str
            Selection of the path

        Examples
        ========
        >>> path = "foo/bar/baz"
        >>> spacepy.datamanager.DataManager.path_slice(path, 1)
        "bar"
        >>> spacepy.datamanager.DataManager.path_slice(path, 0, step=2)
        "foo/baz"
        >>> spacepy.datamanager.DataManager.path_slice(path, 1, 3)
        "bar/baz"
        """
        if stop is None and step is None:
            return RePath.path_split(path)[start]
        else:
            return os.path.join(RePath.path_split(path)[start:stop:step])


def insert_fill(times, data, fillval=numpy.nan, tol=1.5, doTimes=True):
    """
    Populate gaps in data with fill.

    Continuous data are often treated differently from discontinuous data,
    e.g., matplotlib will draw lines connecting data points but break the line
    at fill. Often data will be irregularly sampled but also contain large
    gaps that are not explicitly marked as fill. This function adds a single
    record of explicit fill to each gap, defined as places where the spacing
    between input times is a certain multiple of the median spacing.

    Parameters
    ==========
    times : sequence
        Values representing when the data were taken. Must be one-dimensional,
        i.e., each value must be scalar. Not modified 
    data : sequence
        Input data.

    Other Parameters
    ================
    fillval : 
        Fill value, same type as ``data``. Default is ``numpy.nan``. If scalar,
        will be repeated to match the shape of ``data`` (minus the time axis).

        .. note::
            The default value of ``nan`` will not produce good results with
            integer input.

    tol : float
        Tolerance. A single fill value is inserted between adjacent values
        where the spacing in ``times`` is greater than ``tol`` times the median
        of the spacing across all ``times``. The inserted time for fill is
        halfway between the time on each side. (Default 1.5)
    doTimes : boolean
        If True (default), will return a tuple of the times (with new values
        inserted for the fill records) and the data with new fill values.
        If False, will only return the data -- useful for applying fill to
        multiple arrays of data on the same timebase.

    Raises
    ======
    ValueError : if can't identify the time axis of data
        Try using :func:`numpy.rollaxis` to put the time axis first in both
        ``data`` and ``times``.

    Returns
    =======
    times, data : tuple of sequence
        Copies of input times and data, fill added in gaps (``doTimes`` True)
    data : sequence
        Copy of input data, with fill added in gaps (``doTimes`` False)

    Examples
    ========
    This example shows simple hourly data with a gap, populated with fill.
    Note that only a single fill value is inserted, to break the sequence
    of valid data rather than trying to match the existing cadence.

    >>> import datetime
    >>> import numpy
    >>> import spacepy.datamanager
    >>> t = [datetime.datetime(2012, 1, 1, 0),
             datetime.datetime(2012, 1, 1, 1),
             datetime.datetime(2012, 1, 1, 2),
             datetime.datetime(2012, 1, 1, 5),
             datetime.datetime(2012, 1, 1, 6)]
    >>> temp = [30.0, 28, 27, 32, 35]
    >>> filled_t, filled_temp = spacepy.datamanager.insert_fill(t, temp)
    >>> filled_t
    array([datetime.datetime(2012, 1, 1, 0, 0),
           datetime.datetime(2012, 1, 1, 1, 0),
           datetime.datetime(2012, 1, 1, 2, 0),
           datetime.datetime(2012, 1, 1, 3, 30),
           datetime.datetime(2012, 1, 1, 5, 0),
           datetime.datetime(2012, 1, 1, 6, 0)], dtype=object)
    >>> filled_temp
    array([ 30.,  28.,  27.,  nan,  32.,  35.])

    .. plot::
       :include-source:

       This example plots "gappy" data with and without explicit fill values.

       >>> import matplotlib.pyplot as plt
       >>> import numpy
       >>> import spacepy.datamanager
       >>> x = numpy.append(numpy.arange(0, 6, 0.1), numpy.arange(12, 18, 0.1))
       >>> y = numpy.sin(x)
       >>> xf, yf = spacepy.datamanager.insert_fill(x, y)
       >>> fig = plt.figure()
       >>> ax0 = fig.add_subplot(211)
       >>> ax0.plot(x, y)
       >>> ax1 = fig.add_subplot(212)
       >>> ax1.plot(xf, yf)
       >>> plt.show()

    """
    times = numpy.asanyarray(times)
    data = numpy.asanyarray(data)
    assert(len(times.shape) == 1)
    if len(times) == data.shape[0]:
        timeaxis = 0
    else:
        matches = numpy.nonzero(numpy.asanyarray(data.shape) == len(times))[0]
        if len(matches) != 1:
            raise ValueError(
                "Unable to uniquely match shape of data to count of times.")
        timeaxis = matches[0]
    fillshape = numpy.delete(data.shape, timeaxis) #shape of data w/o time axis
    if numpy.shape(fillval) != fillshape:
        if numpy.shape(fillval) == ():
            fillval = numpy.tile(fillval, fillshape)
        else:
            raise ValueError("Cannot match shape of fill to shape of data")
    diff = numpy.diff(times)
    if hasattr(diff[0], 'seconds'): #datetime
        diff = numpy.vectorize(lambda x: x.days * 3600.0 + x.seconds + x.microseconds / 1.0e6)(diff)
    idx = numpy.nonzero((diff > (numpy.median(diff) * tol)))[0] + 1
    data = numpy.insert(data, idx, numpy.repeat(fillval, len(idx)),
                        axis=timeaxis) #NOOP if no fill
    if not doTimes:
        return data
    try:
        filltimes = (times[idx] + times[idx - 1]) / 2
    except TypeError:
        filltimes = times[idx - 1] + numpy.vectorize(lambda x: datetime.timedelta(seconds=x / 2.0))(diff[idx - 1])
    times = numpy.insert(times, idx, filltimes)
    return times, data
