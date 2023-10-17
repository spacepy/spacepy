#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The datamanager classes and functions are useful for locating the correct
data file for a particular day and manipulating data and subsets in a generic
way.

Authors: Jon Niehof

Institution: University of New Hampshire

Contact: Jonathan.Niehof@unh.edu

Copyright 2015-2020 contributors


About datamanager
-----------------


Examples
--------
Examples go here

.. currentmodule:: spacepy.datamanager

.. rubric:: Classes

.. autosummary::
    :template: clean_class.rst

    DataManager

.. rubric:: Functions

.. autosummary::

    apply_index
    array_interleave
    axis_index
    flatten_idx
    insert_fill
    rebin
    rev_index
    values_to_steps
"""

__all__ = ["DataManager", "apply_index", "array_interleave", "axis_index",
           "flatten_idx", "insert_fill", "rebin", "rev_index",
           "values_to_steps"]

import datetime
import operator
import os.path
import posixpath
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
    ========
    """
    def __init__(self, directories, file_fmt, descend=False, period=None):
        """
        """
        #Duplicate from class docstring. Consider autoclass_content = "both"
        #Convert to system path format
        directories = [os.path.expandvars(os.path.expanduser(
            os.path.normpath(d))) for d in directories]
        file_fmt = posixpath.normpath(file_fmt)
        #The period matching might go into RePath
        if period is None:
            period = next(
                ('1{0}'.format(i) for i in 'sMHdmyY' if i in file_fmt),
                None)
        elif not re.match(r'^\d+[yYmdHMs]$', period):
            period = None #irregular, fun!
        self.directories = directories
        self.descend = descend
        self.period = period
        self.file_fmt = RePath(file_fmt)

    def get_filename(self, dt):
        """
        Returns the filename corresponding to a particular point in time
        """
        if self.period:
            flist = self.files_matching(dt)
        else:
            raise NotImplementedError
        #Now figure out the priority..

    def files_matching(self, dt=None):
        """
        Return all the files matching this file format

        Parameters
        ==========
        dt : datetime
            Optional; if specified, match only files for this date.

        Returns
        =======
        out : generator
            Iterates over every file matching the format specified at
            creation. Note this is specified in native path format!
        """
        #Use os.walk. If descend is False, only continue for matching
        #the re to this point. If True, compare branch to entire re but
        #walk everything
        for d in self.directories:
            native_d = os.path.normpath(d) #Do the walk in native paths
            for (dirpath, dirnames, filenames) in \
                os.walk(native_d, topdown=True, followlinks=True):
                #dirpath is FULL DIRECTORY to this point, make relative
                relpath = dirpath[len(native_d) + 1:]
                #Convert to POSIX for comparisons
                if os.path.sep != posixpath.sep and relpath:
                    relpath = posixpath.join(*RePath.path_split(
                        relpath, native=True))
                if not self.descend:
                    if relpath and not \
                       self.file_fmt.match(relpath, dt, 'start'):
                        continue
                    for i in range(-len(dirnames), 0):
                        dirname = dirnames[i]
                        if os.path.sep != posixpath.sep and dirname:
                            dirname = posixpath.join(*RePath.path_split(
                                dirname, native=True))
                        if not self.file_fmt.match(posixpath.join(
                                relpath, dirname), dt, 'start'):
                            del dirnames[i]
                for f in filenames:
                    if self.file_fmt.match(posixpath.join(relpath, f), dt,
                                           'end' if self.descend else None):
                        yield os.path.join(dirpath, f)


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
        #This should have backrefs where the same format occurs twice
        #(year in both directory and file, for example)
        self.file_re = re.sub(r'(?!(?:%%)+)%([wdmyYHMsjUW])',
                             lambda x: self.fmt_to_re[x.group(1)],
                              expression)
        self.file_fmt_split = self.path_split(self.file_fmt)
        self.file_re_split = self.path_split(self.file_re)

    def match(self, string, dt=None, where=None):
        """
        Matches a string against a path or part thereof, optionally anchored
        at a particular date/time

        Other Parameters
        ================
        dt : datetime.datetime
            The time to specifically match; otherwise matches all files.
        where : str
            Where to match: None to match the string to the entire path
            (default); ``'start'`` to match entire string against the start
            of the path; ``'end'`` to match entire path against end of string.
            Note this matches the last elements of the string, not just
            the last characters, i.e., ``oo/bar`` will not match ``foo/bar``.
            Similarly, ``'start'`` matches the first elements of the path,
            i.e., ``foo/ba`` will not match ``foo/bar``
            ``start`` matches subset of path, i.e., the string is a directory
            that may contain full matches further down the tree. ``end``
            matches a subset of the string, i.e, the string is a
            path to a file that would be a complete match except it has
            additional path elements leading. This is the order that tends to
            be useful.

        Returns
        =======
        out : re.MatchObject
            The result of the match.
        """
        datestr = (datetime.datetime.strftime(dt, self.file_fmt)
                   if dt else self.file_re)
        if where is None:
            pat = datestr
        elif where.lower() == 'end': #Cut down string to match path pattern
            pat = datestr
            string = self.path_slice(string, -len(self.file_re_split), None, 1)
        elif where.lower() == 'start': #Does path pattern start like string?
            pat = self.path_slice(datestr,
                                  0, len(self.path_split(string)))
        else:
            raise(ValueError("where must be 'start', 'stop', or None, not {0}".
                             format(where)))
        return re.match('^' + pat + '$', string)

    @staticmethod
    def path_split(path, native=False):
        """
        Break a path apart into a list for each path element.

        Parameters
        ==========
        path : str
            Path to split
        native : bool
            Is this a native path or UNIX-style? (default False, UNIX)

        Returns
        =======
        out : list of str
            One path element (directory or file) per item
        """
        split = os.path.split if native else posixpath.split
        res = []
        while path:
            base, tail = split(path)
            if base == path: #No further splitting
                res.insert(0, path)
                break
            res.insert(0, tail)
            path = base
        return res

    @staticmethod
    def path_slice(path, start, stop=None, step=None, native=False):
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
            First path element to NOT return, i.e., one past the last.
            If ``stop`` is not specified but ``step`` is, return to end.
        step : int
            Increment between each.
        native : bool
            Is this a native path or UNIX-style? (default False, UNIX)

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
            return RePath.path_split(path, native=native)[start]
        else:
            join = (os.path if native else posixpath).join
            return join(*RePath.path_split(path, native=native)
                        [start:stop:step])


def insert_fill(times, data, fillval=numpy.nan, tol=1.5, absolute=None, doTimes=True):
    """Populate gaps in data with fill.

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
        Tolerance. A single fill value is inserted between adjacent
        values where the spacing in ``times`` is strictly greater than
        ``tol`` times the median of the spacing across all
        ``times``. The inserted time for fill is halfway between the
        time on each side. (Default 1.5)
    absolute :
        An absolute value for maximum spacing, of a type that would result from
        a difference in ``times``. If specified, ``tol`` is ignored and any gap
        strictly larger than ``absolute`` will have fill inserted.
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
    if (numpy.shape(fillval) != fillshape).any():
        if numpy.shape(fillval) == ():
            fillval = numpy.tile(fillval, fillshape)
        else:
            raise ValueError("Cannot match shape of fill to shape of data")
    diff = numpy.diff(times)
    if hasattr(diff[0], 'seconds'): #datetime
        diff = numpy.vectorize(lambda x: x.days * 86400.0 + x.seconds +
                               x.microseconds / 1.0e6)(diff)
        if absolute is not None:
            absolute = absolute.days * 86400.0 + absolute.seconds + \
                       absolute.microseconds / 1.0e6
    if absolute is None:
        idx = numpy.nonzero(diff > (numpy.median(diff) * tol))[0] + 1
    else:
        idx = numpy.nonzero(diff > absolute)[0] + 1
    data = numpy.insert(data, idx, numpy.repeat(fillval, len(idx)),
                        axis=timeaxis) #NOOP if no fill
    if not doTimes:
        return data
    try:
        filltimes = (times[idx] + times[idx - 1]) / 2.0
    except TypeError:
        filltimes = times[idx - 1] + numpy.vectorize(lambda x: datetime.timedelta(seconds=x / 2.0))(diff[idx - 1])
    times = numpy.insert(times, idx, filltimes)
    return times, data


def apply_index(data, idx):
    """Apply an array of indices to data.

    Most useful in dealing with the output from :func:`numpy.argsort`, and
    best explained by the example.

    Parameters
    ==========
    data : array
        Input data, at least two dimensional. The 0th dimension is treated as
        a "time" or "record" dimension.
    idx : sequence
        2D index to apply to the import data. The 0th dimension must be the 
        same size as ``data``'s 0th dimension. Dimension 1 must be the same
        size as one other dimension in data (the first match found is used);
        this is referred to as the "index dimension."

    Raises
    ======
    ValueError : if can't match the shape of data and indices

    Returns
    =======
    data : sequence
        View of ``data``, with index applied. For each index of the 0th
        dimension, the values along the index dimension are obtained by applying
        the value of ``idx`` at the same index in the 0th dimension. This is
        repeated across any other dimensions in ``data``.

        .. warning::
            No guarantee is made whether the returned data is a copy
            of the input data. Modifying values in the input may change
            the values of the input. Call :meth:`~numpy.ndarray.copy` if
            a copy is required.

    Examples
    ========
    Assume ``flux`` is a 3D array of fluxes, with a value for each of
    time, pitch angle, and energy. Assume energy is not necessarily 
    constant in time, nor is ordered in the energy dimension. If
    ``energy`` is a 2D array of the energies as a function of energy
    step for each time, then the following will sort the flux at each
    time and pitch angle in energy order.

    >>> idx = numpy.argsort(energy, axis=1)
    >>> flux_sorted = spacepy.datamanager.apply_index(flux, idx)
    """
    data = numpy.asanyarray(data)
    idx = numpy.asanyarray(idx)
    if len(idx.shape) != 2:
        raise ValueError("idx must have dimensions 2, not {0}".format(
            len(idx.shape)))
    if len(data.shape) < 2:
        raise ValueError("data must have at least dimensions 2")
    if idx.shape[0] != data.shape[0]:
        raise ValueError("data and idx must have same size in "
                         "0th dimension")
    if not idx.shape[1] in data.shape[1:]:
        raise ValueError("Size of idx dimension 1 must match a dimension in "
                         "data")
    idx_dim = data.shape[1:].index(idx.shape[1]) + 1
    return numpy.rollaxis(
        numpy.rollaxis(data, idx_dim, 1) #make time and index dim adjacent
        #get a 2d array where every element matches index of first axis
        [numpy.mgrid[0:idx.shape[0], slice(idx.shape[1])][0],
         idx, #2d array, every element is desired index of second axis
         ...] #and the other axes come along for the ride
        , 1, idx_dim + 1) #and put index dim back in place


def array_interleave(array1, array2, idx):
    """Create an array containing all elements of both array1 and array2

    ``idx`` is an index on the output array which indicates which elements will
    be populated from ``array1``, i.e., ``out[idx] == array1`` (in order.)
    The other elements of ``out`` will be filled, in order, from ``array2``.

    Parameters
    ==========
    array1 : array
        Input data.
    array2 : array
        Input data. Must have same number of dimensions as ``array1``, and all
        dimensions except the zeroth must also have the same length.
    idx : array
        A 1D array of indices on the zeroth dimension of the output array. Must
        have the same length as the zeroth dimension of ``array1``.

    Returns
    =======
    out : array
        All elements from ``array1`` and ``array2``, interleaved according
        to ``idx``.

    Examples
    ========
    >>> import numpy
    >>> import spacepy.datamanager
    >>> a = numpy.array([10, 20, 30])
    >>> b = numpy.array([1, 2])
    >>> idx = numpy.array([1, 2, 4])
    >>> spacepy.datamanager.array_interleave(a, b, idx)
    array([ 1, 10, 20,  2, 30])
    """
    array1 = numpy.asanyarray(array1)
    array2 = numpy.asanyarray(array2)
    idx = numpy.asanyarray(idx)
    assert(len(array1.shape) == len(array2.shape))
    assert(array1.shape[1:] == array2.shape[1:])
    assert(array1.dtype == array2.dtype)
    outarray = numpy.empty(dtype=array1.dtype,
        shape=((array1.shape[0] + array2.shape[0],) + array1.shape[1:]))
    outarray[idx, ...] = array1
    idx_comp = numpy.ones((outarray.shape[0],), dtype=numpy.bool)
    idx_comp[idx] = False
    outarray[idx_comp, ...] = array2
    return outarray


def values_to_steps(array, axis=-1):
    """Transform values along an axis to their order in a unique sequence.

    Useful in, e.g., converting a list of energies to their steps.

    Parameters
    ==========
    array : array
        Input data.

    Other Parameters
    ================
    axis : int
        Axis along which to find the steps.

    Returns
    =======
    steps : array
        An array, the same size as ``array``, with values along ``axis``
        corresponding to the position of the value in ``array`` in a unique,
        sorted, set of the values in ``array`` along that axis. Differs from
        :func:`~numpy.argsort` in that identical values will have identical
        step numbers in the output.

    Examples
    ========
    >>> import numpy
    >>> import spacepy.datamanager
    >>> data = [[10., 12., 11., 9., 10., 12., 11., 9.],
                [10., 12., 11., 9., 14., 16., 15., 13.]]
    >>> spacepy.datamanager.values_to_steps(data)
    array([[1, 3, 2, 0, 1, 3, 2, 0],
           [1, 3, 2, 0, 5, 7, 6, 4]])
    """
    array = numpy.asanyarray(array)
    sortidx = array.argsort(axis=axis)
    steps = rev_index(sortidx, axis=axis)
    d = numpy.diff(apply_index(array, sortidx), axis=axis)
    #Everywhere in SORTED array that the VALUE is same as one before
    same = numpy.insert(d==0, 0, False, axis=axis)
    #Number of duplicate indices BEFORE current index, in SORTED array
    delta = numpy.cumsum(same, axis=-1)
    #Get delta into the unsorted frame, and correct for uniqueness
    return steps - delta.ravel()[
        flatten_idx(rev_index(sortidx, axis=axis), axis=axis)]\
        .reshape(steps.shape)


def flatten_idx(idx, axis=-1):
    """Convert multidimensional index into index on flattened array.

    Convert a multidimensional index, that is values along a particular axis,
    so that it can derefence the flattened array properly. Note this is not
    the same as :func:`~numpy.ravel_multi_index`.

    Parameters
    ==========
    idx : array
        Input index, i.e. a list of elements along a particular axis,
        in the style of :func:`~numpy.argsort`.

    Other Parameters
    ================
    axis : int
        Axis along which ``idx`` operates, defaults to the last axis.

    Returns
    =======
    flat : array
        A 1D array of indices suitable for indexing the flat version of the
        array 

    See Also
    ========
    apply_index

    Examples
    ========
    >>> import numpy
    >>> import spacepy.datamanager
    >>> data = numpy.array([[3, 1, 2], [3, 2, 1]])
    >>> idx = numpy.argsort(data, -1)
    >>> idx_flat = spacepy.datamanager.flatten_idx(idx)
    >>> data.ravel() #flat array
    array([3, 1, 2, 3, 2, 1])
    >>> idx_flat #indices into the flat array
    array([1, 2, 0, 5, 4, 3])
    >>> data.ravel()[idx_flat] #index applied to the flat array
    array([1, 2, 3, 1, 2, 3])
    """
    idx = numpy.asanyarray(idx)
    if not idx.dtype.kind in ('i', 'u'):
        idx = idx.astype(int)
    preshape = idx.shape[:axis]
    postshape = idx.shape[axis:]
    stride = int(numpy.prod(postshape[1:])) #1 if applied to empty
    #The index on this axis moves stride elements in flat
    outidx = idx.flatten() * stride #makes a copy
    #First add the offsets to get us to [..., idx @ axis = 0, 0...)
    outidx += numpy.repeat(
        numpy.arange(0, len(outidx), int(numpy.prod(postshape)),
                     dtype=idx.dtype),
        numpy.prod(postshape))
    #Now offsets for non-zero on the trailing axes [0, 0, ... 0@axis, ...]
    outidx += numpy.tile(numpy.arange(0, stride, dtype=idx.dtype),
                           int(numpy.prod(preshape)) * idx.shape[axis])
    return outidx


def axis_index(shape, axis=-1):
    """Returns array of indices along axis, for all other axes

    Parameters
    ==========
    shape : tuple
        Shape of the output array

    Other Parameters
    ================
    axis : int
        Axis along which to return indices, defaults to the last axis.
       
    Returns
    =======
    idx : array
        An array of indices. The value of each element is that element's index
        along ``axis``.

    See Also
    ========
    numpy.mgrid : This function is a special case

    Examples
    ========
    For a shape of ``(i, j, k, l)`` and ``axis`` = -1,
    ``idx[i, j, k, :] = range(l)`` for all ``i``, ``j``, ``k``.

    Similarly, for the same shape and ``axis = 1``,
    ``idx[i, :, k, l] = range(j)`` for all ``i``, ``k``, ``l``.

    >>> import numpy
    >>> import spacepy.datamanager
    >>> spacepy.datamanager.axis_index((5, 3))
    array([[0, 1, 2],
           [0, 1, 2],
           [0, 1, 2],
           [0, 1, 2],
           [0, 1, 2]])
    >>> spacepy.datamanager.axis_index((5, 3), 0)
        array([[0, 0, 0],
               [1, 1, 1],
               [2, 2, 2],
               [3, 3, 3],
               [4, 4, 4]])
"""
    return operator.getitem(numpy.mgrid, [slice(i) for i in shape])[axis]


def rev_index(idx, axis=-1):
    """From an index, return an index that reverses the action of that index

    Essentially, ``a[idx][rev_index(idx)] == a``

    .. note::
        This becomes more complicated in multiple dimensions, due to the
        vagaries of applying a multidimensional index.

    Parameters
    ==========
    idx : array
        Indices onto an array, often the output of :func:`~numpy.argsort`.

    Other Parameters
    ================
    axis : int
        Axis along which to return indices, defaults to the last axis.

    Returns
    =======
    rev_idx : array
        Indices that, when applied to an array after ``idx``, will return
        the original array (before the application of ``idx``).
 
    See Also
    ========
    apply_index

    Examples
    ========
    >>> import numpy
    >>> import spacepy.datamanager
    >>> data = numpy.array([7, 2, 4, 6, 3])
    >>> idx = numpy.argsort(data)
    >>> data[idx] #sorted
    array([2, 3, 4, 6, 7])
    >>> data[idx][spacepy.datamanager.rev_index(idx)] #original
    array([7, 2, 4, 6, 3])
"""
    #Want an idx2 such that x[idx][idx2] == x
    #idx is position to value map
    #Populate every POSITION in idx2 with the POSITION in idx that
    #has the VALUE of the idx2 position
    #searchsorted on range?
    idx_out = numpy.empty_like(idx).ravel()
    idx_out[flatten_idx(idx, axis)] = axis_index(idx.shape, axis).ravel()
    return idx_out.reshape(idx.shape)


def _find_shape(bigshape, littleshape):
    """Increase dimensionality of small array to match larger

    Add degenerate dimensions to a lesser-dimensioned array to match a
    larger-dimensioned array.

    Parameters
    ==========
    bigshape : tuple
        The shape of the larger-dimmed array
    littleshape : tuple
        The shape of the smaller-dimmed array

    Returns
    =======
    tuple
        A shape of the same length/dimensions of ``bigshape`` with
        dimension sizes taken from ``littleshape`` or of size 1. Suitable
        for calling :func:`numpy.reshape` on array of shape ``littleshape``.
    """
    newshape = []
    littleidx = 0
    for i, size in enumerate(bigshape):
        if littleidx < len(littleshape): # Still have "little" to consume
            if size == littleshape[littleidx]: # Match, just move along
                newshape.append(size)
                littleidx += 1
                continue
            if littleshape[littleidx] == 1: # Consume degenerate dim
                littleidx += 1
        newshape.append(1) # Add a degenerate dim
    if littleidx < len(littleshape): # Didn't consume them all
        raise ValueError('Unable to make shape {} compatible with shape {}'
                         .format(littleshape, bigshape))
    return tuple(newshape)


def rebin(data, bindata, bins, axis=-1, bintype='mean',
          weights=None, clip=False, bindatadelta=None):
    """Rebin one axis of input data based on values of another array

    This is clearest with an example. Consider a flux as a function of time,
    energy, and the look direction of a detector (could be multiple detectors,
    or spin sectors.) The flux is then 3-D, dimensioned Nt x Ne x Nl. Now
    consider that each look direction has an associated pitch angle that is
    a function of time and thus stored in an array Nt x Nl. Then this function
    will redimension the flux into pitch angle bins (rather than tags.)

    So consider the PA bins to have dimension Np + 1 (because it represents
    the edges, the number of bins is one less than the dimension.) Then
    the output will be dimensioned Nt x Ne x Np.

    ``bindata`` must be same or lesser dimensionality than ``data``. Any
    axes which are present must be either of size 1, or the same size as
    ``data``. So for ``data`` 100x5x20, ``bindata`` may be 100x5x20, or
    100, or 100x1x20, but not 100x20x5. This function will insert axes of
    size 1 as needed to match dimensionality.

    Parameters
    ==========
    data : :class:`~numpy.ndarray`
        N-dimensional array of data to be rebinned. :class:`~numpy.nan`
        are ignored.
    bindata : :class:`~numpy.ndarray`
        M-dimensional (M<=N) array of values to be compared to the bins.
    bins : :class:`~numpy.ndarray`
        1-D array of bin edges. Output dimension will be this size
        minus 1. Any values in ``bindata`` that don't fall in the
        bins will be omitted from the output. (See ``clip`` to change
        this behavior).

    Other Parameters
    ================
    axis : int
        Axis of `data` to rebin. This axis will disappear in the
        output and be replaced with an axis of the size of
        ``bins`` less one. (Default -1, last axis)
    bintype : str
        Type of rebinning to perform:
            mean
                Return the mean of all values in the bin (default)
            unc
                Return the quadrature mean of all values in the bin,
                for propagating uncertainty
            count
                Return the count of values that fall in each bin.
    weights : :class:`~numpy.ndarray`
        Relative weight of each sample in ``bindata``. Must be same
        shape as ``bindata``. Purely relative, i.e. the output is
        only affected based on the total of ``weights`` if ``bintype``
        is ``count``. Note if ``weights`` is specified, ``count`` returns
        the sum of the weights, not the count of individual samples.
    clip : boolean
        Clip data to the bins. If true, all input data will be assigned
        a bin and data outside the range of the bin edges will be assigned
        to the extreme bins. If false (default), input data outside the bin
        ranges will be ignored.
    bindatadelta : :class:`~numpy.ndarray`
        By default, the ``bindata`` are treated as point values. If
        ``bindatadelta`` is specified, it is treated as the half-width of
        the ``bindata``, allowing a single input value to be split between
        output bins. Must be scalar, or same shape as `bindata`. Note that
        input values are not weighted by the bin width, but by number of
        input values or by ``weights``. (Combining ``weights`` with
        ``bindatadelta`` is not comprehensively tested.)

    Returns
    =======
    :class:`~numpy.ndarray`
        ``data`` with one axis redimensioned, from its original
        dimension to the bin dimension.

    Examples
    ========
    Consider a particle flux distribution that's a function of energy and
    pitch angle. For simplicity, assume that the energy dependence is a
    simple power law and the pitch angle dependence is Gaussian, with a peak
    whose position oscillates in time over a period of about one hour. This is
    fairly non-physical but illustrative.

    First making the relevant imports::

        >>> import matplotlib.pyplot
        >>> import numpy
        >>> import spacepy.datamanager
        >>> import spacepy.plot

    The functional form of the flux is then::

        >>> def j(e, t, a):
        ...     return e ** -2 * (1 / (90 * numpy.sqrt(2 * numpy.pi))) \\
        ...     * numpy.exp(
        ...         -0.5 * ((a - 90 + 90 * numpy.sin(t / 573.)) / 90.) ** 2)

    Illustrating the flux at one energy as a function of pitch angle::

        >>> times = numpy.arange(0., 7200, 5)
        >>> alpha = numpy.arange(0, 181., 2)
        # Add a dimension so the flux is a 2D array
        >>> flux = j(1., numpy.expand_dims(times, 1),
        ...          numpy.expand_dims(alpha, 0))
        >>> spacepy.plot.simpleSpectrogram(times, alpha, flux, cb=False,
        ...                                ylog=False)
        >>> matplotlib.pyplot.ylabel('Pitch angle (deg)')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux at 1 MeV')

    .. plot:: pyplots/datamanager_rebin_1.py

    Or the flux at one pitch angle as a function of energy::

        >>> energies = numpy.logspace(0, 3, 50)
        >>> flux = j(numpy.expand_dims(energies, 0),
        ...          numpy.expand_dims(times, 1), 90.)
        >>> spacepy.plot.simpleSpectrogram(times, energies, flux, cb=False)
        >>> matplotlib.pyplot.ylabel('Energy (MeV)')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux at 90 degrees')

    .. plot:: pyplots/datamanager_rebin_2.py

    The measurement is usually not aligned with a pitch angle grid, and
    the detector pointing in pitch angle space usually varies with time.
    Taking a very simple case of eight detectors that sweep through pitch
    angle space in an organized fashion at ten degrees per minute, the
    measured pitch angle as a function of detector and time is::

        >>> def pa(d, t):
        ...     return (d * 22.5 + t * (2 * (d % 2) - 1)) % 180
        >>> lines = matplotlib.pyplot.plot(
        ...     times, pa(numpy.arange(8).reshape(1, -1), times.reshape(-1, 1)),
        ...     marker='o', ms=1, linestyle='')
        >>> matplotlib.pyplot.legend(lines,
        ...     ['Detector {}'.format(i) for i in range(4)], loc='best')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.ylabel('Pitch angle (deg)')
        >>> matplotlib.pyplot.title('Measured pitch angle by detector')

    .. plot:: pyplots/datamanager_rebin_3.py

    Assuming a coarser measurement in time and energy than used to illustrate
    the distribution above, the measured flux as a function of time, detector,
    and energy is constructed::

        >>> times = numpy.arange(0., 7200, 300) #5 min cadence
        >>> alpha = pa(numpy.arange(8).reshape(1, -1), times.reshape(-1, 1))
        >>> energies = numpy.logspace(0, 3, 10) #10 energy channels (3/decade)
        # Every dimension (t, detector, e) gets its own numpy axis
        >>> flux = j(numpy.reshape(energies, (1, 1, -1)),
                     numpy.reshape(times, (-1, 1, 1)),
                     numpy.expand_dims(alpha, -1))
        >>> flux.shape
        (24, 8, 10)

    The flux at an energy as a function of detector is not very useful::

        >>> spacepy.plot.simpleSpectrogram(times, numpy.arange(8),
        ...                                flux[..., 0], cb=False, ylog=False)
        >>> matplotlib.pyplot.ylabel('Detector')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux at 1 MeV')

    .. plot:: pyplots/datamanager_rebin_4.py

    As a function of energy for one detector, the energy dependence is
    apparent but time and pitch angle effects are confounded::

        >>> spacepy.plot.simpleSpectrogram(times, energies, flux[:, 0, :],
        ...                                cb=False)
        >>> matplotlib.pyplot.ylabel('Energy (MeV)')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux in detector 0')

    .. plot:: pyplots/datamanager_rebin_5.py

    What is needed is to recover the array of flux dimensioned by
    time, pitch angle, and energy, with appropriate pitch angle bins.
    The assumption is that the pitch angle as a function of time and
    detector is measured and thus the ``alpha`` array is available.
    Using that array, ``rebin`` can change flux from time, detector,
    energy bins to time, pitch angle, energy bins. The axis 1 changes
    from a detector dimension to pitch angle::

        >>> pa_bins = numpy.arange(0, 181, 36)
        >>> flux_by_pa = spacepy.datamanager.rebin(
        ...     flux, alpha, pa_bins, axis=1)
        >>> flux_by_pa.shape
        (24, 6, 10)

    This can then be visualized. The pitch angle coverage is not perfect,
    but the original shape of the distribution is apparent, and further
    analysis can be performed on the regular pitch angle grid::

        >>> spacepy.plot.simpleSpectrogram(times, pa_bins, flux_by_pa[..., 0],
        ...                                cb=False, ylog=False)
        >>> matplotlib.pyplot.ylabel('Pitch angle (deg)')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux at 1MeV')

    .. plot:: pyplots/datamanager_rebin_6.py

    Or by energy::

        >>> spacepy.plot.simpleSpectrogram(times, energies, flux_by_pa[:, 2, :],
        ...                                cb=False)
        >>> matplotlib.pyplot.ylabel('Energy (MeV)')
        >>> matplotlib.pyplot.xlabel('Time (sec)')
        >>> matplotlib.pyplot.title('Flux at 90 degrees')

    .. plot:: pyplots/datamanager_rebin_7.py    

    ``rebin`` can be used for higher dimension data, if the pitch angle
    itself depends on energy (e.g. if an energy sweep takes substantial
    time), and to propagate uncertainties through the rebinning. It can
    also be used to rebin on the time axis, e.g. for transforming time
    base.
    """
    makefloat = lambda x: x if isinstance(x, numpy.ndarray)\
                and issubclass(x.dtype.type, numpy.floating)\
                else numpy.require(x, dtype=numpy.float_)
    bintype = bintype.lower()
    assert bintype in ('mean', 'count', 'unc')
    data = makefloat(data)
    bindata = makefloat(bindata)
    if axis < 0:
        axis = len(data.shape) + axis
    binnedshape = _find_shape(data.shape, bindata.shape)
    if weights is not None:
        weights = makefloat(weights)
        assert weights.shape == bindata.shape
        weights = numpy.reshape(weights, binnedshape)
    if bindatadelta is not None:
        if not numpy.isscalar(bindatadelta):
            bindatadelta = makefloat(bindatadelta)
            assert bindata.shape == bindatadelta.shape
            bindatadelta = numpy.reshape(bindatadelta, binnedshape)
            bindatadelta = numpy.rollaxis(
                bindatadelta, axis=axis, start=len(data.shape))
    # Add axes to match shapes. Move the axis to rebin to the end of the line.
    bindata = numpy.reshape(bindata, binnedshape)
    indata = data # Holding a reference without transformations
    data = numpy.rollaxis(data, axis=axis, start=len(data.shape))
    bindata = numpy.rollaxis(bindata, axis=axis, start=len(data.shape))
    nbins = len(bins) - 1
    # Get an array, last two axes of which are a matrix of
    # (newbin, oldbin) that is the fraction of oldbin (the input bin) that
    # falls within newbin
    # Special case if bindata is point: 1 where oldbin is in newbin, else 0
    # This involves adding a newbin axis to the old bins, and an oldbin axis
    # to the (end of) the bins
    outbin_shape = (1,) * (len(data.shape) - 1) + (-1, 1)
    if bindatadelta is None:
        # What bin is every data point of the binned data in
        whichbin = numpy.digitize(bindata, bins) - 1
        if clip:
            whichbin[whichbin < 0] = 0
            whichbin[whichbin >= nbins] = (nbins - 1)
        select = numpy.require(
            numpy.expand_dims(whichbin, axis=-2)
            == numpy.reshape(
                numpy.arange(nbins, dtype=int), outbin_shape),
            dtype=int)
    else:
        bindatamin = bindata - bindatadelta
        bindatamax = bindata + bindatadelta
        if clip:
            bindatamin[bindatamin < bins[0]] = bins[0]
            bindatamax[bindatamax > bins[-1]] = bins[-1]
        # The edges of the output and input bins, on a common shape
        # that ends with (outputbins, inputbins)
        bindatamin = numpy.expand_dims(bindatamin, axis=-2)
        bindatamax = numpy.expand_dims(bindatamax, axis=-2)
        binmin = numpy.reshape(bins[:-1], outbin_shape)
        binmax = numpy.reshape(bins[1:], outbin_shape)
        # Overlap of input/output bins, if it exists,  ends when either does
        overlap_top = numpy.minimum(bindatamax, binmax)
        # And it start at the higher of the two
        overlap_bottom = numpy.maximum(bindatamin, binmin)
        # No overlap if top < bottom
        overlap = numpy.maximum(overlap_top - overlap_bottom , 0)
        # Normalize the overlap to fraction of the input bin
        select = overlap / (bindatamax - bindatamin)
    if weights is not None: # Apply weights to the inputs
        select = select * numpy.expand_dims(weights, axis=-2)
    # Add a degenerate dimension to the end of data, so now
    # the data end in dims (oldbin, newbin)
    data = numpy.expand_dims(data, axis=-1)
    idx = numpy.isnan(data)
    counts = numpy.matmul(select, ~idx)[..., 0]
    if bintype == 'count':
        return numpy.rollaxis(counts, axis=-1, start=axis)
    if bintype == 'unc': # Square what we're summing over
        data = data ** 2
        select = select ** 2
    if numpy.may_share_memory(indata, data): # Don't overwrite input
        data = data.copy()
    data[idx] = 0
    data_sum = numpy.matmul(select, data)[..., 0]
    data_sum[counts == 0] = numpy.nan
    counts[counts == 0] = 1 # Suppress error when nan/0
    if bintype == 'unc':
        data_sum = numpy.sqrt(data_sum)
    avg = data_sum / counts
    # Put the binned axis back in place
    avg = numpy.rollaxis(avg, axis=-1, start=axis)
    return avg
