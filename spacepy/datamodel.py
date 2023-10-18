#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The datamodel classes constitute a data model implementation
meant to mirror the functionality of the data model output from pycdf, though
implemented slightly differently.

This contains the following classes:
 * :py:class:`dmarray` - numpy arrays that support .attrs for information about the data
 * :py:class:`SpaceData` - base class that extends dict, to be extended by others

Authors: Steve Morley and Brian Larsen

Additional Contributors: Charles Kiyanda and Miles Engel

Institution: Los Alamos National Laboratory

Contact: smorley@lanl.gov; balarsen@lanl.gov

Copyright 2010-2016 Los Alamos National Security, LLC.


About datamodel
---------------

The SpacePy datamodel module implements classes that are designed to make implementing a standard
data model easy. The concepts are very similar to those used in standards like HDF5, netCDF and
NASA CDF.

The basic container type is analogous to a folder (on a filesystem; HDF5 calls this a
group): Here we implement this as a dictionary-like object, a datamodel.SpaceData object, which
also carries attributes. These attributes can be considered to be global, i.e. relevant for the
entire folder. The next container type is for storing data and is based on a numpy array, this
class is datamodel.dmarray and also carries attributes. The dmarray class is analogous to an
HDF5 dataset.

In fact, HDF5 can be loaded directly into a SpacePy datamodel, carrying across all attributes,
using the function fromHDF5:

>>> import spacepy.datamodel as dm
>>> data = dm.fromHDF5('test.h5')

Functions are also available to directly load data and metadata into a
SpacePy datamodel from NASA CDF as well as JSON-headed ASCII. Writers also
exist to output a SpacePy datamodel directly to HDF5 or JSON-headed ASCII.
See `datamodel.fromCDF`, `datamodel.readJSONheadedASCII`,
`datamodel.toHDF5`, and `datamodel.toJSONheadedASCII` for more details.


Examples
--------

Imagine representing some satellite data within the global attributes might be
the mission name and the instrument PI, the variables might be the
instrument counts [n-dimensional array], timestamps[1-dimensional array and an orbit number [scalar].
Each variable will have one attribute (for this example).

>>> import spacepy.datamodel as dm
>>> mydata = dm.SpaceData(attrs={'MissionName': 'BigSat1'})
>>> mydata['Counts'] = dm.dmarray([[42, 69, 77], [100, 200, 250]], attrs={'Units': 'cnts/s'})
>>> mydata['Epoch'] = dm.dmarray([1, 2, 3], attrs={'units': 'minutes'})
>>> mydata['OrbitNumber'] = dm.dmarray(16, attrs={'StartsFrom': 1})
>>> mydata.attrs['PI'] 'Prof. Big Shot'

This has now populated a structure that can map directly to a NASA CDF, HDF5 or JSON-headed ASCII file.
To visualize our datamodel, we can use tree method (which can be applied to any dictionary-like object
using :func:`~spacepy.toolbox.dictree`).

>>> mydata.tree(attrs=True)

::

    +
    :|____MissionName
    :|____PI
    |____Counts
         :|____Units
    |____Epoch
         :|____units
    |____OrbitNumber
         :|____StartsFrom


Guide for NASA CDF users
------------------------
By definition, a NASA CDF only has a single 'layer'. That is, a CDF contains a series of records
(stored variables of various types) and a set of attributes that are either global or local in
scope. Thus to use SpacePy's datamodel to capture the functionality of CDF the two basic data types
are all that is required, and the main constraint is that datamodel.SpaceData objects cannot be
nested (more on this later, if conversion from a nested datamodel to a flat datamodel is required).


Opening a CDF and working directly with the contents can be easily done using the PyCDF module, however,
if you wish to load the entire contents of a CDF directly into a datamodel (complete with attributes)
the following will make life easier:

>>> import spacepy.datamodel as dm
>>> data = dm.fromCDF('inFile.cdf')


A quick guide to JSON-headed ASCII
----------------------------------
In many cases it is preferred to have a human-readable ASCII file, rather than a binary file like CDF
or HDF5. To make it easier to carry all the same metadata that is available in HDF5 or CDF we have
developed an ASCII data storage format that encodes the metadata using JSON (JavaScript Object Notation).
This notation supports two basic datatypes: key/value collections (like a SpaceData) and ordered lists
(which can represent arrays). JSON is human-readable, but if large arrays are stored in metadata is quickly
becomes difficult to read. For this reason we use JSON to encode the metadata (usually smaller datasets)
and store the data in a standard flat-ASCII format. The metadata is provided as a header that describes
the contents of the file.


To use JSON for storing only metadata associated with the data to be written to an ASCII file a minimal
metadata standard must be implemented. We use the following attribute names: DIMENSION and START_COLUMN.
We also recommend using the NASA ISTP metadata standard to assign attribute names. The biggest limitation
of flat ASCII is that sensibly formatting datasets of more than 2-dimensions (i.e. ranks greater than 2)
is not possible. For this reason if you have datasets of rank 3 or greater then we recommend using HDF5.
If text is absolutely required then it is possible to encode multi-dimensional arrays in the JSON metadata,
but this is not recommended.


This format is best understood by illustration. The following example builds a toy SpacePy datamodel and
writes it to a JSON-headed ASCII file. The contents of the file are then shown.

>>> import spacepy.datamodel as dm
>>> data = dm.SpaceData()
>>> data.attrs['Global'] = 'A global attribute'
>>> data['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
>>> data['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
>>> data['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
>>> dm.toJSONheadedASCII('outFile.txt', data, depend0='Var1', order=['Var1'])
#Note that not all field names are required, those not given will be listed
#alphabetically after those that are specified

The file looks like:

.. code-block:: none

    #{
    #    "MVar": {
    #        "Note": "Metadata",
    #        "VALUES": [7.8]
    #    },
    #    "Global": "A global attribute",
    #    "Var1": {
    #        "Local1": "A local attribute",
    #        "DIMENSION": [1],
    #        "START_COLUMN": 0
    #    },
    #    "Var2": {
    #        "DIMENSION": [2],
    #        "START_COLUMN": 2
    #    }
    #}
    1 8 9
    2 9 1
    3 3 4
    4 8 9
    5 7 8

"""

import collections.abc
import copy
import datetime
import gzip
import io
import itertools
import json
from functools import partial
import os
import re
import warnings

import numpy
# from . import toolbox # handled in functions that use it


__contact__ = 'Steve Morley, smorley@lanl.gov'

str_classes = (str, bytes)

class DMWarning(Warning):
    """
    Warnings class for datamodel, subclassed so it can be set to always
    """
    pass
warnings.simplefilter('always', DMWarning)

class MetaMixin(object):
    """Mixin class that adds a 'meta' attribute that acts like 'attrs'

    Recommendation from the Python Heliophysics community is to allow
    access to metadata via either an ``attrs`` attribute or ``meta``.
    This mixin class supports that recommendation.
    """

    @property
    def meta(self):
        """Equivalent to ``attrs``

        Some APIs use ``attrs`` for metadata; some use ``meta``. This
        is a convenience property to make it easier for those familiar
        with the ``meta`` convention.
        """
        return self.attrs

    @meta.setter
    def meta(self, v):
        """Set meta as with attrs"""
        self.attrs = v

    @meta.deleter
    def meta(self):
        """Remove meta (and thus attrs)

        This isn't a good idea but you can do it with attrs,
        so might as well support it in meta.

        This still leaves the meta property hanging around, but
        cannot delete a property from an instance (have to delete
        from the class).
        """
        del self.attrs


class ISTPArray:
    """Mixin class for array using ISTP metadata.

    Array types like `dmarray` provide all these methods; they assume
    attributes of the array use the
    `ISTP metadata standard <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_
    and are unlikely to give good results if that is not the case.

    Note that some operations that may seem to relate to an array (e.g.
    uncertainties) may require the use of other arrays in a container;
    these are in `ISTPContainer`.

    .. versionadded:: 0.5.0

    .. autosummary::
        ~ISTPArray.plot_as_line
        ~ISTPArray.replace_invalid
    .. automethod:: plot_as_line
    .. automethod:: replace_invalid
    """
    attrs: collections.abc.Mapping

    def replace_invalid(self):
        """Return data from array with invalid values replaced by `~numpy.nan`.

        Makes a copy of the data and, for any values equal to the
        ``FILLLVAL`` attribute, greater than ``VALIDMAX``, or less than
        ``VALIDMIN``, replace with NaN.

        Returns
        -------
        `~numpy.ndarray`
            Transformed data

        See Also
        --------
        .pycdf.istp.nanfill : an in-place variant

        Notes
        -----
        .. versionadded:: 0.5.0

        Comparisons with ``FILLVAL`` are done using `~numpy.isclose` and
        so may replace values that are near, but not identical, to fill.
        """
        data = numpy.array(self)
        idx = numpy.zeros_like(data, dtype=bool)
        if self.attrs.get('FILLVAL') is not None:
            idx |= numpy.isclose(data, self.attrs['FILLVAL'])
        if self.attrs.get('VALIDMIN') is not None:
            idx |= data < self.attrs['VALIDMIN']
        if self.attrs.get('VALIDMAX') is not None:
            idx |= data > self.attrs['VALIDMAX']
        data[idx] = numpy.nan
        return data

    def plot_as_line(self):
        """Determines if this array is better plotted as a lineplot or spectrogram.

        Uses array shape and the ``DISPLAY_TYPE`` attribute to determine
        if should be plotted as a lineplot (potentially stacked) or spectrogram.

        Returns
        -------
        `bool`
            ``True`` if should be a lineplot, ``False`` if should be a
            spectrogram

        Notes
        -----
        .. versionadded:: 0.5.0
        """
        if 'DISPLAY_TYPE' in self.attrs:
            return self.attrs['DISPLAY_TYPE'] == 'time_series'
        dims = len(self.shape)
        if dims == 1:
            return True
        if dims > 2:
            return True
        # Reasonable dividing line is probably 4 stacked line plots
        return self.shape[-1] < 5


class ISTPContainer(collections.abc.Mapping):
    """Mixin class for containers using ISTP metadata.

    Container types like `SpaceData` provide all these methods; they assume
    attributes of the container and the arrays it contains use the
    `ISTP metadata standard <https://spdf.gsfc.nasa.gov/sp_use_of_cdf.html>`_
    and are unlikely to give good results if that is not the case.

    .. versionadded:: 0.5.0

    .. autosummary::
        ~ISTPContainer.lineplot
        ~ISTPContainer.main_vars
        ~ISTPContainer.plot
        ~ISTPContainer.spectrogram
    .. automethod:: lineplot
    .. automethod:: main_vars
    .. automethod:: plot
    .. automethod:: spectrogram
    """
    attrs:  collections.abc.Mapping

    def lineplot(self, vname, target=None):
        """Line plot of a value (array) from this container

        Parameters
        ----------
        vname : `str`
            The key into this container of the value to plot (i.e.,
            the name of the variable).

        target : `matplotlib.axes.Axes` or `matplotlib.figure.Figure`, optional
            Where to draw the plot. Default is to create a new figure with
            a single subplot. If ``Axes``, will draw into that subplot (and
            will not draw a legend or figure title); if ``Figure``, will
            make a single subplot (and not set figure title). Handled by
            `~.plot.utils.set_target`.

        Returns
        -------
        ax : `matplotlib.axes.Axes`
            The subplot on which the variable was plotted

        Notes
        -----
        .. versionadded:: 0.5.0
        """
        import spacepy.plot.utils
        v = self[vname]
        fig, ax = spacepy.plot.utils.set_target(target)
        x = self[v.attrs['DEPEND_0']]
        data = v.replace_invalid()
        labels = None
        if v.attrs.get('LABL_PTR_1'):
            labels = self[v.attrs['LABL_PTR_1']]
        deltas = self.get_deltas(vname)
        plot_kwargs = {}
        if len(data.shape) == 1:
            data = data[..., None]
        if deltas and len(deltas[0].shape) == 1:
            deltas = tuple([d[..., None] for d in deltas])
        for dim in range(data.shape[-1]):
            if labels is not None:
                plot_kwargs['label'] = labels[dim]
            if deltas:
                if len(deltas) == 1:
                    yerr = deltas[0][:, dim]
                else:
                    yerr = numpy.stack((deltas[0][:, dim], deltas[1][:, dim]))
                ax.errorbar(numpy.array(x), data[:, dim], yerr=yerr, **plot_kwargs)
            else:
                ax.plot(numpy.array(x), data[:, dim], **plot_kwargs)
        ylabel = v.attrs.get('LABLAXIS', '')
        if v.attrs.get('UNITS'):
            ylabel = '{}{}({})'.format(
                ylabel, ' ' if ylabel else '', v.attrs['UNITS'])
        if ylabel:
            ax.set_ylabel(ylabel)
        if x.attrs.get('LABLAXIS'):
            ax.set_xlabel(x.attrs['LABLAXIS'])
        if labels is not None and target is not ax:
            ax.legend(loc='best')
        if target is None and v.attrs.get('CATDESC'):
            fig.suptitle(v.attrs['CATDESC'])
        spacepy.plot.utils.applySmartTimeTicks(ax, x)
        return ax

    def main_vars(self):
        """Return names of the 'main' variables in this container.

        These are variables that are likely to be of direct interest, rather
        than dependencies and support data. They are chosen primarily by
        not being dependencies of other variables, but if the ``VAR_TYPE``
        attribute is present it must be ``data``.

        Returns
        -------
        `list` of `str`

        Notes
        -----
        .. versionadded:: 0.5.0
        """
        referenced = set()
        for k, v in self.items():
            referenced.update([v.attrs[a] for a in v.attrs
                               if a.startswith(('DEPEND_', 'LABL_PTR_', 'DELTA_'))])
        main = sorted(set(self).difference(referenced))
        if any(('VAR_TYPE' in v.attrs for v in self.values())):
            main = [m for m in main if self[m].attrs.get('VAR_TYPE', '') == 'data']
        return main

    def plot(self, vnames=None, fig=None):
        """Plot one or more values (arrays) from this container

        Parameters
        ----------
        vnames : `list` of `str`, optional.
            The key into this container of the value(s) to plot (i.e.,
            the name of the variable). If not specified, plots all
            'main' variables which are not dependencies of others;
            see `main_vars`.

        fig : `matplotlib.figure.Figure`, optional
            Where to draw the plot. Default is to create a new figure. If
            given, subplots will be added to this figure (it should start
            empty).

        Returns
        -------
        fig : `matplotlib.figure.Figure`
            The figure on which the variables were plotted

        See Also
        --------
        lineplot : to line plot a single variable
        spectrogram : to make a spectrogram of a single variable

        Notes
        -----
        .. versionadded:: 0.5.0

        Examples
        --------

        >>> import spacepy.datamodel
        # https://rbsp-ect.newmexicoconsortium.org/data_pub/rbspa/ECT/level2/
        >>> data = spacepy.datamodel.fromCDF(
        ...     'rbspa_ect-elec-L2_20140115_v2.1.0.cdf')
        >>> fig = data.plot(['FESA', 'Position'])
        >>> fig.show()  # if needed

        >>> import spacepy.pycdf
        # https://rbsp-ect.newmexicoconsortium.org/data_pub/rbspa/hope/level2/spinaverage/
        >>> with spacepy.pycdf.CDF('rbspa_rel04_ect-hope-sci-L2SA_20140108_v6.1.0.cdf') as f:
        ...     data = f.copy()
        >>> fig = data.plot(['FESA', 'FPSA'])
        >>> fig.show()  # if needed

        >>> import spacepy.pycdf
        # https://spp-isois.sr.unh.edu/data_public/ISOIS/level2/
        >>> with spacepy.pycdf.CDF('psp_isois_l2-summary_20201130_v13.cdf') as f:
        ...     data = f.copy()
        >>> fig = data.plot(['A_H_Rate_TS', 'H_CountRate_ChanP_SP'])
        >>> fig.show()  # if needed
        """
        if fig is None:
            import matplotlib.pyplot
            fig = matplotlib.pyplot.figure()
        if isinstance(vnames, collections.abc.Hashable) and vnames in self:
            vnames = [vnames]
        if vnames is None:
            vnames = self.main_vars()
        n_plots = len(vnames)
        for i, k in enumerate(vnames):
            ax = fig.add_subplot(n_plots, 1, i + 1)
            if self[k].plot_as_line():
                self.lineplot(k, target=ax)
                h, l = ax.get_legend_handles_labels()
                if l:
                    ax.legend(h, l, loc='best')
            else:
                self.spectrogram(k, target=ax)
        return fig

    def spectrogram(self, vname, target=None):
        """Spectrogram plot of a value (array) from this container

        Parameters
        ----------
        vname : `str`
            The key into this container of the value to plot (i.e.,
            the name of the variable).

        target : `matplotlib.axes.Axes` or `matplotlib.figure.Figure`, optional
            Where to draw the plot. Default is to create a new figure with
            a single subplot. If ``Axes``, will draw into that subplot (and
            will not set figure title); if ``Figure``, will make a single
            subplot (and not set figure title). Handled by
            `~.plot.utils.set_target`.

        Returns
        -------
        ax : `matplotlib.axes.Axes`
            The subplot on which the variable was plotted

        Notes
        -----
        .. versionadded:: 0.5.0
        """
        import matplotlib.cm
        import spacepy.plot.utils
        v = self[vname]
        fig, ax = spacepy.plot.utils.set_target(target)
        x = self[v.attrs['DEPEND_0']]
        data = v.replace_invalid()
        x = self[v.attrs['DEPEND_0']]
        y = self[v.attrs['DEPEND_1']]
        zlabel = v.attrs.get('LABLAXIS', '')
        if v.attrs.get('UNITS'):
            zlabel = '{}{}({})'.format(
                zlabel, ' ' if zlabel else '', v.attrs['UNITS'])
        zlabel = zlabel if zlabel else None
        try:  # mpl >=3.7
            cmap = matplotlib.colormaps.get_cmap(None)
        except AttributeError:
            cmap = matplotlib.cm.get_cmap()
        cmap = copy.copy(cmap)
        if cmap(-1.)[:3] == cmap(0.)[:3]:  # Underflow to black if not specified
            cmap.set_under('k')
        # Fill to grey or white
        if cmap(numpy.nan)[:3] == cmap(0.)[:3] and cmap(numpy.nan)[-1] > 0.:
            cmap.set_bad((.5, .5, .5, 0.) if cmap(1.)[:3] == (1., 1., 1.)
                         else (1., 1., 1., 0.))
        ax = spacepy.plot.simpleSpectrogram(numpy.array(x), numpy.array(y), data, cbtitle=zlabel,
                                            ax=ax, zero_valid=True, cmap=cmap)
        ylabel = y.attrs.get('LABLAXIS', '')
        if y.attrs.get('UNITS'):
            ylabel = '{}{}({})'.format(
                ylabel, ' ' if ylabel else '', y.attrs['UNITS'])
        if ylabel:
            ax.set_ylabel(ylabel)
        if x.attrs.get('LABLAXIS'):
            ax.set_xlabel(x.attrs['LABLAXIS'])
        if target is None and v.attrs.get('CATDESC'):
            fig.suptitle(v.attrs['CATDESC'])
        spacepy.plot.utils.applySmartTimeTicks(ax, x)
        return ax

    def get_deltas(self, vname):
        """Return deltas for an array

        Returns ISTP delta values. These may be uncertainties or may
        be e.g. bin widths; interpretation is undefined.

        Invalid values are replaced with `~numpy.nan`.

        Parameters
        ----------
        vname : `str`
            The key into this container of the value to get delta
            (i.e.,  the name of the variable).

        Returns
        -------
        deltas : `tuple` of `~numpy.ndarray`
            Deltas for ``vname``. Empty if no deltas available;
            one-element if symmetric; two-element if not symmetric.

        Notes
        -----
        .. versionadded:: 0.5.0
        """
        v = self[vname]
        asymmetric_msg = 'Only one of DELTA_(MINUS|PLUS)_VAR specified.'
        if 'DELTA_PLUS_VAR' not in v.attrs:
            if 'DELTA_MINUS_VAR' in v.attrs:
                raise ValueError(asymmetric_msg)
            return ()
        elif 'DELTA_MINUS_VAR' not in v.attrs:
            raise ValueError(asymmetric_msg)
        dp = self[v.attrs['DELTA_PLUS_VAR']].replace_invalid()
        if v.attrs['DELTA_PLUS_VAR'] == v.attrs['DELTA_MINUS_VAR']:
            return(dp,)
        return(self[v.attrs['DELTA_MINUS_VAR']].replace_invalid(), dp)


class dmarray(numpy.ndarray, MetaMixin, ISTPArray):
    """
    Container for data within a SpaceData object

    Although the format of attributes is not enforced, using ISTP metadata
    enables the use of methods from `ISTPArray`.

    Raises
    ------
    NameError
        raised is the request name was not added to the allowed attributes list

    Examples
    --------
    >>> import spacepy.datamodel as datamodel
    >>> position = datamodel.dmarray([1,2,3], attrs={'coord_system':'GSM'})
    >>> position
    dmarray([1, 2, 3])
    >>> position.attrs
    {'coord_system': 'GSM'}a

    The dmarray, like a numpy ndarray, is versatile and can store
    any datatype; dmarrays are not just for arrays.

    >>> name = datamodel.dmarray('TestName')
    dmarray('TestName')

    To extract the string (or scalar quantity), use the tolist method

    >>> name.tolist()
    'TestName'

    See methods of `ISTPArray` if attributes are ISTP-compliant.

    .. currentmodule:: spacepy.datamodel
    .. autosummary::
        ~dmarray.addAttribute
    .. automethod:: addAttribute
    """
    Allowed_Attributes = ['attrs']

    def __new__(cls, input_array, attrs=None, dtype=None):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        if not dtype:
            obj = numpy.asarray(input_array).view(cls)
        else:
            obj = numpy.asarray(input_array).view(cls).astype(dtype)
        # add the new attribute to the created instance
        if attrs != None:
            obj.attrs = attrs
        else:
            obj.attrs = {}
        # Finally, return the newly created object:
        return obj

    def __array_finalize__(self, obj):
       # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        for val in self.Allowed_Attributes:
            self.__setattr__(val, copy.deepcopy(getattr(obj, val, {})))

    def __array_wrap__(self, out_arr, context=None):
        #check for zero-dims (numpy bug means subclass behaviour isn't consistent with ndarray
        #this traps most of the bad behaviour ( std() and var() still problems)
        if out_arr.ndim > 0:
            return numpy.ndarray.__array_wrap__(self, out_arr, context)
        return numpy.ndarray.__array_wrap__(self, out_arr, context).tolist()

    def __reduce__(self):
        """This is called when pickling, see:
        http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html
        for this particular example.
        Only the attributes in Allowed_Attributes can exist
        """
        object_state = list(numpy.ndarray.__reduce__(self))
        subclass_state = tuple([tuple([val, self.__getattribute__(val)]) for val in self.Allowed_Attributes])
        object_state[2] = (object_state[2], subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Used for unpickling after __reduce__ the self.attrs is recovered from
        the way it was saved and reset.
        """
        nd_state, own_state = state
        numpy.ndarray.__setstate__(self, nd_state)
        for i, val in enumerate(own_state):
            if not val[0] in self.Allowed_Attributes: # this is attrs
                self.Allowed_Attributes.append(own_state[i][0])
            self.__setattr__(own_state[i][0], own_state[i][1])

    def __setattr__(self, name, value):
        """Make sure that .attrs is the only attribute that we are allowing
        dmarray_ne took 15.324803 s
        dmarray_eq took 15.665865 s
        dmarray_assert took 16.025478 s
        It looks like != is the fastest, but not by much over 10000000 __setattr__
        """
        #meta is special-handled because it should NOT be pickled
        if name in ('Allowed_Attributes', 'meta'):
            pass
        elif not name in self.Allowed_Attributes:
            raise TypeError("Only attribute listed in Allowed_Attributes can be set")
        super(dmarray, self).__setattr__(name, value)

    def addAttribute(self, name, value=None):
        """Method to add an attribute to a dmarray
        equivalent to
        a = datamodel.dmarray([1,2,3])
        a.Allowed_Attributes = a.Allowed_Attributes + ['blabla']
        """
        if name in self.Allowed_Attributes:
            raise NameError('{0} is already an attribute cannot add again'.format(name))
        self.Allowed_Attributes.append(name)
        self.__setattr__(name, value)

    def count(self, srchval):
        """
        Equivalent to count method on list

        """
        mask = self == srchval
        return int(mask.sum())
    def _saveAttrs(self):
        Allowed_Attributes = self.Allowed_Attributes
        backup = []
        for atr in Allowed_Attributes:
            backup.append((atr, dmcopy(self.__getattribute__(atr))))
        return backup

    @classmethod
    def _replaceAttrs(cls, arr, backup):
        for key, val in backup:
            if key != 'attrs':
                try:
                    arr.addAttribute(key)
                except NameError:
                    pass
            arr.__setattr__(key, val)
        return arr

    @classmethod
    def append(cls, one, other):
        """
        append data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.append(one, other))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def vstack(cls, one, other):
        """
        vstack data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.vstack((one, other)))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def hstack(cls, one, other):
        """
        hstack data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.hstack((one, other)))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def dstack(cls, one, other):
        """
        dstack data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.dstack((one, other)))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def concatenate(cls, one, other, axis=0):
        """
        concatenate data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.concatenate((one, other), axis=axis))
        return cls._replaceAttrs(outarr, backup)

def dmfilled(shape, fillval=0, dtype=None, order='C', attrs=None):
    """
    Return a new dmarray of given shape and type, filled with a specified value (default=0).

    See Also
    --------
    numpy.ones

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> dm.dmfilled(5, attrs={'units': 'nT'})
    dmarray([ 0.,  0.,  0.,  0.,  0.])

    >>> dm.dmfilled((5,), fillval=1, dtype=np.int)
    dmarray([1, 1, 1, 1, 1])

    >>> dm.dmfilled((2, 1), fillval=np.nan)
    dmarray([[ nan],
           [ nan]])

    >>> a = dm.dmfilled((2, 1), np.nan, attrs={'units': 'nT'})
    >>> a
    dmarray([[ nan],
           [ nan]])
    >>> a.attrs
    {'units': 'nT'}
    """
    a = dmarray(numpy.empty(shape, dtype, order), attrs=attrs)
    a.fill(fillval)
    return a


class SpaceData(dict, MetaMixin, ISTPContainer):
    """
    Datamodel class extending dict by adding attributes.

    Although the format of attributes is not enforced, using ISTP metadata
    enables the use of methods from `ISTPContainer`.

    .. currentmodule:: spacepy.datamodel
    .. autosummary::
        ~SpaceData.flatten
        ~SpaceData.tree
        ~SpaceData.toCDF
        ~SpaceData.toHDF5
        ~SpaceData.toJSONheadedASCII
    .. automethod:: flatten
    .. automethod:: tree
    .. automethod:: toCDF
    .. automethod:: toHDF5
    .. automethod:: toJSONheadedASCII
    """
    def __getitem__(self, key):
        """
        This allows one to make a SpaceData indexed with an iterable of
        keys to return a new spacedata made of the subset of keys
        """
        try:
            return super(SpaceData, self).__getitem__(key)
        except (KeyError, TypeError):
            if isinstance(key, (tuple, list)):
                # make a new SpaceData from these keys
                out = SpaceData()
                out.attrs = self.attrs
                for k in key:
                    out[k] = self[k]
                return out
            else:
                raise KeyError('{0}'.format(key))

    def __init__(self, *args, **kwargs):
        """
        Base class for "Data Model" representation data

        Abstract method, reimplement

        Attributes
        ----------
        attrs : dict
            dictionary of the attributes of the SpaceData object

        """
        #raise(ValueError("Abstract method called, reimplement __init__"))
        self.attrs = {}
        if 'attrs' in kwargs:
            if hasattr(kwargs['attrs'], '__getitem__'):
                self.attrs = kwargs['attrs']
            del kwargs['attrs']

        super(SpaceData, self).__init__(*args, **kwargs)
        self.toCDF = partial(toCDF, SDobject=self)
        self.toCDF.__doc__ = toCDF.__doc__
        self.toHDF5 = partial(toHDF5, SDobject=self)
        self.toHDF5.__doc__ = toHDF5.__doc__
        self.toJSONheadedASCII = partial(toJSONheadedASCII, insd=self)
        self.toJSONheadedASCII.__doc__ = toJSONheadedASCII.__doc__

## To enable string output of repr, instead of just printing, uncomment his block
#    def __repr__(self):
#        #redirect stdout to StringIO
#        import io, sys
#        dum = io.StringIO()
#        sys_stdout_save = sys.stdout
#        sys.stdout = dum
#        self.tree(verbose=True)
#        sys.stdout = sys_stdout_save
#        dum.seek(0)
#        return ''.join(dum.readlines())

    def tree(self, **kwargs):
        '''Print the contents of the SpaceData object in a visual tree

        Other Parameters
        ----------------
        verbose : bool, default False
            print more info
        spaces : str (optional)
            string will added for every line
        levels : int (optional)
            number of levels to recurse through (True, the default,  means all)
        attrs : bool, default False
            display information for attributes
        print_out : bool, default True

                .. versionadded:: 0.5.0

            Print output (original behavior); if ``False``, return the output.

        Examples
        --------
        >>> import spacepy.datamodel as dm
        >>> import spacepy.toolbox as tb
        >>> a = dm.SpaceData()
        >>> a['1'] = dm.SpaceData(dog = 5)
        >>> a['4'] = dm.SpaceData(cat = 'kitty')
        >>> a['5'] = 4
        >>> a.tree()
        +
        |____1
             |____dog
        |____4
             |____cat
        |____5

        See Also
        --------
        toolbox.dictree
        '''
        from . import toolbox
        return toolbox.dictree(self, **kwargs)

    def flatten(self):
        '''
        Method to collapse datamodel to one level deep

        Examples
        --------
        >>> import spacepy.datamodel as dm
        >>> import spacepy.toolbox as tb
        >>> a = dm.SpaceData()
        >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        >>> a['4'] = dm.SpaceData(cat = 'kitty')
        >>> a['5'] = 4
        >>> a.tree()
        +
        |____1
             |____dog
             |____pig
                  |____fish
                       |____a
                       |____b
        |____4
             |____cat
        |____5

        >>> b = dm.flatten(a)
        >>> b.tree()
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5

        >>> a.flatten()
        >>> a.tree()
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5

        '''
        flatobj = flatten(self)
        remkeys = [key for key in self]
        for key in remkeys:
            del self[key]
        for key in flatobj:
            self[key] = copy.copy(flatobj[key])

    # Stubs of partialed-in functions for docs; actual versions populated
    # when class instantiated

    def toCDF(fname, **kwargs):
        """Create CDF file from this SpaceData.

        See `toCDF`; this object is provided for ``SDobject``."""

    def toHDF5(fname, **kwargs):
        """Create HDF5 file from this SpaceData.

        See `toHDF5`; this object is provided for ``SDObject``. """

    def toJSONheadedASCII(fname, **kwargs):
        """Create JSON-headed ASCII file from this SpaceData.

        See `toJSONheadedASCII`; this object is provided for ``insd``."""


def convertKeysToStr(SDobject):
    if isinstance(SDobject, SpaceData):
        newSDobject = SpaceData()
        newSDobject.attrs = SDobject.attrs
    else:
        newSDobject = {}
    for key in SDobject:
        if not isinstance(key, str_classes):
            if isinstance(SDobject[key], dict):
                newSDobject[str(key)] = convertKeysToStr(SDobject[key])
            else:
                newSDobject[str(key)] = SDobject[key]
        else:
            if isinstance(SDobject[key], dict):
                newSDobject[key] = convertKeysToStr(SDobject[key])
            else:
                newSDobject[key] = SDobject[key]

    return newSDobject


def flatten(dobj):
    '''Collapse datamodel to one level deep

    Examples
    --------

    >>> import spacepy.datamodel as dm
    >>> import spacepy.toolbox as tb
    >>> a = dm.SpaceData()
    >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
    >>> a['4'] = dm.SpaceData(cat = 'kitty')
    >>> a['5'] = 4
    >>> a.tree()
    +
    |____1
         |____dog
         |____pig
              |____fish
                   |____a
                   |____b
    |____4
         |____cat
    |____5

    >>> b = dm.flatten(a)
    >>> b.tree()
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5

    >>> a.flatten()
    >>> a.tree()
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5


    See Also
    --------
    unflatten
    SpaceData.flatten

    '''

    try:
        addme = dobj.__class__()
    except TypeError:
        addme = SpaceData()
    remlist = []
    for key in dobj: #iterate over keys in SpaceData
        if isinstance(dobj[key], dict):
            remlist.append(key)
            newname = str(key) + '<--'
            for levkey in dobj[key]:
                if hasattr(dobj[key][levkey], 'keys'):
                    retdict = flatten(dobj[key][levkey])
                    for key2 in retdict:
                        addme[newname+levkey+'<--'+key2] = retdict[key2]
                else:
                    addme[newname+levkey] = copy.copy(dobj[key][levkey])
        else:
            addme[key] = copy.copy(dobj[key])
    return addme

def unflatten(dobj, marker='<--'):
    '''Collapse datamodel to one level deep

    Examples
    --------

    >>> import spacepy.datamodel as dm
    >>> import spacepy.toolbox as tb
    >>> a = dm.SpaceData()
    >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
    >>> a['4'] = dm.SpaceData(cat = 'kitty')
    >>> a['5'] = 4
    >>> a.tree()
    +
    |____1
         |____dog
         |____pig
              |____fish
                   |____a
                   |____b
    |____4
         |____cat
    |____5

    >>> b = dm.flatten(a)
    >>> b.tree()
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5

    >>> c = dm.unflatten(b)
    >>> c.tree()
    +
    |____1
         |____dog
         |____pig
              |____fish
                   |____a
                   |____b
    |____4
         |____cat
    |____5


    '''
    #set up a new object for return
    try:
        addme = dobj.__class__()
    except TypeError:
        addme = SpaceData()
    #the input is assumed to be single level (i.e. it is flat)

    #find all keys that have at least one marker,
    #then unpack. Recurse over these until no more markers are found.
    keydict = {}
    for key in dobj:
        if isinstance(dobj[key], dict):
            raise TypeError('Flat datamodel should not contain dict-likes')
        try:
            if marker in key:
                #get 'group'
                group = key.split(marker)[0]
                if not group in keydict:
                    keydict[group] = {key: ''}
                else:
                    keydict[group][key] = ''
            else: #not nested, just copy key
                addme[key] = dmcopy(dobj[key])
        except:
            addme[key] = dmcopy(dobj[key])
    #now we have all the groups at this level
    #move members of groups into new SpaceDatas
    for grp in keydict:
        addme[grp] = SpaceData()
        for key in keydict[grp]:
            newkey = marker.join(key.split(marker)[1:])
            addme[grp][newkey] = dmcopy(dobj[key])
        # recurse to make sure everything inside is unpacked
        addme[grp] = unflatten(addme[grp], marker=marker)
    return addme


def fromCDF(fname):
    '''
    Create a SpacePy datamodel representation of a NASA CDF file

    Parameters
    ----------
    file : string
        the name of the cdf file to be loaded into a datamodel

    Returns
    -------
    out : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> data = dm.fromCDF('test.cdf')

    See Also
    --------
    .pycdf.CDF.copy
    .pycdf.istp.VarBundle
    '''
    #TODO: add unflatten keyword and restore flattened variables
    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and SpacePy's pyCDF")

    with pycdf.CDF(fname) as cdfdata:
        return cdfdata.copy()

def toCDF(fname, SDobject, skeleton='', flatten=False, overwrite=False,
          autoNRV=False, backward=None, TT2000=None, verbose=False):
    '''
    Create a CDF file from a SpacePy datamodel representation

    Parameters
    ----------
    fname : str
        Filename to write to

    SDobject : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    skeleton : str (optional)
        create new CDF from a skeleton file (default '')

    flatten : bool (optional)
        flatten incoming datamodel - if SpaceData objects are nested (default False)

    overwrite : bool (optional)
        allow overwrite of an existing target file (default False)

    autoNRV : bool (optional)
        attempt automatic identification of non-record varying entries in CDF

    backward : bool (optional)
        ``False`` to create CDF in backward-compatible format; ``True``
         to force v3+ compatibility only. (Default: do not change current
         state, see :meth:`~.pycdf.Library.set_backward`).

         .. versionchanged:: 0.5.0
            Now supports specifying backward compatible or no change;
            previous versions always wrote v3+ CDFs (even if ``False``).

    TT2000 : bool (optional)
        Specify type for variables with names beginning 'Epoch'. Default
        CDF_EPOCH for backward-compatible CDF (``backward`` True) and
        CDF_TT20000 otherwise (``backward`` False or unspecified).

        .. versionchanged:: 0.5.0
           Current handling introduced.

        .. versionchanged:: 0.3.0
           Always write TT2000 variables (due to change in :mod:`~.pycdf`).

    verbose : bool (optional)
        verbosity flag

    Returns
    -------
    None

    Notes
    -----

    .. versionchanged:: 0.5.0
       Invalid keyword arguments now raise :exc:`TypeError` rather than being ignored.
    '''
    if flatten:
        SDobject = SDobject.flatten()
    if overwrite:
        raise NotImplementedError('Overwriting CDFs is not currently enabled '
                                  '- please remove the file manually')
    if TT2000 and backward:
        raise ValueError('Cannot use TT2000 in backward-compatible CDF.')
    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and"
                          " SpacePy's pyCDF")
    if backward is None:
        former_backward = None
    else:
        former_backward = pycdf.lib.set_backward(backward)
    force_epoch = not backward and TT2000 is False  # backward defaults falsey
    with pycdf.CDF(fname, skeleton) as outdata:
        if hasattr(SDobject, 'attrs'):
            for akey in SDobject.attrs:
                outdata.attrs[akey] = dmcopy(SDobject.attrs[akey])
        varLengths = [len(SDobject[var]) for var in SDobject]
        modeLength = next(itertools.groupby((reversed(sorted(varLengths)))))[0]
        for key, val in SDobject.items():
            if isinstance(val, dict):
                raise TypeError('This data structure appears to be nested,'
                                ' please try spacepy.datamodel.flatten')
            if not skeleton:
                if not val.shape:
                    shape_tup = -1
                else:
                    shape_tup = val.shape
                if 'Epoch' not in SDobject:
                    NRVtest = modeLength
                else:
                    NRVtest = len(SDobject['Epoch'])
                if shape_tup[0] != NRVtest: #naive check for 'should-be' NRV
                    try:
                        v = outdata.new(key, val[...], recVary=False)
                        if verbose:
                            print('{0} is being made NRV'.format(key))
                        v.attrs = dmcopy(val.attrs)
                    except ValueError:
                        v = outdata.new(key, val.tolist, recVary=False)
                        v.attrs = dmcopy(val.attrs)
                if force_epoch and 'Epoch' in key:
                    outdata.new(key, val[...], type=pycdf.const.CDF_EPOCH)
                else:
                    try:
                        outdata[key] = val
                    except ValueError:
                        try:
                            outdata[key] = dmarray(
                                [val.tolist()], attrs=dmcopy(val.attrs)).squeeze()
                        except UnicodeEncodeError:
                            tmpAttrs = dmcopy(val.attrs)
                            for akey in tmpAttrs:
                                try: #strings
                                    if hasattr(tmpAttrs[akey], 'encode'):
                                        tmpAttrs[akey] = tmpAttrs[akey].encode('utf-8')
                                    else:
                                        tmpAttrs[akey] = tmpAttrs[akey]
                                except AttributeError: #probably a list of strings
                                    for id, el in enumerate(tmpAttrs[akey]):
                                        tmpAttrs[akey][id] = el.encode('utf-8')

            else:
                outdata[key][...] = val[...]
                for akey in outdata[key].attrs:
                    try:
                        outdata[key].attrs[akey] = dmcopy(val.attrs[akey])
                    except ValueError:
                        outdata[key][...] = dmarray([val.tolist()], attrs=dmcopy(val.attrs))
                    except KeyError:
                        pass
    if former_backward is not None:
        pycdf.lib.set_backward(former_backward)


def fromHDF5(fname, **kwargs):
    '''
    Create a SpacePy datamodel representation of an HDF5 file or netCDF4 file which is HDF5 compliant

    Parameters
    ----------
    file : string
        the name of the HDF5/netCDF4 file to be loaded into a datamodel

    Returns
    -------
    out : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> data = dm.fromHDF5('test.hdf')

    Notes
    -----
    Zero-sized datasets will break in h5py. This is kluged by returning a
    dmarray containing a None.

    This function is expected to work with any HDF5-compliant files, including
    netCDF4 (not netCDF3) and MatLab save files from v7.3 or later, but some
    datatypes are not supported, e.g., non-string vlen datatypes, and will
    raise a warning.
    '''
    def hdfcarryattrs(SDobject, hfile, path):
        if hasattr(hfile[path], 'attrs'):
            #for key, value in hfile[path].attrs.iteritems():
            for key in hfile[path].attrs:
                try:
                    value = hfile[path].attrs[key]
                except TypeError:
                    warnings.warn('Unsupported datatype in dataset {}.attrs[{}]'.format(path, key))
                    continue
                try:
                    SDobject.attrs[key] = value
                except:
                    warnings.warn('The following key:value pair is not permitted\n' +
                                  'key = {0} ({1})\n'.format(key, type(key)) +
                                  'value = {0} ({1})'.format(value, type(value)), DMWarning)

    try:
        import h5py
    except ImportError:
        raise ImportError('HDF5 converter requires h5py')

    if isinstance(fname, str_classes):
        hfile = h5py.File(fname, mode='r')
    else:
        hfile = fname
        #should test here for HDF file object
    path = kwargs.get('path', '/')

    SDobject = SpaceData()
    allowed_elems = [h5py.Group, h5py.Dataset]
    ##carry over the attributes
    hdfcarryattrs(SDobject, hfile, path)
    ##carry over the groups and datasets
    for key, value in hfile[path].items():
        if isinstance(value, allowed_elems[0]):  # if a group
            SDobject[key] = fromHDF5(hfile, path=path+'/'+key)
        elif isinstance(value, allowed_elems[1]):  # if a dataset
            isuni = h5py.check_vlen_dtype(value.dtype) is str
            try:
                if isuni:
                    if hasattr(value, 'asstr'):  # h5py 3+
                        value = value.asstr()
                    value = numpy.require(value[...], dtype=str)
                SDobject[key] = dmarray(value)
            except (TypeError, ZeroDivisionError): #ZeroDivisionError catches zero-sized DataSets
                SDobject[key] = dmarray(None)
            hdfcarryattrs(SDobject[key], hfile, path+'/'+key)
    if path == '/':
        hfile.close()
    return SDobject

def toHDF5(fname, SDobject, **kwargs):
    '''
    Create an HDF5 file from a SpacePy datamodel representation

    Parameters
    ----------
    fname : str
        Filename to write to

    SDobject : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    overwrite : bool (optional)
        allow overwrite of an existing target file (default True)
    mode : str (optional)
        HDF5 file open mode (a, w, r) (default 'a')
    compression : str (optional)
        compress all non-scalar variables using this method (default None)
        (gzip, shuffle, fletcher32, szip, lzf)

        .. versionchanged:: 0.4.0
            No longer compresses scalars (which usually fails).

    compression_opts : str (optional)
        options to the compression, see h5py documentation for more details

    Returns
    -------
    None

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> a = dm.SpaceData()
    >>> a['data'] = dm.dmarray(range(100000), dtype=float)
    >>> dm.toHDF5('test_gzip.h5', a, overwrite=True, compression='gzip')
    >>> dm.toHDF5('test.h5', a, overwrite=True)
    >>> # test_gzip.h5 was 118k, test.h5 was 785k
    '''
    def SDcarryattrs(SDobject, hfile, path, allowed_attrs):
        if hasattr(SDobject, 'attrs'):
            for key, value in SDobject.attrs.items():
                dumval, dumkey = copy.copy(value), copy.copy(key)
                if isinstance(value, allowed_attrs):
                    #test for datetimes in iterables
                    if hasattr(value, '__iter__') and not isinstance(value, str_classes):
                        dumval = [b.isoformat() if isinstance(b, datetime.datetime) else b for b in value]
                    truth = False
                    try:
                        if value.nbytes: truth = True #empty arrays of any dimension are nbytes=0
                    except AttributeError: #not an array
                        if value or value == 0: truth = True

                    if truth:
                        uni = False #No special unicode handling
                        dumval = numpy.asanyarray(dumval)
                        if dumval.size and dumval.dtype.kind == 'U':
                            uni = True #Unicode list, special handling
                        try:
                            if uni:
                                #Tell hdf5 this is unicode. Numpy is UCS-4, HDF5 is UTF-8
                                hfile[path].attrs.create(
                                    dumkey, dumval, dtype=h5py.string_dtype(encoding='utf-8'))
                            else:
                                hfile[path].attrs[dumkey] = dumval
                        except TypeError:
                            hfile[path].attrs[dumkey] = str(dumval)
                            warnings.warn(
                                'The following value is not permitted\n' +
                                'key, value, type = {0}, {1}, {2})\n'.format(
                                    key, value, type(value)) +
                                'value has been converted to a string for output',
                                DMWarning)
                    else:
                        hfile[path].attrs[dumkey] = ''
                elif isinstance(value, datetime.datetime):
                    dumval = value.isoformat()
                    hfile[path].attrs[dumkey] = dumval
                else:
                    #TODO: add support for arrays(?) in attrs (convert to isoformat)
                    warnings.warn('The following key:value pair is not permitted\n' +
                                  'key = {0} ({1})\n'.format(key, type(key)) +
                                  'value type {0} is not in the allowed attribute list'.format(type(value)),
                                  DMWarning)

    try:
        import h5py
    except ImportError:
        raise ImportError('h5py is required to use HDF5 files')

    if not isinstance(SDobject, SpaceData):
        raise ValueError("Input data is not of type SpaceData, check usage:"
                         " toHDF5(fname, datamodel)")
    #mash these into a defaults dict...
    wr_mo = kwargs.get('mode', 'a')
    h5_compr_type = kwargs.get('compression', None)
    if h5_compr_type not in ['gzip', 'szip', 'lzf', 'shuffle', 'fletcher32', None]:
        raise NotImplementedError('Specified compression type not supported')
    h5_compr_opts = None if h5_compr_type == 'lzf'\
                    else kwargs.get('compression_opts', None)

    if 'overwrite' not in kwargs:
        kwargs['overwrite'] = True
    if isinstance(fname, str_classes):
        if os.path.isfile(fname):
            if kwargs['overwrite']:
                os.remove(fname)
            else:
                raise IOError('Cannot write HDF5, file exists (see overwrite) "{!s}"'.format(fname))
        hfile = h5py.File(fname, mode=wr_mo)
        must_close = True
    else:
        hfile = fname
        #should test here for HDF file object
        must_close = False
    path = kwargs.get('path', '/')

    allowed_attrs = [int, float, bytes, str, numpy.ndarray, list, tuple, numpy.string_]
    for v in numpy.typecodes['AllInteger']:
        allowed_attrs.append(numpy.sctypeDict[v])
    for v in numpy.typecodes['AllFloat']:
        allowed_attrs.append(numpy.sctypeDict[v])
    allowed_attrs = tuple(allowed_attrs)

    allowed_elems = (SpaceData, dmarray)

    #first convert non-string keys to str
    SDobject = convertKeysToStr(SDobject)
    SDcarryattrs(SDobject, hfile, path, allowed_attrs)

    try:
        for key, value in SDobject.items():
            if isinstance(value, allowed_elems[0]):
                hfile[path].create_group(key)
                toHDF5(
                    hfile, value, path=path + '/' + key,
                    compression=h5_compr_type, compression_opts=h5_compr_opts)
            elif isinstance(value, allowed_elems[1]):
                comptype, compopts = (None, None) if value.shape == ()\
                                     else (h5_compr_type, h5_compr_opts)
                try:
                    hfile[path].create_dataset(key, data=value,
                                               compression=comptype, compression_opts=compopts)
                except:
                    dumval = numpy.asanyarray(value.copy())
                    dtype = None
                    if dumval.dtype.kind == 'U':
                        dumval = numpy.char.encode(dumval, 'utf-8')
                        dtype = h5py.string_dtype(encoding='utf-8')
                    elif isinstance(value[0], datetime.datetime):
                        for i, val in enumerate(value):
                            dumval[i] = val.isoformat()
                        dumval = dumval.astype('|S35')
                    else:
                        dumval = dumval.atsype('|S35')
                    hfile[path].create_dataset(key, data=dumval, compression=comptype,
                                               compression_opts=compopts, dtype=dtype)
                    #else:
                    #    hfile[path].create_dataset(key, data=value.astype(float))
                SDcarryattrs(SDobject[key], hfile, path+'/'+key, allowed_attrs)
            else:
                warnings.warn('The following data is not being written as is not of an allowed type\n' +
                              'key = {0} ({1})\n'.format(key, type(key)) +
                              'value type {} is not in the allowed data type list'.format(
                                  type(value)), DMWarning)
    finally:
        if must_close:
            hfile.close()


def fromNC3(fname):
    try:
        from scipy.io import netcdf as nc
    except ImportError:
        raise ImportError('SciPy is required to import netcdf3')

    ncfile = nc.netcdf_file(fname, mode='r', mmap=False)

    SDobject = SpaceData(attrs=dmcopy(ncfile._attributes))

    ##carry over the groups and datasets
    for key, value in ncfile.variables.items():
        #try:
            SDobject[key] = dmarray(dmcopy(value.data), attrs=dmcopy(value._attributes))
        #except (TypeError, ZeroDivisionError): #ZeroDivisionError catches zero-sized DataSets
        #    SDobject[key] = dmarray(None)
    ncfile.close()
    return SDobject



def toHTML(fname, SDobject, attrs=(),
           varLinks=False, echo=False, tableTag='<table border="1">'):
    """
    Create an HTML dump of the structure of a spacedata

    Parameters
    ----------
    fname : str
        Filename to write to

    SDobject : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    overwrite : bool (optional)
        allow overwrite of an existing target file (default True)
    mode : str (optional)
        HDF5 file open mode (a, w, r) (default 'a')
    echo : bool
        echo the html to the screen
    varLinks : bool
        make the variable name a link to a stub page

    """
    output = io.StringIO() # put the output into a StringIO
    keys = sorted(SDobject.keys())

    output.write(tableTag)
    output.write('\n')
    output.write('<tr><th>{0}</th>'.format('Variable'))
    for attr in attrs:
        output.write('<th>{0}</th>'.format(attr))
    output.write('</tr>')

    for ii, key in enumerate(keys):
        if ii % 2 == 0:
            output.write('<tr>')
        else:
            output.write('<tr class="alt">')
        output.write('<td>')
        if varLinks:
            output.write('<a href="{0}.html">'.format(key))
        output.write('{0}'.format(key))
        if varLinks:
            output.write('</a>')
        output.write('</td>')

        for attr in attrs:
            try:
                if not isinstance(SDobject[key].attrs[attr], str):
                    tmp = str(SDobject[key].attrs[attr])
                    output.write('<td>{0}</td>'.format(_idl2html(tmp)))
                else:
                    output.write('<td>{0}</td>'.format(_idl2html(SDobject[key].attrs[attr])))
            except KeyError:
                output.write('<td></td>')
        output.write('</tr>\n')
    output.write('</table>\n')
    with open(fname, 'w') as fp:
        fp.write(output.getvalue())
    if echo:
        print(output.getvalue())
    output.close()

def _idl2html(idl):
    """
    given an idl format string for text change it to html

    Parameters
    ==========
    idl : str
        idl formated string

    Returns
    =======
    out : str
        html formatted string
    """
    html = idl
    conv = {'!!': '!',
            '!E': '<sup>',
            '!I': '<sub>'}
    while True: # hate it but for now
        ind = html.find('!')
        if ind == -1:
            break
        code = html[ind:ind+2]
        html = html.replace(code, conv[code])
        if code == '!I':
            if '!N' in html:
                html = html.replace('!N', '</sub>', 1) # just replace 1
            else:
                html = html + '</sub>'
        elif code == '!E':
            if '!N' in html:
                html = html.replace('!N', '</sup>', 1) # just replace 1
            else:
                html = html + '</sup>'
    return html

def readJSONMetadata(fname, **kwargs):
    '''Read JSON metadata from an ASCII data file

    Parameters
    ----------
    fname : str
        Filename to read metadata from

        .. versionchanged:: 0.5.0
                Filename can now be a .gz to indicate the file is gzipped

    Other Parameters
    ----------------
    verbose : bool (optional)
        set verbose output so metadata tree prints on read (default False)

    Returns
    -------
    mdata: spacepy.datamodel.SpaceData
        SpaceData with the metadata from the file
    '''
    if hasattr(fname, 'read'):
        lines = fname.read()
    else:
        if fname.endswith('.gz'):
            with gzip.open(filename=fname, mode='rt', encoding='latin=1') as gzh:
                lines = gzh.read()
        else:
            with open(fname, 'r') as f:
                lines = f.read()

    # isolate header
    p_srch = re.compile(r"^#(.*)$", re.M)
    hreg = re.findall(p_srch, lines)
    header = "".join(hreg)

    # isolate JSON field
    srch = re.search(r'\{\s*(.*)\s*\}', header)
    if isinstance(srch, type(None)):
        raise IOError(
            'The input file has no valid JSON header. Must be valid JSON bounded by braces "{ }".')
    js = srch.group(1)
    inx = js.rfind('end JSON')

    if inx == -1:
        js = ' '.join(('{', js, '}'))
        mdatadict = json.loads(js)
    else:
        js = ' '.join(('{', js[:inx]))
        mdatadict = json.loads(js)

    mdata = SpaceData()
    for key in mdatadict:
        if not hasattr(mdatadict[key], 'keys'): # not dict-like, must be global attrs
            mdata.attrs[key] = mdatadict[key]
        elif 'START_COLUMN' in mdatadict[key]: # is a variable
            mdata[key] = SpaceData(attrs=mdatadict[key])
        elif 'VALUES' in mdatadict[key]: # is global metadata
            dum = mdatadict[key].pop('VALUES')
            mdata[key] = dmarray(dum, attrs=mdatadict[key])
        else: # don't know how to deal with this, store as global attrs
            mdata.attrs[key] = mdatadict[key]

    if 'verbose' in kwargs:
        if kwargs['verbose']:
            mdata.tree(verbose=True, attrs=True)
    return mdata

def readJSONheadedASCII(fname, mdata=None, comment='#', convert=False, restrict=None):
    """read JSON-headed ASCII data files into a SpacePy datamodel

    Parameters
    ----------
    fname : str or list
        Filename(s) to read data from

            .. versionchanged:: 0.5.0
                Filename can now be a .gz to indicate the file is gzipped

    Other Parameters
    ----------------
    mdata : spacepy.datamodel.SpaceData (optional)
        supply metadata object, otherwise is read from fname (default None)
    comment: str (optional)
        comment string in file to be read; lines starting with comment are
        ignored (default '#')
    convert: bool or dict-like (optional)
        If True, uses common names to try conversion from string. If a dict-
        like then uses the functions specified as the dict values to convert
        each element of 'key' to a non-string
    restrict: list of strings (optional)
        If present, restrict the variables stored to only those on this list

    Returns
    -------
    mdata: spacepy.datamodel.SpaceData
        SpaceData with the data and metadata from the file
    """
    import dateutil.parser as dup
    filelike = False
    if isinstance(fname, str_classes):
        fname = [fname]
    elif hasattr(fname, 'readlines'):
        fname = [fname]
        filelike = True
    if not mdata:
        mdata = readJSONMetadata(fname[0])
    if restrict:
        delkeys = [kk for kk in mdata.keys() if kk not in restrict]
        for val in delkeys:
            del mdata[val] #remove undesired keys
    mdata_copy = dmcopy(mdata)
    def innerloop(fh, mdata, mdata_copy):
        line = fh.readline()
        line = line.decode('latin1')
        while (line and line[0] == comment):
            line = fh.readline()
            line = line.decode('latin1')
        fh.seek(-len(line), os.SEEK_CUR) # fixes the missing first data bug
        alldata = fh.readlines()
        if not alldata:
            return mdata
        alldata = [d.decode('latin1') for d in alldata]
        ncols = len(alldata[0].rstrip().split())
        # fixes None in the data from empty lines at the end
        for row in range(len(alldata)): # reverse order
            if not alldata[-1].rstrip(): # blank line (or al white space)
                alldata.pop(-1)
            else:
                break
        nrows = len(alldata)
        data = numpy.empty((nrows, ncols), dtype=object)
        for ridx, line in enumerate(alldata):
            for cidx, el in enumerate(line.rstrip().split()):
                data[ridx, cidx] = el
        for key in mdata_copy.keys():
            if 'START_COLUMN' in mdata_copy[key].attrs:
                st = mdata_copy[key].attrs['START_COLUMN']
                if 'DIMENSION' in mdata_copy[key].attrs:
                    varDims = numpy.array(mdata_copy[key].attrs['DIMENSION'])
                    if not varDims.shape:
                        varDims = numpy.array([varDims])
                    singleDim = True
                    if len(varDims) > 1 or varDims[0] > 1:
                        singleDim = False
                if ('DIMENSION' in mdata_copy[key].attrs) and not singleDim:
                    en = int(mdata_copy[key].attrs['DIMENSION'][0]) + int(st)
                    try:
                        assert mdata[key] == {}
                        mdata[key] = data[:, int(st):int(en)]
                    except (AssertionError, ValueError):
                        mdata[key] = numpy.vstack((mdata[key], data[:, int(st):int(en)]))
                else:
                    try:
                        assert mdata[key] == {}
                        mdata[key] = data[:, int(st)]
                    except (AssertionError, ValueError):
                        mdata[key] = numpy.hstack((mdata[key], data[:, int(st)]))
        return mdata
    for fn in fname:
        if not filelike:
            if fn.endswith('.gz'):
                with gzip.open(filename=fn) as gzh:
                    mdata = innerloop(gzh, mdata, mdata_copy)
            else:
                with open(fn, 'rb') as fh: # fixes windows bug with seek()
                    mdata = innerloop(fh, mdata, mdata_copy)
        else:
            mdata = innerloop(fh, mdata, mdata_copy)
    #now add the attributres to the variables
    keys = list(mdata_copy.keys())
    for key in keys:
        if isinstance(mdata[key], SpaceData):
            mdata[key] = dmarray(None, attrs=mdata_copy[key].attrs)
        else:
            mdata[key] = dmarray(mdata[key], attrs=mdata_copy[key].attrs)

    if convert:
        if isinstance(convert, dict):
            conversions = convert
        else:
            conversions = {'DateTime': lambda x: dup.parse(x, ignoretz=True),
                           'ExtModel': str}
        for conkey in conversions:
            try:
                name = keys.pop(keys.index(conkey)) #remove from keylist
            except ValueError:
                warnings.warn('Key {} for conversion not found in file'.format(conkey), UserWarning)
                continue
            for i, element in numpy.ndenumerate(mdata[name]):
                mdata[name][i] = conversions[name](element)

    for remkey in keys:
        try:
            mdata[remkey] = numpy.asanyarray(mdata[remkey], dtype=float)
        except ValueError:
            pass #this will skip any unspecified string fields
    return mdata

def writeJSONMetadata(fname, insd, depend0=None, order=None, verbose=False, returnString=False):
    '''Scrape metadata from SpaceData object and make a JSON header

    Parameters
    ----------
    fname : str
        Filename to write to (can also use a file-like object)
        None can be given in conjunction with the returnString keyword to skip writing output

    insd : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    depend0 : str (optional)
        variable name to use to indicate parameter on which other data depend (e.g. Time)
    order : list (optional)
        list of key names in order of start column in output JSON file
    verbose: bool (optional)
        verbose output
    returnString: bool (optional)
        return JSON header as string instead of returning None

    Returns
    -------
    None  (unless returnString keyword is True)
    '''
    js_out = {}

    def stripNL(text):
        out = text.group().replace('\n', '').replace('  ', '')
        return out

    #if required, identify depend0 for deciding what's data/metadata
    if depend0 is None:
        #search for DEPEND_0 in metadata
        for key in insd:
            if not hasattr(insd[key], 'attrs'):
                insd[key] = dmarray(insd[key])
            if 'DEPEND_0' in insd[key].attrs:
                depend0 = insd[key].attrs['DEPEND_0']
                if not isinstance(depend0, str_classes):
                    #assume it's a singleton list
                    depend0 = depend0[0]
                if not isinstance(depend0, str_classes):
                    depend0 = None #Failed to get a depend0
                else:
                    break #we're done here
        if depend0 is None:
            #fall back to most common var length
            tmp, keylist = [], list(insd.keys())
            for key in keylist:
                tmp.append(len(insd[key]))
            depend0 = keylist[tmp.index(numpy.bincount(tmp).argmax())]
            #TODO Set using Time, or Epoch, or similar...
    elif not depend0 in insd:
        raise KeyError('Invalid key supplied for ordering metadata on write')
    datalen = len(insd[depend0])

    #start with global attrs
    #TODO: check for datetime objs in attributes
    if insd.attrs:
        glattr = _dateToISO(insd.attrs)
        for key in glattr:
            js_out[key] = dmcopy(glattr[key])
            #TODO Mark these as global somehow (by omission of some metadata?)
            try:
                js_out[key] = js_out[key].tolist()
            except:
                pass
    #collect keys and put in order for output
    #TODO first check for extant START_COLUMN
    #then check dimensionality so that start column and dims can be added, if not present
    if hasattr(order, '__iter__'):
        keylist = order
        #now make sure that all missing keys are added to end
        for key in sorted(insd.keys()):
            if key not in order:
                keylist.append(key)
    else:
        ##TODO do we want to have DEPEND0 first in order by default?
        keylist = sorted(insd.keys())

    idx = 0
    for key in keylist:
        js_out[key] = dmcopy(_dateToISO(insd[key].attrs))
        if len(insd[key]) == datalen: #is data
            if verbose:
                print('data: {0}'.format(key))
            try:
                js_out[key]['DIMENSION'] = list(insd[key].shape[1:])
                if not js_out[key]['DIMENSION']:
                    js_out[key]['DIMENSION'] = [1]
                js_out[key]['START_COLUMN'] = idx
                dims = js_out[key]['DIMENSION']
                idx += int(dims[0])
                if len(dims) > 1:
                    l1 = 'The data cannot be properly represented in JSON-headed ASCII'\
                         ' as it has too high a rank\n'
                    l2 = 'key = {0} ({1})\n'.format(key, insd[key].shape)
                    l3 = 'Maximum allowed number of dimensions is 2\n'
                    warnings.warn(''.join([l1, l2, l3]), DMWarning)
            except AttributeError: #AttrErr if just metadata
                #js_out[key]['DIMENSION'] = insd[key].attrs['DIMENSION']
                pass
        else: #is metadata
            if verbose:
                print('metadata: {0}'.format(key))
            js_out[key]['VALUES'] = dmcopy(_dateToISO(insd[key]))
            js_out[key]['DIMENSION'] = [len(js_out[key]['VALUES'])]
        for kk in js_out[key]:
            try:
                js_out[key][kk] = js_out[key][kk].tolist()
            except:
                pass
    json_str = json.dumps(js_out, indent=4, sort_keys=True)
    reob = re.compile(r'\[.*?\]', re.DOTALL)
    json_str = re.sub(reob, stripNL, json_str) #put lists back onto one line
    #add comment field for header
    json_str = ''.join(['#', json_str])
    json_str = '\n#'.join(json_str.split('\n'))
    json_str = ''.join([json_str, '\n'])

    if isinstance(fname, str_classes):
        with open(fname, 'w') as fh:
            fh.writelines(json_str)
    elif hasattr(fname, 'writelines'):
        fname.writelines(json_str)
    elif (fname is None) and (returnString):
        return json_str

    if returnString:
        return json_str


def _dateToISO(indict):
    """
    covert datetimes to iso strings inside of datamodel attributes
    """
    retdict = dmcopy(indict)
    if isinstance(indict, dict):
        for key in indict:
            if isinstance(indict[key], datetime.datetime):
                retdict[key] = retdict[key].isoformat()
            elif hasattr(indict[key], '__iter__'):
                for idx, el in enumerate(indict[key]):
                    if isinstance(el, datetime.datetime):
                        retdict[key][idx] = el.isoformat()
    else:
        if isinstance(indict, datetime.datetime):
            retdict = retdict.isoformat()
        elif hasattr(indict, '__iter__'):
            retdict = numpy.asanyarray(retdict)
            for idx, el in numpy.ndenumerate(indict):
                if isinstance(el, datetime.datetime):
                    retdict[idx] = el.isoformat()
    return retdict


def toJSONheadedASCII(fname, insd, metadata=None, depend0=None, order=None, **kwargs):
    '''Write JSON-headed ASCII file of data with metadata from SpaceData object

    Parameters
    ----------
    fname : str
        Filename to write to (can also use a file-like object)
        None can be given in conjunction with the returnString keyword to skip writing output

    insd : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    depend0 : str (optional)
        variable name to use to indicate parameter on which other data depend (e.g. Time)
    order : list (optional)
        list of key names in order of start column in output JSON file
    metadata: str or file-like (optional)
        filename with JSON header to use (or file-like with JSON metadata)
    delimiter: str
        delimiter to use in ASCII output (default is whitespace), for tab, use '\t'

    Returns
    -------
    None

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> data = dm.SpaceData()
    >>> data.attrs['Global'] = 'A global attribute'
    >>> data['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
    >>> data['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
    >>> data['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
    >>> dm.toJSONheadedASCII('outFile.txt', data, depend0='Var1', order=['Var1'])
    #Note that not all field names are required, those not given will be listed
    #alphabetically after those that are specified
    '''
    kwarg_dict = {'delimiter': ' '}
    for key in kwarg_dict:
        if key in kwargs:
            kwarg_dict[key] = kwargs[key]
    if not metadata:
        metadata = io.StringIO()
        writeJSONMetadata(metadata, insd, depend0=depend0, order=order)
        metadata.seek(0) #rewind StringIO object to start
    hdr = readJSONMetadata(metadata)

    datlist = []
    for key in hdr:
        if 'START_COLUMN' in hdr[key].attrs:
            #add to list of (start_col, keyname) pairs
            datlist.append((hdr[key].attrs['START_COLUMN'], key, hdr[key].attrs['DIMENSION'][0]))
            #also use for data length
            datlen = len(insd[key])
            if datlen == 0:
                raise ValueError('No data present to write: Use writeJSONmetadata')
                #TODO: Set this to just default to writing the header out and raise a warning
    datlist.sort()
    ncols = datlist[-1][0]+datlist[-1][-1]

    #now open file (file-like) and for each line in len(data)
    #write the line using start_column, name, dimension
    data = numpy.zeros([datlen, ncols], dtype=object)
    for stcol, name, dim in datlist:
        if dim == 1:
            data[:, stcol] = _dateToISO(insd[name])
        else:
            data[:, stcol:stcol+dim] = _dateToISO(insd[name])
    hdstr = writeJSONMetadata(None, hdr, depend0=depend0, order=order, returnString=True)
    with open(fname, 'w') as fh:
        fh.writelines(hdstr)
        for line in data:
            prline = kwarg_dict['delimiter'].join([str(el) for el in line])
            fh.write(''.join([prline, '\n']))


def fromRecArray(recarr):
    '''Takes a numpy recarray and returns each field as a dmarray in a SpaceData container

    Parameters
    ----------
    recarr : numpy record array
        object to parse into SpaceData container

    Returns
    -------
    sd: spacepy.datamodel.SpaceData
        dict-like containing arrays of named records in recarr

    Examples
    --------
    >>> import numpy as np
    >>> import spacepy.datamodel as dm
    >>> x = np.array([(1.0, 2), (3.0, 4)], dtype=[('x', float), ('y', int)])
    >>> print(x, x.dtype)
    array([(1.0, 2), (3.0, 4)], dtype=[('x', '<f8'), ('y', '<i4')])
    >>> sd = dm.fromRecArray(x)
    >>> sd.tree(verbose=1)
    +
    |____x (spacepy.datamodel.dmarray (2,))
    |____y (spacepy.datamodel.dmarray (2,))
    '''
    sd = SpaceData()
    for key in recarr.dtype.fields.keys():
        sd[key] = dmarray(recarr[key])
    return sd

def toRecArray(sdo):
    '''Takes a SpaceData and creates a numpy recarray

    Parameters
    ----------
    sdo : SpaceData
        SpaceData to change to a numpy recarray

    Returns
    -------
    recarr: numpy record array
        numpy.recarray object with the same values (attributes are lost)

    Examples
    --------
    >>> import numpy as np
    >>> import spacepy.datamodel as dm
    >>> sd = dm.SpaceData()
    >>> sd['x'] = dm.dmarray([1.0, 2.0])
    >>> sd['y'] = dm.dmarray([2,4])
    >>> sd.tree(verbose=1)
    +
    |____x (spacepy.datamodel.dmarray (2,))
    |____y (spacepy.datamodel.dmarray (2,))
    >>> ra = dm.toRecArray(sd)
    >>> print(ra, ra.dtype)
    [(2, 1.0) (4, 2.0)] (numpy.record, [('y', '<i8'), ('x', '<f8')])
    '''
    nametype = numpy.dtype([(k, sdo[k].dtype.str) for k in sdo])
    recarr = numpy.rec.fromarrays([sdo[k] for k in sdo], dtype=nametype)
    return recarr


def dmcopy(dobj):
    '''Generic copy utility to return a copy of a (datamodel) object

    Parameters
    ----------
    dobj : object
        object to return a copy of

    Returns
    -------
    copy_obj: object (same type as input)
        copy of input oibject

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> dat = dm.dmarray([2,3], attrs={'units': 'T'})
    >>> dat1 = dm.dmcopy(dat)
    >>> dat1.attrs['copy': True]
    >>> dat is dat1
    False
    >>> dat1.attrs
    {'copy': True, 'units': 'T'}
    >>> dat.attrs
    {'units': 'T'}
    '''
    if isinstance(dobj, (SpaceData, dmarray)):
        return copy.deepcopy(dobj)
    if isinstance(dobj, numpy.ndarray):
        return numpy.copy(dobj)
    return copy.copy(dobj)

def createISTPattrs(datatype, ndims=1, vartype=None, units=' ', NRV=False):
    '''Return set of unpopulated attributes for ISTP compliant variable

    Parameters
    ----------
    datatype : {'data', 'support_data', 'metadata'}
        datatype of variable to create metadata for.
    ndims : int
        number of dimensions, default=1
    vartype : {'float', 'char', 'int', 'epoch', 'tt2000'}
        The type of the variable, default=float
    units : str
        The units of the variable, default=' '
    NRV : bool
        Is the variable NRV (non-record varying), default=False

    Returns
    -------
    attrs : dict
        dictionary of attributes for the variable

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> dm.createISTPattrs('data', ndims=2, vartype='float', units='MeV')
    {'CATDESC': '',
     'DISPLAY_TYPE': 'spectrogram',
     'FIELDNAM': '',
     'FILLVAL': -1e+31,
     'FORMAT': 'F18.6',
     'LABLAXIS': '',
     'SI_CONVERSION': ' > ',
     'UNITS': 'MeV',
     'VALIDMIN': '',
     'VALIDMAX': '',
     'VAR_TYPE': 'data',
     'DEPEND_0': 'Epoch',
     'DEPEND_1': ''}
    '''
    fillvals = {'float': -1e31,
                'char': '',
                'int': numpy.array(-2147483648).astype(numpy.int32),
                'epoch': -1.0E31, #datetime.datetime(9999,12,31,23,59,59,999)}
                'tt2000': numpy.array(-9223372036854775808).astype(numpy.int64)}
    formats = {'float': 'F18.6',
               'char': 'A30',
               'int': 'I11',
               'epoch': '',
               'tt2000': 'I21'}
    disp = {1: 'time_series',
            2: 'spectrogram',
            3: 'spectrogram',
            4: 'spectrogram'}
    if vartype not in fillvals:
        fill = -1e31
        form = 'F15.6'
    else:
        fill = fillvals[vartype]
        form = formats[vartype]

    unit = units

    if datatype == 'data':
        attrs = {
            'CATDESC': '',
            'DISPLAY_TYPE': disp[ndims],
            'FIELDNAM': '',
            'FILLVAL': fill,
            'FORMAT': form,
            'LABLAXIS': '',
            'SI_CONVERSION': ' > ',
            'UNITS': unit,
            'VALIDMIN': '',
            'VALIDMAX': '',
            'VAR_TYPE': 'data'
            }
        for dim in range(ndims):
            attrs['DEPEND_{0}'.format(dim)] = ''
        attrs['DEPEND_0'] = 'Epoch'
    elif datatype == 'support_data':
        attrs = {
            'CATDESC': '',
            'FIELDNAM': '',
            'FORMAT': form,
            'UNITS': unit,
            'VAR_TYPE': 'support_data'
            }
        for dim in range(ndims):
            attrs['DEPEND_{0}'.format(dim)] = ''
        if not NRV:
            attrs['VALIDMIN'] = ''
            attrs['VALIDMAX'] = ''
            attrs['FILLVAL'] = fill
            attrs['DEPEND_0'] = 'Epoch'
        else:
            del attrs['DEPEND_0']
    elif datatype == 'metadata':
        attrs = {
            'CATDESC': '',
            'FIELDNAM': '',
            'FORMAT': form,
            'UNITS': unit,
            'VAR_TYPE': 'metadata'
            }
        for dim in range(ndims):
            attrs['DEPEND_{0}'.format(dim)] = ''
        if not NRV:
            attrs['FILLVAL'] = fill
            attrs['DEPEND_0'] = 'Epoch'
        else:
            del attrs['DEPEND_0']
    else:
        raise ValueError("Invalid datatype (data|support_data|metadata)")

    return attrs


def _getVarLengths(data):
    """
    get the length of all the variables

    Parameters
    ----------
    data : SpaceData
        SpaceData object to return the length of the variables

    Returns
    -------
    data : dict
        dict of the names and lengths of a SpaceData
    """
    ans = {}
    for k, v in data.items():
        ans[k] = len(v)
    return ans

def resample(data, time=None, winsize=0, overlap=0, st_time=None, outtimename='Epoch'):
    """
    resample a SpaceData to a new time interval

    Parameters
    ----------
    data : SpaceData or dmarray
        SpaceData with data to resample or dmarray with data to resample,
        variables can only be 1d or 2d, if time is specified only variables
        the same length as time are resampled, otherwise only variables
        with length equal to the longest length are resampled

    time : array-like
        dmarray of times the correspond to the data

    winsize : datetime.timedelta
        Time frame to average the data over

    overlap : datetime.timedelta
        Overlap in the moving average

    st_time : datetime.datetime
        Starting time for the resample, if not specified the time of the first
        data point is used (see spacepy.toolbox.windowMean)

    Returns
    -------
    ans : SpaceData
        Resampled data, included keys are in the input keys (with the data caveats above)
        and Epoch which contains the output time

    Examples
    --------
    >>> import datetime
    >>> import spacepy.datamodel as dm
    >>> a = dm.SpaceData()
    >>> a.attrs['foo'] = 'bar'
    >>> a['a'] = dm.dmarray(range(10*2)).reshape(10,2)
    >>> a['b'] = dm.dmarray(range(10)) + 4
    >>> a['c'] = dm.dmarray(range(3)) + 10
    >>> times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i) for i in range(10)]
    >>> out = dm.resample(a, times, winsize=datetime.timedelta(hours=2), overlap=datetime.timedelta(hours=0))
    >>> out.tree(verbose=1, attrs=1)
    # +
    # :|____foo (str [3])
    # |____Epoch (spacepy.datamodel.dmarray (4,))
    # |____a (spacepy.datamodel.dmarray (4, 2))
    # :|____DEPEND_0 (str [5])
    #
    # Things to note:
    #    - attributes are preserved
    #    - the output variables have their DEPEND_0 changed to Epoch (or outtimename)
    #    - each dimension of a 2d array is resampled individually
    """
    from . import toolbox
    # check for SpaceData or dmarray input before going to a bunch of work
    if not isinstance(data, (SpaceData, dmarray)):
        raise TypeError('Input must be a SpaceData or dmarray object')
    if time is None:
        time = []
    # can only resample variables that have the same length as time,
    #    if time is default then use all the vals that are the same
    #    as the longest var
    lent = len(time)
    if lent == 0:
        lent = len(data[max(data, key=lambda k: len(data[k]))])
    keys = [k for k in data if len(data[k]) == lent]
    # what time are we starting at?
    try:
        t_int = time.UTC
    except AttributeError:
        t_int = dmarray(time)
    if t_int.any() and ((st_time is None) and isinstance(t_int[0], datetime.datetime)):
        st_time = t_int[0].replace(hour=0, minute=0, second=0, microsecond=0)

    ans = SpaceData()
    ans.attrs = data.attrs

    for k in keys:
        if len(data[k].shape) > 1:
            if len(data[k].shape) > 2:
                raise IndexError("Variables can only be 1d or 2d")
            for i in range(data[k].shape[1]):
                d, t = toolbox.windowMean(data[k][:, i], time=t_int, winsize=winsize,
                                          overlap=overlap, st_time=st_time)
                if k not in ans:
                    ans[k] = dmarray(d)
                else:
                    ans[k] = dmarray.vstack(ans[k], d)
            ans[k] = ans[k].T
        else:
            d, t = toolbox.windowMean(data[k], time=t_int, winsize=winsize,
                                      overlap=overlap, st_time=st_time)
            ans[k] = dmarray(d)
        try:
            ans[k].attrs = data[k].attrs
        except AttributeError: # was not a dmarray
            pass
        ans[k].attrs['DEPEND_0'] = outtimename
    ans[outtimename] = dmarray(t)

    return ans
