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

Functions are also available to directly load data and metadata into a SpacePy datamodel from
NASA CDF as well as JSON-headed ASCII. Writers also exist to output a SpacePy datamodel directly
to HDF5 or JSON-headed ASCII. See :py:func:`datamodel.fromCDF`, :py:func:`datamodel.readJSONheadedASCII`,
:py:func:`datamodel.toHDF5`, and :py:func:`datamodel.toJSONheadedASCII` for more details.


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

from __future__ import division
import copy
import datetime
import itertools
import json
from functools import partial
import os
import re
import warnings

try:
    import StringIO # can't use cStringIO as we might have unicode
except ImportError:
    import io as StringIO

import numpy
from . import toolbox
from . import time as spt


__contact__ = 'Steve Morley, smorley@lanl.gov'

try:
    str_classes = (str, bytes, unicode)
except NameError:
    str_classes = (str, bytes)
    unicode = str

class DMWarning(Warning):
    """
    Warnings class for datamodel, subclassed so it can be set to always
    """
    pass
warnings.simplefilter('always', DMWarning)

class dmarray(numpy.ndarray):
    """
    Container for data within a SpaceData object

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
        else:
            return numpy.ndarray.__array_wrap__(self, out_arr, context).tolist()

    def __reduce__(self):
        """This is called when pickling, see:
        http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html
        for this particular example.
        Only the attributes in Allowed_Attributes can exist
        """
        object_state = list(numpy.ndarray.__reduce__(self))
        subclass_state = tuple([tuple([val, self.__getattribute__(val)]) for val in self.Allowed_Attributes])
        object_state[2] = (object_state[2],subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Used for unpickling after __reduce__ the self.attrs is recovered from
        the way it was saved and reset.
        """
        nd_state, own_state = state
        numpy.ndarray.__setstate__(self,nd_state)
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
        if name == 'Allowed_Attributes':
            pass
        elif not name in self.Allowed_Attributes:
            raise(TypeError("Only attribute listed in Allowed_Attributes can be set"))
        super(dmarray, self).__setattr__(name, value)

    def addAttribute(self, name, value=None):
        """Method to add an attribute to a dmarray
        equivalent to
        a = datamodel.dmarray([1,2,3])
        a.Allowed_Attributes = a.Allowed_Attributes + ['blabla']
        """
        if name in self.Allowed_Attributes:
            raise(NameError('{0} is already an attribute cannot add again'.format(name)))
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
            backup.append( (atr, dmcopy(self.__getattribute__(atr)) ) )
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
        outarr = dmarray(numpy.vstack( (one, other) ))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def hstack(cls, one, other):
        """
        hstack data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.hstack( (one, other) ))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def dstack(cls, one, other):
        """
        dstack data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.dstack( (one, other) ))
        return cls._replaceAttrs(outarr, backup)

    @classmethod
    def concatenate(cls, one, other, axis=0):
        """
        concatenate data to an existing dmarray
        """
        backup = one._saveAttrs()
        outarr = dmarray(numpy.concatenate( (one, other) , axis=axis ))
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
    try:
        a.fill(fillval)
    except TypeError:
        obj = numpy.core.numeric._maketup(dtype, fillval)
        a.fill(obj)
    return a


class SpaceData(dict):
    """
    Datamodel class extending dict by adding attributes.

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
        This allows one to make a SpaceData indexed with an iterable of keys to return a new spacedata
        made of the subset of keys
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
                raise(KeyError('{0}'.format(key)))

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
        self.toCDF = partial(toCDF, SDobject=self, *args, **kwargs)
        self.toCDF.__doc__ = toCDF.__doc__
        self.toHDF5 = partial(toHDF5, SDobject=self, *args, **kwargs)
        self.toHDF5.__doc__ = toHDF5.__doc__
        self.toJSONheadedASCII = partial(toJSONheadedASCII, insd=self, *args, **kwargs)
        self.toJSONheadedASCII.__doc__ = toJSONheadedASCII.__doc__

    #Need to remove the partials on copy, http://bugs.python.org/issue4380
    #They will be recreated by __init__
    def __getstate__(self):
        d = copy.copy(self.__dict__)
        del d['toCDF']
        del d['toHDF5']
        del d['toJSONheadedASCII']
        return d

## To enable string output of repr, instead of just printing, uncomment his block
#    def __repr__(self):
#        #redirect stdout to StringIO
#        import StringIO, sys
#        dum = StringIO.StringIO()
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
        verbose : boolean (optional)
            print more info
        spaces : string (optional)
            string will added for every line
        levels : integer (optional)
            number of levels to recurse through (True means all)
        attrs : boolean (optional)
            display information for attributes

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
        toolbox.dictree(self, **kwargs)

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

            

def convertKeysToStr(SDobject):
    if isinstance(SDobject, SpaceData):
        newSDobject = SpaceData()
        newSDobject.attrs = SDobject.attrs
    else:
        newSDobject = {}
    for key in SDobject:
        if not isinstance(key, str):
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
    except (TypeError):
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
    except (TypeError):
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
        addme[grp] = unflatten(addme[grp], marker=marker) #recurse to make sure everything inside is unpacked
    return addme


def fromCDF(fname, **kwargs):
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
    spacepy.pycdf.CDF.copy
    '''
    #TODO: add unflatten keyword and restore flattened variables
    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and SpacePy's pyCDF")

    with pycdf.CDF(fname) as cdfdata:
        return cdfdata.copy()

def toCDF(fname, SDobject, **kwargs):
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
        create CDF in backward-compatible format (default is v3+ compatibility only)

    TT2000 : bool (optional)
        write variables beginning with 'Epoch' as datatype CDF_TT2000 (default is automatic selection of EPOCH or EPOCH16)

    verbose : bool (optional)
        verbosity flag

    Returns
    -------
    None
    '''
    defaults = {'skeleton': '',
                'flatten': False,
                'overwrite': False,
                'compress': False,
                'autoNRV': False,
                'backward': False,
                'TT2000': False,
                'verbose': False}
    for key in kwargs:
        if key in defaults:
            defaults[key] = kwargs[key]

    if defaults['flatten']:
        SDobject = SDobject.flatten()
    if defaults['overwrite']:
        raise NotImplementedError('Overwriting CDFs is not currently enabled - please remove the file manually')

    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and SpacePy's pyCDF")
    pycdf.lib.set_backward(False)
    with pycdf.CDF(fname, defaults['skeleton']) as outdata:
        if hasattr(SDobject, 'attrs'):
            for akey in SDobject.attrs:
                outdata.attrs[akey] = dmcopy(SDobject.attrs[akey])
        varLengths = [len(SDobject[var]) for var in SDobject]
        modeLength = next(itertools.groupby((reversed(sorted(varLengths)))))[0]
        for key in SDobject:
            if isinstance(SDobject[key], dict):
                raise TypeError('This data structure appears to be nested, please try spacepy.datamodel.flatten')
            if not defaults['skeleton']:
                if not SDobject[key].shape:
                    shape_tup=-1
                else:
                    shape_tup = SDobject[key].shape
                if 'Epoch' not in SDobject:
                    NRVtest = modeLength
                else:
                    NRVtest = len(SDobject['Epoch'])
                if shape_tup[0] != NRVtest: #naive check for 'should-be' NRV
                    try:
                        foo = outdata.new(key, SDobject[key][...], recVary=False)
                        if defaults['verbose']: print('{0} is being made NRV'.format(key))
                        outdata[key].attrs = dmcopy(SDobject[key].attrs)
                    except ValueError:
                        foo = outdata.new(key, SDobject[key].tolist, recVary=False)
                        outdata[key].attrs = dmcopy(SDobject[key].attrs)
                if defaults['TT2000'] and 'Epoch' in key:
                    foo = outdata.new(key, SDobject[key][...], type=pycdf.const.CDF_TIME_TT2000)
                else:
                    try:
                        outdata[key] = SDobject[key]
                    except ValueError:
                        try:
                            outdata[key] = dmarray([SDobject[key].tolist()], attrs=dmcopy(SDobject[key].attrs)).squeeze()
                        except UnicodeEncodeError:
                            tmpAttrs = dmcopy(SDobject[key].attrs)
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
                outdata[key][...] = SDobject[key][...]
                for akey in outdata[key].attrs:
                    try:
                        outdata[key].attrs[akey] = dmcopy(SDobject[key].attrs[akey])
                    except ValueError:
                        outdata[key][...] = dmarray([SDobject[key].tolist()], attrs=dmcopy(SDobject[key].attrs))
                    except KeyError:
                        pass
    return None


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
        if hasattr(hfile[path],'attrs'):
            #for key, value in hfile[path].attrs.iteritems():
            for key in hfile[path].attrs:
                try:
                    value = hfile[path].attrs[key]
                except TypeError:
                    warnings.warn('Unsupported datatype in dataset {0}.attrs[{1}]'.format(path,key))
                    continue
                try:
                    SDobject.attrs[key] = value
                except:
                    warnings.warn('The following key:value pair is not permitted\n' +
                                    'key = {0} ({1})\n'.format(key, type(key)) +
                                    'value = {0} ({1})'.format(value, type(value)), DMWarning)

    try:
        import h5py as hdf
    except ImportError:
        raise ImportError('HDF5 converter requires h5py')

    if type(fname) in str_classes:
        hfile = hdf.File(fname, mode='r')
    else:
        hfile = fname
        #should test here for HDF file object

    if 'path' not in kwargs:
        path = '/'
    else:
        path = kwargs['path']

    SDobject = SpaceData()
    allowed_elems = [hdf.Group, hdf.Dataset]
    ##carry over the attributes
    hdfcarryattrs(SDobject, hfile, path)
    ##carry over the groups and datasets
    for key, value in hfile[path].items():
        #try:
            if type(value) is allowed_elems[0]: #if a group
                SDobject[key] = SpaceData()
                SDobject[key] = fromHDF5(hfile, path=path+'/'+key)
            elif type(value) is allowed_elems[1]: #if a dataset
                try:
                    SDobject[key] = dmarray(value)
                except (TypeError, ZeroDivisionError): #ZeroDivisionError catches zero-sized DataSets
                    SDobject[key] = dmarray(None)
                hdfcarryattrs(SDobject[key], hfile, path+'/'+key)
        #except:
        #    raise ValueError('HDF5 file contains type other than Group or Dataset')
    if path=='/': hfile.close()
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
        compress all the variables using this method (default None) (gzip, shuffle, fletcher32, szip, lzf)
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
                if type(value) in allowed_attrs:
                    #test for datetimes in iterables
                    if hasattr(value, '__iter__') and not isinstance(value, str_classes):
                        dumval = [b.isoformat() if isinstance(b, datetime.datetime) else b for b in value]
                    truth = False
                    try:
                        if value.nbytes: truth = True #empty arrays of any dimension are nbytes=0
                    except AttributeError: #not an array
                        if value or value is 0: truth = True

                    if truth:
                        if bytes is str:
                            if type(key) is unicode:
                                dumkey = key.encode('utf-8')
                            if type(value) is unicode:
                                dumval = value.encode('utf-8')
                        uni = False #No special unicode handling
                        if not bytes is str: #Python 3
                            dumval = numpy.asanyarray(dumval)
                            if dumval.size and dumval.dtype.kind == 'U':
                                uni = True #Unicode list, special handling
                        try:
                            if uni:
                                #Tell hdf5 this is unicode. Numpy is UCS-4, HDF5 is UTF-8
                                hfile[path].attrs.create(dumkey, dumval,
                                    dtype=hdf.special_dtype(vlen=unicode))
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
                    if bytes is str and type(key) is unicode:
                        dumkey = str(key)
                    hfile[path].attrs[dumkey] = dumval
                else:
                    #TODO: add support for arrays(?) in attrs (convert to isoformat)
                    warnings.warn('The following key:value pair is not permitted\n' +
                                    'key = {0} ({1})\n'.format(key, type(key)) +
                                    'value type {0} is not in the allowed attribute list'.format(type(value)),
                                        DMWarning)

    try:
        import h5py as hdf
    except ImportError:
        raise ImportError('h5py is required to use HDF5 files')
    
    try:
        assert isinstance(SDobject, SpaceData)
    except AssertionError:
        raise ValueError("Input data is not of type SpaceData, check usage: toHDF5(fname, datamodel)")
    #mash these into a defaults dict...
    if 'mode' not in kwargs:
        wr_mo = 'a'
    else:
        wr_mo = kwargs['mode']
    if 'compression' not in kwargs:
        h5_compr_type = None
    else:
        h5_compr_type = kwargs['compression']
        if h5_compr_type not in ['gzip', 'szip', 'lzf', 'shuffle', 'fletcher32', None]:
            raise NotImplementedError('Specified compression type not supported')
    if ('compression_opts' not in kwargs) or (h5_compr_type == 'lzf'):
        h5_compr_opts = None
    else:
        h5_compr_opts = kwargs['compression_opts']

    if 'overwrite' not in kwargs: kwargs['overwrite'] = True
    if type(fname) == str:
        if os.path.isfile(fname) and not kwargs['overwrite']:
            raise(IOError('Cannot write HDF5, file exists (see overwrite) "{0!s}"'.format(fname)))
        if os.path.isfile(fname) and kwargs['overwrite']:
            os.remove(fname)
        hfile = hdf.File(fname, mode=wr_mo)
        must_close = True
    else:
        hfile = fname
        #should test here for HDF file object
        must_close = False
    if 'path' in kwargs:
        path = kwargs['path']
    else:
        path = '/'

    try:
        allowed_attrs = [int, long, float, str, unicode, numpy.ndarray, list, tuple, numpy.string_]
    except NameError:
        allowed_attrs = [int,       float, bytes, str, numpy.ndarray, list, tuple, numpy.string_]
    for v in numpy.typecodes['AllInteger']:
        allowed_attrs.append(numpy.typeDict[v])
    for v in numpy.typecodes['AllFloat']:
        allowed_attrs.append(numpy.typeDict[v])

    allowed_elems = [SpaceData, dmarray]

    #first convert non-string keys to str
    SDobject = convertKeysToStr(SDobject)
    SDcarryattrs(SDobject,hfile,path,allowed_attrs)

    try:
        for key, value in SDobject.items():
            if isinstance(value, allowed_elems[0]):
                hfile[path].create_group(key)
                toHDF5(hfile, SDobject[key], path=path+'/'+key, compression=h5_compr_type, compression_opts=h5_compr_opts)
            elif isinstance(value, allowed_elems[1]):
                try:
                    hfile[path].create_dataset(key, data=value, compression=h5_compr_type, compression_opts=h5_compr_opts)
                except:
                    dumval = value.copy()
                    if isinstance(value[0], datetime.datetime):
                        for i, val in enumerate(value): dumval[i] = val.isoformat()
                    hfile[path].create_dataset(key, data=dumval.astype('|S35'), compression=h5_compr_type, compression_opts=h5_compr_opts)
                    #else:
                    #    hfile[path].create_dataset(key, data=value.astype(float))
                SDcarryattrs(SDobject[key], hfile, path+'/'+key, allowed_attrs)
            else:
                warnings.warn('The following data is not being written as is not of an allowed type\n' +
                               'key = {0} ({1})\n'.format(key, type(key)) +
                                  'value type {0} is not in the allowed data type list'.format(type(value)),
                                      DMWarning)
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
           varLinks=False, linkFormat=None, echo=False, tableTag='<table border="1">'):
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
    output = StringIO.StringIO() # put the output into a StringIO
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
                if not isinstance(SDobject[key].attrs[attr], (str, unicode)):
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
                html = html.replace('!N', '</sub>', 1 ) # just replace 1
            else:
                html = html + '</sub>'
        elif code == '!E':
            if '!N' in html:
                html = html.replace('!N', '</sup>', 1 ) # just replace 1
            else:
                html = html + '</sup>'
    return html

def readJSONMetadata(fname, **kwargs):
    '''Read JSON metadata from an ASCII data file

    Parameters
    ----------
    fname : str
        Filename to read metadata from

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
        with open(fname, 'r') as f:
            lines = f.read()

    # isolate header
    p_srch = re.compile(r"^#(.*)$", re.M)
    hreg = re.findall(p_srch, lines)
    header = "".join(hreg)

    # isolate JSON field
    srch = re.search( r'\{\s*(.*)\s*\}', header )
    if isinstance(srch, type(None)):
        raise IOError('The input file has no valid JSON header. Must be valid JSON bounded by braces "{ }".')
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
       if 'START_COLUMN' in mdatadict[key]:
           mdata[key] = SpaceData(attrs=mdatadict[key])
       elif 'VALUES' in mdatadict[key]:
           dum = mdatadict[key].pop('VALUES')
           mdata[key] = dmarray(dum, attrs=mdatadict[key])
       else:
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
    try:
        if isinstance(fname, (str, unicode)):
            fname=[fname]
        elif hasattr(fname, 'readlines'):
            fname = [fname]
            filelike = True
    except NameError: # for Py3
        if isinstance(fname, str):
            fname=[fname]
    if not mdata:
        mdata = readJSONMetadata(fname[0])
    if restrict:
        delkeys = [kk for kk in mdata.keys() if kk not in restrict]
        for val in delkeys:
            del mdata[val] #remove undesired keys
    mdata_copy = dmcopy(mdata)
    def innerloop(fh, mdata, mdata_copy):
        line = fh.readline()
        if not str is bytes:
            line = line.decode('latin1')
        while (line and line[0]==comment):
            line = fh.readline()
            if not str is bytes:
                line = line.decode('latin1')
        fh.seek(-len(line), os.SEEK_CUR) # fixes the missing first data bug
        alldata = fh.readlines()
        if not alldata:
            return mdata
        if not str is bytes:
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
                    if not varDims.shape: varDims = numpy.array([varDims])
                    singleDim = True
                    if len(varDims)>1 or varDims[0]>1:
                        singleDim = False
                if ('DIMENSION' in mdata_copy[key].attrs) and not singleDim:
                    en = int(mdata_copy[key].attrs['DIMENSION'][0]) + int(st)
                    try:
                        assert mdata[key]=={}
                        mdata[key] = data[:,int(st):int(en)]
                    except (AssertionError, ValueError):
                        mdata[key] = numpy.vstack((mdata[key], data[:,int(st):int(en)]))
                else:
                    try:
                        assert mdata[key]=={}
                        mdata[key] = data[:,int(st)]
                    except (AssertionError, ValueError):
                        mdata[key] = numpy.hstack((mdata[key], data[:,int(st)]))
        return mdata
    for fn in fname:
        if not filelike:
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
            conversions=convert
        else:
            conversions = {'DateTime': lambda x: dup.parse(x, ignoretz=True),
                           'ExtModel': lambda x: str(x)}
        for conkey in conversions:
            try:
                name = keys.pop(keys.index(conkey)) #remove from keylist
            except ValueError:
                warnings.warn('Key {0} for conversion not found in file'.format(conkey), UserWarning)
                continue
            #https://github.com/numpy/numpy/issues/2160
            if len(mdata[name].shape) == 1:
                for i,element in numpy.ndenumerate(mdata[name]):
                    mdata[name][i[0]] = conversions[name](element)
            else:
                for i,element in numpy.ndenumerate(mdata[name]):
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
        out = text.group().replace('\n','').replace('  ','')
        return out

    #if required, identify depend0 for deciding what's data/metadata
    if depend0 is None:
        #search for DEPEND_0 in metadata
        for key in insd:
            if not hasattr(insd[key], 'attrs'):
                insd[key] = dmarray(insd[key])
            if 'DEPEND_0' in insd[key].attrs:
                depend0 = insd[key].attrs['DEPEND_0']
                if not isinstance(depend0, str):
                    #assume it's a singleton list
                    depend0 = depend0[0]
                if not isinstance(depend0, str):
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
    else:
        if not depend0 in insd: raise KeyError('Invalid key supplied for ordering metadata on write')
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
            if key not in order: keylist.append(key)
    else:
        ##TODO do we want to have DEPEND0 first in order by default?
        keylist = sorted(insd.keys())

    idx = 0
    for key in keylist:
        js_out[key] = dmcopy(_dateToISO(insd[key].attrs))
        if len(insd[key]) == datalen: #is data
            if verbose: print('data: {0}'.format(key))
            try:
                js_out[key]['DIMENSION'] = list(insd[key].shape[1:])
                if not js_out[key]['DIMENSION']: js_out[key]['DIMENSION'] = [1]
                js_out[key]['START_COLUMN'] = idx
                dims = js_out[key]['DIMENSION']
                idx += int(dims[0])
                if len(dims)>1:
                    l1 = 'The data cannot be properly represented in JSON-headed ASCII as it has too high a rank\n'
                    l2 = 'key = {0} ({1})\n'.format(key, insd[key].shape)
                    l3 = 'Maximum allowed number of dimensions is 2\n'
                    warnings.warn(''.join([l1, l2, l3]), DMWarning)
            except AttributeError: #AttrErr if just metadata
                #js_out[key]['DIMENSION'] = insd[key].attrs['DIMENSION']
                pass
        else: #is metadata
            if verbose: print('metadata: {0}'.format(key))
            js_out[key]['VALUES'] = dmcopy(_dateToISO(insd[key]))
            js_out[key]['DIMENSION'] = [len(js_out[key]['VALUES'])]
        for kk in js_out[key]:
            try:
                js_out[key][kk] = js_out[key][kk].tolist()
            except:
                pass
    json_str = json.dumps(js_out, indent=4, sort_keys=True)
    reob = re.compile('\[.*?\]', re.DOTALL)
    json_str = re.sub(reob, stripNL, json_str) #put lists back onto one line
    #add comment field for header
    json_str = ''.join(['#', json_str])
    json_str = '\n#'.join(json_str.split('\n'))
    json_str = ''.join([json_str,'\n'])

    if isinstance(fname, str):
        with open(fname,'w') as fh:
            fh.writelines(json_str)
    elif hasattr(fname, 'writelines'):
        fname.writelines(json_str)
    elif (fname is None) and (returnString):
        return json_str

    if returnString: return json_str


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
            retdict = numpy.asanyarray(indict)
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
        metadata = StringIO.StringIO()
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
            if datlen==0:
                raise ValueError('No data present to write: Use writeJSONmetadata')
                #TODO: Set this to just default to writing the header out and raise a warning
    datlist.sort()
    ncols = datlist[-1][0]+datlist[-1][-1]

    #now open file (file-like) and for each line in len(data)
    #write the line using start_column, name, dimension
    data = numpy.zeros([datlen,ncols], dtype=object)
    for stcol, name, dim in datlist:
        if dim==1:
            data[:, stcol] = _dateToISO(insd[name])
        else:
            data[:, stcol:stcol+dim] = _dateToISO(insd[name])
    hdstr = writeJSONMetadata(None, hdr, depend0=depend0, order=order, returnString=True)
    with open(fname, 'w') as fh:
        fh.writelines(hdstr)
        for line in data:
            prline = kwarg_dict['delimiter'].join([str(el) for el in line])
            fh.write(''.join([prline,'\n']))


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
    keys = list(sdo.keys())
    recarr = numpy.rec.fromarrays( [sdo[k] for k in sdo], names=keys)
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
    elif isinstance(dobj, numpy.ndarray):
        return numpy.copy(dobj)
    else:
        return copy.copy(dobj)

def createISTPattrs(datatype, ndims=1, vartype=None, units=None, NRV=False):
    '''Return set of unpopulated attributes for ISTP compliant variable
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
    if units:
        unit = units
    else:
        unit = ''
    if datatype == 'data':
        attrs = {'CATDESC': '',
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
        attrs = {'CATDESC': '',
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
        attrs = {'CATDESC': '',
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

def resample(data, time=[], winsize=0, overlap=0, st_time=None, outtimename='Epoch'):
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
    # check for SpaceData or dmarray input before going to a bunch of work
    if not isinstance(data, (SpaceData, dmarray)):
        raise(TypeError('Input must be a SpaceData or dmarray object'))

    # can only resample variables that have the same length as time,
    #    if time is default then use all the vals that are the same
    #    as the longest var
    lent = len(time)
    if lent == 0: lent = len(data[max(data, key=lambda k: len(data[k]))])
    keys = [k for k in data if len(data[k]) == lent]
    d2 = data[keys]
    # what time are we starting at?
    if isinstance(time, spt.Ticktock):
        t_int = time.UTC
    else:
        t_int = dmarray(time)
    if t_int.any() and ((st_time is None) and isinstance(t_int[0], datetime.datetime)):
        st_time = t_int[0].replace(hour=0, minute=0, second=0, microsecond=0)

    ans = SpaceData()
    ans.attrs = data.attrs
    
    for k in keys:
        if len(data[k].shape) > 1:
            if len(data[k].shape) > 2:
                raise(IndexError("Variables can only be 1d or 2d"))
            for i in range(data[k].shape[1]):
                d, t = toolbox.windowMean(data[k][:,i], time=t_int, winsize=winsize, overlap=overlap, st_time=st_time)
                if k not in ans:
                    ans[k] = dmarray(d)
                else:
                    ans[k] = dmarray.vstack(ans[k], d)
            ans[k] = ans[k].T
        else:
            d, t = toolbox.windowMean(data[k], time=t_int, winsize=winsize, overlap=overlap, st_time=st_time)
            ans[k] = dmarray(d)
        try:
            ans[k].attrs = data[k].attrs
        except AttributeError: # was not a dmarray
            pass
        ans[k].attrs['DEPEND_0'] = outtimename
    ans[outtimename] = dmarray(t)
        
    return ans


    
