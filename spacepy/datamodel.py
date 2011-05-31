#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The datamodel classes consitute a data model implementation
meant to mirror the functionality of the data model output from pycdf, though
implemented slightlydifferently.

This contains the following classes:
 * dmarray - numpy arrays that support .attrs for information about the data
 * SpaceData - base class that extends dict, to be extended by others
    Currently used in GPScode and other projects

Authors: Steve Morley and Brian Larsen
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov; balarsen@lanl.gov

About datamodel
---------------

The SpacePy datamodel module implents classes that are designed to make implementing a standard
data model easy. The concepts are very similar to those used in standards like HDF5, netCDF and
NASA CDF.

The basic container type is analagous to a folder (on a filesystem; HDF5 calls this a
group): Here we implement this as a dictionary-like object, a datamodel.SpaceData object, which
also carries attributes. These attributes can be considered to be global, i.e. relevant for the
entire folder. The next container type is for storing data and is based on a numpy array, this
class is datamodel.dmarray and also carries attributes. The dmarray class is analagous to an
HDF5 dataset.


Guide for NASA CDF users
------------------------

By definition, a NASA CDF only has a single `layer'. That is, a CDF contains a series of records 
(stored variables of various types) and a set of attributes that are either global or local in
scope. Thus to use SpacePy's datamodel to capture the functionality of CDF the two basic data types
are all that is required, and the main constraint is that datamodel.SpaceData objects cannot be
nested (more on this later, if conversion from a nested datamodel to a flat datamodel is required).

This is best illustrated with an example. Imagine representing some satellite data within a CDF -- 
the global attributes might be the mission name and the instrument PI, the variables might be the
instrument counts [n-dimensional array], timestamps[1-dimensional array and an orbit number [scalar].
Each variable will have one attribute (for this example).

    >>> import spacepy.datamodel as dm
    >>> mydata = dm.SpaceData(attrs={'MissionName': 'BigSat1'})
    >>> mydata['Counts'] = dm.dmarray([[42, 69, 77], [100, 200, 250]], attrs={'Units': 'cnts/s'})
    >>> mydata['Epoch'] = dm.dmarray([1, 2, 3], attrs={'units': 'minutes'})
    >>> mydata['OrbitNumber'] = dm.dmarray(16, attrs={'StartsFrom': 1})
    >>> mydata.attrs['PI'] 'Prof. Big Shot'

This has now populated a structure that can map directly to a NASA CDF. To visualize our datamodel, 
we can use the toolbox function dictree (which works for any dictionary-like object, including PyCDF
file objects).

    >>> import spacepy.toolbox as tb
    >>> tb.dictree(mydata, attrs=True)
    +
    :|____MissionName
    :|____PI
    |____Counts
         :|____Units
    |____Epoch
         :|____units
    |____OrbitNumber
         :|____StartsFrom


Copyright ©2010 Los Alamos National Security, LLC.
"""

from __future__ import division
import numpy, copy

class dmarray(numpy.ndarray):
    """
    Container for data within a SpaceData object

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

    """
    Allowed_Attributes = ['attrs']

    def __new__(cls, input_array, attrs=None):
       # Input array is an already formed ndarray instance
       # We first cast to be our class type
       obj = numpy.asarray(input_array).view(cls)
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
            self.__setattr__(val, getattr(obj, val, {}))

    def __reduce__(self):
        """This is called when pickling, see:
        http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html
        for this particular examnple.
        Only the attributes in Allowed_Attributes can exist
        """
        object_state = list(numpy.ndarray.__reduce__(self))
        subclass_state = tuple([tuple([val, self.__getattribute__(val)]) for val in self.Allowed_Attributes])
        object_state[2] = (object_state[2],subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Used for unpickling after __reduce__ the self.attrs is recoved from
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
        equilivant to
        a = datamodel.dmarray([1,2,3])
        a.Allowed_Attributes = a.Allowed_Attributes + ['blabla']
        """
        if name in self.Allowed_Attributes:
            raise(NameError('{0} is aleady an attribute cannot add again'.format(name)))
        self.Allowed_Attributes.append(name)
        self.__setattr__(name, value)


class SpaceData(dict):
    """
    Datamodel class extending dict

    """

    def __init__(self, *args, **kwargs):
        """
        Base class for "Data Model" representation data

        Abstract method, reimplement
        """
        #raise(ValueError("Abstract method called, reimplement __init__"))
        self.attrs = {}
        if 'attrs' in kwargs:
            if hasattr(kwargs['attrs'], '__getitem__'):
                self.attrs = kwargs['attrs']
            del kwargs['attrs']

        super(SpaceData, self).__init__(*args, **kwargs)

    def flatten(self):
        '''Method to collapse datamodel to one level deep

        Examples
        --------

        >>> import spacepy.datamodel as dm
        >>> import spacepy.toolbox as tb
        >>> a = dm.SpaceData()
        >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        >>> a['4'] = dm.SpaceData(cat = 'kitty')
        >>> a['5'] = 4
        >>> tb.dictree(a)
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
        >>> tb.dictree(b)
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5

        >>> a.flatten()
        >>> tb.dictree(a)
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5


        See Also
        --------
        spacepy.datamodel.flatten

        '''


        flatobj = flatten(self)
        remkeys = [key for key in self]
        for key in remkeys:
            del self[key]
        for key in flatobj:
            self[key] = copy.copy(flatobj[key])

def flatten(dobj):
    '''Function to collapse datamodel to one level deep

    Examples
    --------

    >>> import spacepy.datamodel as dm
    >>> import spacepy.toolbox as tb
    >>> a = dm.SpaceData()
    >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
    >>> a['4'] = dm.SpaceData(cat = 'kitty')
    >>> a['5'] = 4
    >>> tb.dictree(a)
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
    >>> tb.dictree(b)
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5

    >>> a.flatten()
    >>> tb.dictree(a)
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5


    See Also
    --------
    spacepy.datamodel.SpaceData.flatten

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
