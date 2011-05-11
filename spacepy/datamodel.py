# -*- coding: utf-8 -*-
"""
The datamodel classes consitute a data model implementation
meant to mirror the functionality of the data model output from pycdf, though
implemented slightlydifferently.

This contains the following classes:
 * dmarray - numpy arrays that support .attrs for information about the data
 * SpaceData - base class that extends dict, to be extended by others
    Currently used in GPScode and other projects

 Example:
 >>> import spacepy.datamodel as datamodel
 >>> position = datamodel.dmarray([1,2,3], attrs={'coord_system':'GSM'})
 >>> position
 dmarray([1, 2, 3])
 >>> position.attrs
 {'coord_system': 'GSM'}
"""

from __future__ import division
import numpy, copy

class dmarray(numpy.ndarray):
    """
    Container for data within a SpaceData object

    @author: Brian Larsen, Steve Morley
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov

    @version: V1: 01-Mar-2011 Based on GPSarray from GPScode codebase

    Example:
    >>> import spacepy.datamodel as datamodel
    >>> position = datamodel.dmarray([1,2,3], attrs={'coord_system':'GSM'})
    >>> position
    dmarray([1, 2, 3])
    >>> position.attrs
    {'coord_system': 'GSM'}
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
    Base datamodel class extending dict

    Currently has just method stubs, no real functionality

    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov

    @version: V1: 01-Mar-2011 Based on GPSarray from GPScode codebase
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

    #def __repr__(self):
        #"""
        #Abstract method, reimplement
        #"""
        #super(SpaceData, self).__repr__()
        ##raise(ValueError("Abstract method called, reimplement __repr__"))


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

    addme = dobj.__class__()
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
