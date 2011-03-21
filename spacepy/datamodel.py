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
import numpy

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
    def __new__(cls, input_array, attrs={}):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = numpy.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.attrs = attrs
        # Finally, return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        self.attrs = getattr(obj, 'attrs', {})

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
