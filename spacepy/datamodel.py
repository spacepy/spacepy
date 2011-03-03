# -*- coding: utf-8 -*-
from __future__ import division

import os, glob, re
import datetime
import mmap, numbers
import itertools
import warnings

import numpy

"""
datamodel class is a data model implementation meant to mirror the functionality
of the data model output from pycdf even if implemented differently
 * dmarray - numpy arrays that support .attrs for information about the data
 * SpaceData - base class that extends dict so be extended by others
    Currently used in GPScode and other projects

 Example:
 >>> import spacepy.datamodel as datamodel
 >>> position = datamodel.dmarray([1,2,3], attrs={'coord_system':'GSM'})
 >>> position
 dmarray([1, 2, 3])
 >>> position.attrs
 {'coord_system': 'GSM'}
"""

class dmarray(numpy.ndarray):
    """
    Container for data within a SpaceData object

    @author: Brian Larsen, Steve Morley
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov

    @version: V1: 01-Mar-2011 Based on GPSarray from GPScode codebase
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
    Base datamodel class
    """

    def __init__(self, *args, **kwargs):
        """
        Base class for "Data Model" representation data
        Abstract method, reimplement
        """
        raise(ValueError("Abstract method called, reimplement __init__"))

    def __repr__(self):
        """
        Abstract method, reimplement
        """
        raise(ValueError("Abstract method called, reimplement __repr__"))
