#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
code to create and plot generic spectrograms for space science - these are different
than the standard spectrogram in other fields in that there is no fft used


Authors: Brian Larsen and Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov, smorley@lanl.gov
Los Alamos National Laboratory

Copyright 2011 Los Alamos National Security, LLC.
"""

import spacepy.datamodel as dm

class spectrogram(dm.SpaceData):
    """
    This class generates and then contains the data binned intot he spectrogram

    It is meant to be used on arbitary data series.  The first series "x" is
    plotted on the abscissa and secnd series "y" is plotted on the ordinate and
    the thre series "z" is plotted in color.

    The series are not passed in independently but instead inside a
    spacepy.datamodel.SpaceData container.  Helper routines are provided to
    facilitate the creation of the SpaceData container if the data are not in the format.
    """

    ## NOTE this will need to set the sphonx var autoclass_content to "both"

    def __init__(self, data, **kwargs):
        """
        Parameters
        ==========
        data : spacepy.datamodel.SpaceData
            The data for the spectrogram, the variables to be used default to
            "Epoch" for x, "Energy" for y, and "Flux" for z.  Other names are
            specified using the 'variables' keyword.

        Other Parameters
        ================
        variables : list
            keyword containing the names of the variables to use for the spectrogram
            the list is a list of the SpaceData keys in x, y, z, order
        """
        # is variables in kwargs?
        if not "variables" in kwargs:
            variables = ['Epoch', 'Energy', 'Flux']
        # check to see if the variables are in the spacedata
        for var in kwargs['variables']:
            if not var in data:  # TODO could check other capitilizations
                raise(ValueError(str(var) + ' not found in the input data' ))
