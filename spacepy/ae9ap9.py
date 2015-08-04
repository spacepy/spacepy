#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for reading and dealing with AE9AP9 data files.
See https://www.vdl.afrl.af.mil/programs/ae9ap9/ to download the model.
This is not a AE9AP9 runner.

Authors: Brian Larsen
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov

Copyright 2015 Los Alamos National Security, LLC.
"""

__contact__ = 'Brian Larsen, balarsen@lanl.gov'

import os
import re

import numpy as np

def readEphem(fname):
    """
    read a model generated ephemeris file into a datamodel.SpaceData object

    Parameters
    ==========
    fname : str
        filename of the Ephemeris file

    Returns
    =======
    out : spacepy.datamodel.SpaceData
        Data contained in the Ephemeris file

    Examples
    ========
    """
    # Ae9Ap9 Command-line Utility orbit data file: Run/ephem_sat.dat
    # generated from specified elements, using propagator: Kepler-J2
    #
    # Time format:              Modified Julian Date
    # Coordinate system:        GEI (Geocentric Equatorial Inertial) Cartesian in Earth radii
    # Data Delimiter:           comma
    #
    # datetime(mjd),posx(Re),posy(Re),posz(Re)

def _parseInfo(header):
    """
    given a header parse and return the common information in all headers

    
    # Time format:              Modified Julian Date
    # Coordinate system:        GEI (Geocentric Equatorial Inertial) Cartesian in Earth radii
    # Data Delimiter:           comma
    """
    # get the time format, it is in all the files
    ans = {}
    for ind, val in enumerate(header):
        if "Time format:" in val:
            if "Modified Julian Date" in val:
                ans['time_format'] = 'MJD'
            else:
                raise(NotImplementedError("Sorry can't read that time format yet"))
        elif "Coordinate system:" in val:
            if "GEI (Geocentric Equatorial Inertial) Cartesian in Earth radii" in val:
                ans['coord_system'] = ('GEI', 'Re')
            else:
                raise(NotImplementedError("Sorry can't read that coordinate system yet"))
        elif "Data Delimiter:" in val:
            if "comma" in val:
                ans['delimiter'] = ','
            else:
                raise(NotImplementedError("Sorry can't read that delimiter yet"))
        elif "Ae9Ap9 Software Version:" in val:
            match = re.match(r'^Ae9Ap9 Software Version:.*(\d\.\d\d\.\d\d\d)$', val)
            version = match.group(1).strip()
            ans['software_version'] = version
        elif "Model type:" in val:
            match = re.match(r'^Model type:(.*)$', val)
            mtype = match.group(1).strip()
            ans['model_type'] = mtype           
        elif "Flux type:" in val:
            match = re.match(r'^Flux type:(.*)$', val)
            ftype = match.group(1).strip()
            ans['flux_type'] = ftype
        elif "Flux direction:" in val:
            match = re.match(r'^Flux direction:(.*)$', val)
            fdir = match.group(1).strip()
            ans['flux_direction'] = fdir
        elif "Energy levels" in val:
            match = re.match(r'^Energy levels.*\((.*)\):(.*)$', val)
            ans['energy'] = (np.asarray(match.group(2).strip().split()).astype(float), match.group(1).strip())
        elif "generated from specified elements" in val:
            match = re.search(r'^generated from specified elements.*:\ (.*)$', val)
            ans['propigator'] = match.group(1).strip()

               
    return ans
                    

def _readHeader(fname):
    """
    read only the header from a Ae9Ap9 file
    """
    dat = []
    with open(fname, 'r') as fp:
        while True:
            tmp = fp.readline()
            if tmp[0] != '#':
                break
            dat.append(tmp.strip())
    dat = [v[1:].strip() for v in dat]
    return dat

def parseHeader(fname):
    """
    given an AE9AP9 output test file parse the header and return the information in a
    dictionary

    """
    ans = {}
    ans['filename'] = os.path.abspath(fname)
           
    header = _readHeader(fname)
    # get the information
    ans.update(_parseInfo(header))
            
    columns = header[-1].split(ans['delimiter'])
    columns_unq = _unique_elements_order(columns)
    len_columns_unq = len(columns)
            
    # pull the units off the column headings
    cols = []
    while columns_unq:
        tmp = columns_unq.pop(0)
        umatch = re.match(r'(.*)\((.*)\).*', tmp)
        if umatch:
            cols.append((umatch.group(1), umatch.group(2)))
        else:
            cols.append(tmp)
    # now go through and see if any non-tuple entry has the same name as a tuple entry
    to_remove = []
    for ii, c1 in enumerate(cols):
        if not isinstance(c1, tuple):
            for jj, c2 in enumerate(cols[ii:]):
                if c2[0] == c1:
                    to_remove.append(ii)
    if to_remove:
        for v in to_remove:
            del cols[v]
    # and now replicate the last entry enough times to make it the right length
    cols.extend([cols[-1]] * (len_columns_unq - len(cols)))
            
    ans['columns'] = cols

    if 'pctile' in os.path.basename(fname): #capture the percentile
        match = re.match(r'.*_(\d\d).txt', os.path.basename(fname))
        ans['percentile'] = float(match.group(1))
    return ans

def _unique_elements_order(seq):
    """
    see http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
    """
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]
