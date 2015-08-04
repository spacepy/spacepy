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

import datetime
import functools
import os
import re

from dateutil import relativedelta
import numpy as np
import spacepy.datamodel as dm
import spacepy.time as spt

def readFile(fname, comments='#'):
    """
    read a model generated file into a datamodel.SpaceData object

    Parameters
    ==========
    fname : str
        filename of the file

    Other Parameters
    ================
    comments : str (optional)
        String that is the comments in the data file

    Returns
    =======
    out : :mod:`~spacepy.datamodel.SpaceData`
        Data contained in the  file

    Examples
    ========
    from spacepy import ae9ap9
    ae9ap9.readFile('ephem_sat.dat').tree(verbose=1)
    # +
    # |____Epoch (spacepy.datamodel.dmarray (121,))
    # |____GSE (spacepy.datamodel.dmarray (121, 3))
    # |____MJD (spacepy.datamodel.dmarray (121,))
    # |____posComp (spacepy.datamodel.dmarray (3,))
    """    
    # get the header information first
    header = parseHeader(fname)
    # and read in all the data
    data = np.loadtxt(fname, delimiter=header['delimiter'], comments=comments)
    # now parse the data    
    ans = dm.SpaceData()
    #parse the datetime if it is there (it is always first)
    if 'datetime' in header['columns'][0]:
        if header['time_format'] == 'eDOY': # have to massage the data first
            year = data[:,0].astype(int)
            frac = data[:,1]
            time = spt.Ticktock([datetime.datetime(y,1,1) + relativedelta.relativedelta(days=v)
                               for y, v in zip(year, frac)], 'UTC')
            ans[header['time_format']] = dm.dmarray(data[:,0:2])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            ans[header['time_format']].attrs['LABL_PTR_1'] = 'timeComp'
            ans['timeComp'] = dm.dmarray(['Year', 'eDOY'])
            ans['timeComp'].attrs['VAR_TYPE'] = 'metadata'
            del header['columns'][0:2]
            data = data[:,2:]
        elif header['time_format'] == 'UTC': # have to massage the data first
            tm = data[:,0:6].astype(int)
            time = spt.Ticktock([datetime.datetime(*v) for v in tm], 'UTC' )
            ans[header['time_format']] = dm.dmarray(data[:,0:6])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            ans[header['time_format']].attrs['LABL_PTR_1'] = 'timeComp'
            ans['timeComp'] = dm.dmarray(['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second'])
            ans['timeComp'].attrs['VAR_TYPE'] = 'metadata'
            del header['columns'][0:6]
            data = data[:,6:]
        else:
            time = spt.Ticktock(data[:,0], header['time_format'])
            ans[header['time_format']] = dm.dmarray(data[:,0])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            del header['columns'][0]
            data = data[:,1:]

        ans['Epoch'] = dm.dmarray(time.UTC)
        ans['Epoch'].attrs['VAR_TYPE'] = 'support_data'
    # parse the position, it is always next
    if 'posx' in header['columns'][0]:
        varname = header['coord_system'][0]
        pos = dm.dmarray(data[:,0:3])
        ans[varname] = pos
        ans[varname].attrs['UNITS'] =  header['coord_system'][1]
        ans[varname].attrs['FIELDNAM'] = varname
        ans[varname].attrs['LABLAXIS'] = 'Position'
        ans[varname].attrs['DEPEND_0'] = 'Epoch'
        ans[varname].attrs['LABL_PTR_1'] = 'posComp'
        ans[varname].attrs['VAR_TYPE'] = 'data'
        ans['posComp'] = dm.dmarray(['X', 'Y', 'Z'])
        ans['posComp'].attrs['VAR_TYPE'] = 'metadata'
        del header['columns'][0:3]
        data = data[:,3:]
    
    # if there are PA they are now
    if header['columns'] and ('pitch' in header['columns'][0]):
        raise(NotImplementedError("Sorry have not implemented pitch angle resolved files yet"))
    ## This is an issue and the files have repeadedlines at the same time for each pitch
    ##  angle. They would all need to be read in them combines into a multi-d array
    
    # now parse fluence, flux, or doserate
    # do all the rest of the column headers match?
    col_arr = np.asarray(header['columns'])
    if header['columns'] and (col_arr == col_arr[0]).all():
        varname = col_arr[0][0].title()
        ans[varname] = dm.dmarray(data)
        ans[varname].attrs['UNITS'] = str(col_arr[0][1])
        ans[varname].attrs['FIELDNAM'] = varname
        ans[varname].attrs['LABLAXIS'] = varname        
        ans[varname].attrs['VAR_TYPE'] = 'data'        
        ans[varname].attrs['SCALETYP'] = 'log'        
        header['varname'] = varname
        ans[varname].attrs['Description'] = '{flux_direction} {varname} based on {flux_type} from the {model_type} model'.format(**header)
        ans[varname].attrs['DEPEND_0'] = 'Epoch'
        if ans[varname].shape[1] == header['energy'][0].shape[0]:
            ans[varname].attrs['DEPEND_1'] = 'Energy'
        if 'percentile' in header:
            ans[varname].attrs['TITLE'] = '{0} percentile'.format(header['percentile'])

    # create the energy variable
    if 'energy' in header:
        ans['Energy'] = dm.dmarray(header['energy'][0])
        ans['Energy'].attrs['UNITS'] = header['energy'][1]
        ans['Energy'].attrs['FIELDNAM'] = 'Energy'
        ans['Energy'].attrs['LABLAXIS'] = 'Energy'
        ans['Energy'].attrs['VAR_TYPE'] = 'data'
        ans['Energy'].attrs['SCALETYP'] = 'log'
        ans['Energy'].attrs['VAR_TYPE'] = 'support_data'

    skip_header = ['columns', 'coord_system', 'delimiter', 'energy', 'flux_direction']
    for k in header:
        if k not in skip_header:
            ans.attrs[k] = header[k]   
    return ans
    
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
            elif "Year, day_of_year.frac" in val:
                ans['time_format'] = 'eDOY'
            elif "Year, Month, Day, Hour, Minute, Seconds" in val:
                ans['time_format'] = 'UTC'
            else:
                raise(NotImplementedError("Sorry can't read that time format yet: {0}".format(val)))
        elif "Coordinate system:" in val:
            coord_sys = val.split(':')[1].strip().split()[0]
            if "in Earth radii" in val:
                ans['coord_system'] = (coord_sys, 'Re')
            else:
                ans['coord_system'] = (coord_sys, '')
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
            if '=' not in val:
                match = re.match(r'^Flux direction:(.*)$', val)
                fdir = match.group(1).strip()
            else:
                pa = val.split('=')[1].strip()
                fdir = np.asarray(pa.split())                
            ans['flux_direction'] = fdir
        elif "Energy levels" in val:
            match = re.match(r'^Energy levels.*\((.*)\):(.*)$', val)
            ans['energy'] = (np.asarray(match.group(2).strip().split()).astype(float), match.group(1).strip())
        elif "generated from specified elements" in val:
            match = re.search(r'^generated from specified elements.*:\ (.*)$', val)
            ans['propagator'] = match.group(1).strip()
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

    Parameters
    ==========
    fname : str
        filename of the file

    Returns
    =======
    out : dict
        Dictionary of the header information in the file
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

def combinePercentiles(files, timeframe='all'):
    """
    combine files at different percentiles into one file with the spectra at different
    percentiles for easy plotting and analysis

    NOTE: Requires pandas for any timeframe other than 'all'

    Parameters
    ==========
    files : str
        Filenames of the percentile files

    Other Parameters
    ================
    timeframe : str
        Timeframe to average itn eh input spectra over (either 'all' or a pandas understoop resample() time
    
    Returns
    =======
    out : :mod:`~spacepy.datamodel.SpaceData`
        Combined output spectra of the file
    """
    if not files:
        raise(ValueError("Must input files"))
    data = {}
    for fname in files:
        tmp = readFile(fname)
        if 'percentile' not in tmp.attrs:
            raise(ValueError("File {0} does not have a percentile key".format(fname)))
        data[tmp.attrs['percentile']] = tmp

    varnames = ['Flux', 'Fluence', 'Dose', None]
    for varname in varnames:
        if varname in data[list(data.keys())[0]]:
            break # this leaves it where it was found
    if varname is None:
        raise(ValueError("Did not understand file type could not find one of {0}".format(varnames)))

    # now that they are all read in, we collapse them down to different time bases
    if timeframe == 'all':
        # this is just a full collapse
        # find the var to collapse
        d_comp = np.asarray([])
        ps = sorted(data.keys())[::-1]
        for p in ps:
            if data[p][varname].attrs['DEPEND_0'] == 'Epoch':
                d_comp = np.concatenate((d_comp, data[p][varname].mean(axis=0))) # full time collapse
            else:
                d_comp = np.concatenate((d_comp, data[p][varname])) # full time collapse
        d_comp = d_comp.reshape(len(ps), -1).T
    else:
        try:
            import pandas as pd
        except ImportError:
            raise(ImportError("panads is required for timeframe other than 'all'"))
        raise(NotImplementedError("Timeframes other than 'all' not yet implemented"))
        # make a dataframe
        # resample to timeframe
        # join everything
        # meet up at same place for output
    # now that we have all the data prep it to a spacedata
    ans = dm.SpaceData()
    ans[varname] = dm.dmarray(d_comp)
    if timeframe == 'all':
        ans[varname].attrs['DEPEND_0'] = 'Energy'
        ans[varname].attrs['DEPEND_1'] = 'Percentile'
        ans[varname].attrs['DISPLAY_TYPE'] = 'time_series'
    for k in data[p][varname].attrs:
        if 'DEPEND' in k:
            continue
        else:
            ans[varname].attrs[k] = data[p][varname].attrs[k]
    ans['Energy'] = data[p]['Energy']
    ans['Percentile'] = dm.dmarray(ps)
    ans['Percentile'].attrs['FIELDNAM'] = 'Percentile'
    ans['Percentile'].attrs['LABLAXIS'] = 'Percentile'
    ans['Percentile'].attrs['VALIDMIN'] = 0.0
    ans['Percentile'].attrs['VALIDMAX'] = 100.0
    ans['Percentile'].attrs['SCALETYP'] = 'support_data'

    return ans

def _getData(fnames):
    """
    helper routine for dm.toCDF and dm.toHDF5 and dm.toJSONheadedASCII since the data
    prep is all the same
    """
    if hasattr(fnames, 'strip'): # it is a string
        return readFile(fnames)
    else:
        return combinePercentiles(fnames)

def _toFile(fnames, outname, ftype):
    """
    generic toXXX() routine to change data files to other formats
    """
    dat = _getData(fnames)
    if ftype == 'CDF':
        dm.toCDF(outname, dat, autoNRV=True)
    elif ftype == 'HDF5':
        dm.toHDF5(outname, dat, compression='gzip', compression_opts=7)
    elif ftype == 'ASCII':
        dm.toJSONheadedASCII(outname, dat)
        
toCDF = functools.partial(_toFile, ftype='CDF')
"""
Convert the input file(s) fnames to CDF. If fnames is a single file then call readFile()
if it is an iterable of filenames call combinePercentiles()

Parameters
==========
fnames : str
    Filename(s) to change to CDF

See Also
========
spacepy.datamodel.toCDF
"""

toHDF5 = functools.partial(_toFile, ftype='HDF5')
"""
Convert the input file(s) fnames to HDF5. If fnames is a single file then call readFile()
if it is an iterable of filenames call combinePercentiles()

Parameters
==========
fnames : str
    Filename(s) to change to HDF5

See Also
========
spacepy.datamodel.toHDF5
"""

toJSONheadedASCII = functools.partial(_toFile, ftype='ASCII')
"""
Convert the input file(s) fnames to JSONheadedASCII. If fnames is a single file then call readFile()
if it is an iterable of filenames call combinePercentiles()

Parameters
==========
fnames : str
    Filename(s) to change to toJSONheadedASCII

See Also
========
spacepy.datamodel.toJSONheadedASCII
"""


