#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Tools to read and process omni data (Qin-Denton, etc.)

Authors: Steve Morley, Josef Koller
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov

Copyright 2010-2014 Los Alamos National Security, LLC.


About omni
----------

The omni module primarily manages the hourly OMNI2 and Qin-Denton data, which
are sourced from the Virtual Radiation Belt Observatory (`ViRBO <http://virbo.org>`_), 
who maintain these data sources. The data can be kept up-to-date in SpacePy 
using the :func:`~spacepy.toolbox.update` function in the :mod:`spacepy.toolbox` module.

The OMNI2 data combines data from a variety of satellites that sample the solar
wind (notably ACE and Wind), and propagates the data to Earth's bow shock nose.
The `Qin-Denton <http://virbo.org/QinDenton>`_ data is derived from the OMNI2 
data and is designed for providing input to the Tsyganenko magnetic field 
models. The later Tsyganenko magnetic field models require subsidiary parameters
(G- and W-parameters) that are pre-calculated in the Qin-Denton data. Further,
the Qin-Denton data contains no data gaps -- all gaps are filled (for details on
the gap filling, see the paper by `Qin et al. <http://dx.doi.org/10.1029/2006SW000296>`_.)


Advanced features
-----------------

Higher resolution data, or custom data sources, can also be managed/accessed 
with this module, although this is considered an advanced use for this module.
This is achieved using custom names for the dbase keyword in get_omni, which
must be defined in the SpacePy configuration file (for a user-install on linux,
this is ~/.spacepy/spacepy.rc; see :doc:`configuration`).
An example of the formatting required is

qd1min: /usr/somedir/QinDenton/YYYY/QinDenton_YYYYMMDD_1min.txt

In this example the custom data source name is qd1min. Wildcard substitutions 
can be made for the year (YYYY), month (MM) and day (DD). Future updates will
give more flexibility in data storage model, but currently we assume that all
custom data sources follow a convention in which the data files are daily, and
the files are organized into folders by year. The year, month and day must all
be specified in the filename.

Currently there are some restrictions on the data format for custom data 
sources. The stored data must currently be stored as JSON-headed ASCII.
If data conversions are required, then a valid dictionary of conversion
functions must be supplied via the convert keyword argument. See 
:func:`~spacepy.datamodel.readJSONheadedASCII` for details.
Additionally, by default this will interpolate the data to the requested time 
ticks. To return only the actual recorded data values for the specified time 
range set the keyword argument interp to False.

"""
import bisect, re, os
import sys
import warnings

import numpy as np
from spacepy.datamodel import SpaceData, dmarray, dmcopy, unflatten, readJSONheadedASCII, dmfilled, fromHDF5
from spacepy.toolbox import tOverlapHalf, indsFromXrange
import spacepy.time as spt

__contact__ = 'Steve Morley, smorley@lanl.gov'

def get_omni(ticks, dbase='QDhourly', **kwargs):
    '''
    Returns Qin-Denton OMNI values, interpolated to any time-base from a default hourly resolution

    The update function in toolbox retrieves all available hourly Qin-Denton data, 
    and this function accesses that and interpolates to the given times,
    returning the OMNI values as a SpaceData (dict-like) with
    Kp, Dst, dens, velo, Pdyn, ByIMF, BzIMF, G1, G2, G3, etc.
    (see also http://www.dartmouth.edu/~rdenton/magpar/index.html and
    http://www.agu.org/pubs/crossref/2007/2006SW000296.shtml )

    Parameters
    ==========
    ticks : Ticktock class or array-like of datetimes
        time values for desired output

    dbase : str (optional)
        Select data source, options are 'QDhourly', 'OMNI2hourly', 'Mergedhourly'
        Note - Custom data sources can be specified in the spacepy config file
        as described in the module documentation.

    Returns
    =======
    out : spacepy.datamodel.SpaceData
        containing all Qin-Denton values at times given by ticks

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> ticks = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> d = om.get_omni(ticks)
    >>> d.tree(levels=1)
    +
    |____ByIMF
    |____Bz1
    |____Bz2
    |____Bz3
    |____Bz4
    |____Bz5
    |____Bz6
    |____BzIMF
    |____DOY
    |____Dst
    |____G1
    |____G2
    |____G3
    |____Hr
    |____Kp
    |____Pdyn
    |____Qbits
    |____RDT
    |____UTC
    |____W1
    |____W2
    |____W3
    |____W4
    |____W5
    |____W6
    |____Year
    |____akp3
    |____dens
    |____ticks
    |____velo


    Notes
    =====
    Note about Qbits: If the status variable is 2, the quantity you are using is fairly well
    determined. If it is 1, the value has some connection to measured values, but is not directly
    measured. These values are still better than just using an average value, but not as good
    as those with the status variable equal to 2. If the status variable is 0, the quantity is
    based on average quantities, and the values listed are no better than an average value. The
    lower the status variable, the less confident you should be in the value.

    '''
    dbase_options = {'QDhourly'    : 1,
                     'OMNI2hourly' : 2,
                     'Mergedhourly': 3,
                     'Test'        : -9,
                     }

    if not isinstance(ticks, spt.Ticktock):
        try:
            ticks = spt.Ticktock(ticks, 'UTC')
        except:
            raise TypeError('get_omni: Input times must be a Ticktock object or a list of datetime objects')

    if not dbase in dbase_options:
        from spacepy import config
        if dbase in config:
            #If a dbase is specified that isn't a default, then it MUST be in the spacepy config
            qdpath = os.path.split(os.path.split(config[dbase])[0])[0]
            if not os.path.isdir(qdpath): raise IOError('Specified dbase ({0}) does not have a valid location ({1})'.format(dbase, config[dbase]))
            days = list(set([tt.date() for tt in ticks.UTC]))
            flist = ['']*len(days)
            fnpath, fnformat = os.path.split(config[dbase])
            for idx, day in enumerate(days):
                dp = fnpath.replace('YYYY', '{0}'.format(day.year))
                df = fnformat.replace('YYYY', '{0}'.format(day.year))
                df = df.replace('MM', '{0:02d}'.format(day.month))
                df = df.replace('DD', '{0:02d}'.format(day.day))
                flist[idx] = os.path.join(dp, df)
            if 'convert' in kwargs:
                convdict = kwargs['convert']
            else:
                convdict = True #set to True as default?
            if 'interp' not in kwargs:
                kwargs['interp'] = True
            data = readJSONheadedASCII(sorted(flist), convert=convdict)
            omniout = SpaceData()

            time_var = [var for var in ['DateTime', 'Time', 'Epoch', 'UTC'] if var in data]
            if time_var:
                use_t_var = time_var[0]
            else:
                #no obvious time variable in input files ... can't continue
                raise ValueError('No clear time variable in file')
            
            if kwargs['interp'] is True:    
                data['RDT'] = spt.Ticktock(data[use_t_var]).RDT
                keylist = sorted(data.keys())
                dum = keylist.pop(keylist.index(use_t_var))
                for key in keylist:
                    try:
                        omniout[key] = dmarray(np.interp(ticks.RDT, data['RDT'], data[key], left=np.NaN, right=np.NaN))
                        omniout[key].attrs = dmcopy(data[key].attrs)
                    except:
                        try:
                            omniout[key] = dmfilled([len(ticks.RDT), data[key].shape[1]], fillval=np.NaN, attrs=dmcopy(data[key].attrs))
                            for col in range(data[key].shape[1]):
                                omniout[key][:,col] = np.interp(ticks.RDT, data['RDT'], data[key][:,col], left=np.NaN, right=np.NaN)
                        except ValueError:
                            print('Failed to interpolate {0} to new time base, skipping variable'.format(key))
                        except IndexError:
                            print('Variable {0} appears to be non-record varying, skipping interpolation'.format(key))
                            omniout[key] = data[key]
                omniout['UTC'] = ticks.UTC 
            else:
                #Trim to specified times
                inds = tOverlapHalf([ticks[0].RDT, ticks[-1].RDT], spt.Ticktock(data['DateTime']).RDT)
                for key in data:
                    if len(inds) == len(data[key]):
                        omniout[key] = data[key][inds]
                    else: #is ancillary data
                        omniout[key] = data[key]
                #TODO: convert to same format as OMNI/QD read (or vice versa)
                omniout['UTC'] = omniout[use_t_var]
            return omniout
        else:
            raise IOError('Specified dbase ({0}) must be specified in spacepy.config'.format(dbase))

    def getattrs(hf, key):
        out = {}
        if hasattr(hf[key],'attrs'):
            for kk, value in hf[key].attrs.items():
                try:
                    out[kk] = value
                except:
                    pass
        return out

    def HrFromDT(indt):
        hour = indt.hour
        minute = indt.minute
        second = indt.second
        musecond = indt.microsecond
        return hour+(minute/60.0)+(second/3600.0)+(musecond/3600.0e3)

    import h5py as h5
    fname, QDkeylist, O2keylist = '', [], []
    omnivals = SpaceData()
    dbase_select = dbase_options[dbase]
    if dbase_select in [1, 3, -9]:
        if dbase_select > 0:
            ldb = 'QDhourly'
            fln = omnifln
        else:
            ldb = 'Test'
            fln = testfln
        with h5.File(fln, mode='r') as hfile:
            QDkeylist = [kk for kk in hfile if kk not in ['Qbits', 'UTC']]
            st, en = ticks[0].RDT, ticks[-1].RDT
            ##check that requested requested times are within range of data
            enval, stval = omnirange(dbase=ldb)[1], omnirange(dbase=ldb)[0]
            if (ticks.UTC[0]>enval) or (ticks[-1]<stval):
                raise ValueError('Requested dates are outside data range')
            if (ticks.UTC[-1]>enval) or (ticks[0]<stval):
                print('Warning: Some requested dates are outside data range ({0})'.format(ldb))
            inds = tOverlapHalf([st, en], hfile['RDT'], presort=True) #returns an xrange
            inds = indsFromXrange(inds)
            if inds[0] < 1: inds[0] = 1
            sl_op = slice(inds[0]-1, inds[-1]+2)
    
            fname = ','.join([fname,hfile.filename])
            omnivals.attrs = getattrs(hfile, '/')
            for key in QDkeylist:
                omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
                omnivals[key].attrs = getattrs(hfile, key)
            for key in hfile['Qbits']:
                omnivals['Qbits<--{0}'.format(key)] = dmarray(hfile['/Qbits/{0}'.format(key)][sl_op])
                omnivals['Qbits<--{0}'.format(key)].attrs = getattrs(hfile, '/Qbits/{0}'.format(key))
                QDkeylist.append('Qbits<--{0}'.format(key))

    if dbase_options[dbase] == 2 or dbase_options[dbase] == 3:
        ldb = 'OMNI2hourly'
        with h5.File(omni2fln, mode='r') as hfile:
            O2keylist = [kk for kk in hfile if kk not in ['Epoch','RDT']]
            st, en = ticks[0].RDT, ticks[-1].RDT
            ##check that requested requested times are within range of data
            enval, stval = omnirange(dbase=ldb)[1], omnirange(dbase=ldb)[0]
            if (ticks[0].UTC>enval) or (ticks[-1]<stval):
                raise ValueError('Requested dates are outside data range')
            if (ticks[-1].UTC>enval) or (ticks[0]<stval):
                print('Warning: Some requested dates are outside data range ({0})'.format(ldb))
            inds = tOverlapHalf([st, en], hfile['RDT'], presort=True) #returns an xrange
            inds = indsFromXrange(inds)
            if inds[0] < 1: inds[0] = 1
            sl_op = slice(inds[0]-1, inds[-1]+2)
        
            fname = ','.join([fname,hfile.filename])
            omnivals.attrs = getattrs(hfile, '/') #TODO: This overwrites the previous set on Merged load... Fix!
            omnivals['RDT_OMNI'] = dmarray(hfile['RDT'][sl_op])
            for key in O2keylist:
                omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
                omnivals[key].attrs = getattrs(hfile, key)

    if dbase_options[dbase] == 3:
        #prune "merged" SpaceData
        sigmas = [key for key in omnivals if 'sigma' in key]
        for sk in sigmas: del omnivals[sk]
        bees = [key for key in omnivals if re.search('B._', key)]
        for bs in bees: del omnivals[bs]
        aves = [key for key in omnivals if ('_ave' in key) or ('ave_' in key)]
        for av in aves: del omnivals[av]

    omniout = SpaceData(attrs=dmcopy(omnivals.attrs))
    omniout.attrs['filename'] = fname[1:]
    ###print('QDkeys: {0}\n\nO2keys: {1}'.format(QDkeylist, O2keylist))
    for key in sorted(omnivals.keys()):
        if key in O2keylist:
            omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT_OMNI'], omnivals[key], left=np.NaN, right=np.NaN))
            #set metadata -- assume this has been set properly in d/l'd file to match ECT-SOC files
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
        elif key in QDkeylist:
            omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT'], omnivals[key], left=np.NaN, right=np.NaN))
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
        if key == 'G3': #then we have all the Gs
            omniout['G'] = dmarray(np.vstack([omniout['G1'], omniout['G2'], omniout['G3']]).T)
            omniout['G'].attrs = dmcopy(omnivals['G1'].attrs)
            for i in range(1,4): del omniout['G{0}'.format(i)]
        if key == 'W6':
            omniout['W'] = dmarray(np.vstack([omniout['W1'], omniout['W2'], omniout['W3'], omniout['W4'], omniout['W5'], omniout['W6']]).T)
            omniout['W'].attrs = dmcopy(omnivals['W1'].attrs)
            for i in range(1,7): del omniout['W{0}'.format(i)]
        if 'Qbits' in key:
            #Qbits are integer vals, higher is better, so floor to get best representation of interpolated val
            omniout[key] = np.floor(omnivals[key]) 
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
            if 'G3' in key: #then we have all the Gs
                omniout['Qbits<--G'] = dmarray(np.vstack([omniout['Qbits<--G1'], omniout['Qbits<--G2'], omniout['Qbits<--G3']]).T)
                for i in range(1,4): del omniout['Qbits<--G{0}'.format(i)]
            if 'W6' in key:
                omniout['Qbits<--W'] = dmarray(np.vstack([omniout['Qbits<--W1'], omniout['Qbits<--W2'], omniout['Qbits<--W3'], omniout['Qbits<--W4'], omniout['Qbits<--W5'], omniout['Qbits<--W6']]).T)
                for i in range(1,7): del omniout['Qbits<--W{0}'.format(i)]

    omniout['ticks'] = ticks
    omniout['UTC'] = ticks.UTC
    omniout['Hr'] = dmarray([HrFromDT(val) for val in omniout['UTC']])
    omniout['Year'] = dmarray([val.year for val in omniout['UTC']])
    omniout = unflatten(omniout)

    return omniout


def omnirange(dbase='QDhourly'):
    '''Returns datetimes giving start and end times in the OMNI/Qin-Denton data

    The update function in toolbox retrieves all available hourly Qin-Denton data, 
    and this function accesses that and looks up the start and end times,
    returning them as datetime objects.

    Parameters
    ==========
    dbase : string (optional)
        name of omni database to check. Currently 'QDhourly' and 'OMNI2hourly'

    Returns
    =======
    omnirange : tuple
        containing two datetimes giving the start and end times of the available data

    Examples
    ========
    >>> import spacepy.omni as om
    >>> om.omnirange()
    (datetime.datetime(1963, 1, 1, 0, 0), datetime.datetime(2011, 11, 30, 23, 0))
    >>> om.omnirange(dbase='OMNI2hourly')
    (datetime.datetime(1963, 1, 1, 0, 0), datetime.datetime(2011, 11, 30, 23, 0))

    '''
    
    import h5py as h5
    infile = {'QDhourly': omnifln,
              'OMNI2hourly': omni2fln,
              'Test': testfln}
    if dbase not in infile:
        raise NotImplementedError('')
    # Possible time variables in the HDF file and their ticktock dtype
    timeinfo = [('UTC', 'UTC'), ('Epoch', 'ISO'), ('RDT', 'RDT')]
    with h5.File(infile[dbase], mode='r') as hfile:
        for varname, dtype in timeinfo:
            if varname in hfile:
                tt = spt.Ticktock([hfile[varname][0], hfile[varname][-1]],
                                  dtype=dtype)
                break
        else:
            raise ValueError('Cannot find time variable in {}'
                             .format(infile[dbase]))
    start, end = tt.UTC
    
    return start, end



#-----------------------------------------------

# check for omni file during import
import os, datetime
from spacepy import DOT_FLN, help
import h5py

omnifln = os.path.join(DOT_FLN,'data','omnidata.h5')
omni2fln = os.path.join(DOT_FLN,'data','omni2data.h5')
# Test data is stored relative to the test script
try:
    import spacepy_testing
    testfln = os.path.join(spacepy_testing.datadir, 'OMNItest.h5')
except ImportError:
    testfln = None
if testfln is None or not os.path.isfile(testfln):
    # Hope it's relative to current!
    testfln = os.path.join(os.path.abspath('data'), 'OMNItest.h5')

presentQD = h5py.is_hdf5(omnifln)
presentO2 = h5py.is_hdf5(omni2fln)
if not (presentQD and presentO2):
    warnings.warn(
        "Qin-Denton/OMNI2 data not found in current format."
        " This module has limited functionality."
        " Run spacepy.toolbox.update(QDomni=True) to download data.")
