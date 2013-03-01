#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools to read and process omni data

Authors: Josef Koller, Steve Morley
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov, smorley@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.


"""
import bisect
import numpy as np
from spacepy.datamodel import SpaceData, dmarray, dmcopy, unflatten
import spacepy.time as spt

__contact__ = 'Steve Morley, smorley@lanl.gov'

#new get_omni
def get_omni(ticks, dbase='QDhourly'):
    '''
    Returns Qin-Denton OMNI values, interpolated to any time-base from 1-hourly resolution

    The update function in toolbox retrieves all available hourly Qin-Denton data, 
    and this function accesses that and interpolates to the given times,
    returning the OMNI values as a SpaceData (dict-like) with
    Kp, Dst, dens, velo, Pdyn, ByIMF, BzIMF, G1, G2, G3, etc.
    (see also http://www.dartmouth.edu/~rdenton/magpar/index.html and
    http://www.agu.org/pubs/crossref/2007/2006SW000296.shtml )

    Parameters
    ==========
    ticks : Ticktock class
        time values for desired output

    Returns
    =======
    out : spacepy.datamodel.SpaceData
        containing all omni values at times given by ticks

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> ticks = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> d = om.get_omni(ticks)
    >>> print(d.keys())
    ['velo', 'Bz6', 'Bz5', 'Bz4', 'Bz3', 'Bz2', 'Bz1', 'RDT', 'Dst',
    'akp3', 'DOY', 'Qbits', 'G3', 'G2', 'G1', 'Hr', 'ticks', 'BzIMF',
    'UTC', 'Kp', 'Pdyn', 'dens', 'ByIMF', 'W6', 'W5', 'W4', 'W3', 'W2',
    'W1', 'Year']

    Notes
    =====
    Note about Qbits: If the status variable is 2, the quantity you are using is fairly well
    determined. If it is 1, the value has some connection to measured values, but is not directly
    measured. These values are still better than just using an average value, but not as good
    as those with the status variable equal to 2. If the status variable is 0, the quantity is
    based on average quantities, and the values listed are no better than an average value. The
    lower the status variable, the less confident you should be in the value.

    '''
    if not dbase=='QDhourly': raise NotImplementedError

    import h5py as h5
    from spacepy.toolbox import tOverlapHalf
    if not isinstance(ticks, spt.Ticktock):
        try:
            ticks = spt.Ticktock(ticks, 'UTC')
        except:
            raise TypeError('get_omni: Input times must be a Ticktock object or a list of datetime objects')

    def getattrs(hf, key):
        out = {}
        if hasattr(hf[key],'attrs'):
            for kk, value in hf[key].attrs.iteritems():
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

    with h5.File(omnifln) as hfile:
        keylist = [kk for kk in hfile.keys() if kk not in ['Qbits', 'UTC']]
        st, en = ticks[0].RDT, ticks[-1].RDT
        ##check that requested requested times are within range of data
        if (ticks[0].UTC>omnirange(dbase=dbase)[1]) or (ticks[-1]<omnirange(dbase=dbase)[0]):
            raise ValueError('Requested dates are outside data range')
        inds = tOverlapHalf([st, en], hfile['RDT'], presort=True) #returns an xrange
        sl_op = slice(inds[0]-1, inds[-1]+2)

        omnivals = SpaceData(attrs=getattrs(hfile,'/'))
        fname = hfile.filename
        for key in keylist:
            omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
            omnivals[key].attrs = getattrs(hfile, key)
        for key in hfile['Qbits'].keys():
            omnivals['Qbits<--{0}'.format(key)] = dmarray(hfile['/Qbits/{0}'.format(key)][sl_op])
            omnivals['Qbits<--{0}'.format(key)].attrs = getattrs(hfile, '/Qbits/{0}'.format(key))

    omniout = SpaceData(attrs=dmcopy(omnivals.attrs))
    omniout.attrs['filename'] = fname
    for key in omnivals.keys():
        omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT'], omnivals[key], left=np.NaN, right=np.NaN))
        if 'Qbits' in key:
            #Qbits are integer vals, higher is better, so floor to get best representation of interpolated val
            omniout[key] = np.floor(omniout[key]) 
        #set metadata -- assume this has been set properly in d/l'd file to match ECT-SOC files
        omniout[key].attrs = dmcopy(omnivals[key].attrs)
    omniout['ticks'] = ticks
    omniout['UTC'] = ticks.UTC
    omniout['Hr'] = dmarray([HrFromDT(val) for val in omniout['UTC']])
    omniout['Year'] = dmarray([val.year for val in omniout['UTC']])
    omniout = unflatten(omniout)

    # return warning if values outside of omni data range
    if np.any(np.isnan(omniout['Kp'])):
        print('Warning: some times are outside of the omni data range\nRange: {0}'.format(omnirange(dbase=dbase)))

    return omniout


def omnirange(dbase='QDhourly'):
    '''Returns datetimes giving start and end times in the OMNI/Qin-Denton data

    The update function in toolbox retrieves all available hourly Qin-Denton data, 
    and this function accesses that and looks up the start and end times,
    returning them as datetime objects.

    Parameters
    ==========
    dbase : string (optional)
        name of omni database to check. Currently only 'QDhourly'

    Returns
    =======
    omnirange : tuple
        containing two datetimes giving the start and end times of the available data

    Examples
    ========
    >>> import spacepy.omni as om
    >>> om.omnirange()
    (datetime.datetime(1963, 1, 1, 0, 0), datetime.datetime(2011, 11, 30, 23, 0))

    '''
    
    import h5py as h5
    infile = {'QDhourly': omnifln}
    if dbase not in infile:
        raise NotImplementedError('')
    with h5.File(omnifln) as hfile:
        start, end = hfile['RDT'][0], hfile['RDT'][-1]
        start = spt.Ticktock(start, 'RDT').UTC[0]
        end = spt.Ticktock(end, 'RDT').UTC[0]
    
    return start, end



#-----------------------------------------------

# check for omni file during import
import os, datetime
from spacepy import DOT_FLN, help
from spacepy.toolbox import loadpickle
try:
    import h5py
    _QDext = '.h5'
except ImportError:
    _QDext = '.pkl'

#dotfln = os.environ['HOME']+'/.spacepy'
omnifln = os.path.join(DOT_FLN,'data','omnidata{0}'.format(_QDext))
omni2fln = os.path.join(DOT_FLN,'data','omni2data.pkl')

try: 
    if _QDext=='.h5':
        present = h5py.is_hdf5(omnifln)
    else:
        present = os.path.isfile(omnifln)
except:
    if present:
        print("OMNI data not found in current format. This module has limited functionality.")
    else:
        print("No OMNI data found. This module has limited functionality.")
    print("Run spacepy.toolbox.update(omni=True) to download OMNI data")
