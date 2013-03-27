#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools to read and process omni data

Authors: Josef Koller, Steve Morley
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov, smorley@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.


"""
import bisect, re
import numpy as np
from spacepy.datamodel import SpaceData, dmarray, dmcopy, unflatten
from spacepy.toolbox import tOverlapHalf, indsFromXrange
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

    dbase : str (optional)
        Select data source, options are 'QDhourly', 'OMNI2', 'Mergedhourly'

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
                     'Mergedhourly': 3
                     }
    if not dbase in dbase_options: raise NotImplementedError

    import h5py as h5
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

    fname, QDkeylist, O2keylist = '', [], []
    omnivals = SpaceData()
    if dbase_options[dbase] == 1 or dbase_options[dbase] == 3:
        ldb = 'QDhourly'
        with h5.File(omnifln) as hfile:
            QDkeylist = [kk for kk in hfile.keys() if kk not in ['Qbits', 'UTC']]
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
            omnivals.attrs = getattrs(hfile, '/')
            for key in QDkeylist:
                omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
                omnivals[key].attrs = getattrs(hfile, key)
            for key in hfile['Qbits'].keys():
                omnivals['Qbits<--{0}'.format(key)] = dmarray(hfile['/Qbits/{0}'.format(key)][sl_op])
                omnivals['Qbits<--{0}'.format(key)].attrs = getattrs(hfile, '/Qbits/{0}'.format(key))
                QDkeylist.append('Qbits<--{0}'.format(key))

    if dbase_options[dbase] == 2 or dbase_options[dbase] == 3:
        ldb = 'OMNI2hourly'
        with h5.File(omni2fln) as hfile:
            O2keylist = [kk for kk in hfile.keys() if kk not in ['Epoch','RDT']]
            st, en = ticks[0].RDT, ticks[-1].RDT
            ##check that requested requested times are within range of data
            if (ticks[0].UTC>omnirange(dbase=ldb)[1]) or (ticks[-1]<omnirange(dbase=ldb)[0]):
                raise ValueError('Requested dates are outside data range')
            if (ticks[-1].UTC>enval) or (ticks[0]<stval):
                print('Warning: Some requested dates are outside data range ({0})'.format(ldb))
            inds = tOverlapHalf([st, en], hfile['RDT'], presort=True) #returns an xrange
            inds = indsFromXrange(inds)
            if inds[0] < 1: inds[0] = 1
            sl_op = slice(inds[0]-1, inds[-1]+2)
        
            fname = ','.join([fname,hfile.filename])
            omnivals.attrs = getattrs(hfile, '/')
            omnivals['RDT_OMNI'] = dmarray(hfile['RDT'][sl_op])
            for key in O2keylist:
                omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
                omnivals[key].attrs = getattrs(hfile, key)
        

    omniout = SpaceData(attrs=dmcopy(omnivals.attrs))
    omniout.attrs['filename'] = fname[1:]
    ###print('QDkeys: {0}\n\nO2keys: {1}'.format(QDkeylist, O2keylist))
    for key in omnivals.keys():
        if key in O2keylist:
            omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT_OMNI'], omnivals[key], left=np.NaN, right=np.NaN))
            #set metadata -- assume this has been set properly in d/l'd file to match ECT-SOC files
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
        elif key in QDkeylist:
            omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT'], omnivals[key], left=np.NaN, right=np.NaN))
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
        if 'Qbits' in key:
            #Qbits are integer vals, higher is better, so floor to get best representation of interpolated val
            omniout[key] = np.floor(omniout[key]) 
            omniout[key].attrs = dmcopy(omnivals[key].attrs)
    omniout['ticks'] = ticks
    omniout['UTC'] = ticks.UTC
    omniout['Hr'] = dmarray([HrFromDT(val) for val in omniout['UTC']])
    omniout['Year'] = dmarray([val.year for val in omniout['UTC']])
    omniout = unflatten(omniout)

    ## return warning if values outside of omni data range
    #if np.any(np.isnan(omniout['Kp'])):
    #    print('Warning: some times are outside of the omni data range\nRange: {0}'.format(omnirange(dbase=dbase)))

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
              'OMNI2hourly': omni2fln}
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
    _ext = '.h5'
except ImportError:
    _ext = '.pkl'

#dotfln = os.environ['HOME']+'/.spacepy'
omnifln = os.path.join(DOT_FLN,'data','omnidata{0}'.format(_ext))
omni2fln = os.path.join(DOT_FLN,'data','omni2data{0}'.format(_ext))

if _ext=='.h5':
    presentQD = h5py.is_hdf5(omnifln)
    presentO2 = h5py.is_hdf5(omni2fln)
    if not (presentQD and presentO2):
        print("Qin-Denton/OMNI2 data not found in current format. This module has limited functionality.")
        print("Run spacepy.toolbox.update(QDomni=True) to download data")
else:
    presentQD = os.path.isfile(omnifln)
    presentO2 = os.path.isfile(omni2fln)
    if not (presentQD and presentO2):
        print("No Qin-Denton/OMNI2 data found. This module has limited functionality.")
        print("Run spacepy.toolbox.update(QDomni=True) to download data")
    else:
        print("Qin-Denton/OMNI2 data not found in current format. This module has limited functionality.")
        print("Run spacepy.toolbox.update(QDomni=True) to download data")
