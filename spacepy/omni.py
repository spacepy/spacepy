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
from spacepy.datamodel import SpaceData, dmarray, dmcopy
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

    with h5.File(omnifln) as hfile:
        keylist = [kk for kk in hfile.keys() if kk not in ['Qbits', 'UTC']]
        st, en = ticks[0].RDT, ticks[-1].RDT
        inds = tOverlapHalf([st, en], hfile['RDT'], presort=True) #returns an xrange
        sl_op = slice(inds[0], inds[-1]+1)

        omnivals = SpaceData(attrs=getattrs(hfile,'/'))
        fname = hfile.filename
        for key in keylist:
            if 'Qbits' in key: continue #skip these for now
            omnivals[key] = dmarray(hfile[key][sl_op]) #TODO: add attrs from h5
            omnivals[key].attrs = getattrs(hfile, key)

    omniout = SpaceData(attrs=dmcopy(omnivals.attrs))
    omniout.attrs['filename'] = fname
    for key in keylist:
        if 'Qbits' in key: continue #skip these for now, also fix flattened Qbits and kill load on import
        omniout[key] = dmarray(np.interp(ticks.RDT, omnivals['RDT'], omnivals[key], left=np.NaN, right=np.NaN))
        #set metadata -- assume this has been set properly in d/l'd file to match ECT-SOC files
        omniout[key].attrs = dmcopy(omnivals[key].attrs)
    omniout['UTC'] = ticks.UTC
    omniout['ticks'] = ticks

    #TODO: Before returning, interpolate Qbits to the requested timebase.
    return omniout


# # -----------------------------------------------
# def get_omni(ticks):
#     if not isinstance(ticks, spt.Ticktock):
#         try:
#             ticks = spt.Ticktock(ticks, 'UTC')
#         except:
#             raise TypeError('get_omni: Input times must be a Ticktock object or a list of datetime objects')
# 
#     # extract RTD from ticks
#     RDTvals = ticks.RDT
#     nRDT = len(ticks)
# 
#     omnikeys = omnidata.keys()
#     omnikeys.remove('Qbits') # remove this item because it is a string (cannot interpolate)
#     omnikeys.remove('UTC')
#     omnikeys.remove('RDT')
#     omnikeys.remove('ticks')
# 
#     omnival = SpaceData(attrs={'URL': 'http://virbo.org/QinDenton'})
#     for key in omnikeys:
#         omnival[key] = dmarray(np.interp(RDTvals, omnidata['RDT'], omnidata[key], left=np.NaN, right=np.NaN))
#         #set attributes/units
#         if ('B' in key) or ('Dst' in key):
#             omnival[key].attrs['Units'] = 'nT'
#         elif 'dens' in key:
#             omnival[key].attrs['Units'] = 'cm^{-3}'
#         elif 'Pdyn' in key:
#             omnival[key].attrs['Units'] = 'nPa'
# 
# 
#     # interpolate in Quality bits as well
#     omnival['Qbits'] = SpaceData()
#     for key in omnidata['Qbits'].keys():
#         omnival['Qbits'][key] = dmarray(np.interp(RDTvals, omnidata['RDT'], omnidata['Qbits'][key], \
#             left=np.NaN, right=np.NaN))
#         #floor interpolation values
#         omnival['Qbits'][key] = omnival['Qbits'][key].astype('int')
# 
#     # add time information back in
#     omnival['UTC'] = ticks.UTC
#     omnival['RDT'] = ticks.RDT #TODO:do we need this, since it's already in ticks?
#     omnival['ticks'] = ticks
# 
#     # return warning if values outside of omni data range
#     if np.any(np.isnan(omnival['Kp'])): print("Warning: time is outside of omni data range")
# 
#     # add original omni2 to output
#     try:
#         omni2keys = omni2data.keys()
#     except:
#         return omnival   #return omnival only and don't add omni2
#     
#     # remove duplicate keys
#     delkeys = ['Epoch', 'RDT', 'Hour', 'Year', 'ticks']
#     for key in delkeys: omni2keys.remove(key)
#     for key in omni2keys:
#         omnival[key] = dmarray(np.interp(RDTvals, omni2data['RDT'], omni2data[key], left=np.NaN, right=np.NaN))
#         #set attributes/units
#         
#         if ('B' in key) or ('Dst' in key):
#             omnival[key].attrs['Units'] = 'nT'
#         elif 'dens' in key:
#             omnival[key].attrs['Units'] = 'cm^{-3}'
#         elif 'Pdyn' in key:
#             omnival[key].attrs['Units'] = 'nPa'
# 
# 
#     return omnival



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
    print("No OMNI data found. This module has limited functionality.")
    print("Run spacepy.toolbox.update(omni=True) to download OMNI data")
# else:
#     if not 'ticks' in omnidata:
#         omnidata['ticks'] = spt.Ticktock(omnidata['RDT'], 'RDT')
#     if not isinstance(omnidata['UTC'][0], datetime.datetime):
#         omnidata['UTC'] = omnidata['ticks'].UTC
#     if not 'Hr' in omnidata:
#         omnidata['Hr'] = np.fromiter((dt.hour for dt in omnidata['UTC']),
#                                      dtype='int16', count=len(omnidata['UTC']))
#     if not 'Year' in omnidata:
#         omnidata['Year'] = np.fromiter((dt.year for dt in omnidata['UTC']),
#                                        dtype='int16',
#                                        count=len(omnidata['UTC']))

# try:
#     omni2data = loadpickle(omni2fln)
# except IOError:
#     print("The full OMNI2 dataset has not been not found. This may limit some functionality.")
#     print("HINT: Run spacepy.toolbox.update(all=False, omni2=True) to download the full OMNI2 data.")
# except ImportError:
#     print("The OMNI2 dataset file format has changed.")
#     print("Run spacepy.toolbox.update(all=False, omni2=True) to redownload the full OMNI2 data.")
# else:
#     if not 'ticks' in omni2data:
#         omni2data['ticks'] = spt.Ticktock(omni2data['Epoch'])
#     if not 'RDT' in omni2data:
#         omni2data['RDT'] = omni2data['ticks'].RDT
