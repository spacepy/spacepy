#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools to read and process omni data

Authors: Josef Koller
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.


"""
import numpy as np
from spacepy.datamodel import SpaceData, dmarray
import spacepy.time as st

__contact__ = 'Josef Koller, jkoller@lanl.gov'

# -----------------------------------------------
def get_omni(ticks):
    """
    Returns Qin-Denton OMNI values, interpolated to any time-base from 1-hourly resolution

    Importing the OMNI module will load the pickled Qin-Denton OMNI file, 
    and this function accesses that and interpolates to the given
    and return the omni values as a SpaceData (dictionary-like) with
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

    """
    if not isinstance(ticks, st.Ticktock):
        try:
            ticks = st.Ticktock(ticks, 'UTC')
        except:
            raise TypeError('get_omni: Input times must be a Ticktock object or a list of datetime objects')

    # extract RTD from ticks
    RDTvals = ticks.RDT
    nRDT = len(ticks)

    omnikeys = omnidata.keys()
    omnikeys.remove('Qbits') # remove this item because it is a string (cannot interpolate)
    omnikeys.remove('UTC')
    omnikeys.remove('RDT')
    omnikeys.remove('ticks')

    omnival = SpaceData(attrs={'URL': 'http://virbo.org/QinDenton'})
    for key in omnikeys:
        omnival[key] = dmarray(np.interp(RDTvals, omnidata['RDT'], omnidata[key], left=np.NaN, right=np.NaN))
        #set attributes/units
        if ('B' in key) or ('Dst' in key):
            omnival[key].attrs['Units'] = 'nT'
        elif 'dens' in key:
            omnival[key].attrs['Units'] = 'cm^{-3}'
        elif 'Pdyn' in key:
            omnival[key].attrs['Units'] = 'nPa'


    # interpolate in Quality bits as well
    omnival['Qbits'] = SpaceData()
    for key in omnidata['Qbits'].keys():
        omnival['Qbits'][key] = dmarray(np.interp(RDTvals, omnidata['RDT'], omnidata['Qbits'][key], \
            left=np.NaN, right=np.NaN))
        #floor interpolation values
        omnival['Qbits'][key] = omnival['Qbits'][key].astype('int')

    # add time information back in
    omnival['UTC'] = ticks.UTC
    omnival['RDT'] = ticks.RDT #TODO:do we need this, since it's already in ticks?
    omnival['ticks'] = ticks

    # return warning if values outside of omni data range
    if np.any(np.isnan(omnival['Kp'])): print("Warning: time is outside of omni data range")

    # add original omni2 to output
    try:
        omni2keys = omni2data.keys()
    except:
        return omnival   #return omnival only and don't add omni2
    
    # remove duplicate keys
    delkeys = ['Epoch', 'RDT', 'Hour', 'Year', 'ticks']
    for key in delkeys: omni2keys.remove(key)
    for key in omni2keys:
        omnival[key] = dmarray(np.interp(RDTvals, omni2data['RDT'], omni2data[key], left=np.NaN, right=np.NaN))
        #set attributes/units
        
        if ('B' in key) or ('Dst' in key):
            omnival[key].attrs['Units'] = 'nT'
        elif 'dens' in key:
            omnival[key].attrs['Units'] = 'cm^{-3}'
        elif 'Pdyn' in key:
            omnival[key].attrs['Units'] = 'nPa'


    return omnival



#-----------------------------------------------

# load omni file during import
from spacepy import DOT_FLN, help, time
from spacepy.toolbox import loadpickle
import os
#dotfln = os.environ['HOME']+'/.spacepy'
omnifln = os.path.join(DOT_FLN,'data','omnidata.pkl')
omni2fln = os.path.join(DOT_FLN,'data','omni2data.pkl')
try:
    omnidata = loadpickle(omnifln)
except:
    print("No OMNI data found. This module has limited functionality.")
    print("Run spacepy.toolbox.update(omni=True) to download OMNI data")
else:
    if not 'ticks' in omnidata:
        omnidata['ticks'] = time.Ticktock(omnidata['UTC'], 'UTC')
    if not 'Hr' in omnidata:
        omnidata['Hr'] = np.fromiter((dt.hour for dt in omnidata['UTC']),
                                     dtype='int16', count=len(omnidata['UTC']))
    if not 'Year' in omnidata:
        omnidata['Year'] = np.fromiter((dt.year for dt in omnidata['UTC']),
                                       dtype='int16',
                                       count=len(omnidata['UTC']))
                                       
try:
    omni2data = loadpickle(omni2fln)
except IOError:
    print("The full OMNI2 dataset has not been not found. This may limit some functionality.")
    print("HINT: Run spacepy.toolbox.update(all=False, omni2=True) to download the full OMNI2 data.")
except ImportError:
    print("ImportError: You may not have a working CDF library installed.")
    print("HINT: Run spacepy.toolbox.update(all=False, omni2=True) to download the full OMNI2 data.")
