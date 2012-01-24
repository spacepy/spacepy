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
    will load the pickled omni file, interpolate to the given ticks time
    and return the omni values as dictionary with
    Kp, Dst, dens, velo, Pdyn, ByIMF, BzIMF, G1, G2, G3, etc.
    (see also http://www.dartmouth.edu/~rdenton/magpar/index.html and
    http://www.agu.org/pubs/crossref/2007/2006SW000296.shtml )

    Note carefully about Qbits: If the status variable is 2, the quantity you are using is fairly well
    determined. If it is 1, the value has some connection to measured values, but is not directly
    measured. These values are still better than just using an average value, but not as good
    as those with the status variable equal to 2. If the status variable is 0, the quantity is
    based on average quantities, and the values listed are no better than an average value. The
    lower the status variable, the less confident you should be in the value.

    Parameters
    ==========
    ticks : Ticktock class
        containing time information

    Returns
    =======
    out : dict
        containing all omni values as a dictionary

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> ticks = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> d = om.get_omni(ticks)
    ['velo', 'Bz6', 'Bz5', 'Bz4', 'Bz3', 'Bz2', 'Bz1', 'RDT', 'Dst',
    'akp3', 'DOY', 'Qbits', 'G3', 'G2', 'G1', 'Hr', 'ticks', 'BzIMF',
    'UTC', 'Kp', 'Pdyn', 'dens', 'ByIMF', 'W6', 'W5', 'W4', 'W3', 'W2',
    'W1', 'Year']


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


    return omnival


# -----------------------------------------------
def get_G123(TAI, omnidata):

    """
    get specific G1, G2, G3 for this TAI
    """

    # define some numbers
    n = 12 # measurement points for averaging
    TAIstart = TAI - 3600
    TAIgrid = linspace(TAIstart, TAI, n)

    # get interpolated values from previous hour
    velogrid = omni.getval('velo', TAIgrid, omnidata)
    Bygrid = omni.getval('ByIMF', TAIgrid, omnidata)
    Bzgrid = omni.getval('BzIMF', TAIgrid, omnidata)
    densgrid = omni.getval('dens', TAIgrid, omnidata)

    # calc. clock angle, Bperp etc
    theta = arctan2(-Bygrid, -Bzgrid)+pi
    Bperp = sqrt(Bygrid**2 + Bzgrid**2)
    hperp = (Bperp/40)**2/(1+Bperp/40)
    a = 0.005
    Bsouth = Bzgrid
    Bsouth[where(Bsouth>0)] = 0
    Bsouth = abs(Bsouth)

    G1 = sum( hperp*velogrid*sin(theta/2)**3 )/n
    G2 = a*sum(velogrid*Bsouth)/n
    G3 = sum(velogrid*densgrid*Bsouth)/n/2000


    return G1, G2, G3


#-----------------------------------------------

# load omni file during import
from spacepy import DOT_FLN, loadpickle, help, time
import os
#dotfln = os.environ['HOME']+'/.spacepy'
omnifln = DOT_FLN+'/data/omnidata.pkl'
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
