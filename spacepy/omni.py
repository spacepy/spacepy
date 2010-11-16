#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools to read and process omni data
"""
__version__ = "$Revision: 1.15 $, $Date: 2010/11/16 23:58:48 $"
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'

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

    Input:
    ======
        - ticks (Ticktock class) : containing time information
        
    Returns:
    ========
        - omnival (dictionary) : containing all omni values as a dictionary

    Example:
    ========
    >>> tick = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> d = get_omni(tick)
    
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 26-Jan-2010 (JK)
    V1.1: 11-Mar-2010: fixed bug in get_omni; will now return the correct 6_status, 8_status (JK)
    V1.2: 11-Jun-2010: rewrote status information and put into 'Qbits' (JK)
    """

    import numpy as n
    import spacepy.time as st

    # extract RTD from ticks
    RDTvals = ticks.RDT
    nRDT = len(ticks)

    omnikeys = omnidata.keys()
    omnikeys.remove('Qbits') # remove this item because it is a string (cannot interpolate)
    omnikeys.remove('UTC')
    omnikeys.remove('RDT')
    omnikeys.remove('ticks')

    omnival = {}
    for key in omnikeys:
        omnival[key] = n.interp(RDTvals, omnidata['RDT'], omnidata[key], left=n.NaN, right=n.NaN)
    
    # interpolate in Quality bits as well
    omnival['Qbits'] = {}
    for key in omnidata['Qbits'].keys():
        omnival['Qbits'][key] = n.interp(RDTvals, omnidata['RDT'], omnidata['Qbits'][key], \
            left=n.NaN, right=n.NaN)
        #floor interpolation values
        omnival['Qbits'][key] = omnival['Qbits'][key].astype('int')
        
    # add time information back in
    omnival['UTC'] = ticks.UTC
    omnival['RDT'] = ticks.RDT
    omnival['ticks'] = ticks
    
    # return warning if values outside of omni data range
    if n.any(n.isnan(omnival['Kp'])): print("Warning: time is outside of omni data range")
    
    
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
from spacepy import DOT_FLN, loadpickle, help
import os
#dotfln = os.environ['HOME']+'/.spacepy'
omnifln = DOT_FLN+'/data/omnidata.pkl'
try:
    omnidata = loadpickle(omnifln)
except:
    print("No OMNI data found. This module has limited functionality.")
    print("Run spacepy.toolbox.update(omni=True) to download OMNI data")
