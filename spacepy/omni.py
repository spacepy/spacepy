#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
tools to read and process omni data
"""
__version__ = "$Revision: 1.9 $, $Date: 2010/06/10 17:25:12 $"
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'


# load omni file during import
import os
dotfln = os.environ['HOME']+'/.spacepy'
from spacepy import loadpickle
omnifln = dotfln+'/data/omnidata.pkl'
try:
    omnidata = loadpickle(omnifln)
    __nodataflag = False
except:
    print "No omni data found. This module has limited functionality."
    __nodataflag = True

# -----------------------------------------------
def pickleomni(fln='', overwrite=True, data=None):
    """
    read in fln='omni_intp.dat'  file from Qin and Denton
    and save as pickle. The pickle will replace the current
    data file in the package, or provide data as strings

    Input:
    ======
        - fln (string) : filename of the original ASCII file
        - overwrite (optional boolean) : If true -> overwrite data file in package
            if false -> save in current working directory
        - data (string) : container with data read-in from file

    Example:
    ========
    >>> pickleomni('omnidata.dat')
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov 
    """

    import re
    import numpy as n
    import datetime
    import spacepy
    import spacepy.time as st

    # open and read the complete file
    if data == None:
        fh = open(fln)
        A = fh.readlines()
        fh.close()
    else:
        A = data

    # through away header
    k = 0
    while re.search(' 1963', A[k]) == None: k += 1

    # create a keylist
    keys = A[k-1].split()
    keys.remove('8')
    keys.remove('6')
    keys[keys.index('status')] = '8_status'
    keys[keys.index('stat')] = '6_status'
    keys[keys.index('dst')] = 'Dst'
    keys[keys.index('kp')] = 'Kp'
    #keys[keys.index('Hr')] = 'Hr'
    keys[keys.index('V_SW')] = 'velo'
    keys[keys.index('Den_P')] = 'dens'
    keys[keys.index('Day')] = 'DOY'
    keys[keys.index('Year')] = 'Year'
    
    # put it into a 2D table
    tab = n.zeros((len(A)-k,len(keys)))
    stat8 = ['']*(len(A)-k)
    stat6 = ['']*(len(A)-k)
    for i in n.arange(len(A)-k):
        tab[i,:] = A[k+i].split()
        stat8[i] = A[k+i].split()[11]
        stat6[i] = A[k+i].split()[27]

    tab = n.reshape(tab, n.shape(tab))
    # take out everything less than 1950 (some lines had a zero)
    idx = n.where(tab[:,0]>1950)[0]
    tab = tab[idx,:]
    stat8 = n.array(stat8)[idx]
    stat6 = n.array(stat6)[idx]
    
    omnidata = {} 
    # sort through and make an omni dictionary
    # extract keys from line above
    for ikey, i  in zip(keys,range(len(keys))):
        omnidata[ikey] = tab[:,i]
    
    omnidata['6_status'] = stat6
    omnidata['8_status'] = stat8

    # add TAI to omnidata
    nTAI = len(omnidata['DOY'])
    omnidata['UTC'] = ['']*nTAI
    omnidata['RDT'] = n.zeros(nTAI)
    
    # add time information to omni pickle (long loop)
    for i in range(nTAI):
        year = int(omnidata['Year'][i])
        doy = int(omnidata['DOY'][i])
        month, day = st.doy2date(year,doy)
        UT_hr = omnidata['Hr'][i]
        hour, minute = divmod(UT_hr*60., 60)
        minute, second = divmod(minute*60., 60)  
        omnidata['UTC'][i] = datetime.datetime(year, month, day, int(hour), int(minute), int(second))
    
    omnidata['ticktock'] = st.Ticktock(omnidata['UTC'], 'UTC')
    omnidata['RDT'] = omnidata['ticktock'].RDT
    
    # save as pickle
    if overwrite == True:
        spacepy.savepickle(omnifln, omnidata)
    else:
        # save in current working directory
        spacepy.savepickle('omnidata.pkl', omnidata)

    return 

# -----------------------------------------------
def get_omni(ticktock):
    """
    will load the pickled omni file, interpolate to the given ticktock time
    and return the omni values as dictionary with 
    Kp, Dst, dens, velo, Pdyn, ByIMF, BzIMF, G1, G2, G3, etc.
    (see also http://www.dartmouth.edu/~rdenton/magpar/index.html and
    http://www.agu.org/pubs/crossref/2007/2006SW000296.shtml )
    
    Note carefully: If the status variable is 2, the quantity you are using is fairly well 
    determined. If it is 1, the value has some connection to measured values, but is not directly 
    measured. These values are still better than just using an average value, but not as good 
    as those with the status variable equal to 2. If the status variable is 0, the quantity is 
    based on average quantities, and the values listed are no better than an average value. The 
    lower the status variable, the less confident you should be in the value.

    Input:
    ======
        - ticktock (Ticktock class) : containing time information
        
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
    V1.2: 10-Jun-2010: took out 6_status etc b/c it's too slow (JK)
    """

    import numpy as n

    # extract RTD from ticktock
    RDTvals = ticktock.RDT
    nRDT = len(ticktock)

    omnikeys = omnidata.keys()
    omnikeys.remove('6_status') # remove this item because it is a string (cannot interpolate)
    omnikeys.remove('8_status')
    omnikeys.remove('UTC')
    omnikeys.remove('ticktock')

    omnival = {}
    for key in omnikeys:
        omnival[key] = n.zeros(nRDT)
        omnival[key] = n.interp(RDTvals, omnidata['RDT'], omnidata[key], left=n.NaN, right=n.NaN)
        
    # add time information back in
    omnival['UTC'] = ticktock.UTC
    # add interpolation parameters back in
    #for key in ['6_status','8_status']:
    #    omnival[key] = ['']*nRDT
    #    for iRDT, RDT in zip( n.arange(nRDT), RDTvals):
    #        idx = n.argmin( abs(omnidata['RDT']-RDT) )
    #        omnival[key][iRDT] = omnidata[key][idx]

    # return warning if values outside of omni data range
    if n.any(n.isnan(omnival['Kp'])): print "Warning: time is outside of omni data range"
    
    
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
#Test whether data file exists in correct location, if not, offer to fetch.
if __nodataflag:
    import spacepy.toolbox as tb
    ans = tb.query_yes_no('\nUpdate OMNI from ViRBO now? (Internet connection required) ', default="yes")
    if ans=='yes':
        omni_fname_zip = dotfln+'/data/WGhour-latest.d.zip'
        omni_fname_dat = dotfln+'/data/omnidata.pkl'
        tb.update(all=False, omni=True, callfromOMNI=True)
        pickleomni(fln=omni_fname_dat)
        # delete left-overs
        os.remove(omni_fname_zip)
        omnidata = loadpickle(omnifln)
    else:
        #raise ImportError("No OMNI data downloaded. This module has limited functionality.")
        print "\nNo OMNI data downloaded. This module has limited functionality."
   
