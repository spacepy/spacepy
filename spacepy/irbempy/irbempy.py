#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
module wrapper for irbem_lib
Reference for this library
https://sourceforge.net/projects/irbem/
D. Boscher, S. Bourdarie, T. P. O'Brien, T. Guild, IRBEM library V4.3, 2004-2008


Authors
-------
Josef Koller, Steve Morley 

Copyright 2010 Los Alamos National Security, LLC.
"""

import itertools, numbers
try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable
import sys, os
import warnings

import numpy as np

from spacepy import help, config
import spacepy.coordinates as spc
import spacepy.datamodel as dm
try:
    from . import irbempylib as oplib
except ImportError:
    if 'sphinx' in sys.argv[0]:
        warnings.warn('Could not import irbempylib. '
                      'You appear to be building docs, so ignoring this error.')
    else:
        raise
import spacepy.toolbox as tb

#check whether TS07_DATA_PATH is set, if not then set to spacepy's installed data directory
if 'TS07_DATA_PATH' not in os.environ:
    import pkg_resources
    dataTS07D = os.path.join('data', 'TS07D')
    spdatapath = pkg_resources.resource_filename('spacepy', dataTS07D)
    os.environ['TS07_DATA_PATH'] = spdatapath #set environment variable here


SYSAXES_TYPES = {'GDZ': {'sph': 0, 'car': None},
    'GEO': {'sph': None, 'car': 1}, 'GSM': {'sph': None, 'car': 2},
    'GSE': {'sph': None, 'car': 3}, 'SM': {'sph': None, 'car': 4},
    'GEI': {'sph': None, 'car': 5}, 'MAG': {'sph': None, 'car': 6},
    'SPH': {'sph': 7, 'car': None}, 'RLL': {'sph': 8, 'car': None},
    'TOD': {'sph': None, 'car': 12}, 'J2000': {'sph': None, 'car': 13}}

# -----------------------------------------------
def updateTS07Coeffs(path=None, force=False, verbose=False, **kwargs):
    '''Update coefficients for TS07 magnetic field model

    Parameters
    ----------
    path : str
        Base path for TS07 dynamic coefficients. If None (default)
        then the TS07_DATA_PATH environment variable is used.
    force : boolean
        If True (default is False) then missing paths will be created.
    verbose : boolean
        Print verbose output. Default is False.
    start : datetime.datetime or string
        Start time for archive retrieval. If start is a string then it should
        ISO8601 format (YYYY-MM-DDTHH:mm:SS). Defaults to 2007-01-01.
    end : datetime.datetime or string
        End time for archive retrieval. Required format same as start.
        Defaults to time of query (i.e. latest available).
    '''
    import glob, tarfile
    import spacepy.time as spt
    dt = spt.datetime
    if sys.version_info[0]<3:
        import urllib as u
    else:
        import urllib.request as u

    if 'user_agent' in config and config['user_agent']:
        class AppURLopener(u.FancyURLopener):
            version = config['user_agent']
        u._urlopener = AppURLopener()

    if not path:
        #test for TS07_DATA_PATH
        if 'TS07_DATA_PATH' in os.environ:
            path = os.environ['TS07_DATA_PATH']
        else:
            err_str = 'updateTS07Coeffs: Directory for TS07 data must be specified by \n' + \
                  'TS07_DATA_PATH environment variable or path keyword.'
            raise ValueError(err_str)

    #test that path exists, unless force keyword is set
    if not os.path.isdir(path):
        if force:
            os.makedirs(path, mode=0o777)
        else:
            raise IOError('updateTS07Coeffs: Specified path for TS07 is not valid - to force creation use "force" keyword')
    
    if 'start' not in kwargs:
        firstDate = dt.datetime(2007,1,1)
    else:
        firstDate = spt.Ticktock(kwargs['start']).UTC[0] #make sure it's a datetime
    if 'end' not in kwargs:
        lastDate = dt.datetime.now()
    else:
        lastDate = spt.Ticktock(kwargs['end']).UTC[0]
        
    baseURL = 'http://rbspgway.jhuapl.edu/models/magneticfieldmodeling/ts07d/coeffs_v02/'

    #make inventory of current TS07 data
    current_days = sorted(glob.glob(os.path.join(path,'Coeffs','*')))
    current_days = [os.path.split(dd)[-1] for dd in current_days if 'tgz' not in dd]
    current_dates = []
    for vals in current_days:
        try:
            adddate = spt.doy2date(int(vals.split('_')[0]),int(vals.split('_')[1]), dtobj=True)
        except ValueError: #not a valid date format
            continue
        current_dates.append(adddate)

    request_days = spt.tickrange(firstDate, lastDate, 1)
    timewant = set(request_days.UTC)
    timegot = set(current_dates)
    stillwant = sorted(list(timewant.difference(timegot)))

    print("Retrieving {2} requested files [{0} - {1}]".format(firstDate.date().isoformat(), lastDate.date().isoformat(), len(stillwant)))
    for day in stillwant:
        doy = spt.Ticktock(day).DOY[0]
        name = '{0}_{1:03d}'.format(day.year, doy)
        target = os.path.join(path, 'Coeffs', name)
        source = baseURL+'/{0}/'.format(day.year)+'{0}.tgz'.format(name)
        targetfile = os.path.join('{0}.tgz'.format(name))
        if verbose:
            print("\nFetching archive for {0}".format(day.date().isoformat()))
            tmp, report = u.urlretrieve(source, targetfile, reporthook=tb.progressbar)
        else:
            tmp, report = u.urlretrieve(source, targetfile)
        #if the file existed on the server then we won't get an HTML message back
        if 'html' not in report['content-type']:
            tar = tarfile.open(tmp, 'r:gz')
            tar.extractall(target)
            tar.close()
        #now remove tgz file
        os.unlink(tmp)
        

# -----------------------------------------------
def get_Bfield(ticks, loci, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
    """
    call get_bfield in irbem lib and return a dictionary with the B-field vector and 
    strenght.

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table
        - (see Lstar documentation for further explanation)

    Returns
    =======
        - results (dictionary) : containing keys: Bvec, and Blocal

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.coordinates as spc
    >>> import spacepy.irbempy as ib
    >>> t = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = spc.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> ib.get_Bfield(t,y)
    {'Blocal': array([  976.42565251,  3396.25991675]),
       'Bvec': array([[ -5.01738885e-01,  -1.65104338e+02,   9.62365503e+02],
       [  3.33497974e+02,  -5.42111173e+02,   3.33608693e+03]])}

    See Also
    ========
    get_Lstar, find_Bmirror, find_magequator

    """

    # prepare input values for irbem call
    d = prep_irbem(ticks, loci, alpha=[], extMag=extMag, options=options, omnivals=omnivals)
    nTAI = len(ticks)
    badval = d['badval']
    kext = d['kext']
    sysaxes = d['sysaxes']
    iyearsat = d['iyearsat']
    idoysat = d['idoysat']
    secs = d['utsat']
    utsat = d['utsat']
    xin1 = d['xin1']
    xin2 = d['xin2']
    xin3 = d['xin3']
    magin = d['magin']

    results = dm.SpaceData()
    results['Blocal'] = np.zeros(nTAI)
    results['Bvec'] = np.zeros((nTAI,3))
    for i in np.arange(nTAI):
        BxyzGEO, Blocal = oplib.get_field1(kext,options,sysaxes,iyearsat[i],idoysat[i],secs[i], \
            xin1[i],xin2[i],xin3[i], magin[:,i])

        # take out all the odd 'bad values' and turn them into NaN
        if np.isclose(Blocal,badval): Blocal = np.NaN
        BxyzGEO[np.where( np.isclose(BxyzGEO, badval)) ] = np.NaN

        results['Blocal'][i] = Blocal
        results['Bvec'][i,:] = BxyzGEO

    return results

# -----------------------------------------------
def find_Bmirror(ticks, loci, alpha, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
    """
    call find_mirror_point from irbem library and return a dictionary with values for 
    Blocal, Bmirr and the GEO (cartesian) coordinates of the mirror point

    Parameters
    ==========
    ticks : Ticktock class
        containing time information
    loci : Coords class
        containing spatial information
    alpha : array-like
        containing the pitch angles
    extMag : str
        optional; will choose the external magnetic field model 
        possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
        OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
        'T05', 'ALEX', 'TS07']
    options : array-like (optional)
        length=5 : explained in Lstar
    omnivals : dict (optional)
        if not provided, will use lookup table 
        (see get_Lstar documentation for further explanation)

    Returns
    =======
    results : dictionary
        containing keys: Blocal, Bmirr, GEOcar

    Examples
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> ib.find_Bmirror(t,y,[90,80,60,10])
    {'Blocal': array([ 0.,  0.]),
     'Bmirr': array([ 0.,  0.]),
     'loci': Coords( [[ NaN  NaN  NaN]
     [ NaN  NaN  NaN]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']}

    See Also
    ========
    get_Lstar, get_Bfield, find_magequator
    """

    # prepare input values for irbem call
    d = prep_irbem(ticks, loci, alpha, extMag=extMag, options=options, omnivals=omnivals)
    nTAI = len(ticks)
    badval = d['badval']
    kext = d['kext']
    sysaxes = d['sysaxes']
    iyearsat = d['iyearsat']
    idoysat = d['idoysat']
    secs = d['utsat']
    xin1 = d['xin1']
    xin2 = d['xin2']
    xin3 = d['xin3']
    magin = d['magin']
    #degalpha = d['degalpha']
    #nalp_max = d['nalp_max']
    #nalpha = d['nalpha']

    results = dm.SpaceData()
    results['Blocal'] = np.zeros(nTAI)
    results['Bmirr'] = np.zeros(nTAI)
    results['loci'] = ['']*nTAI
    for i in np.arange(nTAI):
        blocal, bmirr, GEOcoord = oplib.find_mirror_point1(kext,options,sysaxes, \
            iyearsat[i],idoysat[i],secs[i], xin1[i],xin2[i],xin3[i], alpha, magin[:,i])

        # take out all the odd 'bad values' and turn them into NaN
        if np.isclose(blocal,badval): blocal = np.NaN
        if np.isclose(bmirr,badval) : bmirr  = np.NaN
        GEOcoord[np.where( np.isclose(GEOcoord,badval)) ] = np.NaN

        results['Blocal'][i] = blocal
        results['Bmirr'][i] = bmirr	
        results['loci'][i] = GEOcoord

    results['loci'] = spc.Coords(results['loci'], 'GEO', 'car')

    return results


# -----------------------------------------------
def find_magequator(ticks, loci, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
    """
    call find_magequator from irbem library and return a dictionary with values for
    Bmin and the GEO (cartesian) coordinates of the magnetic equator

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - (see Lstar documentation for further explanation)

    Returns
    =======
        - results (dictionary) : containing keys: Bmin, Coords instance with GEO coordinates of 
            the magnetic equator

    Examples
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> op.find_magequator(t,y)
    {'Bmin': array([  945.63652413,  3373.64496167]),
        'loci': Coords( [[ 2.99938371  0.00534151 -0.03213603]
        [ 2.00298822 -0.0073077   0.04584859]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']}

    See Also
    ========
    get_Lstar, get_Bfield, find_Bmirr
    """
    # prepare input values for irbem call
    d = prep_irbem(ticks, loci, alpha=[], extMag=extMag, options=options, omnivals=omnivals)
    nTAI = len(ticks)
    badval = d['badval']
    kext = d['kext']
    sysaxes = d['sysaxes']
    iyearsat = d['iyearsat']
    idoysat = d['idoysat']
    secs = d['utsat']
    xin1 = d['xin1']
    xin2 = d['xin2']
    xin3 = d['xin3']
    magin = d['magin']

    results = {}
    results['Bmin'] = np.zeros(nTAI)
    results['loci'] = ['']*nTAI
    for i in np.arange(nTAI):
        bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
            iyearsat[i],idoysat[i],secs[i], xin1[i],xin2[i],xin3[i],magin[:,i])

        # take out all the odd 'bad values' and turn them into NaN
        if np.isclose(bmin,badval): bmin = np.NaN
        GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN

        results['Bmin'][i] = bmin
        results['loci'][i] = GEOcoord

    results['loci'] = spc.Coords(results['loci'], 'GEO', 'car')

    return results


# -----------------------------------------------
def find_LCDS(ticks, alpha, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None, tol=0.05, bracket=[3,12], mlt=0, **kwargs):
    """
    Find the last closed drift shell (LCDS) for a given equatorial pitch angle.

    Uses the IRBEM library to determine L* and searches via bisection to find LCDS
    to a given tolerance in R (radial distance along GSM equator at local midnight).

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information **for a single time**
        - alpha (numeric) : equatorial pitch angle for search
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - (see Lstar documentation for further explanation)
        - tol (float) : tolerance for search in radial distance
        - bracket (list): X-GSM coordinates to bracket bisection search

    Returns
    =======
        results (SpaceData, dictionary-like): contains keys LCDS, K, AlphaEq and UTC

    Examples
    ========
    >>> t = spacepy.time.Ticktock(['2002-02-02T12:00:00'], 'ISO')
    >>> spacepy.irbempy.find_LCDS(t, 90, extMag='T89')

    See Also
    ========
    get_Lstar, get_Bfield, find_Bmirr
    """
    # prepare input values for irbem call
    nTAI = len(ticks)

    #First set inner bracket (default to R of 3)
    try:
        assert len(bracket)==2
    except:
        raise ValueError('Specified initial bracket is invalid')

    if not isinstance(alpha, Iterable): alpha = [alpha]
    results = dm.SpaceData()
    if len(alpha)!=1:
        results['LCDS'] = dm.dmfilled([nTAI,len(alpha)])
        results['K'] = dm.dmfilled([nTAI,len(alpha)])
    else:
        results['LCDS'] = dm.dmfilled([nTAI,])
        results['K'] = dm.dmfilled([nTAI,])
    results['LCDS'].attrs['DESCRIPTION'] = "Last closed drift shell calculated with SpacePy's irbempy module"
    results['LCDS'].attrs['UNITS'] = "dimensionless"
    results['LCDS'].attrs['DEPEND_0'] = "UTC"
    results['LCDS'].attrs['DEPEND_1'] = "K"
    results['LCDS'].attrs['MODEL'] = "{0} (IRBEM implementation)".format(extMag)
    results['K'].attrs['DESCRIPTION'] = "Modified 2nd adiabatic invariant at LCDS"
    results['K'].attrs['UNITS'] = "R_E.G^{1/2}"
    results['K'].attrs['DEPEND_0'] = "UTC"
    results['K'].attrs['DEPEND_1'] = "AlphaEq"
    results['K'].attrs['MODEL'] = "{0} (IRBEM implementation)".format(extMag)
    results['UTC'] = ticks.UTC
    results['AlphaEq'] = dm.dmarray(alpha, 
                                    attrs={'DESCRIPTION': 'Equatorial pitch angle for LCDS calculation',
                                           'UNITS': 'degrees'})
    
    mlt *= 15 #hours to degrees
    mlt = np.deg2rad(mlt)

    for idxt, tt in enumerate(ticks):
        if not omnivals:
            #prep_irbem will get omni if not specified, but to save on repeated calls, do it once here
            import spacepy.omni as omni
            omnivals = omni.get_omni(tt)
        for idxa, pa in enumerate(alpha):
            b1x = -1.0*bracket[0]*np.cos(mlt)
            b1y = -1.0*bracket[0]*np.sin(mlt)
            loci_brac1 = spc.Coords([b1x,b1y,0], 'GSM', 'car')
            if 'verbose' in kwargs: print('Initial inner bracket: {0}'.format(loci_brac1))

            d = prep_irbem(tt, loci_brac1, alpha=[pa], extMag=extMag, options=options, omnivals=omnivals)
            badval = d['badval']
            kext = d['kext']
            sysaxes = d['sysaxes']
            iyearsat = d['iyearsat']
            idoysat = d['idoysat']
            secs = d['utsat']
            xin1 = d['xin1']
            xin2 = d['xin2']
            xin3 = d['xin3']
            magin = d['magin']
            nTtoG = 1.0e-5

            bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])

            # take out all the odd 'bad values' and turn them into NaN
            if np.isclose(bmin,badval): bmin = np.NaN
            GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
            #Now get Lstar at this location...
            if GEOcoord[0]!=np.NaN:
                pos1 = spc.Coords(GEOcoord, 'GEO', 'car')
                LS1 = get_Lstar(tt, pos1, pa, extMag=extMag, options=options, omnivals=omnivals)
                if np.isnan(LS1['Lstar']).any():
                    raise ValueError('Specified inner bracket ({0}) is on an open drift shell'.format(loci_brac1))
                LS1['K'] = LS1['Xj']*np.sqrt(LS1['Bmirr']*nTtoG)
            else:
                raise ValueError('Specified inner bracket ({0}) is on an open drift shell'.format(loci_brac1))
            #print('L* at inner bracket: {0}'.format(LS1['Lstar']))
            LCDS, LCDS_K = LS1['Lstar'][0], LS1['K'][0]
                
            #Set outer bracket (default to R of 12)
            b2x = -1.0*bracket[1]*np.cos(mlt)
            b2y = -1.0*bracket[1]*np.sin(mlt)
            loci_brac2 = spc.Coords([b2x,b2y,0], 'GSM', 'car')
            if 'verbose' in kwargs: print('Initial outer bracket: {0}'.format(loci_brac2))

            d2 = prep_irbem(tt, loci_brac2, alpha=[pa], extMag=extMag, options=options, omnivals=omnivals)
            badval = d2['badval']
            kext = d2['kext']
            sysaxes = d2['sysaxes']
            iyearsat = d2['iyearsat']
            idoysat = d2['idoysat']
            secs = d2['utsat']
            xin1 = d2['xin1']
            xin2 = d2['xin2']
            xin3 = d2['xin3']
            magin = d2['magin']

            bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])

            # take out all the odd 'bad values' and turn them into NaN
            if np.isclose(bmin,badval): bmin = np.NaN
            GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
            #Now get Lstar at this location...
            if GEOcoord[0]!=np.NaN:
                pos2 = spc.Coords(GEOcoord, 'GEO', 'car')
                LS2 = get_Lstar(tt, pos2, pa, extMag=extMag, options=options, omnivals=omnivals)
                if not np.isnan(LS2['Lstar']).any():
                    raise ValueError('Specified outer bracket is on a closed drift shell')
            else:
                LS2 = {'Lstar': np.NaN}
            #print('L* at outer bracket: {0}; Xgsm = {1}'.format(LS2['Lstar'], loci_brac2.x))
            
            #now search by bisection
            while (tb.hypot(loci_brac2.x, loci_brac2.y) - tb.hypot(loci_brac1.x, loci_brac1.y) > tol):
                newdist = (tb.hypot(loci_brac2.x, loci_brac2.y) + tb.hypot(loci_brac1.x, loci_brac1.y))/2.0
                newx = -1.0*newdist*np.cos(mlt)
                newy = -1.0*newdist*np.sin(mlt)
                pos_test = spc.Coords([newx, newy, 0], 'GSM', 'car')

                dtest = prep_irbem(tt, pos_test, alpha=[pa], extMag=extMag, options=options, omnivals=omnivals)
                badval = dtest['badval']
                kext = dtest['kext']
                sysaxes = dtest['sysaxes']
                iyearsat = dtest['iyearsat']
                idoysat = dtest['idoysat']
                secs = dtest['utsat']
                xin1 = dtest['xin1']
                xin2 = dtest['xin2']
                xin3 = dtest['xin3']
                magin = dtest['magin']

                bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                    iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])
                # take out all the odd 'bad values' and turn them into NaN
                if np.isclose(bmin,badval): bmin = np.NaN
                GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
                #print('bmin, GEOcoord = {0},{1}'.format(bmin, GEOcoord))
                if not (np.isnan(bmin)):
                    #Now get Lstar at this location...
                    postry = spc.Coords(GEOcoord, 'GEO', 'car')
                    LStry = get_Lstar(tt, postry, pa, extMag=extMag, options=options, omnivals=omnivals)
                    LStry['K'] = LStry['Xj']*np.sqrt(LStry['Bmirr']*nTtoG)
                else:
                    LStry = {'Lstar': np.NaN, 'K': np.NaN}
                if 'verbose' in kwargs: print('L* at test point: {0}; Xgsm = {1}'.format(LStry['Lstar'], pos_test))

                if np.isnan(LStry['Lstar']).any():
                    loci_brac2 = pos_test
                else:
                    loci_brac1 = pos_test
                    LCDS, LCDS_K = LStry['Lstar'][0], LStry['K'][0]
            try:
                results['LCDS'][idxt, idxa] = LCDS
                results['K'][idxt, idxa] = LCDS_K
            except IndexError:
                results['LCDS'][idxt] = LCDS
                results['K'][idxt] = LCDS_K

    results['LCDS'] = results['LCDS']
    results['K'] = results['K']

    return results

# -----------------------------------------------
def find_LCDS_K(ticks, K, extMag='T01STORM', options=[1,1,3,0,0], omnivals=None, tol=0.05, bracket=[3,12], mlt=0, **kwargs):
    """
    Find the last closed drift shell (LCDS) for a given equatorial pitch angle.

    Uses the IRBEM library to determine L* and searches via bisection to find LCDS
    to a given tolerance in R (radial distance along GSM equator at local midnight).

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information **for a single time**
        - alpha (numeric) : equatorial pitch angle for search
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - (see Lstar documentation for further explanation)
        - tol (float) : tolerance for search in radial distance
        - bracket (list): X-GSM coordinates to bracket bisection search

    Returns
    =======
        results (SpaceData, dictionary-like): contains keys LCDS, K, AlphaEq and UTC

    Examples
    ========
    >>> t = spacepy.time.Ticktock(['2002-02-02T12:00:00'], 'ISO')
    >>> spacepy.irbempy.find_LCDS(t, 90, extMag='T89')

    See Also
    ========
    get_Lstar, get_Bfield, find_Bmirr
    """
    # prepare input values for irbem call
    nTAI = len(ticks)
    Aopt = [0]
    Aopt.extend(options[1:])

    #First set inner bracket (default to R of 3)
    try:
        assert len(bracket)==2
    except:
        raise ValueError('Specified initial bracket is invalid')

    if not isinstance(K, Iterable): K = [K]
    results = dm.SpaceData()
    if len(K)!=1:
        results['LCDS'] = dm.dmfilled([nTAI,len(K)])
        results['AlphaEq'] = dm.dmfilled([nTAI,len(K)])
        results['Success'] = dm.dmfilled([nTAI,len(K)], fillval='Success', dtype='|S24')
    else:
        results['LCDS'] = dm.dmfilled([nTAI,])
        results['AlphaEq'] = dm.dmfilled([nTAI,])
    results['LCDS'].attrs['DESCRIPTION'] = "Last closed drift shell calculated with SpacePy's irbempy module"
    results['LCDS'].attrs['UNITS'] = "dimensionless"
    results['LCDS'].attrs['DEPEND_0'] = "UTC"
    results['LCDS'].attrs['DEPEND_1'] = "K"
    results['LCDS'].attrs['MODEL'] = "{0} (IRBEM implementation)".format(extMag)
    results['AlphaEq'].attrs['DESCRIPTION'] = "Modified 2nd adiabatic invariant at LCDS"
    results['AlphaEq'].attrs['UNITS'] = "R_E.G^{1/2}"
    results['AlphaEq'].attrs['MODEL'] = "{0} (IRBEM implementation)".format(extMag)
    results['UTC'] = ticks.UTC
    results['K'] = dm.dmarray(K, attrs={'DESCRIPTION': 'Equatorial pitch angle for LCDS calculation',
                                           'UNITS': 'degrees'})
    
    mlt *= 15 #hours to degrees
    mlt = np.deg2rad(mlt)

    for idxt, tt in enumerate(ticks):
        if not omnivals:
            #prep_irbem will get omni if not specified, but to save on repeated calls, do it once here
            import spacepy.omni as omni
            omnivals = omni.get_omni(tt)
        for idxk, k_i in enumerate(K):
            b1x = -1.0*bracket[0]*np.cos(mlt)
            b1y = -1.0*bracket[0]*np.sin(mlt)
            loci_brac1 = spc.Coords([b1x,b1y,0], 'GSM', 'car')
            if 'verbose' in kwargs: print('Initial inner bracket: {0}'.format(loci_brac1))

            d = prep_irbem(tt, loci_brac1, extMag=extMag, options=options, omnivals=omnivals)
            badval = d['badval']
            kext = d['kext']
            sysaxes = d['sysaxes']
            iyearsat = d['iyearsat']
            idoysat = d['idoysat']
            secs = d['utsat']
            xin1 = d['xin1']
            xin2 = d['xin2']
            xin3 = d['xin3']
            magin = d['magin']
            nTtoG = 1.0e-5
            pa = np.NaN

            bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])

            # take out all the odd 'bad values' and turn them into NaN
            if np.isclose(bmin,badval): bmin = np.NaN
            GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
            #Now get Lstar at this location...
            if np.isfinite(GEOcoord[0]):
                pos1 = spc.Coords(GEOcoord, 'GEO', 'car')
                pa = AlphaOfK(tt, pos1, k_i, extMag=extMag, options=Aopt, omnivals=omnivals)
                LS1 = get_Lstar(tt, pos1, pa, extMag=extMag, options=options, omnivals=omnivals)
                if np.isnan(LS1['Lstar']).any():
                    try:
                        results['LCDS'][idxt, idxk] = np.NaN
                        results['AlphaEq'][idxt, idxk] = pa
                        results['Success'][idxt, idxk] = 'Invalid inner bracket'
                    except IndexError:
                        results['LCDS'][idxt] = np.NaN
                        results['AlphaEq'][idxt] = pa
                        results['Success'][idxt] = 'Invalid inner bracket'
                    continue
                LS1['K'] = LS1['Xj']*np.sqrt(LS1['Bmirr']*nTtoG)
            else:
                try:
                    results['LCDS'][idxt, idxk] = np.NaN
                    results['AlphaEq'][idxt, idxk] = pa
                    results['Success'][idxt, idxk] = 'Invalid inner bracket'
                except IndexError:
                    results['LCDS'][idxt] = np.NaN
                    results['AlphaEq'][idxt] = pa
                    results['Success'][idxt] = 'Invalid inner bracket'
                continue

            #print('L* at inner bracket: {0}'.format(LS1['Lstar']))
            LCDS, LCDS_PA = LS1['Lstar'][0], pa
                
            #Set outer bracket (default to R of 12)
            b2x = -1.0*bracket[1]*np.cos(mlt)
            b2y = -1.0*bracket[1]*np.sin(mlt)
            loci_brac2 = spc.Coords([b2x,b2y,0], 'GSM', 'car')
            if 'verbose' in kwargs: print('Initial outer bracket: {0}'.format(loci_brac2))

            d2 = prep_irbem(tt, loci_brac2, extMag=extMag, options=options, omnivals=omnivals)
            badval = d2['badval']
            kext = d2['kext']
            sysaxes = d2['sysaxes']
            iyearsat = d2['iyearsat']
            idoysat = d2['idoysat']
            secs = d2['utsat']
            xin1 = d2['xin1']
            xin2 = d2['xin2']
            xin3 = d2['xin3']
            magin = d2['magin']

            bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])

            # take out all the odd 'bad values' and turn them into NaN
            if np.isclose(bmin,badval): bmin = np.NaN
            GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
            #Now get Lstar at this location...
            if np.isfinite(GEOcoord[0]):
                pos2 = spc.Coords(GEOcoord, 'GEO', 'car')
                pa = AlphaOfK(tt, pos2, k_i, extMag=extMag, options=Aopt, omnivals=omnivals)
                LS2 = get_Lstar(tt, pos2, pa, extMag=extMag, options=options, omnivals=omnivals)
                if not np.isnan(LS2['Lstar']).any():
                    try:
                        results['LCDS'][idxt, idxk] = np.NaN
                        results['AlphaEq'][idxt, idxk] = pa
                        results['Success'][idxt, idxk] = 'Invalid outer bracket'
                    except IndexError:
                        results['LCDS'][idxt] = np.NaN
                        results['AlphaEq'][idxt] = pa
                        results['Success'][idxt] = 'Invalid outer bracket'
                    continue
            else:
                LS2 = {'Lstar': np.NaN}
            #print('L* at outer bracket: {0}; Xgsm = {1}'.format(LS2['Lstar'], loci_brac2.x))
            
            #now search by bisection
            while (tb.hypot(loci_brac2.x, loci_brac2.y) - tb.hypot(loci_brac1.x, loci_brac1.y) > tol):
                newdist = (tb.hypot(loci_brac2.x, loci_brac2.y) + tb.hypot(loci_brac1.x, loci_brac1.y))/2.0
                newx = -1.0*newdist*np.cos(mlt)
                newy = -1.0*newdist*np.sin(mlt)
                pos_test = spc.Coords([newx, newy, 0], 'GSM', 'car')

                dtest = prep_irbem(tt, pos_test, extMag=extMag, options=options, omnivals=omnivals)
                badval = dtest['badval']
                kext = dtest['kext']
                sysaxes = dtest['sysaxes']
                iyearsat = dtest['iyearsat']
                idoysat = dtest['idoysat']
                secs = dtest['utsat']
                xin1 = dtest['xin1']
                xin2 = dtest['xin2']
                xin3 = dtest['xin3']
                magin = dtest['magin']

                bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                    iyearsat[0],idoysat[0],secs[0], xin1[0],xin2[0],xin3[0],magin[:,0])
                # take out all the odd 'bad values' and turn them into NaN
                if np.isclose(bmin,badval): bmin = np.NaN
                GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN
                #print('bmin, GEOcoord = {0},{1}'.format(bmin, GEOcoord))
                if not (np.isnan(bmin)):
                    #Now get Lstar at this location...
                    postry = spc.Coords(GEOcoord, 'GEO', 'car')
                    pa = AlphaOfK(tt, postry, k_i, extMag=extMag, options=Aopt, omnivals=omnivals)
                    LStry = get_Lstar(tt, postry, pa, extMag=extMag, options=options, omnivals=omnivals)
                    LStry['K'] = LStry['Xj']*np.sqrt(LStry['Bmirr']*nTtoG)
                else:
                    LStry = {'Lstar': np.NaN, 'K': np.NaN}
                if 'verbose' in kwargs: print('L* at test point: {0}; Xgsm = {1}'.format(LStry['Lstar'], pos_test))

                if np.isnan(LStry['Lstar']).any():
                    loci_brac2 = pos_test
                else:
                    loci_brac1 = pos_test
                    LCDS, LCDS_PA = LStry['Lstar'][0], pa
            try:
                results['LCDS'][idxt, idxk] = LCDS
                results['AlphaEq'][idxt, idxk] = LCDS_PA
            except IndexError:
                results['LCDS'][idxt] = LCDS
                results['AlphaEq'][idxt] = LCDS_PA

    results['LCDS'] = results['LCDS']
    results['AlphaEq'] = results['AlphaEq']

    return results

# -----------------------------------------------
def AlphaOfK(ticks, loci, K, extMag='T01STORM', options=[0,0,3,0,0], omnivals=None):
    """
    Find the equatorial pitch angle corresponding to a given second invariant K.

    Uses the IRBEM library to determine K and searches via bisection to find Alpha(K).

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information **for a single time**
        - loci (Coords class) : containing position information **for a single point**
        - K (numeric) : value of second invariant (K) for search
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - (see Lstar documentation for further explanation)
        - tol (float) : tolerance for search in radial distance
        - bracket (list): X-GSM coordinates to bracket bisection search

    Returns
    =======
        AlphaEq : equatorial pitch angle corresponding to K

    Examples
    ========
    >>> t = spacepy.time.Ticktock(['2002-09-01T04:00:00'], 'ISO')
    >>> loci = spacepy.coordinates.Coords([-4,0,0], 'GSM', 'car')
    >>> spacepy.irbempy.AlphaOfK(t, loci, 0.11, extMag='T89')
    48.984375

    See Also
    ========
    get_Lstar, get_Bfield, find_Bmirr, find_LCDS
    """
    # prepare input values for irbem call
    d = prep_irbem(ticks, loci, extMag=extMag, options=options, omnivals=omnivals)
    nTAI = len(ticks)
    badval = d['badval']
    kext = d['kext']
    sysaxes = d['sysaxes']
    iyearsat = d['iyearsat']
    idoysat = d['idoysat']
    secs = d['utsat']
    xin1 = d['xin1']
    xin2 = d['xin2']
    xin3 = d['xin3']
    magin = d['magin']
    nTtoG = 1.0e-5

    outvals = np.zeros(nTAI)
    outvals.fill(np.NaN)
    for i in np.arange(nTAI):
        bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
                         iyearsat[i],idoysat[i],secs[i], xin1[i],xin2[i],xin3[i],magin[:,i])

        # take out all the odd 'bad values' and turn them into NaN
        if np.isclose(bmin,badval): bmin = np.NaN
        GEOcoord[np.where( np.isclose(GEOcoord, badval)) ] = np.NaN

        pa0 = 90 #start with equatorially mirroring
        #Now get K for initial alpha at this location...
        if np.isfinite(GEOcoord[0]):
            pos1 = spc.Coords(GEOcoord, 'GEO', 'car')
            LS1 = get_Lstar(ticks[i], pos1, pa0, extMag=extMag, options=options, omnivals=omnivals)
            if np.isnan(LS1['Xj']).any():
                return np.NaN
            LS1['K'] = LS1['Xj']*np.sqrt(LS1['Bmirr']*nTtoG)
        else:
            #print('i = {0}, skipping because of bad position'.format(i))
            continue
        K0 = LS1['K'][0]

        if np.abs(K0-K) < 1e-3:
            outvals[i] = pa0
            #print('i= {0}, found alpha = {1}'.format(i, pa0))
            continue

        Done, count = False, 0
        pa_upper, pa_lower, pa_test = pa0, 0, 30
        while not Done:
            #print('Testing PA={0}'.format(pa_test))
            #Now get K for initial alpha at this location...
            LS1 = get_Lstar(ticks[i], pos1, pa_test, extMag=extMag, options=options, omnivals=omnivals)
            if np.isnan(LS1['Xj']).any():
                break
            LS1['K'] = LS1['Xj']*np.sqrt(LS1['Bmirr']*nTtoG)
            #print('L* at inner bracket: {0}'.format(LS1['Lstar']))
            K_test = LS1['K'][0]

            if np.abs(K_test-K) < 1e-3:
                #print('***Found K={0} (Req. {1}) at Alpha={2}'.format(K_test, K, pa_test))
                outvals[i] = pa_test
                break
            #print('Not Done: Kfound = {0}, Kwant = {1}'.format(K_test, K))
            if K_test > K: # K is too large, hence alpha is too small, increase.
                pa_lower = pa_test
                #print('Reset lower bound. U,L = {0},{1}'.format(pa_upper, pa_lower))
            else: # K is too small, hence alpha is too large, reduce
                pa_upper = pa_test
                #print('Reset upper bound. U,L = {0},{1}'.format(pa_upper, pa_lower))
            pa_test = pa_lower+np.abs((pa_upper-pa_lower))/2.0
            #print('Change alpha to {0}'.format(pa_test))
            count += 1
            if count>20: break
    return outvals
        

# -----------------------------------------------
def find_footpoint(ticks, loci, extMag='T01STORM', options=[1,0,3,0,0], hemi='same', alt=100, omnivals=None):
    """
    call find_foot_point1 from irbem library and return a dictionary with values for
    Bmin and the GEO (cartesian) coordinates of the magnetic equator

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained in Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - (see Lstar documentation for further explanation)
        - hemi (string) : optional (valid cases are 'same', 'other', 'north' or 'south')
                          will set the target hemisphere for tracing the footpoint
        - alt (numeric) : optional keyword to set stop height [km] of fieldline trace (default 100km)

    Returns
    =======
        - results (spacepy.datamodel.SpaceData) : containing keys
                   Bfoot    - Magnitude of B-field at footpoint [nT]
                   loci     - Coords instance with GDZ coordinates of the magnetic footpoint [alt, lat, lon]
                   Bfootvec - Components of B-field at footpoint in cartesian GEO coordinates [nT]

    Examples
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[3,0,0]], 'GEO', 'car')
    >>> spacepy.irbempy.find_footpoint(t, y)
    {'Bfoot': array([ 47559.04643444,  47542.84688657]),
     'Bfootvec': array([[-38428.07217246,   4497.31549786, -27657.19291928],
           [-38419.08514332,   4501.45390964, -27641.14866517]]),
     'loci': Coords( [[ 99.31443778  55.71415787 -10.21888955]
     [ 99.99397026  55.70716296 -10.22797462]] ), dtype=GDZ,sph, units=['km', 'deg', 'deg']}

    See Also
    ========
    get_Lstar, get_Bfield, find_Bmirr, find_magequator
    """
    # prepare input values for irbem call
    d = prep_irbem(ticks, loci, extMag=extMag, options=options, omnivals=omnivals)
    nTAI = len(ticks)
    badval = d['badval']
    kext = d['kext']
    sysaxes = d['sysaxes']
    iyearsat = d['iyearsat']
    idoysat = d['idoysat']
    secs = d['utsat']
    xin1 = d['xin1']
    xin2 = d['xin2']
    xin3 = d['xin3']
    magin = d['magin']

    hemi_dict = {'same': 0, 'north': 1, 'south': -1, 'other': 2}
    if hemi.lower() in ['same', 'other', 'north', 'south']:
        hemi_flag = hemi_dict[hemi.lower()]
    else:
        raise ValueError('Option for hemisphere to trace to ({0}) is invalid.\n'.format(hemi) +
                         'Valid options are: same, other, north and south')

    results = dm.SpaceData()
    results['Bfoot'] = np.zeros(nTAI)
    results['Bfootvec'] = np.zeros([nTAI, 3])
    results['loci'] = ['']*nTAI
    for i in np.arange(nTAI):
        xfoot, bfoot, bfootmag = oplib.find_foot_point1(kext, options, sysaxes,\
            iyearsat[i], idoysat[i], secs[i], xin1[i], xin2[i], xin3[i], alt, hemi_flag, magin[:,i])

        # take out all the odd 'bad values' and turn them into NaN
        if np.isclose(bfootmag,badval): bmin = np.NaN
        xfoot[np.where( np.isclose(xfoot, badval)) ] = np.NaN
        bfoot[np.where( np.isclose(bfoot, badval)) ] = np.NaN

        results['Bfoot'][i] = bfootmag
        results['loci'][i] = xfoot
        results['Bfootvec'][i] = bfoot

    results['loci'] = spc.Coords(results['loci'], 'GDZ', 'sph')
    results.attrs
    return results


# -----------------------------------------------
def coord_trans(loci, returntype, returncarsph ):
    """
    thin layer to call coor_trans1 from irbem lib
    this will convert between systems GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL

    Parameters
    ==========
        - loci (Coords instance) : containing coordinate information, can contain n points
        - returntype (str) : describing system as GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
        - returncarsph (str) : cartesian or spherical units 'car', 'sph'
        
    Returns
    =======
        - xout (ndarray) : values after transformation in (n,3) dimensions
        
    Examples
    ========
    >>> c = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> c.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> coord_trans(c, 'GSM', 'car')
    array([[ 2.8639301 , -0.01848784,  0.89306361],
    [ 1.9124434 ,  0.07209424,  0.58082929]])

    See Also
    ========
    sph2car, car2sph

    """
    sysaxesin = get_sysaxes( loci.dtype, loci.carsph )
    sysaxesout = get_sysaxes( returntype, returncarsph )

    # swap carsph if sysaxesout is None
    if (sysaxesout is None) and (returncarsph == 'sph'):
        sysaxesout = get_sysaxes( returntype, 'car' )
        aflag = True
    elif (sysaxesout is None) and (returncarsph == 'car'):
        sysaxesout = get_sysaxes( returntype, 'sph' )
        aflag = True
    else:
        aflag = False

    xout = np.zeros(np.shape(loci.data))
    for i in np.arange(len(loci)):
        iyear = loci.ticks.UTC[i].year
        idoy = loci.ticks.DOY[i]
        secs = loci.ticks.UTC[i].hour*3600. + loci.ticks.UTC[i].minute*60. + \
            loci.ticks.UTC[i].second
            
        xout[i,:] = oplib.coord_trans1(sysaxesin, sysaxesout, \
            iyear, idoy, secs, loci.data[i])
        
    # add  sph to car or v/v convertion if initial sysaxesout was None
    if  aflag == True:
        if returncarsph == 'sph':
            xout = car2sph(xout)
        else: # 'car' needs to be returned
            xout = sph2car(xout)

    return xout

# -----------------------------------------------
def car2sph(CARin):
    """
    coordinate transformation from cartesian to spherical

    Parameters
    ==========
        - CARin (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
            units of [Re, Re, Re] = [x,y,z]
        
    Returns
    =======
        - results (ndarray) : values after conversion to spherical coordinates in
            radius, latitude, longitude in units of [Re, deg, deg]
        
    Examples
    ========
    >>> sph2car([1,45,0])
    array([ 0.70710678,  0.        ,  0.70710678])

    See Also
    ========
    sph2car
    """	

    if isinstance(CARin[0], numbers.Number):
        CAR = np.array([CARin])
    else:
        CAR = np.asanyarray(CARin)

    res = np.zeros(np.shape(CAR))
    for i in np.arange(len(CAR)):
        x, y, z = CAR[i,0], CAR[i,1], CAR[i,2]
        r = np.sqrt(x*x+y*y+z*z)
        sq = np.sqrt(x*x+y*y)
        if (x == 0) & (y == 0): # on the poles
            longi = 0.
            if z < 0:
                lati = -90.
            else:
                lati = 90.0
        else:
            longi = np.arctan2(y,x)*180./np.pi
            lati = 90. - np.arctan2(sq, z)*180./np.pi
        res[i,:] = [r, lati, longi]
        
    if isinstance(CARin[0], numbers.Number):
        return res[0]
    else:
        return res

# -----------------------------------------------
def sph2car(SPHin):
    """
    coordinate transformation from spherical to cartesian

    Parameters
    ==========
        - SPHin (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
            units of [Re, deg, deg] = [r, latitude, longitude]
        
    Returns
    =======
        - results (ndarray) : values after conversion to cartesian coordinates x,y,z
        
    Examples
    ========
    >>> sph2car([1,45,45])
    array([ 0.5       ,  0.5       ,  0.70710678])

    See Also
    ========
    car2sph

    """

    if isinstance(SPHin[0], numbers.Number):
        SPH = np.array([SPHin])
    else:
        SPH = np.asanyarray(SPHin)

    res = np.zeros(np.shape(SPH))
    for i in np.arange(len(SPH)):
        r,lati,longi = SPH[i,0], SPH[i,1], SPH[i,2]
        colat = np.pi/2. - lati*np.pi/180.
        x = r*np.sin(colat)*np.cos(longi*np.pi/180.)
        y = r*np.sin(colat)*np.sin(longi*np.pi/180.)
        z = r*np.cos(colat)
        res[i,:] = [x, y, z]


    if isinstance(SPHin[0], numbers.Number):
        return res[0]
    else:
        return res

# -----------------------------------------------    
def get_sysaxes(dtype, carsph):
    """
    will return the sysaxes according to the irbem library

    Parameters
    ==========
    - dtype (str) : coordinate system, possible values: GDZ, GEO, GSM, GSE, SM, 
            GEI, MAG, SPH, RLL
    - carsph (str) : cartesian or spherical, possible values: 'sph', 'car'
        
    Returns
    =======
    - sysaxes (int) : value after oner_desp library from 0-8 (or None if not available)

    Examples
    ========
    >>> get_sysaxes('GSM', 'car')
    2

    See Also
    ========
    get_dtype

    """

    #typedict = {'GDZ': {'sph': 0, 'car': 10},
        #'GEO': {'sph': 1, 'car': 11}, 'GSM': {'sph': 22, 'car': 2},
        #'GSE': {'sph': 23, 'car': 3}, 'SM': {'sph': 24, 'car': 4},
        #'GEI': {'sph': 25, 'car': 5}, 'MAG': {'sph': 26, 'car': 6},
        #'SPH': {'sph': 7, 'car': 17}, 'RLL': {'sph': 8, 'car': 18}}
        
    
    sysaxes = SYSAXES_TYPES[dtype][carsph]

    return sysaxes
    
# -----------------------------------------------    
def get_dtype(sysaxes):
    """
    will return the coordinate system type as string

    Parameters
    ==========
        - sysaxes (int) : number according to the irbem, possible values: 0-8
        
    Returns
    =======
        - dtype (str) : coordinate system GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
        - carsph (str) : cartesian or spherical 'car', 'sph'

    Examples
    ========
    >>> get_dtype(3)
    ('GSE', 'car')

    See Also
    ========
    get_sysaxes

    """

    for key in SYSAXES_TYPES:
        for subkey in SYSAXES_TYPES[key]:
            if SYSAXES_TYPES[key][subkey] == sysaxes:
                dtype = key
                carsph = subkey

    return dtype, carsph

# -----------------------------------------------    
def get_AEP8(energy, loci, model='min', fluxtype='diff', particles='e'):
    """
    will return the flux from the AE8-AP8 model

    Parameters
    ==========
        - energy (float) : center energy in MeV; if fluxtype=RANGE, then needs to be a list [Emin, Emax]
        - loci (Coords)  : a Coords instance with the location inside the magnetosphere
             optional      instead of a Coords instance, one can also provide a list with [BBo, L] combination
        - model (str)    : MIN or MAX for solar cycle dependence
        - fluxtype (str) : DIFF, RANGE, INT are possible types
        - particles (str): e or p or electrons or protons
        
    Returns
    =======
        - float : flux from AE8/AP8 model
        
    Examples
    ========
    >>> spacepy.irbempy.get_aep8()

    >>> import spacepy.time as spt
    >>> import spacepy.coordinates as spc
    >>> import spacepy.irbempy as ib
    >>> t = spt.Ticktock(['2017-02-02T12:00:00'], 'ISO')
    >>> y = spc.Coords([3,0,0], 'GEO', 'car')
    >>> y.ticks = t
    >>> energy = 1.0 #MeV
    >>> ib.get_AEP8(energy, y, model='max')
    1932209.4427359989

    """
    # find bad values and dimensions
    
    if particles.lower() == 'e': # then choose electron model
        if model.upper() == 'MIN':
            whichm = 1
        elif model.upper() == 'MAX':
            whichm = 2
        else:
            print('Warning: model='+model+' is not implemented: Choose MIN or MAX')
    elif particles.lower() == 'p': # then choose proton model
        if model.upper() == 'MIN':
            whichm = 3
        elif model.upper() == 'MAX':
            whichm = 4
        else:
            print('Warning: model='+model+' is not implemented: Choose MIN or MAX')
    else:
        print('Warning: particles='+particles+' is not available: choose e or p')
    
    if fluxtype.upper() == 'DIFF':
        whatf = 1
    elif fluxtype.upper() == 'RANGE':
        whatf = 2
    elif fluxtype.upper() == 'INT':
        whatf = 3
    else:
        print('Warning: fluxtype='+fluxtype+' is not implemented: Choose DIFF, INT or RANGE')
    
    # need range if whatf=2
    Nene  = 1
    if whatf == 2: assert len(energy) == 2, 'Need energy range with this choice of fluxtype=RANGE'
    ntmax = 1

    #build dummy OMNI dictionary for prep_irbem (AE/AP8 doesn't need these inputs)
    dum_omni = dict()
    magkeys = ['Kp', 'Dst', 'dens', 'velo', 'Pdyn', 'ByIMF', 'BzIMF',\
                    'G1', 'G2', 'G3', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    for key in magkeys:
        dum_omni[key] = [0]*len(loci)
    
    if isinstance(loci, spc.Coords):
        assert loci.ticks, 'Coords require time information with a Ticktock object'
        d = prep_irbem(ticks=loci.ticks, loci=loci, omnivals=dum_omni)
        E_array = np.zeros((2,d['nalp_max']))
        E_array[:,0] = energy
        # now get the flux
        flux = oplib.fly_in_nasa_aeap1(ntmax, d['sysaxes'], whichm, whatf, Nene, E_array, d['iyearsat'], d['idoysat'], d['utsat'], \
            d['xin1'], d['xin2'], d['xin3'])
    elif isinstance(loci, (list, np.ndarray)):
        BBo, L = loci
        d = prep_irbem(omnivals=dum_omni)
        E_array = np.zeros((2,d['nalp_max']))
        E_array[:,0] = energy
        B_array = np.zeros(d['ntime_max'])
        B_array[0] = BBo
        L_array = np.zeros(d['ntime_max'])
        L_array[0] = L
        # now get the flux
        flux =  oplib.get_ae8_ap8_flux(ntmax, whichm, whatf, Nene, E_array, B_array, L_array)
    
    else:
        print('Warning: coords need to be either a spacepy.coordinates.Coords instance or a list of [BBo, L]')
        
    
    flux[np.where( np.isclose(flux, d['badval'])) ] = np.NaN
    
    return flux[0,0]
    

# -----------------------------------------------
def _get_Lstar(ticks, loci, alpha, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None, landi2lstar=False): 
    """
    This will call make_lstar1 or make_lstar_shell_splitting_1 from the irbem library
    and will lookup omni values for given time if not provided (optional). If pitch angles
    are provided, drift shell splitting will be calculated and "Bmirr" will be returned. If they
    are not provided, then no drift shell splitting is calculated and "Blocal" is returned.

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - alpha (list or ndarray) : pitch angles in degrees; if provided will 
            calculate shell splitting; max 25 values
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained below
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - landi2lstar : if True, will use the faster landi2lstar routine if possible. This 
            routine can only be used with OPQUIET+IGRF magnetic field models.

    Returns
    =======
        - results (dictionary) : containing keys: Lm, Lstar, Bmin, Blocal (or Bmirr), 
            Xj (I - the 2nd invariant), MLT
            if pitch angles provided in "alpha" then drift shells are calculated and "Bmirr" 
            is returned if not provided, then "Blocal" at spacecraft is returned.

    Examples
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> spacepy.irbempy.Lstar(t,y)
    {'Blocal': array([ 1020.40493286,  3446.08845227]),
        'Bmin': array([ 1019.98404311,  3437.63865243]),
        'Lm': array([ 3.08948304,  2.06022102]),
        'Lstar': array([ 2.97684043,  1.97868577]),
        'MLT': array([ 23.5728333 ,  23.57287944]),
        'Xj': array([ 0.00112884,  0.00286955])}


    Notes
    =====
    External Field
        - 0	   : no external field
        - MEAD	: Mead & Fairfield [1975] (uses 0<=Kp<=9 - Valid for rGEO<=17. Re)
        - T87SHORT: Tsyganenko short [1987] (uses 0<=Kp<=9 - Valid for rGEO<=30. Re)
        - T87LONG : Tsyganenko long [1987] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - T89	 : Tsyganenko [1989] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - OPQUIET : Olson & Pfitzer quiet [1977] (default - Valid for rGEO<=15. Re)
        - OPDYN   : Olson & Pfitzer dynamic [1988] (uses 5.<=dens<=50., 300.<=velo<=500., 
            -100.<=Dst<=20. - Valid for rGEO<=60. Re)
        - T96	 : Tsyganenko [1996] (uses -100.<=Dst (nT)<=20., 0.5<=Pdyn (nPa)<10., 
            |ByIMF| (nT)<1=0., |BzIMF| (nT)<=10. - Valid for rGEO<=40. Re)
        - OSTA	: Ostapenko & Maltsev [1997] (uses dst,Pdyn,BzIMF, Kp)
            T01QUIET: Tsyganenko [2002a,b] (uses -50.<Dst (nT)<20., 0.5<Pdyn (nPa)<=5., 
            |ByIMF| (nT)<=5., |BzIMF| (nT)<=5., 0.<=G1<=10., 0.<=G2<=10. - Valid for xGSM>=-15. Re)
        - T01STORM: Tsyganenko, Singer & Kasper [2003] storm  (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 - 
            there is no upper or lower limit for those inputs - Valid for xGSM>=-15. Re)
        - T05	 : Tsyganenko & Sitnov [2005] storm  (uses Dst, Pdyn, ByIMF, BzIMF, 
            W1, W2, W3, W4, W5, W6 - no upper or lower limit for inputs - Valid for xGSM>=-15. Re)
        - TS07   : Tsyganenko and Sitnov [2007] model. Uses specially calculated coefficient files.

    OMNI contents
        - Kp: value of Kp as in OMNI2 files but has to be double instead of integer type
        - Dst: Dst index (nT)
        - dens: Solar Wind density (cm-3)
        - velo: Solar Wind velocity (km/s)
        - Pdyn: Solar Wind dynamic pressure (nPa)
        - ByIMF: GSM y component of IMF mag. field (nT)
        - BzIMF: GSM z component of IMF mag. field (nT)
        - G1:  G1=< Vsw*(Bperp/40)2/(1+Bperp/40)*sin3(theta/2) > where the <> mean an average over the 
            previous 1 hour, Vsw is the solar wind speed, Bperp is the transverse IMF 
            component (GSM) and theta it's clock angle.
        - G2: G2=< a*Vsw*Bs > where the <> mean an average over the previous 1 hour, 
                Vsw is solar wind speed, Bs=|IMF Bz| when IMF Bz < 0 and Bs=0 when IMF Bz > 0, a=0.005.
        - G3:  G3=< Vsw*Dsw*Bs /2000.> where the <> mean an average over the previous 1 hour, 
                Vsw is the solar wind speed, Dsw is the solar wind density, Bs=|IMF Bz| when IMF 
                Bz < 0 and Bs=0 when IMF Bz > 0.
        - W1 - W6: see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
        - AL: the auroral index

    Options
        - 1st element: 0 - don't compute L* or phi ;  1 - compute L*; 2- compute phi
        - 2nd element: 0 - initialize IGRF field once per year (year.5);  
            n - n is the  frequency (in days) starting on January 1st of each year 
            (i.e. if options(2nd element)=15 then IGRF will be updated on the following days of the 
            year: 1, 15, 30, 45 ...)
        - 3rd element: resolution to compute L* (0 to 9) where 0 is the recomended value to ensure a 
            good ratio precision/computation time
            (i.e. an error of ~2% at L=6). The higher the value the better will be the precision, the 
            longer will be the computing time. Generally there is not much improvement for values 
            larger than 4. Note that this parameter defines the integration step (theta) 
            along the field line such as dtheta=(2pi)/(720*[options(3rd element)+1])
        - 4th element: resolution to compute L* (0 to 9). The higher the value the better will be 
            the precision, the longer will be 
            the computing time. It is recommended to use 0 (usually sufficient) unless L* is not 
            computed on a LEO orbit. For LEO orbit higher values are recommended. Note that this 
            parameter defines the integration step (phi) along the drift shell such as 
            dphi=(2pi)/(25*[options(4th element)+1])
        - 5th element: allows to select an internal magnetic field model (default is set to IGRF)
            - 0 = IGRF
            - 1 = Eccentric tilted dipole
            - 2 = Jensen&Cain 1960
            - 3 = GSFC 12/66 updated to 1970
            - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
            - 5 = Centred dipole
            
    """
    nTAI = len(ticks)
    d = prep_irbem(ticks, loci, alpha, extMag, options, omnivals)
    nalpha = d['nalpha']
    if d['kext'] != 5 or d['options'][4] != 0:
        landi2lstar = False
    no_shell_splitting = (nalpha == 0) or (nalpha == 1 and alpha[0] == 90)

    # Arguments that are common to all flavors of L* functions
    args = [nTAI, d['kext'], d['options'],
            d['sysaxes'],d['iyearsat'], d['idoysat'], d['utsat'],
            d['xin1'], d['xin2'], d['xin3'], d['magin']]
    if no_shell_splitting: # no drift shell splitting
        func = oplib.landi2lstar1 if landi2lstar else oplib.make_lstar1
    else: # with drift shell splitting
        # Drift shell splitting requires pitch angle positional args
        args.insert(1, nalpha)
        args.insert(-1, d['degalpha'])
        func = oplib.landi2lstar_shell_splitting1 if landi2lstar \
               else oplib.make_lstar_shell_splitting1
    # For a locally 90-degree particle, bmirr is blocal
    lm, lstar, bmirr, bmin, xj, mlt = func(*args)

    # take out all the odd 'bad values' and turn them into NaN
    lm[np.where( np.isclose(lm,d['badval'])) ] = np.NaN
    lstar[np.where( np.isclose(lstar,d['badval'])) ] = np.NaN
    bmin[np.where( np.isclose(bmin,d['badval'])) ] = np.NaN
    xj[np.where( np.isclose(xj,d['badval'])) ] = np.NaN
    mlt[np.where( np.isclose(mlt,d['badval'])) ] = np.NaN

    results = {}
    if no_shell_splitting:
        results['Lm'] = lm[0:nTAI][:,None]
        results['Lstar'] = lstar[0:nTAI][:,None]
        bmirr[np.where( np.isclose(bmirr,d['badval'])) ] = np.NaN
        results['Blocal'] = bmirr[0:nTAI]
        results['Bmirr'] = results['Blocal'][:,None]
        results['Bmin'] = bmin[0:nTAI]
        results['Xj'] = xj[0:nTAI][:,None]
        results['MLT'] = mlt[0:nTAI]
    else:		
        results['Lm'] = lm[0:nTAI, 0:nalpha]
        results['Lstar'] = lstar[0:nTAI, 0:nalpha]
        bmirr[np.where( np.isclose(bmirr, d['badval'])) ] = np.NaN
        results['Bmirr'] = bmirr[0:nTAI, 0:nalpha]
        results['Bmin'] = bmin[0:nTAI]
        results['Xj'] = xj[0:nTAI, 0:nalpha]
        results['MLT'] = mlt[0:nTAI]
        
    return results

# -----------------------------------------------
def get_Lm(ticks, loci, alpha, extMag='T01STORM', intMag='IGRF', IGRFset=0, omnivals=None):
    """
    Return the MacIlwain L value for a given location, time and model

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - alpha (list or ndarray) : pitch angles in degrees
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - intMag (string) : optional: select the internal field model
                            possible values ['IGRF','EDIP','JC','GSFC','DUN','CDIP']
                            For full details see get_Lstar
        - omni values as dictionary (optional) : if not provided, will use lookup table 

    Returns
    =======

        - results (dictionary) : containing keys: Lm, Bmin, Blocal (or Bmirr), Xj, MLT 
            if pitch angles provided in "alpha" then drift shells are calculated and "Bmirr" 
            is returned if not provided, then "Blocal" at spacecraft is returned. A negative
            value for Lm indicates the field line is closed but particles are lost to the
            atmosphere; the absolute value indicates the L value.

    Examples
    ========


    Notes
    =====



    """
    intMaglookup = {'IGRF': 0, 'EDIP': 1, 'JC': 2, 'GSFC': 3, 'DUN': 4, 'CDIP': 5}
    if intMag not in intMaglookup:
        raise ValueError('Invalid value of intMag: valid values are: {0}'.format(intMaglookup.keys()))
    if IGRFset != 0:
        try:
            assert IGRFset > 0
            ##TODO: test for numeric type
        except AssertionError:
            raise ValueError('IGRFset must be positive-valued and numeric')
        
    opts = [0, IGRFset, 0, 0, intMaglookup[intMag]]

    results = get_Lstar(ticks, loci, alpha, extMag=extMag, options=opts, omnivals=omnivals)
    dum = results.pop('Lstar')
    return results

# -----------------------------------------------
def get_Lstar(ticks, loci, alpha=90, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None, landi2lstar=False):
    """
    This will call make_lstar1 or make_lstar_shell_splitting_1 from the irbem library
    and will lookup omni values for given time if not provided (optional). If pitch angles
    are provided, drift shell splitting will be calculated and "Bmirr" will be returned. If they
    are not provided, then no drift shell splitting is calculated and "Blocal" is returned.

    Parameters
    ==========
        - ticks (Ticktock class) : containing time information
        - loci (Coords class) : containing spatial information
        - alpha (list or ndarray) : optional pitch angles in degrees (default is 90);
            if provided will calculate shell splitting; max 25 values
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX', 'TS07']
        - options (optional list or array of integers length=5) : explained below
        - omni values as dictionary (optional) : if not provided, will use lookup table 
        - landi2lstar : if True, will use the faster landi2lstar routine if possible. This 
            routine can only be used with OPQUIET+IGRF magnetic field models.

    Returns
    =======
        - results (dictionary) : containing keys: Lm, Lstar, Bmin, Blocal (or Bmirr), Xj, MLT 
            if pitch angles provided in "alpha" then drift shells are calculated and "Bmirr" 
            is returned if not provided, then "Blocal" at spacecraft is returned. A negative
            value for Lm indicates the field line is closed but particles are lost to the
            atmosphere; the absolute value indicates the L value. A negative value for Lstar
            indicates the field line is closed but particles are lost to the atmosphere
            before completing a drift orbit; the absolute value indicates the drift shell.

    Examples
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> spacepy.irbempy.Lstar(t,y)
    {'Blocal': array([ 1020.40493286,  3446.08845227]),
        'Bmin': array([ 1019.98404311,  3437.63865243]),
        'Lm': array([ 3.08948304,  2.06022102]),
        'Lstar': array([ 2.97684043,  1.97868577]),
        'MLT': array([ 23.5728333 ,  23.57287944]),
        'Xj': array([ 0.00112884,  0.00286955])}


    Notes
    =====
    External Field
        - 0	   : no external field
        - MEAD	: Mead & Fairfield [1975] (uses 0<=Kp<=9 - Valid for rGEO<=17. Re)
        - T87SHORT: Tsyganenko short [1987] (uses 0<=Kp<=9 - Valid for rGEO<=30. Re)
        - T87LONG : Tsyganenko long [1987] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - T89	 : Tsyganenko [1989] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - OPQUIET : Olson & Pfitzer quiet [1977] (default - Valid for rGEO<=15. Re)
        - OPDYN   : Olson & Pfitzer dynamic [1988] (uses 5.<=dens<=50., 300.<=velo<=500., 
            -100.<=Dst<=20. - Valid for rGEO<=60. Re)
        - T96	 : Tsyganenko [1996] (uses -100.<=Dst (nT)<=20., 0.5<=Pdyn (nPa)<10., 
            abs(ByIMF) (nT)<1=0., abs(BzIMF) (nT)<=10. - Valid for rGEO<=40. Re)
        - OSTA	: Ostapenko & Maltsev [1997] (uses dst,Pdyn,BzIMF, Kp)
            T01QUIET: Tsyganenko [2002a,b] (uses -50.<Dst (nT)<20., 0.5<Pdyn (nPa)<=5., 
            abs(ByIMF) (nT)<=5., abs(BzIMF) (nT)<=5., 0.<=G1<=10., 0.<=G2<=10. - Valid for xGSM>=-15. Re)
        - T01STORM: Tsyganenko, Singer & Kasper [2003] storm  (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 - 
            there is no upper or lower limit for those inputs - Valid for xGSM>=-15. Re)
        - T05	 : Tsyganenko & Sitnov [2005] storm  (uses Dst, Pdyn, ByIMF, BzIMF, 
            W1, W2, W3, W4, W5, W6 - no upper or lower limit for inputs - Valid for xGSM>=-15. Re)

    OMNI contents
        - Kp: value of Kp as in OMNI2 files but has to be double instead of integer type
        - Dst: Dst index (nT)
        - dens: Solar Wind density (cm-3)
        - velo: Solar Wind velocity (km/s)
        - Pdyn: Solar Wind dynamic pressure (nPa)
        - ByIMF: GSM y component of IMF mag. field (nT)
        - BzIMF: GSM z component of IMF mag. field (nT)
        - G1:  G1=< Vsw*(Bperp/40)2/(1+Bperp/40)*sin3(theta/2) > where the <> mean an average over the 
            previous 1 hour, Vsw is the solar wind speed, Bperp is the transverse IMF 
            component (GSM) and theta it's clock angle.
        - G2: G2=< a*Vsw*Bs > where the <> mean an average over the previous 1 hour, 
                Vsw is solar wind speed, Bs=|IMF Bz| when IMF Bz < 0 and Bs=0 when IMF Bz > 0, a=0.005.
        - G3:  G3=< Vsw*Dsw*Bs /2000.> where the <> mean an average over the previous 1 hour, 
                Vsw is the solar wind speed, Dsw is the solar wind density, Bs=|IMF Bz| when IMF 
                Bz < 0 and Bs=0 when IMF Bz > 0.
        - W1 - W6: see definition in (JGR-A, v.110(A3), 2005.) (PDF 1.2MB)
        - AL: the auroral index

    Options
        - 1st element: 0 - don't compute L* or phi ;  1 - compute L*; 2- compute phi
        - 2nd element: 0 - initialize IGRF field once per year (year.5);  
            n - n is the  frequency (in days) starting on January 1st of each year 
            (i.e. if options(2nd element)=15 then IGRF will be updated on the following days of the 
            year: 1, 15, 30, 45 ...)
        - 3rd element: resolution to compute L* (0 to 9) where 0 is the recomended value to ensure a 
            good ratio precision/computation time
            (i.e. an error of ~2% at L=6). The higher the value the better will be the precision, the 
            longer will be the computing time. Generally there is not much improvement for values 
            larger than 4. Note that this parameter defines the integration step (theta) 
            along the field line such as dtheta=(2pi)/(720*[options(3rd element)+1])
        - 4th element: resolution to compute L* (0 to 9). The higher the value the better will be 
            the precision, the longer will be 
            the computing time. It is recommended to use 0 (usually sufficient) unless L* is not 
            computed on a LEO orbit. For LEO orbit higher values are recommended. Note that this 
            parameter defines the integration step (phi) along the drift shell such as 
            dphi=(2pi)/(25*[options(4th element)+1])
        - 5th element: allows to select an internal magnetic field model (default is set to IGRF)
            - 0 = IGRF
            - 1 = Eccentric tilted dipole
            - 2 = Jensen&Cain 1960
            - 3 = GSFC 12/66 updated to 1970
            - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
            - 5 = Centred dipole
            
    """

    from spacepy import config as config_dict

    def get_ov(fullov, stind, enind):
        '''Chop up omni data for multiprocessing'''
        out = dm.SpaceData()
        keylist = list(fullov.keys())
        dum = keylist.pop(keylist.index('Qbits'))
        for key in keylist:
            out[key] = fullov[key][stind:enind]
        for key in fullov['Qbits']:
            out['Qbits'][key] = fullov['Qbits'][key][stind:enind]
        return out
    
    def reassemble(result):
        '''Reassemble the results from the multiprocessing'''
        funcs = {'Bmin': np.hstack,
                 'Blocal': np.hstack,
                 'Bmirr': np.vstack,
                 'Lm': np.vstack,
                 'Lstar': np.vstack,
                 'MLT': np.hstack,
                 'Xj': np.vstack}
        out = dm.SpaceData() 
        for el in result:
            for key in result[0]:
                if key in list(out.keys()): #not first chunk
                    out[key] = funcs[key]([out[key], el[key]])
                else: #first chunk
                    out[key] = el[key].copy()
        return out

    if isinstance(alpha, numbers.Number):
        alpha = [alpha]

    ncpus = config_dict['ncpus']
    ncalc = len(ticks)
    nalpha = len(alpha)

    if ncalc < ncpus * 2: #Don't multiprocess if not worth it
        ncpus = 1
    if ncpus > 1:
        import __main__ as main
        if hasattr(main, '__file__'):
            try:
                from multiprocessing import Pool
                pool = Pool(ncpus)
            except (ImportError, OSError):
                ncpus = 1 #if pool setup fails, e.g. Wine, go single
        else:
            ncpus = 1 #won't multiprocess in interactive mode

    if ncpus > 1:
        nblocks = ncpus
        blocklen = np.floor_divide(ncalc, ncpus)
        tt, cc = [], []
        if omnivals:
            ov = []
        else:
            ov = None
        for block in range(nblocks):
            startind = block*blocklen
            if block != nblocks-1: #not last block
                endind = block*blocklen + blocklen
            else:
                endind = block*blocklen + 2*blocklen #going past the end of an array in a slice is fine
            tt.append(ticks[startind:endind]) #chunk time tags
            cc.append(loci[startind:endind]) #chunk positions
            if omnivals:
                ov.append(get_ov(omnivals, startind, endind))
        inputs = [[tch, cch, alpha, extMag, options, ov, landi2lstar] for tch, cch in zip(tt,cc)]
        result = pool.map(_multi_get_Lstar, inputs)
        pool.close()
        pool.join()
        DALL = reassemble(result)
    else: # single NCPU, no chunking
        DALL = _get_Lstar(ticks, loci, alpha, extMag, options, omnivals, landi2lstar)

    return DALL


def _multi_get_Lstar(inputs):
    '''
    '''
    ticks = inputs[0]
    loci = inputs[1]
    alpha = inputs[2]
    extMag = inputs[3]
    options = inputs[4]
    omnivals = inputs[5]
    landi2lstar = inputs[6]
    DALL = _get_Lstar(ticks, loci, alpha, extMag, options, omnivals, landi2lstar)
    
    return DALL

# -----------------------------------------------
def prep_irbem(ticks=None, loci=None, alpha=[], extMag='T01STORM', options=[1,0,0,0,0], omnivals=None): 
    """
    Prepare inputs for direct IRBEM-LIB calls. Not expected to be called by the user.
    """
    # setup dictionary to return input values for irbem
    d= {}
    d['badval'] = -1e31
    d['nalp_max'] = 25
    d['ntime_max'] = 100000
    d['options'] = options
    badval = d['badval']
    nalp_max = d['nalp_max']
    ntime_max = d['ntime_max']

    if ticks is None:
        return d

    UTC = ticks.UTC
    DOY = ticks.DOY
    eDOY = ticks.eDOY
    nTAI = len(ticks)

    # setup mag array and move omni values
    magin = np.zeros((nalp_max,ntime_max),float)
    magkeys = ['Kp', 'Dst', 'dens', 'velo', 'Pdyn', 'ByIMF', 'BzIMF',\
                    'G1', 'G2', 'G3', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    def fakeOMNI(npts, mk):
        fakeomni = {}
        for kk in magkeys:
            fakeomni[kk] = [1]*npts
        return fakeomni

    # get omni values
    if omnivals is None:
        # No OMNI values provided so look up (requires data files)
        if extMag.upper() not in ['0', 'OPQUIET','TS07']:
            import spacepy.omni as omni
            omnivals = omni.get_omni(ticks)
        # UNLESS we're asking for a static model... then we spoof the input as the numbers are irrelevant
        # OR TS07, which uses special files
        elif extMag.upper() == 'TS07':
            omnivals = fakeOMNI(len(ticks), magkeys) #keep this separate, just in case TS07 params ever get passed in the usual way
        else:
            omnivals = fakeOMNI(len(ticks), magkeys)
    if 'G' in omnivals:
        for n in range(1,4):
            dum = omnivals['G'][...,n-1]
            if dum.ndim == 0: dum = np.array([dum])
            omnivals['G{0}'.format(n)] = dum
        del omnivals['G']
    if 'W' in omnivals:
        for n in range(1,7):
            dum = omnivals['W'][...,n-1]
            if dum.ndim == 0: dum = np.array([dum])
            omnivals['W{0}'.format(n)] = dum
        del omnivals['W']
        
    for iTAI in np.arange(nTAI):
        for ikey, key in enumerate(magkeys):
            magin[ikey, iTAI] = omnivals[key][iTAI]

    # multiply Kp*10 to look like omni database
    # this is what irbem lib is looking for
    magin[0,:] = magin[0,:]*10.
    
    d['magin'] = magin

    # setup time array
    iyearsat = np.zeros(ntime_max, dtype=int)
    idoysat = np.zeros(ntime_max, dtype=int)
    utsat = np.zeros(ntime_max, dtype=float)
    for i in np.arange(nTAI):
        iyearsat[i] = UTC[i].year
        idoysat[i] = int(DOY[i])
        utsat[i] = (eDOY[i]-np.floor(eDOY[i]))*86400.
    d['iyearsat'] = iyearsat
    d['idoysat'] = idoysat
    d['utsat'] = utsat

    # copy coordinates into array
    # prepare coordinates
    if loci.sysaxes is None:
        # System type not supported by IRBEM
        # Convert car -> sph or vice versa as required
        newcarsph = [key for (key, val) in SYSAXES_TYPES[loci.dtype].items()
                     if val is not None][0]
        posi = loci.convert(loci.dtype, newcarsph)
    else:
        posi = loci
    d['sysaxes'] = posi.sysaxes
    xin1 = np.zeros(ntime_max, dtype=float)
    xin2 = np.zeros(ntime_max, dtype=float)
    xin3 = np.zeros(ntime_max, dtype=float) 
    if posi.carsph == 'sph':
        xin1[0:nTAI] = posi.radi[:]
        xin2[0:nTAI] = posi.lati[:]
        xin3[0:nTAI] = posi.long[:]
    else:
        xin1[0:nTAI] = posi.x[:]
        xin2[0:nTAI] = posi.y[:]
        xin3[0:nTAI] = posi.z[:]
    d['xin1'] = xin1
    d['xin2'] = xin2
    d['xin3'] = xin3

    # convert external magnetic field flag
    extkeys = ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 'OPQUIET', 'OPDYN', 'T96', \
                'OSTA', 'T01QUIET', 'T01STORM', 'T05', 'ALEX', 'TS07']
    assert extMag in  extkeys, 'extMag not available: %s' % extMag
    kext = extkeys.index(extMag.upper())
    d['kext'] = kext

    # calc at given pitch angles 'alpha'?
    degalpha = np.zeros(nalp_max, dtype=float)
    if isinstance(alpha, numbers.Number):
        nalpha = 1
        alpha = [alpha]
    nalpha = len(alpha)
    if nalpha > d['nalp_max']:
        raise ValueError('Too many pitch angles requested; {} is maximum.'
                         .format(d['nalp_max']))
    d['nalpha'] = nalpha
    if nalpha > 0:
        degalpha[0:nalpha] = alpha
    d['degalpha'] = degalpha
        
    return d





