#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lstar and Lmax calculation using artificial neural network (ANN) technique.

Authors: Josef Koller, Yiqun Yu
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov, yiqun@lanl.gov

Copyright 2012 Los Alamos National Security, LLC.

"""
import os.path

try:
    import ffnet
except ImportError:
    raise RuntimeError(
        'LANLstar requires ffnet (http://ffnet.sourceforge.net/)')
import numpy as np


def _get_net_path(filename):
    """Gets the full path for a network file given the filename"""
    fspec = os.path.join(
        os.path.split(__file__)[0], 'data', 'LANLstar', filename)
    if os.path.exists(fspec):
        return fspec
    else:
        raise RuntimeError("Could not find neural network file " + filename)

# ------------------------------------------------
def _LANLcommon(inputdict, extMag, domax):
    """
    Shared code between LANLstar and LANLmax

    domax is True for LANLmax, False for LANLstar
    """
    lstar_keylists = {
        'OPDyn': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF',
                  'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
        'OPQuiet': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo','BzIMF',
                    'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
        'T01QUIET': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G1', 'G2','Lm', 'Bmirr', 'PA', 'rGSM',
                     'latGSM', 'lonGSM'],
        'T01STORM': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G2', 'G3', 'Lm', 'Bmirr', 'PA', 'rGSM',
                     'latGSM', 'lonGSM'],
        'T05': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF','BzIMF',
                'W1','W2','W3','W4','W5','W6',
                'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
        'T89': ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF','BzIMF',
                'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
        'T96': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF','BzIMF',
                'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
                 }
    lmax_keylists = {
        'OPDyn': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA'],
        'OPQuiet': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF','PA'],
        'T01QUIET': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G1','G2','PA'],
        'T01STORM': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G2','G3','PA'],
        'T05': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
                'W1','W2','W3','W4','W5','W6','PA'],
        'T89': ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF', 'BzIMF', 'PA'],
        'T96': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'PA'],
                 }
    lstar_nets = { 'OPDyn'   : 'LANLstar_OPDyn.net',
                    'OPQuiet' : 'LANLstar_OPQuiet.net',
                    'T01QUIET': 'LANLstar_T01QUIET.net',
                    'T01STORM': 'LANLstar_T01STORM.net',
                    'T05': 'LANLstar_T05.net',
                    'T89': 'LANLstar_T89.net',
                    'T96': 'LANLstar_T96.net',
                    }
    lmax_nets = {  'OPDyn'   : 'Lmax_OPDyn.net',
                    'OPQuiet' : 'Lmax_OPQuiet.net',
                    'T01QUIET': 'Lmax_T01QUIET.net',
                    'T01STORM': 'Lmax_T01STORM.net',
                    'T05': 'Lmax_T05.net',
                    'T89': 'Lmax_T89.net',
                    'T96': 'Lmax_T96.net',
                    }
    
    npt = len(inputdict['Year'])
    Lstar_out = {} 
    
    if isinstance(extMag, str): extMag = [extMag]

    for modelkey in extMag:
        if domax:
            keylist = lmax_keylists[modelkey]
        else:
            keylist = lstar_keylists[modelkey]
        #T89 checks Kp, everything else Dst
        specialkey = 'Dst' if 'Dst' in keylist else 'Kp'
        if isinstance(inputdict[specialkey], float):
            arrayflag = False
            for key in inputdict.keys():
                inputdict[key] = [inputdict[key]]
        else:
            arrayflag = True
	    ncalc = len(inputdict['Dst'])
	
        ncalc = len(inputdict[specialkey])
        Lstar = np.zeros(ncalc)
        inpar = np.zeros(len(keylist))

        if domax:
            netfile = lmax_nets[modelkey]
        else:
            netfile = lstar_nets[modelkey]
        network = ffnet.loadnet(_get_net_path(netfile))
        for i in range(ncalc):	
            # copy over keylist into inpar
            for ikey, key in enumerate(keylist):
    		inpar[ikey] = inputdict[key][i]
            Lstar[i] = network(inpar)
            
        if arrayflag is False:
            Lstar_out[modelkey] = Lstar[0]
        else:
            Lstar_out[modelkey] = Lstar

    return Lstar_out

# ------------------------------------------------
def LANLstar(inputdict, extMag):
    """
    This will calculate Lstar based on the L* artificial neural network (ANN) trained from 
    different magnetospheric field models.

    Parameters
    ==========
    extMag : list of string(s)
        containing one or more of the following external magnetic field models: 
        'OPDyn', 'OPQuiet', 'T89', 'T96', 'T01QUIET', 'T01STORM', 'T05'
    
    inputdict : dictionary
        containing the following keys, each entry is a list or array. Note the keys for the above models are different.

        -- For OPDyn:  
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']

        -- For OPQuiet:
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']       
        
        -- For T89:
          ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF', 'BzIMF', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']       
       
        -- For T96:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
       
        -- For T01QUIET:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'G1', 'G2', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']

        -- For T01STORM:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'G2', 'G3', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']

        -- For T05:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'W1','W2','W3','W4','W5','W6',
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']

       
       Dictionaries with numpy vectors are allowed.
            
    Returns
    =======
    out : dictionary
        Lstar array for each key which corresponds to the specified magnetic field model.

    Examples
    ========
    >>> import spacepy.LANLstar as LS
    >>>
    >>> inputdict = {}
    >>> inputdict['Kp']     = [2.7      ]            # Kp index
    >>> inputdict['Dst']    = [7.7777   ]            # Dst index (nT)
    >>> inputdict['dens']   = [4.1011   ]            # solar wind density (/cc)
    >>> inputdict['velo']   = [400.1011 ]            # solar wind velocity (km/s)
    >>> inputdict['Pdyn']   = [4.1011   ]            # solar wind dynamic pressure (nPa)
    >>> inputdict['ByIMF']  = [3.7244   ]            # GSM y component of IMF magnetic field (nT)  
    >>> inputdict['BzIMF']  = [-0.1266  ]            # GSM z component of IMF magnetic field (nT)
    >>> inputdict['G1']     = [1.029666 ]            # as defined in Tsganenko 2003 
    >>> inputdict['G2']     = [0.549334 ]
    >>> inputdict['G3']     = [0.813999 ]                
    >>> inputdict['W1']     = [0.122444 ]            # as defined in Tsyganenko and Sitnov 2005
    >>> inputdict['W2']     = [0.2514   ]                                                                                      
    >>> inputdict['W3']     = [0.0892   ]
    >>> inputdict['W4']     = [0.0478   ]
    >>> inputdict['W5']     = [0.2258   ]
    >>> inputdict['W6']     = [1.0461   ]
    >>>
    >>> inputdict['Year']   = [1996     ]
    >>> inputdict['DOY']    = [6        ]
    >>> inputdict['Hr']     = [1.2444   ]
    >>>
    >>> inputdict['Lm']     = [4.9360   ]             # McIllwain L
    >>> inputdict['Bmirr']  = [315.6202 ]             # magnetic field strength at the mirror point
    >>> inputdict['rGSM']   = [4.8341   ]             # radial coordinate in GSM [Re]
    >>> inputdict['lonGSM'] = [-40.2663 ]             # longitude coodrinate in GSM [deg]
    >>> inputdict['latGSM'] = [36.44696 ]             # latitude coordiante in GSM [deg]
    >>> inputdict['PA']     = [57.3874  ]             # pitch angle [deg]
    >>> 
    >>> LS.LANLstar(inputdict, ['OPDyn','OPQuiet','T01QUIET','T01STORM','T89','T96','T05'])
    {'OPDyn': array([4.7171]),
     'OPQuiet': array([4.6673]),
     'T01QUIET': array([4.8427]),
     'T01STORM': array([4.8669]), 
     'T89': array([4.5187]),
     'T96': array([4.6439]),
     'TS05': array([4.7174])}
    """
    return _LANLcommon(inputdict, extMag, False)

def LANLmax(inputdict, extMag):

    """
    This will calculate the last closed drift shell Lmax based on the L* artificial neural network (ANN) trained from 
    different magnetospheric field models.

    Parameters
    ==========
    extMag : list of string(s)
        containing one or more of the following external Magnetic field models: 
        'OPDyn', 'OPQuiet', 'T89', 'T96', 'T01QUIET', 'T01STORM', 'T05'
    
    inputdict : dictionary
        containing the following keys, each entry is a list or array. Note the keys for the above models are different.

        -- For OPDyn:  
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA']

        -- For OPQuiet:
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA']

        -- For T89:
          ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF', 'BzIMF', 'PA']
       
        -- For T96:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF','PA']
       
        -- For T01QUIET:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'G1', 'G2','PA']

        -- For T01STORM:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'G2', 'G3', 'PA']

        -- For T05:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'W1','W2','W3','W4','W5','W6', 'PA']

        Dictionaries with numpy vectors are allowed.
            
    Returns
    =======
    out : dictionary
        Lmax array for each key which corresponds to the specified magnetic field model.

    Examples
    ========
    >>> import spacepy.LANLstar as LS
    >>>
    >>> inputdict = {}
    >>> inputdict['Kp']     = [2.7      ]            # Kp index
    >>> inputdict['Dst']    = [7.7777   ]            # Dst index (nT)
    >>> inputdict['dens']   = [4.1011   ]            # solar wind density (/cc)
    >>> inputdict['velo']   = [400.1011 ]            # solar wind velocity (km/s)
    >>> inputdict['Pdyn']   = [4.1011   ]            # solar wind dynamic pressure (nPa)
    >>> inputdict['ByIMF']  = [3.7244   ]            # GSM y component of IMF magnetic field (nT)  
    >>> inputdict['BzIMF']  = [-0.1266  ]            # GSM z component of IMF magnetic field (nT)
    >>> inputdict['G1']     = [1.029666 ]            # as defined in Tsganenko 2003 
    >>> inputdict['G2']     = [0.549334 ]
    >>> inputdict['G3']     = [0.813999 ]                
    >>> inputdict['W1']     = [0.122444 ]            # as defined in Tsyganenko and Sitnov 2005
    >>> inputdict['W2']     = [0.2514   ]                                                                                      
    >>> inputdict['W3']     = [0.0892   ]
    >>> inputdict['W4']     = [0.0478   ]
    >>> inputdict['W5']     = [0.2258   ]
    >>> inputdict['W6']     = [1.0461   ]
    >>>
    >>> inputdict['Year']   = [1996     ]
    >>> inputdict['DOY']    = [6        ]
    >>> inputdict['Hr']     = [1.2444   ]
    >>>
    >>> inputdict['PA']     = [57.3874  ]             # pitch angle [deg]
    >>> 
    >>> LS.LANLmax(inputdict, ['OPDyn','OPQuiet','T01QUIET','T01STORM','T89','T96','T05'])
    {'OPDyn': array([10.6278]),
     'OPQuiet': array([9.3352]),
     'T01QUIET': array([10.0538]),
     'T03STORM': array([9.9300]), 
     'T89': array([8.2888]),
     'T96': array([9.2410]),
     'T05': array([9.9295])}
     """
    return _LANLcommon(inputdict, extMag, True)
