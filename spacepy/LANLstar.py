#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lstar and Lmax calculation using artificial neural network (ANN) technique.

Authors: Steve Morley, Josef Koller, Yiqun Yu, Aaron Hendry
Contact: smorley@lanl.gov, yiqunyu17@gmail.com

Copyright 2012 Los Alamos National Security, LLC.

.. autosummary::
    :toctree: autosummary

    LANLstar
    LANLmax
"""
import os.path
import sys
import warnings

import numpy as np
from . import toolbox
from . import datamodel as dm


def _get_net_path(filename):
    """Gets the full path for a network file given the filename"""
    fspec = os.path.join(
        os.path.split(__file__)[0], 'data', 'LANLstar', filename)
    if os.path.exists(fspec) or os.path.exists(fspec + '.txt'):
        return fspec
    else:
        raise RuntimeError("Could not find neural network file " + filename)


def _LANLcommon(indict, extmag, lmax=False):
    """
    Shared code between LANLstar and LANLmax

    lmax is True for LANLmax, False for LANLstar
    """
    n = len(indict['Year'])
    
    lstar_keylists = {
        'OPDYN': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF',
                  'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM'],
        'OPQUIET': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo','BzIMF',
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
        'RAMSCB':['Year','DOY','Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                  'PA','SMx', 'SMy', 'SMz'],
         }

    lmax_keylists = {
        'OPDYN': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA'],
        'OPQUIET': ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF','PA'],
        'T01QUIET': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G1','G2','PA'],
        'T01STORM': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF',
                     'G2','G3','PA'],
        'T05': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
                'W1','W2','W3','W4','W5','W6','PA'],
        'T89': ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF', 'BzIMF', 'PA'],
        'T96': ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'PA'],
                 }

    ls = dm.SpaceData()

    if isinstance(extmag, str):
        extmag = [extmag]

    # make sure that if we futz with Gs/Ws we need to not modify the inputs
    inputdict = indict.copy()
    if 'G' in inputdict:
        for x in range(1, 4):
            dum = inputdict['G'][..., x-1]
            if dum.ndim == 0: dum = np.array([dum])
            inputdict['G{0}'.format(x)] = dum
        del inputdict['G']
    if 'W' in inputdict:
        for x in range(1, 7):
            dum = inputdict['W'][..., x-1]
            if dum.ndim == 0: dum = np.array([dum])
            inputdict['W{0}'.format(x)] = dum
        del inputdict['W']

    parlist = lstar_keylists if not lmax else lmax_keylists
    for modname in extmag:
        # Concatenate the input parameters into single array
        parset = [inputdict[kk] for kk in parlist[modname]]
        if n == 1:
            params = np.hstack(parset)
        else:
            params = np.vstack(parset).T
        # Load the ANN parameters from the ported text files
        net = _get_model(modelstr=modname, lmax=lmax)
        # Make output array filled with badvals
        ls[modname] = dm.dmfilled((n), fillval=np.nan)

        # This can/should be done as a single matrix op. Need to check the dimensionality and how it all broadcasts...
        for idx, row in enumerate(np.atleast_2d(params)): 
            # This is the actual neural network calculation
            inlayer = (net['inweights'].T*row + net['inbias'].T)
            hidden = 1.0/(1+np.exp(-(net['ihbias'].T + (inlayer[::-1]*net['ihweights']).sum(axis=1))))
            if net['two_layers']:
                # If we have two hidden layers, we need to include them both
                hidden = 1.0/(1+np.exp(-(net['hhbias'].T + (hidden[:, np.newaxis]*net['hhweights']).sum(axis=0))))
            output = 1.0/(1+np.exp(-(net['hobias'] + (hidden*net['howeights'].T).sum())))
    
            # apply linear output function to get L*, ls[modname] is known 1D
            ls[modname][idx] = (net['outweight']*output + net['outbias']).item()

    return ls


def _get_model(modelstr=None, lmax=False):
    '''Find and load coefficient set defining neural network model

    Other Parameters
    ================
    modelstr : string
        Name of external magnetic field model. Valid choices are [OPDYN, OPQUIET,
        T89, T96, T01QUIET, T01STORM, T05, RAMSCB]. Note that RAMSCB is only defined
        for the LANLstar model and not for the LANLmax model.
    lmax : bool
        If True, load the coefficient set for the LANLmax network. Otherwise, load
        the coefficient set for the LANLstar network. Note that RAMSCB is not defined
        if this is set to True.
    '''
    name = modelstr.upper()
    optl = {'OPDYN': 'LANLstar_OPDyn.txt',
            'OPQUIET': 'LANLstar_OPQuiet.txt',
            'T01QUIET': 'LANLstar_T01QUIET.txt',
            'T01STORM': 'LANLstar_T01STORM.txt',
            'T05': 'LANLstar_T05.txt',
            'T89': 'LANLstar_T89.txt',
            'T96': 'LANLstar_T96.txt',
            'RAMSCB': 'LANLstar_RAMSCB.txt'
            }
    optm = {'OPDYN': 'Lmax_OPDyn.txt',
            'OPQUIET': 'Lmax_OPQuiet.txt',
            'T01QUIET': 'Lmax_T01QUIET.txt',
            'T01STORM': 'Lmax_T01STORM.txt',
            'T05': 'Lmax_T05.txt',
            'T89': 'Lmax_T89.txt',
            'T96': 'Lmax_T96.txt',
            }
    optdict = optl if not lmax else optm
    model = dm.readJSONheadedASCII(_get_net_path(optdict[name]))
    model['ihbias'] = model['ihbias'].T
    model['hobias'] = model['hobias'].T
    model['hhweights'] = model['hhweights'].T
    model['howeights'] = model['howeights'].T
    model['two_layers'] = model['two_layers'][0]
    model['two_layers'] = True if (model['two_layers'] in [1.0, 'True']) else False

    return model


# ------------------------------------------------
def LANLstar(inputdict, extMag):
    """
    Calculate Lstar

    Based on the L* artificial neural network (ANN) trained from 
    different magnetospheric field models.

    Parameters
    ==========
    extMag : list of string(s)
        containing one or more of the following external magnetic field models: 
        'OPDYN', 'OPQUIET', 'T89', 'T96', 'T01QUIET', 'T01STORM', 'T05'
    
    inputdict : dictionary
        containing the following keys, each entry is a list or array. Note the keys for the above models are different.

        -- For OPDYN:  
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 
          'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']

        -- For OPQUIET:
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

        -- For RAMSCB:
          ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
          'PA', 'SMx','SMy','SMz']
       
       Dictionaries with numpy vectors are allowed.
            
    Returns
    =======
    out : dictionary
        Lstar array for each key which corresponds to the specified magnetic field model.

    Examples
    ========
    >>> import spacepy.LANLstar as LS
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
    >>> # now add date
    >>> inputdict['Year']   = [1996     ]
    >>> inputdict['DOY']    = [6        ]
    >>> inputdict['Hr']     = [1.2444   ]
    >>> # and pitch angle, which doesn't come if taking params from OMNI
    >>> inputdict['Lm']     = [4.9360   ]             # McIllwain L
    >>> inputdict['Bmirr']  = [315.6202 ]             # magnetic field strength at the mirror point
    >>> inputdict['rGSM']   = [4.8341   ]             # radial coordinate in GSM [Re]
    >>> inputdict['lonGSM'] = [-40.2663 ]             # longitude coodrinate in GSM [deg]
    >>> inputdict['latGSM'] = [36.44696 ]             # latitude coordiante in GSM [deg]
    >>> inputdict['PA']     = [57.3874  ]             # pitch angle [deg]
    >>> inputdict['SMx']    = [3.9783   ]
    >>> inputdict['SMy']    = [-2.51335 ]
    >>> inputdict['SMz']    = [1.106617 ]
    >>> # and then call the neural network
    >>> LS.LANLstar(inputdict, ['OPDYN','OPQUIET','T01QUIET','T01STORM','T89','T96','T05','RAMSCB'])
    {'OPDYN': array([4.7171]),
     'OPQUIET': array([4.6673]),
     'T01QUIET': array([4.8427]),
     'T01STORM': array([4.8669]), 
     'T89': array([4.5187]),
     'T96': array([4.6439]),
     'TS05': array([4.7174]),
     'RAMSCB','array([5.9609])}
     """
    return _LANLcommon(inputdict, extMag, lmax=False)


def LANLmax(inputdict, extMag):
    """
    Calculate last closed drift shell (Lmax)

    Based on the L* artificial neural network (ANN) trained from 
    different magnetospheric field models.

    Parameters
    ==========
    extMag : list of string(s)
        containing one or more of the following external Magnetic field models: 
        'OPDYN', 'OPQUIET', 'T89', 'T96', 'T01QUIET', 'T01STORM', 'T05'
    
    inputdict : dictionary
        containing the following keys, each entry is a list or array. Note the keys for the above models are different.

        -- For OPDYN:  
          ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA']

        -- For OPQUIET:
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
    >>> # now add date
    >>> inputdict['Year']   = [1996     ]
    >>> inputdict['DOY']    = [6        ]
    >>> inputdict['Hr']     = [1.2444   ]
    >>> # and pitch angle, which doesn't come if taking params from OMNI
    >>> inputdict['PA']     = [57.3874  ]             # pitch angle [deg]
    >>> # and then call the neural network
    >>> LS.LANLmax(inputdict, ['OPDYN','OPQUIET','T01QUIET','T01STORM','T89','T96','T05'])
    {'OPDYN': array([10.6278]),
     'OPQUIET': array([9.3352]),
     'T01QUIET': array([10.0538]),
     'T01STORM': array([9.9300]), 
     'T89': array([8.2888]),
     'T96': array([9.2410]),
     'T05': array([9.9295])}
     """
    return _LANLcommon(inputdict, extMag, lmax=True)

def addPA(indict, PA):
    '''Function to add pitch angle to input dictionary from, e.g., omni module

    Parameters
    ==========
    indict : dictionary-like
        containing keys required for LANLstar and LANLmax
    
    PA : float
        pitch angle

    Returns
    =======
    out : dictionary
        input dictionary with input pitch angle added as 'PA' key having the length of other inputs

    Examples
    ========
    >>> import spacepy.LANLstar as LS
    >>> inputdict = {}
    >>> inputdict['Year']
    >>> # additional keys would be defined here
    >>> LS.addPA(indict, PA)

    '''
    ll = len(indict['Year'])
    indict['PA'] = [PA]*ll
    return indict
