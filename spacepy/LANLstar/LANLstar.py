#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Lstar and Lmax calculation using artificial neural network (ANN) technique.

Authors: Josef Koller, Yiqun Yu
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov, yiqun@lanl.gov

Copyright ?2012 Los Alamos National Security, LLC.

"""

import libLANLstar
import numpy as np

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
    
    Bmodel = {'OPDyn':_LANLstar_OPDyn, 'OPQuiet':_LANLstar_OPQuiet, 'T01QUIET':_LANLstar_T01Quiet,
              'T01STORM':_LANLstar_T01Storm,'T05':_LANLstar_TS05, 'T89':_LANLstar_T89, 'T96':_LANLstar_T96}

    npt = len(inputdict['Year'])
    Lstar_out = {} 
    
    if isinstance(extMag, str): extMag = [extMag]

    for key in extMag:
        Lstar_out[key] = np.zeros(npt)
        Lstar_out[key] = Bmodel[key](inputdict)

    return Lstar_out

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

    Bmodel = {'OPDyn':_LANLmax_OPDyn, 'OPQuiet':_LANLmax_OPQuiet, 'T01QUIET':_LANLmax_T01Quiet,
              'T01STORM':_LANLmax_T01Storm, 'T05':_LANLmax_TS05,'T89':_LANLmax_T89, 'T96':_LANLmax_T96}
    
    npt = len(inputdict['Year'])
    Lmax_out = {} 
    
    if isinstance(extMag, str): extMag = [extMag]
    
    for key in extMag:
        Lmax_out[key] = np.zeros(npt)
        Lmax_out[key] = Bmodel[key](inputdict)

    return Lmax_out

# ------------------------------------------------------
def _LANLstar_OPDyn(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Olson-Pfitzer Dynamic model.
    
    """    
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
	
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo','BzIMF', \
                   'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
    	# copy over keylist into inpar
    	for ikey, key in enumerate(keylist):
    		inpar[ikey] = inputdict[key][i]
    	Lstar[i] = libLANLstar.lanlstar_opdyn(inpar)

    if arrayflag is False:
    	return Lstar[0]
    else:
        return Lstar
    

# ----------------------------------------------------------------
def _LANLstar_OPQuiet(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Olson-Pfitzer Quiet model.
    
    """    	
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
    	
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
    	
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo','BzIMF', \
                   'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
    	
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_opquiet(inpar)
    
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLstar_T01Quiet(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 2001 model.
    
    """
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
        
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
                   'G1', 'G2','Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
        
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_t01quiet(inpar)
    
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLstar_T01Storm(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 2003  model.

    """

    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF','BzIMF', \
                   'G2', 'G3', 'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_t03storm(inpar)
        
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLstar_TS05(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Tsyganenko & Sitnov (2005) model.
 
    """
        
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
        
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF','BzIMF', \
                   'W1','W2','W3','W4','W5','W6', \
                   'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
        
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_ts05(inpar)
        
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLstar_T89(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 1989 model.
    
    """
        
    if isinstance(inputdict['Kp'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Kp'])
    Lstar = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF','BzIMF', \
                   'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_t89(inpar)
    
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLstar_T96(inputdict):
    """
    This will calculate Lstar based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 1996 model.
    
    """        
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lstar = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF','BzIMF', \
                   'Lm', 'Bmirr', 'PA', 'rGSM', 'latGSM', 'lonGSM']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lstar[i] = libLANLstar.lanlstar_t96(inpar)
        
    if arrayflag is False:
        return Lstar[0]
    else:
        return Lstar
    
    
# --------------------------------------------------------------
def _LANLmax_OPDyn(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Olson-Pfitzer Dynamic model.
        
    """    
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF', 'PA']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_opdyn(inpar)
        
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax
    
    
#--------------------------------------------------------------
def _LANLmax_OPQuiet(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Olson-Pfitzer Quiet model.
    
    """
    
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'dens', 'velo', 'BzIMF','PA']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_opquiet(inpar)
        
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax
    
    
# --------------------------------------------------------------
def _LANLmax_T01Quiet(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 2001 model.
        
    """    
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
                   'G1','G2','PA']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_t01quiet(inpar)
    
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax
    
    
#--------------------------------------------------------------
def _LANLmax_T01Storm(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 2003 model.

    """
    
    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
                   'G2','G3','PA']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_t03storm(inpar)
        
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax

# --------------------------------------------------------------
def _LANLmax_TS05(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which 
    was trained with the Tsyganenko & Sitnov (2005) model.
    """

    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
               'W1','W2','W3','W4','W5','W6','PA']
    inpar = np.zeros(len(keylist))
        
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_ts05(inpar)
    
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax
        
# --------------------------------------------------------------
def _LANLmax_T89(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 1989 model.
        
    """
        
    if isinstance(inputdict['Kp'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Kp'])
    Lmax = np.zeros(ncalc)
    
    keylist = ['Year', 'DOY', 'Hr', 'Kp', 'Pdyn', 'ByIMF', 'BzIMF', 'PA']
    inpar = np.zeros(len(keylist))
        
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_t89(inpar)
    
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax
    
    
# --------------------------------------------------------------
def _LANLmax_T96(inputdict):
    """
    This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
    was trained with the Tsyganenko 1996 model.
        
    """

    if isinstance(inputdict['Dst'], float):
        arrayflag = False
        for key in inputdict.keys():
            inputdict[key] = [inputdict[key]]
    else:
        arrayflag = True
        
    ncalc = len(inputdict['Dst'])
    Lmax = np.zeros(ncalc)
        
    keylist = ['Year', 'DOY', 'Hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 'PA']
    inpar = np.zeros(len(keylist))
    
    for i in range(ncalc):	
        # copy over keylist into inpar
        for ikey, key in enumerate(keylist):
            inpar[ikey] = inputdict[key][i]
        Lmax[i] = libLANLstar.lanlmax_t96(inpar)
        
    if arrayflag is False:
        return Lmax[0]
    else:
        return Lmax

    
