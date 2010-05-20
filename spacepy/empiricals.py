#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
empirical models
"""

from spacepy import help
__version__ = "$Revision: 1.1 $, $Date: 2010/05/20 17:19:44 $"
__author__ = ['J. Koller, Los Alamos National Lab (jkoller@lanl.gov)',
'Steve Morley (smorley@lanl.gov/morley_steve@hotmail.com)']


# -----------------------------------------------
# Last Closed Drift Shell - Lmax
# -----------------------------------------------    
def get_Lmax(ticks, Lmax_model):
    """
    calculate a simple empirical model for Lmax
    
    """
    import numpy as n
    import spacepy.omni as om
        
    omni = om.getomni(ticks)
    Dst = omni['Dst']
    
    Lmax = n.zeros(len(Dst))
    
    if Lmax_model is 'JKemp':
        for i, iDst in enumerate(Dst):
            Lmax[i] = 6.7e-5*iDst*iDst + 0.0421*iDst + 8.945
        
    return Lmax

# -----------------------------------------------
# Plasmapause Location
# -----------------------------------------------       
def get_plasma_pause(ticks, Lpp_model):
    """
    Plasmapause location model(s) 
    """
    
    import numpy as np
    import datetime as dt
    import spacepy.omni as om
    import spacepy.time as spt
    
    if Lpp_model == 'CA2000': #Carpenter & Anderson (2000)
        omni = om.getomni(ticks)
        Kp = omni['Kp']    
        Lpp = np.zeros(len(Kp))
        for i, iKp in enumerate(Kp):
            Lpp[i] = 5.6 - 0.48*iKp
    
    if Lpp_model == 'M2002': #Moldwin et al. (2002)
        ticksplus = spt.tickrange(ticks.UTC[0]-dt.timedelta(hours=12), ticks.UTC[0], dt.timedelta(hours=3))
        ticks = spt.Ticktock()
        omni = om.getomni(ticks)
        prev3 = [np.nan, np.nan, np.nan]
        Lpp = np.zeros(len(Kp))
        for i, iKp in enumerate(Kp):
            prev3.append(iKp)
            Kpmax = max(prev3)
            Lpp[i] = 5.99 - 0.382*Kpmax
            prev3 = prev3[1:] 
        
    return Lpp

# -----------------------------------------------
# Magnetopause Location
# -----------------------------------------------    
    
def ShueMP(P,Bz):
    """Calculates the Shue et al. (1997) subsolar magnetopause radius
    
    Ported from Drew Turner's (LASP) MatLab script
    
    Inputs:
    =======
    SW ram pressure [nPa], IMF Bz (GSM) [nT]
    
    Output:
    =======
    Magnetopause (sub-solar point) standoff distance [Re]
    """
    
    import numpy as np
    
    # Initialize r0 and make P and Bz numpy arrays
    r0 = np.zeros((len(P)),dtype=float)
    Bz = np.array(Bz)
    P = np.array(P)
    
    # Find where Bz >= 0 and where it is < 0
    iBzPos = np.where(Bz >= 0)
    iBzNeg = np.where(Bz < 0)
    
    # Calculate r0
    r0[iBzPos] = (11.4 + 0.013*Bz[iBzPos])*P[iBzPos]**(-1/6.6)
    r0[iBzNeg] = (11.4 + 0.140*Bz[iBzNeg])*P[iBzNeg]**(-1/6.6)
    
    return r0