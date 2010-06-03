#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
empirical models
"""

from spacepy import help
__version__ = "$Revision: 1.4 $, $Date: 2010/06/03 17:26:55 $"
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
        
    omni = om.get_omni(ticks)
    Dst = omni['Dst']
    
    Lmax = n.zeros(len(Dst))
    
    if Lmax_model is 'JKemp':
        for i, iDst in enumerate(Dst):
            Lmax[i] = 6.7e-5*iDst*iDst + 0.0421*iDst + 8.945
        
    return Lmax

# -----------------------------------------------
# Plasmapause Location
# -----------------------------------------------       
def get_plasma_pause(ticks, Lpp_model='M2002'):
    """
    Plasmapause location model(s)
    
    Input:
    ======
    - ticks: TickTock object of desired times
    - Lpp_model (kwarg): 'CA1992' or 'M2002' (default)
    CA1992 returns the Carpenter and Anderson model,
    M2002 returns the Moldwin et al. model
    
    Authors:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    Steve Morley, Los Alamos National Lab (smorley@lanl.gov)
    
    """
    
    import numpy as np
    import datetime as dt
    import spacepy.omni as om
    import spacepy.time as spt
    import spacepy.toolbox as tb
    
    if Lpp_model == 'CA1992': #Carpenter & Anderson (1992)
        ticksplus = spt.tickrange(ticks.UTC[0]-dt.timedelta(hours=24), 
            ticks.UTC[-1], dt.timedelta(hours=3), 'UTC')
        omni = om.get_omni(ticksplus)
        Kp = omni['Kp']
        prev6 = [np.nan]*6
        Lpp = np.zeros(len(Kp))
        for i, iKp in enumerate(Kp):
            prev6.append(iKp)
            Kpmax = max(prev6)
            Lpp[i] = 5.6 - 0.46*Kpmax
            prev6 = prev6[1:]
        [tpint, tintp] = tb.tCommon(ticksplus.UTC, ticks.UTC)
    
    if Lpp_model == 'M2002': #Moldwin et al. (2002)
        ticksplus = spt.tickrange(ticks.UTC[0]-dt.timedelta(hours=12), 
            ticks.UTC[-1], dt.timedelta(hours=3), 'UTC')
        omni = om.get_omni(ticksplus)
        Kp = omni['Kp']
        prev3 = [np.nan]*3
        Lpp = np.zeros(len(Kp))
        for i, iKp in enumerate(Kp):
            prev3.append(iKp)
            Kpmax = max(prev3)
            Lpp[i] = 5.99 - 0.382*Kpmax
            prev3 = prev3[1:]
        [tpint, tintp] = tb.tCommon(ticksplus.UTC, ticks.UTC)
        
    return Lpp[tpint]

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