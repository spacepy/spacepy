#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module with some useful empirical models (plasmapause, magnetopause, Lmax)


Authors
-------
Steve Morley, Josef Koller

Copyright ©2010 Los Alamos National Security, LLC.
"""
import datetime

from spacepy import help

import spacepy.omni as om
import spacepy.time as spt
import numpy as np

def getLmax(ticks, model='JKemp'):
    """
    calculate a simple empirical model for Lmax - last closed drift-shell

    @param ticks: Ticktock object of desired times
    @type ticks: spacepy.time.Ticktock

    @keyword model: 'JKemp' (default - empirical model of J. Koller)

    @return: Lmax - L* of last closed drift shell
    @rtype: numpy.ndarray

    """

    omni = om.get_omni(ticks)
    Dst = omni['Dst']

    Lmax = np.zeros(len(Dst))

    if model is 'JKemp':
        for i, iDst in enumerate(Dst):
            Lmax[i] = 6.07e-5*iDst*iDst + 0.0436*iDst + 9.37
    else:
        raise(ValueError('Invalid model selected'))

    return Lmax

def getPlasmaPause(ticks, model='M2002', LT='all'):
    """
    Plasmapause location model(s)

    Input:
    ======
        - ticks: TickTock object of desired times
        - Lpp_model (kwarg): 'CA1992' or 'M2002' (default)
    CA1992 returns the Carpenter and Anderson model,
    M2002 returns the Moldwin et al. model
        - LT: requested local time sector, 'all' is valid option

    Returns:
    ========
        - Lpp: Plasmapause radius in Earth radii

    Example:
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
    >>> emp.getPlasmaPause(ticks)
    array([ 4.9586,  4.9586,  4.9586,  4.9586,  4.9586,  4.9586,  4.9586,
        5.1114,  5.608 ,  5.7226,  5.7226])

    """

    import datetime as dt
    import spacepy.toolbox as tb

    def calcLpp(Kpmax, A, B):
        currLpp = A - B*Kpmax
        return currLpp

    model_list = ['CA1992', 'M2002']

    if model == 'CA1992':
        if LT!='all':
            print('No LT dependence currently supported for this model')
    if model not in model_list:
        raise ValueError("Please specify a valid model:\n'M2002' or 'CA1992'")

    if LT=='all':
        parA = {'CA1992': 5.6, 'M2002': 5.39}
        parB = {'CA1992': 0.46, 'M2002': -0.382}
        priorvals = {'CA1992': datetime.timedelta(hours=24),
                     'M2002': datetime.timedelta(hours=12)}
        A, B = parA[model], parB[model]
        prior = priorvals[model]
    else:
        try:
            float(LT)
        except (ValueError, TypeError):
            raise ValueError("Please specify a valid LT:\n'all' or a numeric type")
        parA = {'CA1992': [5.6]*24,
                'M2002': [5.7]*3+[6.05]*6+[5.2]*6+[4.45]*6+[5.7]*3}
        parB = {'CA1992': [0.46]*24,
                'M2002': [-0.42]*3+[-0.573]*6+[-0.425]*6+[-0.167]*6+[-0.42]*3}
        priorvals = {'CA1992': [datetime.timedelta(hours=24)]*24,
                     'M2002': [datetime.timedelta(hours=12)]*24}
        LThr = long(LT)
        prior = priorvals[model][LThr]
        A, B = parA[model][LThr], parB[model][LThr]

    #TODO: allow calling with ticks as dict of Kp (will also need UT for Kpmax)
    st, en = ticks.UTC[0]-prior, ticks.UTC[-1]
    einds, oinds = tb.tOverlap([st, en], om.omnidata['UTC'])
    utc = np.array(om.omnidata['UTC'])[oinds]
    Kp = np.array(om.omnidata['Kp'])[oinds]
    Lpp = np.zeros(len(ticks))

    for i, t1 in enumerate(ticks.UTC):
        t0 = t1-prior
        iprevday, dum = tb.tOverlap(utc, [t0, t1])
        if iprevday:
            Kpmax = max(Kp[iprevday])
            Lpp[i] = calcLpp(Kpmax, A, B)
        else:
            Lpp[i] = np.nan

    return Lpp

def getMPstandoff(ticks):
    """Calculates the Shue et al. (1997) subsolar magnetopause radius

    Inputs:
    =======
        - ticks: TickTock object of desired times
    (will be interpolated from hourly OMNI data)
    OR
        - dictionary of form {'P': [], 'Bz': []}
    Where P is SW ram pressure [nPa] and Bz is IMF Bz (GSM) [nT]

    Returns:
    ========
    Magnetopause (sub-solar point) standoff distance [Re]

    Example:
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
    >>> emp.ShueMP(ticks)
    array([ 10.57319532,  10.91327759,  10.75086872,  10.77577211,
         9.7818026 ,  11.03744739,  11.4065    ,  11.27555453,
        11.47988569,  11.82025823,  11.23834817])
    >>> data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
    >>> emp.ShueMP(data)
    array([ 10.29156018,   8.96790412])

    """
    if type(ticks) == spt.Ticktock:
        omni = om.get_omni(ticks)
        P, Bz = omni['Pdyn'], omni['BzIMF']
    elif type(ticks) == dict:
        P, Bz = ticks['P'], ticks['Bz']
        try:
            len(P)
        except TypeError:
            P = [P]
        try:
            len(Bz)
        except TypeError:
            Bz = [Bz]
    else:
        raise(TypeError('Invalid Input type'))

    try:
        # Initialize r0 and make P and Bz numpy arrays
        r0 = np.zeros((len(P)), dtype=float)
        Bz = np.array(Bz)
        P = np.array(P)

        # Find where Bz >= 0 and where it is < 0
        iBzPos = np.where(Bz >= 0)
        iBzNeg = np.where(Bz < 0)

        # Calculate r0
        r0[iBzPos] = (11.4 + 0.013*Bz[iBzPos])*P[iBzPos]**(-1/6.6)
        r0[iBzNeg] = (11.4 + 0.140*Bz[iBzNeg])*P[iBzNeg]**(-1/6.6)

        return r0
    except TypeError:
        raise TypeError("Please check for valid input types")


def getDststar(ticks, model='OBrien'):
    """Calculate the pressure-corrected Dst index, Dst*

    Inputs:
    =======
        - ticks: TickTock object of desired times
    (will be interpolated from hourly OMNI data)
    OR
        - dictionary including 'Pdyn' and 'Dst' keys where data are lists or arrays
          and Dst is in [nT], and Pdyn is in [nPa]

    Returns:
    ========
    Dst* - the pressure corrected Dst index from OMNI [nT]

    Use:
    ====
    Coefficients are applied to the standard formulation e.g. Burton et al., 1975
    of Dst* = Dst - b*sqrt(Pdyn) + c
    The default is the O'Brien and McPherron model (2002).
    Other options are Burton et al. (1975) and Borovsky and Denton (2010)

    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> ticks = spt.tickrange('2000-10-16T00:00:00', '2000-10-31T12:00:00', 1/24.)
    >>> dststar = getDststar(ticks)

    User-determined coefficients can also be supplied as a two-element list
    or tuple of the form (b,c), e.g.

    >>> dststar = getDststar(ticks, model=(2,11)) #b is extreme driving from O'Brien

    We have chosen the OBrien model as the default here as this was rigorously
    determined from a very long data set and is pertinent to most conditions.
    It is, however, the most conservative correction. Additionally, Siscoe,
    McPherron and Jordanova (2005) argue that the pressure contribution to Dst diminishes
    during magnetic storms.

    To show the relative differences, run the following example:

    >>> params = [('Burton','k-'), ('OBrien','r-'), ('Borovsky','b-')]
    >>> for model, col in params:
            dststar = getDststar(ticks, model=model)
            plt.plot(ticks.UTC, dststar, col)

    """
    model_params = {'Burton': (15.8, 20),
                    'OBrien': (7.26, 11),
                    'Borovsky': (20.7,27.7)}

    if isinstance(model, str):
        try:
            b, c = model_params[model]
        except KeyError:
            raise ValueError('Invalid pressure correction model selected')
    else:
        try:
            b, c = model[0], model[1]
        except (KeyError, IndexError):
            raise ValueError('Invalid coefficients set: must be of form (b,c)')

    if isinstance(ticks, spt.Ticktock):
        omni = om.get_omni(ticks)
        P, Dst = omni['Pdyn'], omni['Dst']
    elif isinstance(ticks, dict):
        P, Dst = ticks['Pdyn'], ticks['Dst']
        if isinstance(P, list):
            P, Dst = np.array(P), np.array(Dst)

    #get Dst*
    Dststar = Dst - b*P**0.5 + c

    return Dststar


ShueMP = getMPstandoff
get_plasma_pause = getPlasmaPause
get_Lmax = getLmax
