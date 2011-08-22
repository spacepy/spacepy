#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module with some useful empirical models (plasmapause, magnetopause, Lmax)


Authors: Steve Morley, Josef Koller
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov

Copyright ©2010 Los Alamos National Security, LLC.
"""
import datetime

import datetime

import numpy as np

from spacepy import help
import spacepy.toolbox as tb
import spacepy.omni as om
import spacepy.time as spt

def getLmax(ticks, model='JKemp'):
    """
    calculate a simple empirical model for Lmax - last closed drift-shell

    What is the paper this model is from?  Put it here!

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        Ticktock object of desired times
    model : string, optional
        'JKemp' (default - empirical model of J. Koller)

    Returns
    =======
    out : np.ndarray
        Lmax - L* of last closed drift shell

    Examples
    ========
    >>> from spacepy.empiricals import getLmax
    >>> import spacepy.time as st
    >>> import datetime
    >>> ticks = st.tickrange(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 1, 3), deltadays=1)
    array([ 7.4928412,  8.3585632,  8.6463423])


    """
    omni = om.get_omni(ticks)
    Dst = omni['Dst']
    Lmax = np.zeros(len(Dst))
    if model is 'JKemp':
        for i, iDst in enumerate(Dst):
            Lmax[i] = 6.07e-5*iDst*iDst + 0.0436*iDst + 9.37
    else:
        raise ValueError('Invalid model selection')
    return Lmax

def getPlasmaPause(ticks, model='M2002', LT='all'):
    """
    Plasmapause location model(s)

    We need to list the references here!

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        TickTock object of desired times
    Lpp_model : string, optional
        'CA1992' or 'M2002' (default)
        CA1992 returns the Carpenter and Anderson model,
        M2002 returns the Moldwin et al. model
    LT : int, float
        requested local time sector, 'all' is valid option

    Returns
    =======
    out : float
        Plasmapause radius in Earth radii

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
    >>> emp.getPlasmaPause(ticks)
    array([ 6.42140002,  6.42140002,  6.42140002,  6.42140002,  6.42140002,
        6.42140002,  6.42140002,  6.26859998,  5.772     ,  5.6574    ,
        5.6574    ])
    """
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

    Lets put the full reference here

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        TickTock object of desired times (will be interpolated from hourly OMNI data)
        OR dictionary of form {'P': [], 'Bz': []}
        Where P is SW ram pressure [nPa] and Bz is IMF Bz (GSM) [nT]

    Returns
    =======
    out : float
        Magnetopause (sub-solar point) standoff distance [Re]

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
    >>> emp.ShueMP(ticks)
    array([ 10.57319537,  10.91327764,  10.75086873,  10.77577207,
         9.78180261,  11.0374474 ,  11.4065    ,  11.27555451,
        11.47988573,  11.8202582 ,  11.23834814])
    >>> data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
    >>> emp.ShueMP(data)
    array([ 9.96096838,  8.96790412])
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

    We need to add in the references to the models here!

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        TickTock object of desired times (will be interpolated from hourly OMNI data)
        OR dictionary including 'Pdyn' and 'Dst' keys where data are lists or arrays
        and Dst is in [nT], and Pdyn is in [nPa]

    Returns
    =======
    out : float
        Dst* - the pressure corrected Dst index from OMNI [nT]

    Examples
    ========
    Coefficients are applied to the standard formulation e.g. Burton et al., 1975
    of Dst* = Dst - b*sqrt(Pdyn) + c
    The default is the O'Brien and McPherron model (2002).
    Other options are Burton et al. (1975) and Borovsky and Denton (2010)

    >>> import spacepy.time as spt
    >>> import spacepy.omni as om
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2000-10-16T00:00:00', '2000-10-31T12:00:00', 1/24.)
    >>> dststar = emp.getDststar(ticks)
    >>> dststar[0]
    -21.317220132108943

    User-determined coefficients can also be supplied as a two-element list
    or tuple of the form (b,c), e.g.

    >>> dststar = emp.getDststar(ticks, model=(2,11)) #b is extreme driving from O'Brien

    We have chosen the OBrien model as the default here as this was rigorously
    determined from a very long data set and is pertinent to most conditions.
    It is, however, the most conservative correction. Additionally, Siscoe,
    McPherron and Jordanova (2005) argue that the pressure contribution to Dst diminishes
    during magnetic storms.

    To show the relative differences, run the following example:

    >>> import matplotlib.pyplot as plt
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
