#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module with some useful empirical models (plasmapause, magnetopause, Lmax)


Authors: Steve Morley, Josef Koller
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.
"""
import datetime
from functools import partial
import numpy as np

import scipy.integrate as integ
from spacepy import help
import spacepy.toolbox as tb
import spacepy.omni as om
import spacepy.time as spt

__contact__ = 'Steve Morley, smorley@lanl.gov'

def getLmax(ticks, model='JKemp'):
    """
    calculate a simple empirical model for Lmax - last closed drift-shell

    Uses the parametrized Lmax from:
    Koller and Morley (2010)
    'Magnetopause shadowing effects for radiation belt models during high-speed solar wind streams'
    American Geophysical Union, Fall Meeting 2010, abstract #SM13A-1787

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

def getPlasmaPause(ticks, model='M2002', LT='all', omnivals=None):
    """
    Plasmapause location model(s)

    CA1992 -- Carpenter, D. L., and R. R. Anderson, An ISEE/whistler 
    model of equatorial electron density in the magnetosphere, 
    J. Geophys. Res., 97, 1097, 1992.
    M2002 -- Moldwin, M. B., L. Downward, H. K. Rassoul, R. Amin, 
    and R. R. Anderson, A new model of the location of the plasmapause: 
    CRRES results, J. Geophys. Res., 107(A11), 1339, 
    doi:10.1029/2001JA009211, 2002.

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
    omnivals : spacepy.datamodel.SpaceData, dict
        dict-like containing UTC (datetimes) and Kp keys

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
        parB = {'CA1992': 0.46, 'M2002': 0.382}
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
                'M2002': [0.42]*3+[0.573]*6+[0.425]*6+[0.167]*6+[0.42]*3}
        priorvals = {'CA1992': [datetime.timedelta(hours=24)]*24,
                     'M2002': [datetime.timedelta(hours=12)]*24}
        LThr = long(LT)
        prior = priorvals[model][LThr]
        A, B = parA[model][LThr], parB[model][LThr]

    #TODO: allow calling with ticks as dict of Kp (will also need UT for Kpmax)
    st, en = ticks.UTC[0]-prior, ticks.UTC[-1]
    if omnivals is None:
        omdat = om.get_omni(spt.tickrange(st, en, 1.0/24.0), dbase='QDhourly')
    else:
        #now test for sanity of input
        try:
            isinstance(omnivals, dict)
        except:
            raise TypeError('Not a valid input type for omnivals, expected spacepy.datamodel.SpaceData')
        try:
            assert 'UTC' in omnivals
            assert 'Kp' in omnivals
        except:
            raise KeyError('Required data not found in input dict-like (omnivals)')
        omdat = omnivals

    einds, oinds = tb.tOverlap([st, en], omdat['UTC'])
    utc = np.array(omdat['UTC'])[oinds]
    Kp = np.array(omdat['Kp'])[oinds]
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
    elif isinstance(ticks, dict): 
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


def vampolaPA(omniflux, **kwargs):
    '''Pitch angle model of sin^n form

    Parameters
    ==========
    omniflux : arraylike or float
        omnidirectional number flux data

    order : integer or float (optional)
        order of sin^n functional form for distribution (default=2)

    alphas : arraylike (optional)
        pitch angles at which to evaluate the differential number flux
        (default is 5 to 90 degrees in 36 steps)

    Returns
    =======
    dnflux : array
        differential number flux corresponding to pitch angles alphas
    alphas : array
        pitch angles at which the differential number flux was evaluated

    Examples
    ========
    Omnidirectional number flux of [3000, 6000]

    >>> from spacepy.empiricals import pamodel
    >>> pamodel.vampolaPA(3000, alpha=[45, 90])
    (array([  954.92965855,  1909.8593171 ]), [45, 90])
    >>> data, pas = pamodel.vampolaPA([3000, 6000], alpha=[45, 90])
    >>> pas
    [45, 90]
    >>> data
    array([[  954.92965855,  1909.8593171 ],
           [ 1909.8593171 ,  3819.71863421]])
    '''
    defaults = {'order': 2,
                'alpha': tb.linspace(5,90,18)}

    if hasattr(omniflux, '__iter__'):
        omniflux = np.asanyarray(omniflux)
    else:
        omniflux = np.asanyarray([omniflux])

    #substitute defaults
    for key in defaults:
        if key not in kwargs:
            kwargs[key] = defaults[key]

    if hasattr(kwargs['order'], '__iter__'):
        try:
            assert len(kwargs['order'])==len(omniflux)
        except AssertionError:
            raise ValueError('order must be either single-valued or the same length as omniflux')
    else:
        kwargs['order'] = np.asanyarray([kwargs['order']]*len(omniflux))
    normfac = np.empty(len(kwargs['order']), dtype=float)
    
    def sinfunc(x, order=kwargs['order']): #define distribution function
        dum = np.sin(x)
        return dum**order

    for idx, tmporder in enumerate(kwargs['order']):
        #def partial function so that  order is fixed
        sinfunc_o = partial(sinfunc, order=tmporder)
        normfac[idx] = integ.quad(sinfunc_o, 0, np.pi)[0] #quad returns (val, accuracy)

    #now make the differential number flux
    dnflux = np.zeros((len(kwargs['alpha']), len(omniflux))).squeeze()
    for i, a_val in enumerate(np.deg2rad(kwargs['alpha'])):
        dnflux[i] = omniflux * sinfunc(a_val) / normfac

    return dnflux, kwargs['alpha']


def getVampolaOrder(L):
    '''Empirical lookup of power for sin^n pitch angle model from Vampola (1996)

    Vampola, A.L. Outer zone energetic electron environment update, 
    Final Report of ESA/ESTEC/WMA/P.O. 151351, ESA-ESTEC, Noordwijk, 
    The Netherlands, 1996.

    Parameters
    ==========
    L : arraylike or float
        
    Returns
    =======
    order : array
        coefficient for sin^n model corresponding to McIlwain L (computed for OP77?)

    '''
    lmc = np.arange(3,8.00001,0.25)
    vamp_n = [5.38, 5.078, 4.669, 3.916, 3.095, 2.494, 2.151, 1.998, 1.899,
              1.942, 1.974, 1.939, 1.970, 2.136, 1.775, 1.438, 1.254, 1.194,
              1.046, 0.989, 0.852]
    
    if not hasattr(L, '__iter__'): L = [L]
    L = np.asanyarray(L)
    #if outside valid range, use end value
    L[L<=3] = 3
    L[L>=8] = 8
    #interpolate to get order for the given L
    order = np.interp(L, lmc, vamp_n)

    return order


ShueMP = getMPstandoff
get_plasma_pause = getPlasmaPause
get_Lmax = getLmax
