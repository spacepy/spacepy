#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module with some useful empirical models (plasmapause, magnetopause, Lmax)


Authors: Steve Morley, Josef Koller
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.
"""
from __future__ import division
import datetime
import warnings

from functools import partial
import numpy as np
import scipy.integrate as integ

from spacepy import help
import spacepy.datamodel as dm
import spacepy.toolbox as tb
import spacepy.omni as om
import spacepy.time as spt

__contact__ = 'Steve Morley, smorley@lanl.gov'

def getLmax(ticks, model='JKemp', dbase='QDhourly'):
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

    See Also
    ========
    spacepy.LANLstar.LANLmax

    """
    omni = om.get_omni(ticks, dbase=dbase)
    Dst = omni['Dst']
    Lmax = np.zeros(len(Dst))
    if model == 'JKemp':
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
    RT1970 -- Rycroft, M. J., and J. O. Thomas, The magnetospheric
    plasmapause and the electron density trough at the alouette i
    orbit, Planetary and Space Science, 18(1), 65-80, 1970


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

    Warns
    =====
    RuntimeWarning
        If the CA1992 model is called with LT as it is not implemented 

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
    def calcLpp(Kpmax, A, B, power=1):
        currLpp = A - B*Kpmax**power
        return currLpp

    model_list = ['CA1992', 'M2002', 'RT1970']

    if model == 'CA1992':
        if LT != 'all':
            warnings.warn('No LT dependence currently supported for CA1992 model', RuntimeWarning)
    if model not in model_list:
        raise ValueError("Please specify a valid model:\n{0}".format(' or '.join(model_list)))

    if LT=='all':
        parA = {'CA1992': 5.6,  'M2002': 5.39,  'RT1970': 5.64}
        parB = {'CA1992': 0.46, 'M2002': 0.382, 'RT1970': 0.78}
        priorvals = {'CA1992': datetime.timedelta(hours=24),
                     'M2002': datetime.timedelta(hours=12),
                     'RT1970': datetime.timedelta(0)}
        A, B = parA[model], parB[model]
        prior = priorvals[model]
    else:
        try:
            float(LT)
        except (ValueError, TypeError):
            raise ValueError("Please specify a valid LT:\n'all' or a numeric type")
        parA = {'CA1992': [5.6]*24,
                'M2002': [5.7]*3+[6.05]*6+[5.2]*6+[4.45]*6+[5.7]*3,
                'RT1970': [5.64]*24}
        parB = {'CA1992': [0.46]*24,
                'M2002': [0.42]*3+[0.573]*6+[0.425]*6+[0.167]*6+[0.42]*3,
                'RT1970': [0.78]*24}
        priorvals = {'CA1992': [datetime.timedelta(hours=24)]*24,
                     'M2002': [datetime.timedelta(hours=12)]*24,
                     'RT1970': [datetime.timedelta(0)]*24}
        try:
            LThr = long(LT)
        except NameError:
            LThr = int(LT)
        prior = priorvals[model][LThr]
        A, B = parA[model][LThr], parB[model][LThr]

    st, en = ticks.UTC[0]-prior, ticks.UTC[-1]
    if omnivals is None:
        omdat = om.get_omni(spt.tickrange(st, en, 1.0/24.0), dbase='QDhourly')
    else:
        #now test for sanity of input
        try:
            assert isinstance(omnivals, dict)
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

    if model == 'RT1970':
        power = 0.5
    else:
        power = 1
    for i, t1 in enumerate(ticks.UTC):
        t0 = t1-prior
        iprevday, dum = tb.tOverlap(utc, [t0, t1])
        if iprevday:
            Kpmax = max(Kp[iprevday])
            Lpp[i] = calcLpp(Kpmax, A, B, power=power)
        else:
            Lpp[i] = np.nan

    return Lpp


def getMagnetopause(ticks, LTs=None, dbase='QDhourly'):
    '''Calculates the Shue et al. (1997) position in equatorial plane

    Shue et al. (1997), A new functional form to study the solar wind 
    control of the magnetopause size and shape, J. Geophys. Res., 102(A5), 
    9497–9511, doi:10.1029/97JA00196.

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        TickTock object of desired times (will be interpolated from hourly OMNI data)
        OR dictionary of form {'P': [], 'Bz': []}
        Where P is SW ram pressure [nPa] and Bz is IMF Bz (GSM) [nT]
    LTs : array-like
        Array-like of local times for evaluating the magnetopause model. Default is
        6 LT to 18 LT in steps of 20 minutes.

    Returns
    =======
    out : array
        NxMx2 array of magnetopause positions [Re]
        N is number of timesteps, M is number of local times. The 2 positions
        are the X_GSE and Y_GSE positions of the magnetopause

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.Ticktock(['2002-01-01T12:00:00','2002-01-04T00:00:00'])
    >>> localtimes = [13,12,11]
    >>> emp.getMagnetopause(ticks, LTs=localtimes)
    array([[[ 10.27674331,  -2.75364507],
        [ 10.52909163,   0.        ],
        [ 10.27674331,   2.75364507]],
       [[ 10.91791834,  -2.9254474 ],
        [ 11.18712131,   0.        ],
        [ 10.91791834,   2.9254474 ]]])
    >>> emp.getMPstandoff(ticks) #should give same result as getMagnetopause for 12LT
    array([ 10.52909163,  11.18712131])
    
    To plot the magnetopause:
    >>> import numpy as np
    >>> import spacepy.plot as splot
    >>> import matplotlib.pyplot as plt
    >>> localtimes = np.arange(5, 19.1, 0.5)
    >>> mp_pos = emp.getMagnetopause(ticks, localtimes)
    >>> plt.plot(mp_pos[0,:,0], mp_pos[0,:,1])
    >>> ax1 = plt.gca()
    >>> ax1.set_xlim([-5,20])
    >>> ax1.set_xlabel('X$_{GSE}$ [R$_E$]')
    >>> ax1.set_ylabel('Y$_{GSE}$ [R$_E$]')
    >>> splot.dual_half_circle(ax=ax1)
    >>> ax1.axes.set_aspect('equal')

    '''

    alpha = []
    r0 = getMPstandoff(ticks, dbase=dbase, alpha=alpha)
    alpha = np.asarray(alpha)

    if LTs is None:
        LTs = np.arange(6, 18.1, 1/3.0)

    r = np.zeros([len(ticks),len(LTs)])
    out = np.zeros([len(ticks),len(LTs), 2])

    angs = (12.0-np.asanyarray(LTs))*15.0

    costheta = np.cos(np.deg2rad(angs))
    sintheta = np.sin(np.deg2rad(angs))
    for idx, ct in enumerate(costheta):
        r[:,idx] = r0 * (2.0/(1+ct))**alpha
    out[:,:,0] = r*costheta #Xgse
    out[:,:,1] = r*sintheta #Ygse

    return out


def getMPstandoff(ticks, dbase='QDhourly', alpha=[]):
    """Calculates the Shue et al. (1997) subsolar magnetopause radius

    Shue et al. (1997), A new functional form to study the solar wind 
    control of the magnetopause size and shape, J. Geophys. Res., 102(A5), 
    9497–9511, doi:10.1029/97JA00196.

    Parameters
    ==========
    ticks : spacepy.time.Ticktock
        TickTock object of desired times (will be interpolated from hourly OMNI data)
        OR dictionary of form {'P': [], 'Bz': []}
        Where P is SW ram pressure [nPa] and Bz is IMF Bz (GSM) [nT]
    alpha : list
        Used as an optional return value to obtain the flaring angles. To use, assign
        an empty list and pass to this function through the keyword argument. The list 
        will be modified in place, adding the flaring angles for each time step.

    Returns
    =======
    out : float
        Magnetopause (sub-solar point) standoff distance [Re]

    Examples
    ========
    >>> import spacepy.time as spt
    >>> import spacepy.empiricals as emp
    >>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)
    >>> emp.getMPstandoff(ticks)
    array([ 10.57319537,  10.91327764,  10.75086873,  10.77577207,
         9.78180261,  11.0374474 ,  11.4065    ,  11.27555451,
        11.47988573,  11.8202582 ,  11.23834814])
    >>> data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
    >>> emp.getMPstandoff(data)
    array([ 9.96096838,  8.96790412])
    """
    if type(ticks) == spt.Ticktock:
        omni = om.get_omni(ticks, dbase=dbase)
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
        flarang = (0.58 - 0.01*Bz)*(1.0 + 0.01*P)
        alpha.extend(flarang)

        return r0
    except TypeError:
        raise TypeError("Please check for valid input types")

def getDststar(ticks, model='OBrien', dbase='QDhourly'):
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
        omni = om.get_omni(ticks, dbase=dbase)
        P, Dst = omni['Pdyn'], omni['Dst']
    elif isinstance(ticks, dict):
        P, Dst = ticks['Pdyn'], ticks['Dst']
        if isinstance(P, list):
            P, Dst = np.array(P), np.array(Dst)

    #get Dst*
    Dststar = Dst - b*P**0.5 + c

    return Dststar


def getExpectedSWTemp(velo, model='XB15', units='K'):
    '''Return the expected solar wind temperature based on the bulk velocity
    
    The formulations used by this function are those given by,
    L87  -- Lopez, R.E., J. Geophys. Res., 92, 11189-11194, 1987
    BS06 -- Borovsky, J.E. and J.T. Steinberg, Geophysical Monograph Series 167, 59-76, 2006
    XB15 -- Xu, F. and J.E. Borovsky, J. Geophys. Res., 120, 70-100, 2015

    Parameters
    ==========
    velo : array-like
        Array like of solar wind bulk velocity values [km/s]

    model : str [optional]
        Name of model to use. Valid choices are L87, BS06 and XB15. Default is XB15

    units : str [optional]
        Units for output temperature, options are eV or K. Default is Kelvin [K]

    Returns
    =======
    Texp : array-like
        The expected solar wind temperature given the bulk velocity [K] or [eV]

    '''
    v = np.asanyarray(velo)
    def bs06(v):
        '''Borovsky and Steinberg 2006, Geophysical Monograph - median values'''
        Texp = np.empty(len(v))
        Texp.fill(np.nan)
        Texp[v<372]  = 1.28e-8*v[v<372]**3.324
        Texp[v>=372] = 0.0572*v[v>=372] - 16.79
        return Texp

    def xb15(v):
        '''Xu and Borovsky 2015, JGR'''
        Texp = (v/258.0)**3.113
        return Texp

    def l87(v):
        '''Lopez 1987, JGR - Tables 1&2 [T(V)] from IMP-8'''
        Texp = np.empty(len(v))
        Texp.fill(np.nan)
        Texp[v<500]  = (0.031*v[v<500] - 5.1)**2
        Texp[v>=500] = (0.02*v[v>=500] + 0.5)**2
        return Texp*1e3/1.16045221e4 #return in eV

    formulae = {'BS06': bs06, 'XB15': xb15, 'L87': l87}

    try:
        mod = model.upper()
        Texp = formulae[mod](v)
    except KeyError:
        raise KeyError('Invalid model specified for SW temperature')

    if units.lower()=='k':
        return Texp*1.16045221e4
    elif units.lower()=='ev':
        return Texp
    else:
        raise ValueError('Invalid units specified for SW temperature, must be "K" or "eV"')
    

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

    >>> from spacepy.empiricals import vampolaPA
    >>> vampolaPA(3000, alpha=[45, 90])
    (array([ 179.04931098,  358.09862196]), [45, 90])
    >>> data, pas = vampolaPA([3000, 6000], alpha=[45, 90])
    >>> pas
    [45, 90]
    >>> data
    array([[ 179.04931098,  358.09862196],
       [ 358.09862196,  716.19724391]])

    Notes
    =====
    Directional number flux integrated over pitch angle from 0 to 90 degrees
    is a factor of 4*pi lower than omnidirectional number flux.

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
        sinfunc_o = partial(sinfunc, order=tmporder+1)
        normfac[idx] = omniflux[idx]/(2*np.pi*integ.quad(sinfunc_o, 0, np.pi)[0])
    #now make the differential number flux
    dnflux = np.zeros((len(kwargs['alpha']), len(omniflux)))
    for i, a_val in enumerate(np.deg2rad(kwargs['alpha'])):
        dnflux[i] = normfac * sinfunc(a_val)

    return dnflux.squeeze(), kwargs['alpha']


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

    Examples
    ========
    Apply Vampola pitch angle model at L=[4, 6.6]

    >>> from spacepy.empiricals import vampolaPA, getVampolaOrder
    >>> order = getVampolaOrder([4,6.6])
    >>> order
    array([ 3.095 ,  1.6402])
    >>> vampolaPA([3000, 3000], alpha=[45, 90], order=order)
    (array([[ 140.08798878,  192.33572182],
        [ 409.49143136,  339.57417256]]), [45, 90])
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


def omniFromDirectionalFlux(fluxarr, alphas, norm=True):
    '''
    Calculate omnidirectional flux [(s cm^2 kev)^-1] from directional flux [(s sr cm^2 keV)^-1] array

    J = 2.pi integ(j sin(a) da)
    If kwarg 'norm' is True (default), the omni flux is normalized by 4.pi to make it per steradian,
    in line with the PRBEM guidelines

    Parameters
    ==========
    fluxarr : arraylike
        Array of directional flux values
    alphas : arraylike
        Array of pitch angles corresponding to fluxarr
        
    Returns
    =======
    omniflux : float
        Omnidirectional flux value

    Examples
    ========
    Roundtrip from omni flux, to directional flux (Vampola model),
    integrate to get back to omni flux.

    >>> from spacepy.empiricals import vampolaPA, omniFromDirectionalFlux
    >>> dir_flux, pa  = vampolaPA(3000, alpha=range(0,181,2), order=4)
    >>> dir_flux[:10], pa[:10]
    (array([  0.00000000e+00,   6.64032473e-04,   1.05986545e-02,
          5.34380898e-02,   1.67932162e-01,   4.06999226e-01,
          8.36427502e-01,   1.53325140e+00,   2.58383611e+00,
          4.08170975e+00]), [0, 2, 4, 6, 8, 10, 12, 14, 16, 18])
    >>> omniFromDirectionalFlux(dir_flux, pa, norm=False)
    3000.0000008112293
    
    Calculate "spin-averaged" flux, giving answer per steradian

    >>> omniFromDirectionalFlux(dir_flux, pa, norm=True)
    238.73241470239859
    '''
    try:
        flen = len(fluxarr)
    except TypeError:
        fluxarr = np.array([fluxarr])
        flen = 1
    finally:
        if flen==1:
            #assume isotropy
            return fluxarr
    if norm:
       fac = 1
       denomina = integ.quad(np.sin, 0, np.pi)[0]
    else:
       fac = 2*np.pi
       denomina = 1
    alphrad = np.deg2rad(alphas)
    numera = integ.simps(fluxarr*np.sin(alphrad), alphrad)
    omniflux = fac*numera/denomina
    return omniflux


def getSolarRotation(ticks, rtype='carrington', fp=False, reverse=False):
    '''Calculates solar rotation number (Carrington or Bartels) for a given date/time

    Parameters
    ==========
    ticks : spacepy.time.Ticktock or datetime.datetime
        
    Returns
    =======
    rnumber : integer or array
        Carrington (or Bartels) rotation number

    '''
    def total_seconds(dobj):
        return dobj.days*24*3600 + dobj.seconds + dobj.microseconds/1e6
    if rtype.lower() == 'carrington':
        start_date = datetime.datetime(1853,11,9,21,38,44,160000)
        #length = datetime.timedelta(days=27, minutes=396, seconds=25, microseconds=919999)
        length = datetime.timedelta(days=27.2753)
    elif rtype.lower() == 'bartels':
        start_date = datetime.datetime(1832,2,8)
        start_JD = spt.Ticktock(start_date).JD
        length = datetime.timedelta(days=27)
    else:
        raise ValueError('Solar rotation type {0} not recognized: Must be either "carrington" or "Bartels"'.format(rtype))
    if not reverse:
        try:
            nels = len(ticks)
        except TypeError:
            try:
                rotation = total_seconds(ticks-start_date)/total_seconds(length)
                rotation += 1
                if not fp:
                    rotation = int(rotation)
                return rotation
            except:
                raise RuntimeError('Unidentified problem with input time {0} in getSolarRotation'.format(ticks))
        if isinstance(ticks, spt.Ticktock):
            rotation = [total_seconds(tt-start_date)/total_seconds(length) for tt in ticks.UTC]
            rotation = np.array(rotation) + 1
        else:
            try:
                rotation = [total_seconds(tt-start_date)/total_seconds(length) for tt in ticks]
                rotation = np.array(rotation) + 1
            except:
                raise RuntimeError('Unidentified problem with input time {0} in getSolarRotation'.format(ticks))
        if not fp:
            rotation = rotation.astype(int)
        return rotation
    else:
        #for now just assume single input, non-iterable
        elapsed = length*(ticks-1)
        date = elapsed + start_date
        return date

def getSolarProtonSpectra(norm=3.20e7, gamma=-0.96, E0=15.0, Emin=.1, Emax=600, nsteps=100):
    '''Returns a SpaceData with energy and fluence spectra of solar particle events

    The formulation follows that of:
    Ellison and Ramaty ApJ 298: 400-408, 1985
    dJ/dE = K^{-\\gamma}exp(-E/E0)
    
    and the defualt values are the 10/16/2003 SEP event of:
    Mewaldt, R. A., et al. (2005), J. Geophys. Res., 110, A09S18, doi:10.1029/2005JA011038.

    Other Parameters
    ================
    norm : float
        Normilization factor for the intensity of the SEP event
    gamma : float
        Power law index
    E0 : float
        Expoential scaling factor
    Emin : float
        Minimum energy for fit
    Emax : float
        Maximum energy for fit
    nsteps : int
        The number of log spaced energy steps to return
        
    Returns
    =======
    data : dm.SpaceData
        SpaceData with the energy and fluence values
    '''
    E = tb.logspace(Emin, Emax, nsteps)
    fluence = norm*E**(gamma)*np.exp(-E/E0)
    ans = dm.SpaceData()
    ans['Energy'] = dm.dmarray(E)
    ans['Energy'].attrs = {'UNITS':'MeV',
                           'DESCRIPTION':'Particle energy per nucleon'}
    ans['Fluence'] = dm.dmarray(fluence)
    ans['Fluence'].attrs = {'UNITS' : 'cm^{-2} sr^{-1} (MeV/nuc)^{-1}',
                            'DESCRIPTION':'Fluence spectra fir to the model'}
    return ans


      
ShueMP = getMPstandoff
get_plasma_pause = getPlasmaPause
get_Lmax = getLmax
