#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
module wrapper for onera_desp_lib
Reference: D. Boscher1, S. Bourdarie1, P. O'Brien2, 
T. Guild2,(1 ONERA-DESP, Toulouse France; 2 Aerospace Corporation, 
Washington DC, USA), ONERA-DESP library V4.2, Toulouse-France, 2004-2008
"""
from spacepy import help
__version__ = "$Revision: 1.1 $, $Date: 2010/05/20 17:19:44 $"
__author__ = 'Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)'

# -----------------------------------------------
def get_Bfield(ticktock, spaco, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
	"""
	call get_bfield in onera_desp lib and return a dictionary with the B-field vector and 
	strenght.
	
	Input:
	======
		- ticktock (Ticktock class) : containing time information
		- spaco (Spaco class) : containing spatial information
		- extMag (string) : optional; will choose the external magnetic field model 
							possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
							'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
							'T05', 'ALEX']
		- options (optional list or array of integers length=5) : explained in Lstar
		- omni values as dictionary (optional) : if not provided, will use lookup table
		- (see Lstar documentation for further explanation)
	
	Returns:
	========
		- results (dictionary) : containing keys: Bvec, and Blocal

	Example:
	========
	>>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
	>>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
	>>> op.get_Bfield(t,y)
	{'Blocal': array([  945.99989101,  3381.71633205]),
		'Bvec': array([[  8.05001055e-01,  -1.54645026e+02,   9.33273841e+02],
		[  3.36352963e+02,  -5.33658140e+02,   3.32236076e+03]])}


	See Also:
	=========
	Lstar, find_Bmirror, find_magequator

	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

	Version:
	========
	V1: 05-Mar-2010 (JK)
	
	"""
	
	import numpy as n
	import spacepy.onerapylib as oplib
	import spacepy.spacetime as st
	import spacepy.toolbox as tb
	
	# prepare input values for onera_desp call
	d = prep_onera(ticktock, spaco, alpha=[], extMag=extMag, options=options, omnivals=omnivals)
	nTAI = len(ticktock)
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

	results = {}
	results['Blocal'] = n.zeros(nTAI)
	results['Bvec'] = n.zeros((nTAI,3))
	for i in n.arange(nTAI):
		BxyzGEO, Blocal = oplib.get_field1(kext,options,sysaxes,iyearsat[i],idoysat[i],secs[i], \
			xin1[i],xin2[i],xin3[i], magin[:,i])

		# take out all the odd 'bad values' and turn them into NaN
		if tb.feq(Blocal,badval): Blocal = n.NaN
		BxyzGEO[n.where( tb.feq(BxyzGEO, badval)) ] = n.NaN
	
		results['Blocal'][i] = Blocal
		results['Bvec'][i,:] = BxyzGEO
	
	return results

# -----------------------------------------------
def find_Bmirror(ticktock, spaco, alpha, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
	"""
	call find_mirror_point from onera_desp library and return a dictionary with values for 
	Blocal, Bmirr and the GEO (cartesian) coordinates of the mirror point
	
	Input:
	======
		- ticktock (Ticktock class) : containing time information
		- spaco (Spaco class) : containing spatial information
		- alpha (list or ndarray) : containing the pitch angles
		- extMag (string) : optional; will choose the external magnetic field model 
							possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
							'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
							'T05', 'ALEX']
		- options (optional list or array of integers length=5) : explained in Lstar
		- omni values as dictionary (optional) : if not provided, will use lookup table 
		- (see Lstar documentation for further explanation)

	Returns:
	========
		- results (dictionary) : containing keys: Blocal, Bmirr, GEOcar

	Example:
	========
	>>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
	>>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
	>>> op.find_Bmirror(t,y,[90,80,60,10])
	{'Blocal': array([ 0.,  0.]),
	 'Bmirr': array([ 0.,  0.]),
	 'GEOcar': Coords( [[ NaN  NaN  NaN]
	 [ NaN  NaN  NaN]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']}

	See Also:
	=========
	Lstar, get_Bfield, find_magequator

	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

	Version:
	========
	V1: 05-Mar-2010 (JK)
	
	"""
	
	import numpy as n
	import spacepy.onerapylib as oplib
	import spacepy.coordinates as c
	import spacepy.toolbox as tb
	
	# prepare input values for onera_desp call
	d = prep_onera(ticktock, spaco, alpha, extMag=extMag, options=options, omnivals=omnivals)
	nTAI = len(ticktock)
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
	degalpha = d['degalpha']
	nalp_max = d['nalp_max']
	nalpha = len(alpha)

	results = {}
	results['Blocal'] = n.zeros(nTAI)
	results['Bmirr'] = n.zeros(nTAI)
	results['GEOcar'] = ['']*nTAI
	for i in n.arange(nTAI):
		blocal, bmirr, GEOcoord = oplib.find_mirror_point1(kext,options,sysaxes, \
			iyearsat[i],idoysat[i],secs[i], xin1[i],xin2[i],xin3[i], alpha, magin[:,i])

		# take out all the odd 'bad values' and turn them into NaN
		if tb.feq(blocal,badval): blocal = n.NaN
		if tb.feq(bmirr,badval) : bmirr  = n.NaN
		GEOcoord[n.where( tb.feq(GEOcoord,badval)) ] = n.NaN

		results['Blocal'][i] = blocal
		results['Bmirr'][i] = bmirr	
		results['GEOcar'][i] = GEOcoord
	
	results['GEOcar'] = c.Coords(results['GEOcar'], 'GEO', 'car')
	
	return results

# -----------------------------------------------
def find_magequator(ticktock, spaco, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
	"""
	call find_magequator from onera_desp library and return a dictionary with values for
	Bmin and the GEO (cartesian) coordinates of the magnetic equator
	
	Input:
	======
		- ticktock (Ticktock class) : containing time information
		- spaco (Spaco class) : containing spatial information
		- extMag (string) : optional; will choose the external magnetic field model 
							possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
							'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
							'T05', 'ALEX']
		- options (optional list or array of integers length=5) : explained in Lstar
		- omni values as dictionary (optional) : if not provided, will use lookup table 
		- (see Lstar documentation for further explanation)

	Returns:
	========
		- results (dictionary) : containing keys: Bmin, Coords instance with GEO coordinates of 
			the magnetic equator

	Example:
	========
	>>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
	>>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
	>>> op.find_magequator(t,y)
	{'Bmin': array([  945.63652413,  3373.64496167]),
	 'GEOcar': Coords( [[ 2.99938371  0.00534151 -0.03213603]
	 [ 2.00298822 -0.0073077   0.04584859]] ), dtype=GEO,car, units=['Re', 'Re', 'Re']}

	See Also:
	=========
	Lstar, get_Bfield, find_Bmirr

	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

	Version:
	========
	V1: 05-Mar-2010 (JK)
	
	"""
	import numpy as n
	import spacepy.onerapylib as oplib
	import spacepy.toolbox as tb
	import spacepy.spacetime as st
	
	# prepare input values for onera_desp call
	d = prep_onera(ticktock, spaco, alpha=[], extMag=extMag, options=options, omnivals=omnivals)
	nTAI = len(ticktock)
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
	results['Bmin'] = n.zeros(nTAI)
	results['GEOcar'] = ['']*nTAI
	for i in n.arange(nTAI):
		bmin, GEOcoord = oplib.find_magequator1(kext,options,sysaxes,\
			iyearsat[i],idoysat[i],secs[i], xin1[i],xin2[i],xin3[i],magin[:,i])
	
		# take out all the odd 'bad values' and turn them into NaN
		if tb.feq(bmin,badval): bmin = n.NaN
		GEOcoord[n.where( tb.feq(GEOcoord, badval)) ] = n.NaN
	
		results['Bmin'][i] = bmin
		results['GEOcar'][i] = GEOcoord
	
	results['GEOcar'] = c.Coords(results['GEOcar'], 'GEO', 'car')
	
	return results


# -----------------------------------------------
def coord_trans( spaco, returntype, returncarsph ):
	"""
	thin layer to call coor_trans1 from onera_desp lib
	this will convert between systems GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
	
	Input:
	======
		- spaco (Coords instance) : containing coordinate information, can contain n points
		- returntype (str) : describing system as GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
		- returncarsph (str) : cartesian or spherical units 'car', 'sph'
		
	Returns:
	========
		- xout (ndarray) : values after transformation in (n,3) dimensions
		
	Example:
	========
	>>> coord = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
	>>> coord.ticktock = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
	>>> coord_trans(coord, 'GSM', 'car')
	array([[ 2.8639301 , -0.01848784,  0.89306361],
	[ 1.9124434 ,  0.07209424,  0.58082929]])
	
	See Also:
	=========
	sph2car, car2sph
	
	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov
	
	Version:
	========
	V1: 05-Mar-2010 (JK)
	
	"""
	import numpy as n
	import onerapylib
	import onerapy as op
	
	sysaxesin = get_sysaxes( spaco.dtype, spaco.carsph )
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
	
	xout = n.zeros(n.shape(spaco.data))
	for i in n.arange(len(spaco)):
		iyear = spaco.ticktock.UTC[i].year
		idoy = spaco.ticktock.DOY[i]
		secs = spaco.ticktock.UTC[i].hour*3600. + spaco.ticktock.UTC[i].minute*60. + \
			spaco.ticktock.UTC[i].second
			
		xout[i,:] = onerapylib.coord_trans1(sysaxesin, sysaxesout, \
			iyear, idoy, secs, spaco.data[i])
		
	# add  sph to car or v/v convertion if initial sysaxesout was None
	if  aflag == True:
		if returncarsph == 'sph':
			xout = op.car2sph(xout)
		else: # 'car' needs to be returned
			xout = op.sph2car(xout)

	return xout

# -----------------------------------------------
def car2sph(CARin):
	"""
	coordinate transformation from cartesian to spherical
	
	Input:
	======
		- CARin (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
			units of [Re, Re, Re] = [x,y,z]
		
	Returns:
	========
		- results (ndarray) : values after conversion to spherical coordinates in
			radius, latitude, longitude in units of [Re, deg, deg]
		
	Example:
	========
	>>> sph2car([1,45,0])
	array([ 0.70710678,  0.        ,  0.70710678])

	See Also:
	=========
	sph2car
	
	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov
	
	Version:
	========
	V1: 05-Mar-2010 (JK)
	
	"""	
	
	import numpy as n

	if isinstance(CARin[0], (float, int)):
		CAR = n.array([CARin])
	else:
		CAR = n.array(CARin)
	
	res = n.zeros(n.shape(CAR))
	for i in n.arange(len(CAR)):
		x, y, z = CAR[i,0], CAR[i,1], CAR[i,2]
		r = n.sqrt(x*x+y*y+z*z)
		sq = n.sqrt(x*x+y*y)
		if (x == 0) & (y == 0): # on the poles
			longi = 0.
			if z < 0:
				lati = -90.
			else:
				lati = 90.0
		else:
			longi = n.arctan2(y,x)*180./n.pi
			lati = 90. - n.arctan2(sq, z)*180./n.pi
		res[i,:] = [r, lati, longi]
		
	if isinstance(CARin[0], (float, int)):
		return res[0]
	else:
		return res

# -----------------------------------------------
def sph2car(SPHin):
	"""
	coordinate transformation from spherical to cartesian
	
	Input:
	======
		- SPHin (list or ndarray) : coordinate points in (n,3) shape with n coordinate points in
			units of [Re, deg, deg] = [r, latitude, longitude]
		
	Returns:
	========
		- results (ndarray) : values after conversion to cartesian coordinates x,y,z
		
	Example:
	========
	>>> sph2car([1,45,45])
	array([ 0.5       ,  0.5       ,  0.70710678])

	See Also:
	=========
	car2sph
	
	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov
	
	Version:
	========
	V1: 05-Mar-2010 (JK)
	"""
	
	import numpy as n
	
	if isinstance(SPHin[0], (float, int)):
		SPH = n.array([SPHin])
	else:
		SPH = n.array(SPHin)
	
	res = n.zeros(n.shape(SPH))
	for i in n.arange(len(SPH)):
		r,lati,longi = SPH[i,0], SPH[i,1], SPH[i,2]
		colat = n.pi/2. - lati*n.pi/180.
		x = r*n.sin(colat)*n.cos(longi*n.pi/180.)
		y = r*n.sin(colat)*n.sin(longi*n.pi/180.)
		z = r*n.cos(colat)
		res[i,:] = [x, y, z]
	
	
	if isinstance(SPHin[0], (float, int)):
		return res[0]
	else:
		return res
	
# -----------------------------------------------    
def get_sysaxes(dtype, carsph):
    """
    will return the sysaxes according to the onera_desp library

    Input:
    ======
        - dtype (str) : coordinate system, possible values: GDZ, GEO, GSM, GSE, SM, 
                GEI, MAG, SPH, RLL
        - carsph (str) : cartesian or spherical, possible values: 'sph', 'car'
        
    Returns:
    ========
        - sysaxes (int) : value after oner_desp library from 0-8 (or None if not available)
        
    Example:
    ========
    >>> get_sysaxes('GSM', 'car')
    2

    See Also:
    =========
    get_dtype

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov

    Version:
    ========
    V1: 05-Mar-2010 (JK)
    V2: 10-May-2010 (SM)
    """

    typedict = {'GDZ': {'sph': 0, 'car': None},
        'GEO': {'sph': 1, 'car': None}, 'GSM': {'sph': None, 'car': 2},
        'GSE': {'sph': None, 'car': 3}, 'SM': {'sph': None, 'car': 4},
        'GEI': {'sph': None, 'car': 5}, 'MAG': {'sph': None, 'car': 6},
        'SPH': {'sph': 7, 'car': None}, 'RLL': {'sph': 8, 'car': None}}
    
    #typedict = {'GDZ': {'sph': 0, 'car': 10},
        #'GEO': {'sph': 1, 'car': 11}, 'GSM': {'sph': 22, 'car': 2},
        #'GSE': {'sph': 23, 'car': 3}, 'SM': {'sph': 24, 'car': 4},
        #'GEI': {'sph': 25, 'car': 5}, 'MAG': {'sph': 26, 'car': 6},
        #'SPH': {'sph': 7, 'car': 17}, 'RLL': {'sph': 8, 'car': 18}}
        
    
    sysaxes = typedict[dtype][carsph]

    return sysaxes
    
# -----------------------------------------------    
def get_dtype(sysaxes):
    """
    will return the coordinate system type as string

    Input:
    ======
        - sysaxes (int) : number according to the onera_desp_lib, possible values: 0-8
        
    Returns:
    ========
        - dtype (str) : coordinate system GDZ, GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL
        - carsph (str) : cartesian or spherical 'car', 'sph'

    Example:
    ========
    >>> get_dtype(3)
    ('GSE', 'car')

    See Also:
    =========
    get_sysaxes

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 05-Mar-2010 (JK)
    """

    oneralist_sph = ['GDZ', 'SPH', 'RLL']
    oneralist_sphidx = [0, 7, 8]
    oneralist_car = ['GEO', 'GSM', 'GSE', 'SM', 'GEI', 'MAG']
    oneralist_caridx = [1,2,3,4,5,6]
    if sysaxes in oneralist_caridx:
        dtype = oneralist_car[oneralist_caridx.index(sysaxes)]
        carsph = 'car'
    elif sysaxes in onerlist_sphidx:
        dtype = oneralist_sph[oneralist_sphidx.index(sysaxes)]
        carsph = 'sph'
    # or None after onera library
    else:
        print "sysaxes="+str(sysaxes)+" not supported"
        dtype = None
        carsph = None

    return dtype, carsph

# -----------------------------------------------
def _get_Lstar(ticktock, spaco, alpha=[], extMag='T01STORM', options=[1,0,0,0,0], omnivals=None): 
	"""
	This will call make_lstar1 or make_lstar_shell_splitting_1 from the onera library
	and will lookup omni values for given time if not provided (optional). If pitch angles
	are provided, drift shell splitting will be calculated and "Bmirr" will be returned. If they
	are not provided, then no drift shell splitting is calculated and "Blocal" is returned.
	
	Input:
	======
		- ticktock (Ticktock class) : containing time information
		- spaco (Spaco class) : containing spatial information
		- alpha (list or ndarray) : optional pitch angles in degrees; if provided will 
			calculate shell splitting; max 25 values
		- extMag (string) : optional; will choose the external magnetic field model 
							possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
							'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
							'T05', 'ALEX']
		- options (optional list or array of integers length=5) : explained below
		- omni values as dictionary (optional) : if not provided, will use lookup table 

	Returns:
	========
		- results (dictionary) : containing keys: Lm, Lstar, Bmin, Blocal (or Bmirr), Xj, MLT 
			if pitch angles provided in "alpha" then drift shells are calculated and "Bmirr" 
			is returned if not provided, then "Blocal" at spacecraft is returned.

	Example:
	========
	>>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
	>>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
	>>> spacepy.onerapy.Lstar(t,y)
	{'Blocal': array([ 1020.40493286,  3446.08845227]),
	 'Bmin': array([ 1019.98404311,  3437.63865243]),
	 'Lm': array([ 3.08948304,  2.06022102]),
	 'Lstar': array([ 2.97684043,  1.97868577]),
	 'MLT': array([ 23.5728333 ,  23.57287944]),
	 'Xj': array([ 0.00112884,  0.00286955])}


	External Magnetic Field:
	========================
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

	OMNI Values:
	============
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
	
	Options:
	========
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
		   
	Author:
	=======
	Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

	Version:
	========
	V1: 03-Feb-2010 (JK)
	V1.1: 22-Feb-2010: fixed Blocal and Bmirr bug (JK)
	V1.2: 05-Mar-2010: ticktock, spaco support (JK)
	"""
	import spacepy
	import spacepy.toolbox as tb
	import spacepy.onerapylib as oplib
	import numpy as n

	nTAI = len(ticktock)
	nalpha = len(alpha)
	d = prep_onera(ticktock, spaco, alpha, extMag, options, omnivals)
	

	if nalpha == 0: # no drift shell splitting
		lm, lstar, blocal, bmin, xj, mlt = oplib.make_lstar1(nTAI, d['kext'], d['options'], d['sysaxes'],\
					d['iyearsat'], d['idoysat'], d['utsat'], d['xin1'], d['xin2'], d['xin3'], d['magin'])

	elif nalpha > 0 and nalpha <= d['nalp_max']: # with drift shell splitting
		lm, lstar, bmirr, bmin, xj, mlt = oplib.make_lstar_shell_splitting1(nTAI, nalpha, \
			d['kext'], d['options'], d['sysaxes'], d['iyearsat'], d['idoysat'], d['utsat'], \
			d['xin1'], d['xin2'], d['xin3'], d['degalpha'], d['magin'])

	else:
		print 'ERROR: too many pitch angles requested; 25 is maximum'
	
	# take out all the odd 'bad values' and turn them into NaN
	lm[n.where( tb.feq(lm,d['badval'])) ] = n.NaN
	lstar[n.where( tb.feq(lstar,d['badval'])) ] = n.NaN
	bmin[n.where( tb.feq(bmin,d['badval'])) ] = n.NaN
	xj[n.where( tb.feq(xj,d['badval'])) ] = n.NaN
	mlt[n.where( tb.feq(mlt,d['badval'])) ] = n.NaN
	
	results = {}
	if nalpha == 0:
		results['Lm'] = lm[0:nTAI]
		results['Lstar'] = lstar[0:nTAI]
		blocal[n.where( tb.feq(blocal,d['badval'])) ] = n.NaN
		results['Blocal'] = blocal[0:nTAI]
		results['Bmin'] = bmin[0:nTAI]
		results['Xj'] = xj[0:nTAI]
		results['MLT'] = mlt[0:nTAI]
	else:		
		results['Lm'] = lm[0:nTAI, 0:nalpha]
		results['Lstar'] = lstar[0:nTAI, 0:nalpha]
		bmirr[n.where( tb.feq(bmirr, d['badval'])) ] = n.NaN
		results['Bmirr'] = bmirr[0:nTAI, 0:nalpha]
		results['Bmin'] = bmin[0:nTAI]
		results['Xj'] = xj[0:nTAI, 0:nalpha]
		results['MLT'] = mlt[0:nTAI]
		
	return results
	
# -----------------------------------------------
def get_Lstar(ticktock, spaco, alpha, extMag='T01STORM', options=[1,0,0,0,0], omnivals=None):
    """
    This will call make_lstar1 or make_lstar_shell_splitting_1 from the onera library
    and will lookup omni values for given time if not provided (optional). If pitch angles
    are provided, drift shell splitting will be calculated and "Bmirr" will be returned. If they
    are not provided, then no drift shell splitting is calculated and "Blocal" is returned.

    Input:
    ======
        - ticktock (Ticktock class) : containing time information
        - spaco (Spaco class) : containing spatial information
        - alpha (list or ndarray) : optional pitch angles in degrees; if provided will 
            calculate shell splitting; max 25 values
        - extMag (string) : optional; will choose the external magnetic field model 
                            possible values ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 
                            'OPQUIET', 'OPDYN', 'T96', 'OSTA', 'T01QUIET', 'T01STORM', 
                            'T05', 'ALEX']
        - options (optional list or array of integers length=5) : explained below
        - omni values as dictionary (optional) : if not provided, will use lookup table 

    Returns:
    ========
        - results (dictionary) : containing keys: Lm, Lstar, Bmin, Blocal (or Bmirr), Xj, MLT 
            if pitch angles provided in "alpha" then drift shells are calculated and "Bmirr" 
            is returned if not provided, then "Blocal" at spacecraft is returned.

    Example:
    ========
    >>> t = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
    >>> y = Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
    >>> spacepy.onerapy.Lstar(t,y)
    {'Blocal': array([ 1020.40493286,  3446.08845227]),
        'Bmin': array([ 1019.98404311,  3437.63865243]),
        'Lm': array([ 3.08948304,  2.06022102]),
        'Lstar': array([ 2.97684043,  1.97868577]),
        'MLT': array([ 23.5728333 ,  23.57287944]),
        'Xj': array([ 0.00112884,  0.00286955])}


    External Magnetic Field:
    ========================
        - 0    : no external field
        - MEAD  : Mead & Fairfield [1975] (uses 0<=Kp<=9 - Valid for rGEO<=17. Re)
        - T87SHORT: Tsyganenko short [1987] (uses 0<=Kp<=9 - Valid for rGEO<=30. Re)
        - T87LONG : Tsyganenko long [1987] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - T89    : Tsyganenko [1989] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
        - OPQUIET : Olson & Pfitzer quiet [1977] (default - Valid for rGEO<=15. Re)
        - OPDYN   : Olson & Pfitzer dynamic [1988] (uses 5.<=dens<=50., 300.<=velo<=500., 
            -100.<=Dst<=20. - Valid for rGEO<=60. Re)
        - T96    : Tsyganenko [1996] (uses -100.<=Dst (nT)<=20., 0.5<=Pdyn (nPa)<10., 
            |ByIMF| (nT)<1=0., |BzIMF| (nT)<=10. - Valid for rGEO<=40. Re)
        - OSTA  : Ostapenko & Maltsev [1997] (uses dst,Pdyn,BzIMF, Kp)
            T01QUIET: Tsyganenko [2002a,b] (uses -50.<Dst (nT)<20., 0.5<Pdyn (nPa)<=5., 
            |ByIMF| (nT)<=5., |BzIMF| (nT)<=5., 0.<=G1<=10., 0.<=G2<=10. - Valid for xGSM>=-15. Re)
        - T01STORM: Tsyganenko, Singer & Kasper [2003] storm  (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 - 
            there is no upper or lower limit for those inputs - Valid for xGSM>=-15. Re)
        - T05    : Tsyganenko & Sitnov [2005] storm  (uses Dst, Pdyn, ByIMF, BzIMF, 
            W1, W2, W3, W4, W5, W6 - no upper or lower limit for inputs - Valid for xGSM>=-15. Re)

    OMNI Values:
    ============
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

    Options:
    ========
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
            
    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 03-Feb-2010 (JK)
    V1.1: 22-Feb-2010: fixed Blocal and Bmirr bug (JK)
    V1.2: 05-Mar-2010: ticktock, spaco support (JK)
    """

    import spacepy
    import numpy as n

    if isinstance(alpha, (float,int)):
        alpha = [alpha]
    assert len(alpha) is 1, 'len(alpha) needs to be 1'	

    ncpus = spacepy.NCPUS
    ncalc = len(ticktock)
    nalpha = len(alpha)

    if ncpus > 1 and ncalc >= ncpus*2:
        import pp
        server = pp.Server(ncpus)
        ncalc = len(ticktock)
        jobs = []

        for ijob in range(ncpus):
            ppidx = range(ijob, ncalc, ncpus)
            jobs.append(server.submit(_get_Lstar, (ticktock[ppidx], spaco[ppidx], alpha, \
                extMag, options, omnivals), depfuncs=(prep_onera,), modules=() ))
            
        # retrieve results from all jobs
        RES = []
        for job in jobs:
            RES.append(job())

        # setup dictionary
        DALL = {}
        for key in RES[0].keys():
            if len(n.shape(RES[0][key])) is 2:
                DALL[key] = n.zeros((ncalc, nalpha))
            else:
                DALL[key] = n.zeros(ncalc)
            
        for i, d in enumerate(RES):
            ppidx = range(i, ncalc, ncpus)
            for key in DALL.keys():
                if len(n.shape(d[key])) is 2:
                    DALL[key][ppidx,:] = d[key]
                else:
                    DALL[key][ppidx] = d[key]
                    

    # single NCPU
    else:
        DALL = _get_Lstar(ticktock, spaco, alpha, extMag, options, omnivals)

    return DALL
	
# -----------------------------------------------
def prep_onera(ticktock=None, spaco=None, alpha=[], extMag='T01STORM', options=[1,0,0,0,0], omnivals=None): 
	"""
	"""
	import numpy as n
	import spacepy.omni as omni

	# setup dictionary to return input values for onera
	d= {}
	d['badval'] = -1e31
	d['nalp_max'] = 25
	d['ntime_max'] = 100000
	d['options'] = options
	badval = d['badval']
	nalp_max = d['nalp_max']
	ntime_max = d['ntime_max']
	
	if ticktock is None:
		return d
	
	UTC = ticktock.UTC
	DOY = ticktock.DOY
	eDOY = ticktock.eDOY
	nTAI = len(ticktock)
	
	# setup mag array and move omni values
	magin = n.zeros((nalp_max,ntime_max),float)	
	magkeys = ['Kp', 'Dst', 'dens', 'velo', 'Pdyn', 'ByIMF', 'BzIMF',\
					'G1', 'G2', 'G3', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
	# get omni values
	if omnivals is None: # nothing provided so use lookup table
		omnivals = omni.getomni(ticktock)
	for iTAI in n.arange(nTAI):
		for ikey, key in enumerate(magkeys):
			magin[ikey, iTAI] = omnivals[key][iTAI]
	d['magin'] = magin
	
	# setup time array
	iyearsat = n.zeros(ntime_max, dtype=int)
	idoysat = n.zeros(ntime_max, dtype=int)
	utsat = n.zeros(ntime_max, dtype=float)
	for i in n.arange(nTAI):
		iyearsat[i] = UTC[i].year
		idoysat[i] = int(DOY[i])
		utsat[i] = (eDOY[i]-n.floor(eDOY[i]))*86400.
	d['iyearsat'] = iyearsat
	d['idoysat'] = idoysat
	d['utsat'] = utsat
	
	# copy coordinates into array
	# prepare coordinates
	d['sysaxes'] = spaco.sysaxes
	xin1 = n.zeros(ntime_max, dtype=float)
	xin2 = n.zeros(ntime_max, dtype=float)
	xin3 = n.zeros(ntime_max, dtype=float) 
	if spaco.carsph == 'sph':
		xin1[0:nTAI] = spaco.radi[:]
		xin2[0:nTAI] = spaco.lati[:]
		xin3[0:nTAI] = spaco.long[:]
	else:
		xin1[0:nTAI] = spaco.x[:]
		xin2[0:nTAI] = spaco.y[:]
		xin3[0:nTAI] = spaco.z[:]
	d['xin1'] = xin1
	d['xin2'] = xin2
	d['xin3'] = xin3

	# convert external magnetic field flag
	extkeys = ['0', 'MEAD', 'T87SHORT', 'T87LONG', 'T89', 'OPQUIET', 'OPDYN', 'T96', \
			   'OSTA', 'T01QUIET', 'T01STORM', 'T05', 'ALEX']
	assert extMag in  extkeys, 'extMag not available: %s' % extMag
	kext = extkeys.index(extMag.upper())
	d['kext'] = kext
 
	# calc at given pitch angels 'alpha'?
	degalpha = n.zeros(nalp_max, dtype=float)
	if isinstance(alpha, float):
		nalpha = 1
		alpha = [alpha]
	nalpha = len(alpha)
	if nalpha > 0:
		degalpha[0:nalpha] = alpha
	d['degalpha'] = degalpha
	
	return d
	
   
	
	

