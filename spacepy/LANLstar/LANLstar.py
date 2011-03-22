#!/usr/bin/env python
# 
# 
#  -------------------------------
#  wrapper for calling libLANLstar 
#
#  Version: 08-Mar-2011
#  Author: J. Koller
#  Email: jkoller@lanl.gov
#  -------------------------------
#
#

# ----------------------------------------------------------------
def LANLstar(inputdict):
	"""
	This will calculate Lstar based on the artificial neural network LANLstar which
	was trained with the Tsyganenko & Sitnov (2005) model.
	
    Input:
    ======
        - inputdict (dictionary) : containing the following keys
            ['Year', 'DOY', 'UT_hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', 
		    'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'Lm', 'Bmirr', 'PA', 
		    'rGSM', 'latGSM', 'lonGSM']
            Dictionaries with numpy vectors are allowed.
            
    Returns:
    ========
        - Lstar (float) : returns Lstar as a float or numpy array of floats
        

    Example:
    ========
    >>> import spacepy.LANLstar as LS
    >>> inputdict = {}
    >>> inputdict['Dst'] = 7.7777777777777777 # Dst index (nT)
    >>> inputdict['Pdyn'] = 4.1011111111111118  # solar wind dynamic pressure (nPA)
    >>> inputdict['ByIMF'] = 3.7244444444444444 # GSM y component of IMF magnetic field (nT)
    >>> inputdict['BzIMF'] = -0.12666666666666665  # GSM z component of IMF magnetic field (nT)
    >>> inputdict['W1'] = 0.12244444444444445  # as defined in Tsyganenko and Sitnov 2005
    >>> inputdict['W2'] = 0.2514               
    >>> inputdict['W3'] = 0.089266666666666661 
    >>> inputdict['W4'] = 0.047866666666666668 
    >>> inputdict['W5'] = 0.22586666666666666 
    >>> inputdict['W6'] = 1.0461333333333334 
    >>> 
    >>> # time info
    >>> inputdict['Year'] = 1996    # year 
    >>> inputdict['DOY'] = 6        # day of the year
    >>> inputdict['UT_hr'] = 1.2444444444444445   # UT in units of hours
    >>> 
    >>> # Bfield info
    >>> inputdict['Lm'] = 9.141032669310043   # McIllwain L for 90 deg PA
    >>> inputdict['Bmirr'] = 565.14440185399371 # magnetic field strength at the mirror point
    >>> 
    >>> # coordinates
    >>> inputdict['rGSM'] = 4.83415065 # radial coordinate in GSM [Re]
    >>> inputdict['lonGSM'] = -40.26632902  # longitude coodrinate in GSM [deg]
    >>> inputdict['latGSM'] = 36.44696388  # latitude coordiante in GSM [deg]
    >>> 
    >>> # pitch angle
    >>> inputdict['PA'] = 57.387447236040067 # Pitch angle
    >>> 
    >>> LS.LANLstar(inputdict)
    8.4505511125920574

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 22-Mar-2011 (JK)
	"""
	
	import libLANLstar
	import numpy as n
	
	if isinstance(inputdict['Dst'], float):
		arrayflag = False
		for key in inputdict.keys():
			inputdict[key] = [inputdict[key]]
	else:
		arrayflag = True
	
	ncalc = len(inputdict['Dst'])
	Lstar = n.zeros(ncalc)
	
	keylist = ['Year', 'DOY', 'UT_hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
		'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'Lm', 'Bmirr', 'PA', \
		'rGSM', 'latGSM', 'lonGSM']
	inpar = n.zeros(len(keylist))
	
	for i in range(ncalc):	
		# copy over keylist into inpar
		for ikey, key in enumerate(keylist):
			inpar[ikey] = inputdict[key][i]
		Lstar[i] = libLANLstar.lanlstar(inpar)

	if arrayflag is False:
		return Lstar[0]
	else:
		return Lstar


# ----------------------------------------------------------------
def LANLmax(inputdict):
	"""
	This will calculate Lstar_max (last closed drift shell) based on the artificial neural network LANLstar which
	was trained with the Tsyganenko & Sitnov (2005) model.
	
    Input:
    ======
        - inputdict (dictionary) : containing the following keys
            ['Year', 'DOY', 'UT_hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
		    'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'PA']
            Dictionaries with numpy vectors are allowed.
            
    Returns:
    ========
        - Lstar (float) : returns Lstar as a float or numpy array of floats
        

    Example:
    ========
    >>> import spacepy.LANLstar as LS
    >>> inputdict = {}
    >>> inputdict['Dst'] = 7.7777777777777777 # Dst index (nT)
    >>> inputdict['Pdyn'] = 4.1011111111111118  # solar wind dynamic pressure (nPA)
    >>> inputdict['ByIMF'] = 3.7244444444444444 # GSM y component of IMF magnetic field (nT)
    >>> inputdict['BzIMF'] = -0.12666666666666665  # GSM z component of IMF magnetic field (nT)
    >>> inputdict['W1'] = 0.12244444444444445  # as defined in Tsyganenko and Sitnov 2005
    >>> inputdict['W2'] = 0.2514               
    >>> inputdict['W3'] = 0.089266666666666661 
    >>> inputdict['W4'] = 0.047866666666666668 
    >>> inputdict['W5'] = 0.22586666666666666 
    >>> inputdict['W6'] = 1.0461333333333334 
    >>> 
    >>> # pitch angle
    >>> inputdict['PA'] = 57.387447236040067 # Pitch angle
    >>> 
    >>> LS.LANLmax(inputdict)
    9.8981001125908552

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 22-Mar-2011 (JK)
	"""
	
	import libLANLstar
	import numpy as n
	
	if isinstance(inputdict['Dst'], float):
		arrayflag = False
		for key in inputdict.keys():
			inputdict[key] = [inputdict[key]]
	else:
		arrayflag = True
	
	ncalc = len(inputdict['Dst'])
	Lmax = n.zeros(ncalc)
	
	keylist = ['Year', 'DOY', 'UT_hr', 'Dst', 'Pdyn', 'ByIMF', 'BzIMF', \
		'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'PA']
	inpar = n.zeros(len(keylist))
	
	for i in range(ncalc):	
		# copy over keylist into inpar
		for ikey, key in enumerate(keylist):
			inpar[ikey] = inputdict[key][i]
		Lmax[i] = libLANLstar.lanlmax(inpar)

	if arrayflag is False:
		return Lmax[0]
	else:
		return Lmax

