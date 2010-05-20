#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division
from spacepy import help
import numpy as np

##Physical constants (SI units))
c = 2.99792458e8 #speed of light in m/s
eVtoK = 1/11604 #eV to Kelvin
k_b = 1.38065812e-23 #boltzmann constant in J/K
mp = 1.673623110e-27 #mass of proton in kg
me = 9.109389754e-31 #mass of electron in kg
mu0 = 4e-7*np.pi #permeability of free space
q = 1.6021773349e-19 #electronic charge
eps0 = 8.85418782e-12 #permittivity of free space
amu = 1.6605e-27 #atomic mass unit


class Ion():
    def __init__(self, vel=None, temp=None, units=None, Z=1, verbose=False):
        
        self.vel = vel
        self.temp = temp
        self.units = units
        self.Z = Z
        self.verbose = verbose
        try:
            assert self.units
            assert self.vel or self.temp
        except AssertionError:
            raise AssertionError
        
        #check for units for temp
        temp_units = ['eV', 'K']
        vel_units = ['m/s']


def vIon(temp, units='eV', Z=1, verbose=False):
    """Calculate ion velocity [m/s] from temperature
    
    Inputs:
    temp - Temperature in eV or Kelvin
    units (default is 'eV') - 'eV' or 'K'
    Z (default is 1) - ratio of ion mass to proton
    
    Returns:
    v1 - ion velocity in metres per second
    """
    try:
        assert type(units.lower()) == str
    except:
        print 'vIon error: units keyword not valid'
        return None
    if units.lower() == 'ev':
        v1 = np.sqrt(2*temp*k_b/(Z*mp*eVtoK))
        if verbose:
            print 'if E = %g eV, then v = %g km/s   for Z = %d mp' % (temp, v1*1e-3, Z)
    if units.lower() == 'k':
        v1 = np.sqrt(2*temp*k_b/(Z*mp))
        if verbose:
            print 'if E = %g K, then v = %g km/s   for Z = %d mp' % (temp, v1*1e-3, Z)
    return v1

def eIon(vel, units='both', Z=1, verbose=False):
    """Calculate ion energy [eV or K] from velocity
    
    Inputs:
    vel - ion velocity in m/s
    units (for output; default is 'both') - 'eV' or 'K'
    Z (default is 1) - ratio of ion mass to proton
    
    Returns:
    energy - tuple of energy in eV and K.
    Single values, in a tuple, are returned for non-default units
    """
    e1 = 0.5*mp*Z*vel**2/(k_b)
    if verbose:
        print 'if v = %g km/s, then E = %g K ~ %g keV  for Z = %d mp' % (vel*1e-3, e1, e1*eVtoK*1e-3, Z)
        
    if units.lower() == 'both':
        return (e1, e1*eVtoK)
    if units.lower() == 'ev':
        return (e1)
    if units.lower() == 'k':
        return (e1*eVtoK)

def vAlfven(B, n, Z=1, verbose=False):
    """Calculate Alfven velocity [m/s] for input plasma paramters
    
    Inputs:
    B - magnetic field strength [nT]
    n - plasma number density [cm^{-3}]
    Z (default is 1) - ratio of ion mass to proton
    
    Returns:
    vA - Alfven velocity
    """
    rho = Z*mp*n*100**3 #mass density [kg/m^3]
    vA = B/np.sqrt(rho*mu0)
    if verbose:
        print 'vA = %g km/s    for B = %g nT, no = %d ions/cc, Z = %d mp' % (vA*1e-3, B*1e9, n, Z)
    return vA