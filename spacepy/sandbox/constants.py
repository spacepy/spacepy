#!/usr/bin/python
# -*- coding: utf-8 -*-
from spacepy import help
from spacepy.toolbox import dictree
import numpy as np

##Physical constants (SI units))
c = 2.99792458e8 #speed of light in m/s
eVtoK = 1.0/11604.0 #eV to Kelvin
k_b = 1.38065812e-23 #boltzmann constant in J/K
mp = 1.673623110e-27 #mass of proton in kg
me = 9.109389754e-31 #mass of electron in kg
mu0 = 4e-7*np.pi #permeability of free space
q = 1.6021773349e-19 #electronic charge
eps0 = 8.85418782e-12 #permittivity of free space
amu = 1.6605e-27 #atomic mass unit


#Would-be 'convenience class'

class Ion():
    def __init__(self, vel=None, temp=None, units=None, Z=1, verbose=False):
        
        self.vel = vel
        self.temp = temp
        self.units = {'Temp': None, 'Vel': None}
        self.Z = Z
        self.verbose = verbose
        #at least one property must be calculable
        try:
            assert units
            assert self.vel or self.temp
        except AssertionError:
            raise AssertionError
        
        #check for units for temp
        temp_units = ['eV', 'K']
        vel_units = ['m/s']
        
        if units in temp_units:
            self.units['Temp'] = units
            self.vel = self.vIon()
    
    def __str__(self):
        return 'Ion instance: %s' % self.__dict__
        
    __repr__ = __str__


    def vIon(self):
        """Calculate ion velocity [m/s] from temperature
        
        Inputs:
        temp - Temperature in eV or Kelvin
        units (default is 'eV') - 'eV' or 'K'
        Z (default is 1) - ratio of ion mass to proton
        
        Returns:
        v1 - ion velocity in metres per second
        """
        if self.units['Temp'].lower() == 'ev':
            self.vel = np.sqrt(2*self.temp*k_b/(self.Z*mp*eVtoK))
            if self.verbose:
                print('if E = %g eV, then v = %g km/s   for Z = %d mp' % (self.temp, self.vel*1e-3, self.Z))
        if self.units['Temp'].lower() == 'k':
            self.vel = np.sqrt(2*self.temp*k_b/(self.Z*mp))
            if self.verbose:
                print('if E = %g K, then v = %g km/s   for Z = %d mp' % (self.temp, self.vel*1e-3, self.Z))
                
        self.units['Vel'] = 'm/s'
        
        return self.vel

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
            print('if v = %g km/s, then E = %g K ~ %g keV  for Z = %d mp' % (vel*1e-3, e1, e1*eVtoK*1e-3, Z))
            
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
        print('vA = %g km/s    for B = %g nT, no = %d ions/cc, Z = %d mp' % (vA*1e-3, B*1e9, n, Z))
    return vA
