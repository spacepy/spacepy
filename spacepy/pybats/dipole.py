#!/usr/bin/env python
'''
Some functions for the generation of a dipole field.

Copyright 2010 Los Alamos National Security, LLC.
'''

import numpy as np
import math
import pylab as plt

def b_mag(x,y):
    '''
    For a position *x*, *y* in units of planetary radius, return the
    strength of the dipole magnetic field in nanoTesla.
    '''
    r = np.sqrt(x**2 + y**2)
    cos = y/r

    return 30400 * np.sqrt(1+3*cos**2)/r**3
    
def b_hat(x, y):
    ''' 
    For given parameters, return two arrays, x and y, corresponding
    to the x and y components of b_hat for a dipole field.  Plotting these
    two matrices using MatPlotLib's quiver function will create a beautiful
    dipole field for tracing and other stuff.
    '''
    xgrid, ygrid = np.meshgrid(x,y)
    
    r = np.sqrt(xgrid**2 + ygrid**2)
    cos = ygrid/r
    sin = xgrid/r

    denom = np.sqrt(1.0 + 3.0*cos**2)

    b_r     = 2.0 * cos / denom
    b_theta =       sin / denom

    b_x = b_r*sin + b_theta*cos
    b_y = b_r*cos - b_theta*sin

    return(b_x, b_y)

def b_line(x, y, npoints=30):
    '''For a starting X, Y point return x and y vectors that trace the dipole
    field line that passes through the given point.'''
    npoints = float(npoints)
    r = np.sqrt(x**2 + y**2)
    try:
        theta = np.arctan(x/y)
    except ZeroDivisionError:
        theta = math.pi/2.0
    R = r/(np.sin(theta)**2)

    if x<0:
        theta = np.arange(math.pi, 2.0*math.pi, math.pi/npoints)
    else:
        theta = np.arange(0, math.pi, math.pi/npoints)
    r_vec = R * np.sin(theta)**2

    x_out = r_vec * np.sin(theta)
    y_out = r_vec * np.cos(theta)

    return (x_out, y_out)

def test():
    '''
    A quick test of the dipole field functions.
    '''
    x = np.arange(-100.0, 101.0, 5.0)
    y = np.arange(-100.0, 101.0, 5.0)

    x_vec, y_vec = b_hat(x,y)

    fig = plt.figure(figsize=(10,8))
    ax1 = plt.subplot('111')

    ax1.quiver(x,y, x_vec, y_vec)

    for i in range(-120, 121, 10):
        (x,y) = b_line(float(i), 0.0, 100)
        ax1.plot(x, y, 'b')
    for theta in np.arange(math.pi/2.0, 3.0*math.pi/2.0, math.pi/100.0):
        x = np.sin(theta)
        y = np.cos(theta)
        (x,y) = b_line(x, y, 100)
        ax1.plot(x, y,'r')
    ax1.set_xlim( [-100, 100] )
    ax1.set_ylim( [-100, 100] )
    plt.title('Unit vectors for an arbitrary dipole field')
    
    fig.show()

if __name__ == '__main__':
    test()
    
