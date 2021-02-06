#!/usr/bin/env python
'''
A set of routines for fast field line tracing.
"Number crunching" is performed in C for speed.

Copyright 2010-2014 Los Alamos National Security, LLC.
'''

import ctypes
import sys
import warnings

import numpy
from .. import lib

assert(lib.have_libspacepy)


def _trace2d_common(func, fieldx, fieldy, xstart, ystart, gridx, gridy,
                    maxstep=20000, ds=0.01):
    #These are a bit redundant with the checking in the lib.py
    #call dict, but need to make sure they're numpy arrays so can get shape.
    gridx = numpy.require(gridx, ctypes.c_double, 'C')
    gridy = numpy.require(gridy, ctypes.c_double, 'C')
    fieldx = numpy.require(fieldx, ctypes.c_double, 'C')
    fieldy = numpy.require(fieldy, ctypes.c_double, 'C')
    outx = numpy.empty(shape=(maxstep,), dtype=ctypes.c_double)
    outy = numpy.empty(shape=(maxstep,), dtype=ctypes.c_double)
    count = func(gridx.shape[0], gridy.shape[0], maxstep, ds, xstart, ystart,
		 gridx, gridy, fieldx, fieldy, outx, outy)
    return (outx[0:count], outy[0:count])

def trace2d_eul(fieldx, fieldy, xstart, ystart, gridx, gridy,
                maxstep=20000, ds=0.01):
    """
    Given a 2D vector field, trace a streamline from a given point
    to the edge of the vector field.  The field is integrated using
    Euler's method.  While this is faster than rk4, it is less accurate.

    Only valid for regular grid with coordinates gridx, gridy.
    If gridx and gridy are not given, assume that xstart and ystart
    are normalized coordinates (e.g., position in terms of array
    indices.)
    """
    return _trace2d_common(lib.cEuler, fieldx, fieldy, xstart, ystart,
                           gridx, gridy, maxstep, ds)


def trace2d_rk4(fieldx, fieldy, xstart, ystart, gridx, gridy,
                maxstep=20000, ds=0.01):
    """
    Given a 2D vector field, trace a streamline from a given point
    to the edge of the vector field.  The field is integrated using
    Runge Kutta 4.  Slower than Euler, but more accurate.  The
    higher accuracy allows for larger step sizes (ds kwarg).  For
    a demonstration of the improved accuracy, run test_asymtote and
    test_dipole, bouth found in the pybats.trace2d module.
    Only valid for regular grid with coordinates gridx, gridy.
    If gridx and gridy are not given, assume that xstart and ystart
    are normalized coordinates (e.g., position in terms of array
    indices.)
    """
    return _trace2d_common(lib.cRk4, fieldx, fieldy, xstart, ystart,
                           gridx, gridy, maxstep, ds)


###################################################
# TEST SUITE  #
###################################################
def test_asymtote():
    '''
    Test streamline tracing by plotting vectors and associated streamlines
    through a simple velocity field where Vx=x, Vy=-y.
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    # Start by creating a velocity vector field.
    xmax = 200.0
    ymax = 20.0
    x = np.arange(-10.0, xmax + 0.25, 0.25)
    y = np.arange(-10.0, ymax + 0.25, 0.25)
    
    xgrid, ygrid = np.meshgrid(x,y)
    
    vx = xgrid * 1.0
    vy = ygrid * -1.0

    xstart = 1.0
    ystart = 10.0
    x1=[0,0]; y1=[0,0]
    x2=[0,0]; y2=[0,0]
    (x1, y1) = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=0.1)
    (x2, y2) = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=0.5)
    (x3, y3) = trace2d_rk4(vx, vy, xstart, ystart, x, y, ds=1.0)
    (x4, y4) = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=0.1)
    (x5, y5) = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=0.5)
    (x6, y6) = trace2d_eul(vx, vy, xstart, ystart, x, y, ds=1.0)

    # analytical solution (1/const) = x*y:
    const = 1 / (xstart * ystart)
    x_anly = np.arange(-10.0, xmax, 0.001)
    y_anly = 1.0 / (const * x_anly)

    fig = plt.figure(1, figsize=(10,8))
    ax1 = plt.subplot('111')
    
    #ax1.quiver(x, y, vx, vy)
    ax1.plot(x_anly, y_anly, 'k', label='Analytic',linewidth=3.0)
    ax1.plot(x1, y1, 'b',   label='RK4 ds=0.1', linewidth=1.5)
    ax1.plot(x2, y2, 'b--', label='RK4 ds=0.5', linewidth=1.5)
    ax1.plot(x3, y3, 'b:',  label='RK4 ds=1.0', linewidth=1.5)
    ax1.plot(x4, y4, 'g',   label='Euler ds=0.1', linewidth=.75)
    ax1.plot(x5, y5, 'g--', label='Euler ds=0.5', linewidth=.75)
    ax1.plot(x6, y6, 'g:',  label='Euler ds=1.0', linewidth=.75)
    ax1.legend()
    ax1.set_title("Runge Kutta 4 vs Euler's: Asymtotic Field")
    ax1.set_xlabel("Normalized 'X' Coordinate")
    ax1.set_ylabel("Normalized 'Y' Coordinate")
    ax1.set_xlim( [3, 30] )
    ax1.set_ylim( [0.25, 2.5 ] )

    # Annotate plot.
    ax1.annotate("Euler's method diverges strongly\nalong curves except when"+
                 " taking \nvery small steps.  RK4 is\nfar more accurate "+
                 "for all dS.", 
                 xy=(10.5,0.83), xycoords='data', xytext=(9,.3),
                 arrowprops=dict(fc='black',shrink=0.05), 
                 horizontalalignment='center')
    ax1.annotate("Tracing begins at x=1, y=10.",
                 xy=(4.8,2.45), xycoords='data', xytext=(6,1.9),
                 arrowprops=dict(fc='black',shrink=0.05))
    ax1.annotate("The trace continues until \n"+
                 "x=200.  At that point, Euler \n"+
                 "dS=0.1 and RK4 dS=1.0 converge\n"+
                 " at the same point, despite \n"+
                 "the 10X difference in step size.",
                 xy=(29,.5), xycoords='data', xytext=(20,1),
                 arrowprops=dict(fc='black',shrink=0.005), 
                 horizontalalignment='center')
    plt.show()

    return ax1

def test_dipole():
    '''
    Trace field lines through a dipole field to test
    the pybats.trace2d module.
    '''
    import numpy as np
    import matplotlib.pyplot as plt
    from spacepy.pybats.dipole import (b_hat, b_line)

    # Now do dipole magnetic field.
    # Start by creating a field of unit vectors...
    x = np.arange(-100.0, 101.0, 5.0)
    y = np.arange(-100.0, 101.0, 5.0)
    xgrid, ygrid = np.meshgrid(x,y)

    bx, by = b_hat(x,y)

    # New figure.
    fig2 = plt.figure(figsize=(14,8))
    fig2.subplots_adjust(wspace=0.15, left=0.08, right=0.94)
    ax2 = plt.subplot('121')
    ax3 = plt.subplot('322')
    ax4a= plt.subplot('347')
    ax4b= plt.subplot('348')
    ax5 = plt.subplot('326')
    ax2.quiver(x,y, bx, by, units='x', pivot='middle')
    ax3.quiver(x,y, bx, by, units='x', pivot='tip')
    ax4a.quiver(x,y,bx, by, units='x', pivot='tip')
    ax4b.quiver(x,y,bx, by, units='x', pivot='tip')
    ax5.quiver(x,y, bx, by, units='x', pivot='tip')

    # Trace through this field.
    xstart = 10.0
    ystart = 25.0
    ds = 0.1
    for ystart in range(0, 31, 5):
        (x1, y1) = trace2d_rk4(bx, by, xstart, ystart, x, y, ds=ds)
        l1 = ax2.plot(x1,y1,'b')[0]
        ax3.plot(x1,y1,'b'); ax4b.plot(x1,y1,'b')
        ax5.plot(x1,y1,'b'); ax4a.plot(x1,y1,'b')
        (x2, y2) = trace2d_eul(bx, by, xstart, ystart, x, y, ds=ds)
        l2 = ax2.plot(x2,y2,'r')[0]
        ax3.plot(x2,y2,'r'); ax4b.plot(x2,y2,'r')
        ax5.plot(x2,y2,'r'); ax4a.plot(x2,y2,'r')
        (x3, y3) = b_line(xstart, ystart, npoints=300)
        l3 = ax2.plot(x3,y3,'k--')[0]
        ax3.plot(x3,y3,'k--'); ax4b.plot(x3,y3,'k--')
        ax5.plot(x3,y3,'k--'); ax4a.plot(x3,y3,'k--')
    
    ax2.set_xlim([-2,  100])
    ax2.set_ylim([-30, 100])
    ax2.set_title("Full View")
    ax2.set_xlabel("Normalized 'X' Coordinate")
    ax2.set_ylabel("Normalized 'Y' Coordinate")
    ax2.legend( (l1, l2, l3),('RK4', 'Euler', 'Analytical'), 'upper left' )

    ax3.set_title("Zoomed Views")
    ax3.set_xlim([8.5, 17.5])
    ax3.set_ylim([3, 33])
    goodpos = ax3.get_position()

    ax4a.set_xlim([20,30])
    ax4a.set_ylim([-12,12])
    pos = ax4a.get_position()
    pos.x0 = goodpos.x0
    pos.x1 = pos.x0 + (goodpos.x1-goodpos.x0)/2.0 -0.01
    ax4a.set_position(pos)

    ax4b.set_xlim([50,60])
    ax4b.set_ylim([-12,12])
    pos = ax4b.get_position()
    pos.x0 = goodpos.x0 + (goodpos.x1-goodpos.x0)/2.0 +0.01
    ax4b.set_position(pos)
    ax4b.set_yticklabels('', visible=False)

    ax5.set_xlim([1,7])
    ax5.set_ylim([-7,-3])
    ax5.set_xlabel("Normalized 'X' Coordinate")

    fig2.suptitle("RK4 vs Euler's Method: Dipole Field for dS=%.3f" % ds)
    plt.show()
    
if __name__ == '__main__':
    test_asymtote()
    test_dipole()
