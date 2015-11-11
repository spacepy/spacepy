'''
Functions for common math problems such as derivatives, etc.
These should typically be called via other interfaces, but are made available
to all users.
'''

import numpy as np


def d_dx(U, dx):
    '''
    Given a 2D array of U values that is *regularly spaced* and is ordered using
    'C' or 'matplotlib' ordering (such that $x$ progresses along the last 
    index, $y$ along the next-to-last, etc.), 
    take spatial the derivative with respect to $x$.
    A 2D array of dU/dx values are returned.  Uses second order 
    central differences (non-edge values) and second order forward/backward 
    differences (edge values) to obtain first derivative without ghost cells.
    '''
    
    ny,nx=U.shape
    du_dx=np.zeros( (ny,nx) )

    # Central differences for central x locations.
    du_dx[:,1:-1] = (U[:,2:] - U[:,0:-2]) / (2.*dx)
    # Forward differences for minimum x locations.
    du_dx[:,0] = (-3.*U[:,0] + 4.*U[:,1] - U[:,2]) / (2.*dx)
    # Backward differences for maximum x locations.
    du_dx[:,-1] = (3.*U[:,-1] - 4.*U[:,-2] + U[:,-3]) / (2.*dx)

    return du_dx

def d_dy(U, dy):
    '''
    Given a 2D array of U values that is *regularly spaced* and is ordered using
    'C' or 'matplotlib' ordering (such that $x$ progresses along the last 
    index, $y$ along the next-to-last, etc.), 
    take spatial the derivative with respect to $y$.
    A 2D array of dU/dx values are returned.  Uses second order 
    central differences (non-edge values) and second order forward/backward 
    differences (edge values) to obtain first derivative without ghost cells.
    '''
    ny,nx=U.shape
    du_dy=np.zeros( (ny,nx) )

    # Central differences for central x locations.
    du_dy[1:-1,:] = (U[2:,:] - U[0:-2,:]) / (2.*dy)
    # Forward differences for minimum x locations.
    du_dy[0,:]  = (-3.*U[0,:] + 4.*U[1,:]  - U[2,:]) / (2.*dy)
    # Backward differences for maximum x locations.
    du_dy[-1,:] = (3.*U[-1,:] - 4.*U[-2,:] + U[-3,:]) / (2.*dy)

    return du_dy

def interp_2d_reg(x, y, xgrid, ygrid, val, dx=None, dy=None):
    '''
    For a set of points (*x*, *y*) that lie inside of the 2D arrays of x and y
    locations, *xgrid*, *ygrid*, interpolate 2D array of values, *val* to 
    those points using simple bilinear interpolation.  This function will
    extrapolate outside of *xgrid*, *ygrid*, so use with caution.
    '''

    # Change from matrices to vectors:
    if xgrid.ndim>1:
        xgrid = xgrid[0,:]
        ygrid = ygrid[:,0]
        
    # Get resolution if not supplied:
    if (not dx) or (not dy):
        dx = xgrid[1]-xgrid[0]
        dy = ygrid[1]-ygrid[0]
    ySize, xSize = ygrid.size, xgrid.size

    # Normalized coords:
    xNorm = (x-xgrid[0])/dx
    yNorm = (y-ygrid[0])/dy

    # LowerLeft index of four surrounding points.
    # Bind points such that four surrounding always within xgrid, ygrid.
    # It is this binding that forces extrapolation!
    xLL = np.array(np.floor(xNorm), dtype=int)
    yLL = np.array(np.floor(yNorm), dtype=int)
    xLL[xLL>(xSize-2)] = xSize-2; xLL[xLL<0] = 0
    yLL[yLL>(ySize-2)] = ySize-2; yLL[yLL<0] = 0

    # Re-normalize to LL values.
    xNorm=xNorm-xLL
    yNorm=yNorm-yLL

    # Interpolate.
    out = \
        (val[yLL  , xLL  ] * (1.0-xNorm) * (1.0-yNorm) ) + \
        (val[yLL  , xLL+1] * (    xNorm) * (1.0-yNorm) ) + \
        (val[yLL+1, xLL  ] * (1.0-xNorm) * (    yNorm) ) + \
        (val[yLL+1, xLL+1] * (    xNorm) * (    yNorm) )
    
    return out

def interp_bilin_scalar(x, y, z, xMin=0., yMin=0., dx=1., dy=1.):
    '''
    Fast, simple bilinear interpolation from 2d regular grid with starting
    points *xMin*, *yMin* and normalized spacing *dx*, *dy* of values on grid
    (*z*) to new location, *x*, *y*.  Used to quickly set up ghost cells for
    advanced tracing.
    '''
    x = (np.asarray(x)-xMin)/dx
    y = (np.asarray(y)-yMin)/dy

    # Get indices encircling interpolation point:
    x0 = np.floor( x )
    x1 = x0 + 1
    y0 = np.floor( y )
    y1 = y0 + 1

    # Bind location indices to array limits:
    x0 = np.clip(x0, 0, z.shape[1]-2);
    x1 = np.clip(x1, 1, z.shape[1]-1);
    y0 = np.clip(y0, 0, z.shape[0]-2);
    y1 = np.clip(y1, 1, z.shape[0]-1);

    Q00 = z[ y0, x0 ]
    Q10 = z[ y1, x0 ]
    Q01 = z[ y0, x1 ]
    Q11 = z[ y1, x1 ]

    wa = (x1-x) * (y1-y)
    wb = (x1-x) * (y-y0)
    wc = (x-x0) * (y1-y)
    wd = (x-x0) * (y-y0)

    return wa*Q00 + wb*Q10 + wc*Q01 + wd*Q11
