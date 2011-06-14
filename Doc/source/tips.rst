==============================
SpacePy Python Progamming Tips
==============================

One often hears that interpreted langauges are too slow for whatever task someone
needs to do.  In many cases this is exactly the wrong-headed.  As the time spent
programming/debugging in an interpeted language is less than a compiled language
the programmer has time to figure out where code is slow and make it faster.  This
page is didicated to that idea, providing examples of code speedup and best practices.

Lists, for loops, and arrays
============================
This example teaches the lesson that every IDL_ or Matlab_ programmer already
knows; do everything in arrays and never use a for loop.

This bit of code takes in a series of points, computes their magnitude, and drops
the larges 100 of them.

This is how the code started out::

    def SortRemove_HighFluxPts_(Shell_x0_y0_z0, ShellCenter, Num_Pts_Removed):
        #Sort the Shell Points based on radial distance (Flux prop to 1/R^2) and remove Num_Pts_Removed points with the highest flux
        Num_Pts_Removed = np.abs(Num_Pts_Removed)  #make sure the number is positive
        #Generate an array of radial distances of points from origin
        R = []
        for xyz in Shell_x0_y0_z0:
            R.append(1/np.linalg.norm(xyz + ShellCenter)) #Flux prop to 1/r^2, but don't need the ^2
        R = np.asarray(R)
        ARG = np.argsort(R)   # array of sorted indeces based on flux in 1st column
        Shell_x0_y0_z0 = np.take(Shell_x0_y0_z0, ARG, axis = 0)  # sort based on index order
        return Shell_x0_y0_z0[:-Num_Pts_Removed,:]   #remove last points that have the "anomalously" high flux



.. _IDL: http://www.ittvis.com/language/en-us/productsservices/idl.aspx
.. _Matlab: http://www.mathworks.com/products/matlab/


Zip
===
the zip_ function is a great thing but it is really slow, if you find yourself
using it then you probably need to reexamine the algorithm that you are using.

This example generate evenly distributed N points on the unit sphere centered at
(0,0,0) using the "Golded Spiral" method.

The origional code::

    def PointsOnSphere_(N):
    # Generate evenly distributed N points on the unit sphere centered at (0,0,0)
    # Uses "Golden Spiral" method
        x0 = np.array((N,), dtype= float)
        y0 = np.array((N,), dtype= float)
        z0 = np.array((N,), dtype= float)
        phi = (1 + np.sqrt(5)) / 2. # the golden ratio
        long_incr = 2.0*np.pi / phi # how much to increment the longitude
        dz = 2.0 / float(N)    # a unit sphere has diameter 2
        bands = np.arange(0, N, 1) # each band will have one point placed on it
        z0 = bands * dz - 1 + (dz/2)  # the z location of each band/point
        r = np.sqrt(1 - z0*z0)    # the radius can be directly determined from height
        az = bands * long_incr # the azimuth where to place the point
        x0 = r * np.cos(az)
        y0 = r * np.sin(az)
        x0_y0_z0 = np.array(zip(x0,y0,z0))     #combine into 3 column (x,y,z) file
        return (x0_y0_z0)

Profileing this with cProfile one can see a lot of time in zip()::

    Tue Jun 14 09:54:41 2011    PointsOnSphere_.prof

             9 function calls in 8.132 seconds

       Ordered by: cumulative time

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.010    0.010    8.132    8.132 <string>:1(<module>)
            1    0.470    0.470    8.122    8.122 test1.py:192(PointsOnSphere_)
            4    6.993    1.748    6.993    1.748 {numpy.core.multiarray.array}
            1    0.654    0.654    0.654    0.654 {zip}
            1    0.005    0.005    0.005    0.005 {numpy.core.multiarray.arange}
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}

So lets try and do a few simple rewites to make this faster.  Using numpy.vstack
is the first one that came to mind.  The change here is to replace building up
the array from the elemets made by zip to just appending the data we already have
to an array that we already have::

    def PointsOnSphere_(N):
    # Generate evenly distributed N points on the unit sphere centered at (0,0,0)
    # Uses "Golded Spiral" method
        x0 = np.array((N,), dtype= float)
        y0 = np.array((N,), dtype= float)
        z0 = np.array((N,), dtype= float)
        phi = (1 + np.sqrt(5)) / 2. # the golden ratio
        long_incr = 2.0*np.pi / phi # how much to increment the longitude
        dz = 2.0 / float(N)    # a unit sphere has diameter 2
        bands = np.arange(0, N, 1) # each band will have one point placed on it
        z0 = bands * dz - 1 + (dz/2)  # the z location of each band/point
        r = np.sqrt(1 - z0*z0)    # the radius can be directly determined from height
        az = bands * long_incr # the azimuth where to place the point
        x0 = r * np.cos(az)
        y0 = r * np.sin(az)
        x0_y0_z0 = np.vstack((x0, y0, z0)).transpose()
        return (x0_y0_z0)

Profileing this with cProfile one can see that this is now fast enough for me,
no more work to do.  We picked up a 48x speed increase, I;m sure this can still
be made better and let the spacepy team know if you rewite it and it will be
included::

    Tue Jun 14 09:57:41 2011    PointsOnSphere_.prof

             32 function calls in 0.168 seconds

       Ordered by: cumulative time
       List reduced from 13 to 10 due to restriction <10>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.010    0.010    0.168    0.168 <string>:1(<module>)
            1    0.123    0.123    0.159    0.159 test1.py:217(PointsOnSphere_)
            1    0.000    0.000    0.034    0.034 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/shape_base.py:177(vstack)
            1    0.034    0.034    0.034    0.034 {numpy.core.multiarray.concatenate}
            1    0.002    0.002    0.002    0.002 {numpy.core.multiarray.arange}
            1    0.000    0.000    0.000    0.000 {map}
            3    0.000    0.000    0.000    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/shape_base.py:58(atleast_2d)
            6    0.000    0.000    0.000    0.000 {numpy.core.multiarray.array}
            3    0.000    0.000    0.000    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/numeric.py:237(asanyarray)
            1    0.000    0.000    0.000    0.000 {method 'transpose' of 'numpy.ndarray' objects}


.. _zip: http://docs.python.org/library/functions.html#zip

--------------------------

:Release: |version|
:Doc generation date: |today|
