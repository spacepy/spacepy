==============================
SpacePy Python Programming Tips
==============================

One often hears that interpreted languages are too slow for whatever task someone
needs to do.  In many cases this is exactly the wrong-headed.  As the time spent
programming/debugging in an interpreted language is less than a compiled language
the programmer has time to figure out where code is slow and make it faster.  This
page is dedicated to that idea, providing examples of code speedup and best practices.

Lists, for loops, and arrays
============================
This example teaches the lesson that every IDL_ or Matlab_ programmer already
knows; do everything in arrays and never use a for loop.

This bit of code takes in a series of points, computes their magnitude, and drops
the larges 100 of them.

This is how the code started out, Shell_x0_y0_z0 is a Nx3 numpy array,
ShellCenter is a 3 element list or array, and Num_Pts_Removed is the number of
points to drop::

    import numpy as np
    def SortRemove_HighFluxPts_(Shell_x0_y0_z0, ShellCenter, Num_Pts_Removed):
        #Sort the Shell Points based on radial distance (Flux prop to 1/R^2) and remove Num_Pts_Removed points with the highest flux
        Num_Pts_Removed = np.abs(Num_Pts_Removed)  #make sure the number is positive
        #Generate an array of radial distances of points from origin
        R = []
        for xyz in Shell_x0_y0_z0:
            R.append(1/np.linalg.norm(xyz + ShellCenter)) #Flux prop to 1/r^2, but don't need the ^2
        R = np.asarray(R)
        ARG = np.argsort(R)   # array of sorted indies based on flux in 1st column
        Shell_x0_y0_z0 = np.take(Shell_x0_y0_z0, ARG, axis = 0)  # sort based on index order
        return Shell_x0_y0_z0[:-Num_Pts_Removed,:]   #remove last points that have the "anomalously" high flux

A cProfile of this yields a lot of time spend just in the function itself, this
is the for loop (list comprehension is a little faster but not much in this case)::

    Tue Jun 14 10:10:56 2011    SortRemove_HighFluxPts_.prof

             700009 function calls in 4.209 seconds

       Ordered by: cumulative time
       List reduced from 14 to 10 due to restriction <10>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.002    0.002    4.209    4.209 <string>:1(<module>)
            1    2.638    2.638    4.207    4.207 test1.py:235(SortRemove_HighFluxPts_)
       100000    0.952    0.000    1.529    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/linalg/linalg.py:1840(norm)
       100001    0.099    0.000    0.240    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/numeric.py:167(asarray)
       100000    0.229    0.000    0.229    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       100001    0.141    0.000    0.141    0.000 {numpy.core.multiarray.array}
       100000    0.082    0.000    0.082    0.000 {method 'ravel' of 'numpy.ndarray' objects}
       100000    0.042    0.000    0.042    0.000 {method 'conj' of 'numpy.ndarray' objects}
       100000    0.016    0.000    0.016    0.000 {method 'append' of 'list' objects}
            1    0.000    0.000    0.005    0.005 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/fromnumeric.py:45(take)

Simply pulling out the addition inside the for loop makes an amazing difference
(2.3x speedup).  We believe the difference is that pulling out the addition lets
numpy do its thing in C and not in python for each element::

    def SortRemove_HighFluxPts_(Shell_x0_y0_z0, ShellCenter, Num_Pts_Removed):
        #Sort the Shell Points based on radial distance (Flux prop to 1/R^2) and remove Num_Pts_Removed points with the highest flux
        Num_Pts_Removed = np.abs(Num_Pts_Removed)  #make sure the number is positive
        #Generate an array of radial distances of points from origin
        R = []
        Shell_xyz = Shell_x0_y0_z0 + ShellCenter
        for xyz in Shell_xyz:
            R.append(1/np.linalg.norm(xyz)) #Flux prop to 1/r^2, but don't need the ^2
        R = np.asarray(R)
        ARG = np.argsort(R)   # array of sorted indies based on flux in 1st column
        Shell_x0_y0_z0 = np.take(Shell_x0_y0_z0, ARG, axis = 0)  # sort based on index order
        return Shell_x0_y0_z0[:-Num_Pts_Removed,:]   #remove last points that have the "anomalously" high flux

A quick profile::

    Tue Jun 14 10:18:39 2011    SortRemove_HighFluxPts_.prof

             700009 function calls in 1.802 seconds

       Ordered by: cumulative time
       List reduced from 14 to 10 due to restriction <10>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.001    0.001    1.802    1.802 <string>:1(<module>)
            1    0.402    0.402    1.801    1.801 test1.py:235(SortRemove_HighFluxPts_)
       100000    0.862    0.000    1.361    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/linalg/linalg.py:1840(norm)
       100000    0.207    0.000    0.207    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       100001    0.080    0.000    0.199    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/numeric.py:167(asarray)
       100001    0.120    0.000    0.120    0.000 {numpy.core.multiarray.array}
       100000    0.067    0.000    0.067    0.000 {method 'ravel' of 'numpy.ndarray' objects}
       100000    0.041    0.000    0.041    0.000 {method 'conj' of 'numpy.ndarray' objects}
       100000    0.014    0.000    0.014    0.000 {method 'append' of 'list' objects}
            1    0.000    0.000    0.005    0.005 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/fromnumeric.py:45(take)

A closer look here reveals that all of this can be done on the arrays without
the for loop (or list comprehension)::

    def SortRemove_HighFluxPts_(Shell_x0_y0_z0, ShellCenter, Num_Pts_Removed):
        #Sort the Shell Points based on radial distance (Flux prop to 1/R^2) and remove # points with the highest flux
        Num_Pts_Removed = np.abs(Num_Pts_Removed)  #make sure the number is positive
        #Generate an array of radial distances of points from origin
        R = 1 / np.sum((Shell_x0_y0_z0 + ShellCenter) ** 2, 1)
        ARG = np.argsort(R)   # array of sorted indies based on flux in 1st column
        Shell_x0_y0_z0 = np.take(Shell_x0_y0_z0, ARG, axis = 0)  # sort based on index order
        return Shell_x0_y0_z0[:-Num_Pts_Removed,:]   #remove last points that have the "anomalously" high flux

The answer is exactly the same and from where we started this is a 382x speedup::

    Tue Jun 14 10:21:54 2011    SortRemove_HighFluxPts_.prof

             10 function calls in 0.011 seconds

       Ordered by: cumulative time

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.000    0.000    0.011    0.011 <string>:1(<module>)
            1    0.002    0.002    0.011    0.011 test1.py:236(SortRemove_HighFluxPts_)
            1    0.000    0.000    0.004    0.004 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/fromnumeric.py:598(argsort)
            1    0.004    0.004    0.004    0.004 {method 'argsort' of 'numpy.ndarray' objects}
            1    0.000    0.000    0.003    0.003 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/fromnumeric.py:45(take)
            1    0.003    0.003    0.003    0.003 {method 'take' of 'numpy.ndarray' objects}
            1    0.000    0.000    0.002    0.002 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/fromnumeric.py:1379(sum)
            1    0.002    0.002    0.002    0.002 {method 'sum' of 'numpy.ndarray' objects}
            1    0.000    0.000    0.000    0.000 {isinstance}
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}

Overall think really hard before you write a for loop or a list comprehension.

.. _IDL: http://www.ittvis.com/language/en-us/productsservices/idl.aspx
.. _Matlab: http://www.mathworks.com/products/matlab/


Zip
===
the zip_ function is a great thing but it is really slow, if you find yourself
using it then you probably need to reexamine the algorithm that you are using.

This example generate evenly distributed N points on the unit sphere centered at
(0,0,0) using the "Golden Spiral" method.

The original code::

    import numpy as np
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

Profiling this with cProfile one can see a lot of time in zip()::

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

So lets try and do a few simple rewrites to make this faster.  Using numpy.vstack
is the first one that came to mind.  The change here is to replace building up
the array from the elements made by zip to just appending the data we already have
to an array that we already have::

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
        x0_y0_z0 = np.vstack((x0, y0, z0)).transpose()
        return (x0_y0_z0)

Profiling this with cProfile one can see that this is now fast enough for me,
no more work to do.  We picked up a 48x speed increase, I'm sure this can still
be made better and let the spacepy team know if you rewrite it and it will be
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

For additions or fixes to this page contact Brian Larsen at Los Alamos.
