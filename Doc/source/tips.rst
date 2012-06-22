===============================
SpacePy Python Programming Tips
===============================

One often hears that interpreted languages are too slow for whatever task someone
needs to do.  In many cases this belief is unfounded.  As the time spent
programming and debugging in an interpreted language is of far less than for a compiled 
language, the programmer has more time to identify bottlenecks in the code and optimize it
where necessary.  This page is dedicated to that idea, providing examples of code speedup 
and best practices.

One often neglected way to speed up development time is to use established libraries, and 
the time spent finding existing code that does what you want can be more productive than
trying to write and optimize your own algorithms. We recommend exploring the SpacePy 
documentation, as well as taking the time to learn some of the vast functionality already in
numpy_ and the Python standard library.

    * `Basic examples`_
    * `Lists, for loops, and arrays`_
    * `Zip`_
    * `External links`_


Basic examples
==============
Though there are some similarities, Python does not look like (or work like) Matlab or IDL. 
As (most of us) are, or have been, Matlab or IDL programmers, we have to rethink how we do 
things -- what is efficient in one language may not be the most efficient in another.  
One truth that Python shares with these other languages, however, is that if you are using 
a for loop there is likely to be a faster way...

Take the simple case of a large data array where you want to add one to each element.
Here wa show four of the possible ways to do this, and by using a profiling tool, we can 
show that the speeds of the different methods can vary substantially.

Create the data

>>> data = range(10000000)

The most C-like way is just a straight up for loop

>>> for i in range(len(data)):
>>>     data[i] += 1

This is 6 function calls in 2.590 CPU seconds (from :py:mod:`cProfile`)

The next, more Pythonic, way is with a list comprehension

>>> data = [val+1 for val in data]

This is 4 function calls in 1.643 CPU seconds (~1.6x)

Next we introduce numpy_ and change our list to an array then add one

>>> data = np.asarray(data)
>>> data += 1

This is 6 function calls in 1.959 CPU seconds (~1.3x), better than the for loop but worse
than the list comprehension

Next we do this the `right` way and just create it in numpy_ and never leave

>>> data = np.arange(10000000)
>>> data += 1

this is 4 function calls in 0.049 CPU seconds (~53x).

While this is a really simple example it shows the basic premise, if you need to work 
with numpy_, start in numpy_ and stay in numpy_. This will usually be true for 
array-based manipulations.

If in doubt, and speed is not a major consideration, use the most human-readable form.
This will make your code more maintainable and encourage its use by others.



Lists, for loops, and arrays
============================
This example teaches the lesson that most advanced IDL_ or Matlab_ programmers already
know; do everything in arrays and never use a for loop if there is another choice. The 
language has optimized array manipulation and it is unlikely that you will find a faster
way with your own code.

The following bit of code takes in a series of coordinates, computes their displacement, and drops
the largest 100 of them.

This is how the code started out, Shell_x0_y0_z0 is an Nx3 numpy array,
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

A cProfile of this yields a lot of time spent just in the function itself; this
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

Simply moving the addition outside the for-loop gives a factor of 2.3 speedup.  
We believe that the difference arising from moving the addition lets
numpy (which works primarily in C) operate once only. This massively reduces the calling overhead 
as array operations are done as for loops in C, and not in element-wise in python::

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

The answer is exactly the same and comparing to the initial version of this code we have managed a speedup of 382x::

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

In summary, when working on arrays it's worth taking the time to think about whether you can get the results you need without for-loops or list comprehensions. The small amount of development time will likely be recouped very quickly.

.. _IDL: http://www.ittvis.com/language/en-us/productsservices/idl.aspx
.. _Matlab: http://www.mathworks.com/products/matlab/


Zip
===
The :py:func:`zip` function is extremely useful, but it is really slow. If you find yourself
using it on large amounts of data then significant time-savings might be achieved by re-writing your code
to make the :py:func:`zip` operation unnecessary. A good alternative, if you do need the functionality 
of :py:func:`zip`, is in :py:func:`itertools.izip`. This is far more efficient as it builds an interator.

This example generates N points, evenly distributed on the unit sphere centered at
(0,0,0) using the "Golden Spiral" method.

The original code::

    import numpy as np
    def PointsOnSphere(N):
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

Profiling this with :py:mod:`cProfile` shows that a lot of time is spent in :py:func:`zip`::

    Tue Jun 14 09:54:41 2011    PointsOnSphere.prof

             9 function calls in 8.132 seconds

       Ordered by: cumulative time

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.010    0.010    8.132    8.132 <string>:1(<module>)
            1    0.470    0.470    8.122    8.122 test1.py:192(PointsOnSphere)
            4    6.993    1.748    6.993    1.748 {numpy.core.multiarray.array}
            1    0.654    0.654    0.654    0.654 {zip}
            1    0.005    0.005    0.005    0.005 {numpy.core.multiarray.arange}
            1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}

So lets try and do a few simple rewrites to make this faster.  Using numpy.vstack
is the first one that came to mind.  The change here is to replace building up
the array from the elements made by :py:func:`zip` to just concatenating the arrays we already have::

    def PointsOnSphere(N):
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

Profiling this with :py:mod:`cProfile` one can see that this is now fast enough for me,
no more work to do.  We picked up a 48x speed increase, I'm sure this can still
be made better and let the SpacePy team know if you rewrite it and it will be
included::

    Tue Jun 14 09:57:41 2011    PointsOnSphere.prof

             32 function calls in 0.168 seconds

       Ordered by: cumulative time
       List reduced from 13 to 10 due to restriction <10>

       ncalls  tottime  percall  cumtime  percall filename:lineno(function)
            1    0.010    0.010    0.168    0.168 <string>:1(<module>)
            1    0.123    0.123    0.159    0.159 test1.py:217(PointsOnSphere)
            1    0.000    0.000    0.034    0.034 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/shape_base.py:177(vstack)
            1    0.034    0.034    0.034    0.034 {numpy.core.multiarray.concatenate}
            1    0.002    0.002    0.002    0.002 {numpy.core.multiarray.arange}
            1    0.000    0.000    0.000    0.000 {map}
            3    0.000    0.000    0.000    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/shape_base.py:58(atleast_2d)
            6    0.000    0.000    0.000    0.000 {numpy.core.multiarray.array}
            3    0.000    0.000    0.000    0.000 /opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/numeric.py:237(asanyarray)
            1    0.000    0.000    0.000    0.000 {method 'transpose' of 'numpy.ndarray' objects}



External links
==============
To learn about NumPy from a MatLab user's perspective
    * `NumPy for MatLab users`_
    * `Mathesaurus`_

And if you're coming from an IDL, or R, background
    * `Mathesaurus`_

Here is a collection of links that serve as a decent reference for Python and speed
    * `PythonSpeed PerformanceTips`_
    * `scipy array tip sheet`_
    * `Python Tips, Tricks, and Hacks`_



.. _numpy: http://docs.scipy.org/doc/numpy/reference/
.. _`NumPy for MatLab users`: http://www.scipy.org/NumPy_for_Matlab_Users
.. _`Mathesaurus`: http://mathesaurus.sourceforge.net/
.. _`PythonSpeed PerformanceTips`: http://wiki.python.org/moin/PythonSpeed/PerformanceTips
.. _`scipy array tip sheet`: http://pages.physics.cornell.edu/~myers/teaching/ComputationalMethods/python/arrays.html
.. _`Python Tips, Tricks, and Hacks`: http://www.siafoo.net/article/52

--------------------------

:Release: |version|
:Doc generation date: |today|

For additions or fixes to this page contact the SpacePy team at Los Alamos.
