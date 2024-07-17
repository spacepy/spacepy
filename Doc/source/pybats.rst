########################################
PyBats - SWMF & BATS-R-US Analysis Tools
########################################

.. contents:: Table of Contents
    :depth: 2
    :local:

.. currentmodule:: spacepy.pybats
    
PyBats!  An open source Python-based interface for reading, manipulating,
and visualizing BATS-R-US and SWMF output.
For more information on the SWMF, please visit the
`Center for Space Environment Modeling <http://csem.engin.umich.edu>`_.

See also the `full API documentation <spacepy.pybats>`.

*************
PyBats Module
*************

pybats
======

Introduction
------------
At its most fundamental level, PyBats provides access to output files written
by the Space Weather Modeling Framework and the codes contained within.
The key task performed by PyBats is loading simulation data into a Spacepy
data model object so that the user can move on to the important tasks of
analyzing and visualizing the values.  The secondary goal of PyBats is to make
common tasks performed with these data as easy as possible.  The result is that
most SWMF output can be opened and visualized using only a few lines of code.
Many complicated tasks, such as field line integration, is included as well.


Organization
------------
Many output files in the SWMF share a common format.  Objects to handle broad
formats like these are found in the base module.  The base module has classes
to handle SWMF *input* files, as well.

The rest of the classes are organized by code, i.e. classes and functions
specifically relevant to BATS-R-US can be found in
:mod:`spacepy.pybats.bats`.  Whenever a certain code included in the SWMF
requires a independent class or a subclass from the PyBats base module, it
will receive its own submodule.


Conventions and Prefixes
------------------------
Nearly every class in PyBats inherits from :class:`spacepy.datamodel.SpaceData`
so it is important for users to understand how to employ and explore SpaceData
objects.  There are a few exceptions, so always pay close attention to the
docstrings and examples.  Legacy code that does not adhere to this pattern is
slowly being brought up-to-date with each release.

Visualization methods have two prefixes: *plot_* and *add_*.  Whenever a method
begins with *plot_*, a quick-look product will be created that is not highly-
configurable.  These methods are meant to yeild either simple
diagnostic plots or static, often-used products.  There are few methods that
use this prefix.  The prefix *add_* is always followed by *plot_type*; it
indicates a plotting method that is highly configurable and meant to be
combined with other *add_*-like methods and matplotlib commands.

Common calculations, such as calculating Alfven wave speeds of MHD results,
are strewn about PyBats' classes.  They are always given the method prefix
*calc_*, i.e. *calc_alfven*.  Methods called *calc_all* will search for all
class methods with the *calc_* prefix and call them.


Submodules
----------
There are submodules for most models included within the SWMF.  The classes
and methods contained within are code-specific, yielding power and
convenience at the cost of flexibility.  A few of the submodules are helper
modules- they are not code specific, but rather provide functionality not
related to an SWMF-included code. See `pybats API <spacepy.pybats>`.


Top-Level Classes & Functions
-----------------------------
Top-level PyBats classes handle common-format input and output from the SWMF
and are very flexible.  However, they do little beyond open files for the user.

There are several functions found in the top-level module.  These are mostly
convenience functions for customizing plots.


class IdlFile
=============
The IdlFile class of pybats.


Introduction
------------
An object class that reads/parses an IDL-formatted output file from the
SWMF and places it into a :class:`spacepy.pybats.PbData` object.

Usage:
>>>data = spacepy.pybats.IdlFile('binary_file.out')

See :class:`spacepy.pybats.PbData` for information on how to explore
data contained within the returned object.

This class serves as a parent class to SWMF component-specific derivative
classes that do more preprocessing of the data before returning the
object.  Hence, using this class to read binary files is typically not
the most efficient way to proceed.  Look for a PyBats sub module that suits
your specific needs, or use this base object to write your own.

Multi-Frame Files
-----------------
Typically, Idl-formatted data has a single *frame*, or a single snapshot
worth of data (a ``*.out`` file).  However, it is possible to (externally)
combine many of these files together such that a time series of data frames
are contained within a single file (``*.outs`` files). This class can read
these files, but only one data frame can be made available at a time.  This
prevents very large memory requirements.

These files are opened and handled similarly to regular IDL-formatted
files  with some important differences.  Upon instantiation, the user
may select which data frame to open with the *iframe* kwarg.  The default
action is to open the first frame.  The user may learn more about the
number of frames within the file and their associated epoch/iteration
information by examining the top level *attrs* (see below.)
The user may switch to any arbitrary frame using the ``switch_frame(iframe)``
object method.  This will load the relevant data from disk into the
object, overwriting the previous contents.

If the user has created any new variables using the data from a certain
data frame, those new values will not be updated automatically.  An
exception is for any *self.calc_* type methods that are set to update
automatically.

Top Level Attributes
--------------------
Critical file information can be found via the top-level object attributes
(e.g., accessing the dictionary ``self.attrs``).  All IdlFile objects and
child classes begin with at least these attributes:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Attribute Name
      - Description
    * - file
      - Path/Name of file represented by object
    * - iter/time/runtime
      - Iteration/datetime/runtime of the current frame
    * - iters/time/runtimes
      - Lists of all iterations/times of each data frame
    * - range
      - The range of iterations/epochs covered in file
    * - ndim
      - Number of spatial dimensions covered by the data
    * - nframe
      - The total number of data frame within the file
    * - iframe 
      - The current frame loaded (zero-based)
    * - format
      - The format of the file, either binary or ascii
    * - header
      - The raw string header of the file

Notes
-----
PyBats assumes little endian byte ordering because
this is what most machines use.  However, there is an autodetect feature
such that, if PyBats doesn't make sense of the first read (a record length
entry, or RecLen), it will proceed using big endian ordering.  If this
doesn't work, the error will manifest itself through the "struct" package
as an "unpack requires a string of argument length 'X'".

.. versionchanged:: 0.5.0

    Unstructured data are presented as in the files. When reading
    3D magnetosphere files, this preserves the 3D block structure,
    as required for the BATSRUS interpolator in the `Kamodo
    Heliophysics model readers package
    <https://github.com/nasa/kamodo>`_. Before 0.5.0, binary
    unstructured data were sorted in an attempt to put nearby
    positions close to each other in the data arrays. This sorting
    was nondeterministic and has been removed; see
    :meth:`~spacepy.pybats.bats.Bats2d.extract` and
    :class:`~spacepy.pybats.qotree.QTree` for processing
    adjacent cells. (ASCII data were never sorted.)


class ImfInput
==============
A class to read, write, manipulate, and visualize solar wind upstream
input files for SWMF simulations.  More about such files can be found
in the SWMF/BATS-R-US documentation for the \\#SOLARWINDFILE command.

Creating an :class:`ImfInput` object is simple:

>>> from spacepy import pybats
>>> obj=pybats.ImfInput(filename='test.dat', load=True)

Upon instantiation, if *filename* is a valid file AND kwarg *load* is set
to boolean True, the contents of *filename* are loaded into the object
and no other work needs to be done.

If *filename* is False or *load* is False, a blank :class:`ImfInput` file
is created for the user to manipulate.
The user can set the time array and the
associated data values (see *obj.attrs['var']* for a list) to any values
desired and use the method *obj.write()* to dump the contents to an SWMF
formatted input file.  See the documentation for the write method for
more details.

Like most :mod:`~spacepy.pybats` objects, you may interact with
:class:`ImfInput` objects as if they were specialized dictionaries.
Access data like so:

>>> obj.keys()
['bx', 'by', 'bz', 'vx', 'vy', 'vz', 'n', 't']
>>> density=obj['n']

Adding new data entries is equally simple so long as you have the values
and the name for the values:

>>> import numpy as np
>>> u = np.sqrt(obj['ux']**2 + obj['uy']**2 + obj['uz']**2)
>>> obj['u']=u

If new data entries are added as :class:`~spacepy.datamodel.dmarray`
objects, the ``label`` and ``units`` attributes can be set to enhance plotting.

>>> from spacepy import datamodel
>>> u = np.sqrt(obj['ux']**2 + obj['uy']**2 + obj['uz']**2)
>>> obj['u']= datamodel.dmarray(u, {'units': '$km/s$', 'label': '|U|'})

Concerning Variable Naming & Order
----------------------------------
By default, IMF files contain the following variables in this order:

Year, month, day, hour, minute, second, millisecond, bx, by, bz,
vx, vy, vz, n, t

If the variable order changes, or if new state variables are included
(e.g., species-specific densities for multi-ion simulations), the
#VAR entry must be included in the solar wind file.  While the SWMF
documentation refers to density and temperature values as having names
'dens' or 'temp', **only 'n', 't', or MHD state variable names as defined
in the BATS-R-US equation file are accepted.**.  In multi-ion simulations,
'n' is the total number density; all other densities must sum to 'n'.

To illustrate, consider this example of converting a single fluid input
file to a multi-fluid input file where density is split evenly across
two ion species: a solar wind fluid and an ionosphere fluid that has
minimal density upstream of the Earth. Each fluid has an explicit state
variable name defined in the BATS-R-US equation module that configures the
multi-fluid configuration.

>>> import numpy as np
>>> from spacepy import pybats, datamodel
>>> imf = pybats.ImfInput('tests/data/pybats_test/imf_single.dat')
>>> imf['IonoRho'] = datamodel.dmarray(np.zeros(imf['n'].shape) + 0.01,
... {'units': '$cm^{-3}$', 'label': '$\\rho_{Iono}$'})
>>> imf['SwRho'] = imf['n'] - imf['IonoRho']
>>> imf['SwRho'].attrs['label'] = '$\\rho_{Sw}$'
>>> imf.quicklook( ['bz', ['n', 'SwRho', 'IonoRho'], 'pram'])

=========== ============================================================
Kwarg       Description
----------- ------------------------------------------------------------
filename    Set the input/output file name.
load        Read file upon instantiation?  Defaults to **True**
npoints     For empty data sets, sets number of points (default is 0)
=========== ============================================================

.. versionchanged:: 0.5.0
    Default variable names for temperature and density are now 't' and 'n'.


**************
Bats Submodule
**************

Bats2d
======
A child class of :class:`~spacepy.pybats.IdlFile` tailored to 2D BATS-R-US output.

Calculations
------------
New values can be added via the addition of new keys.  For example,
a user could add radial distance to an equatorial Bats2d object as follows:

>>> import numpy as np
>>> from spacepy.pybats import bats
>>> mhd = bats.Bats2d('tests/data/pybats_test/z=0_sine.out')
>>> mhd['rad'] = np.sqrt( mhd['x']**2 + mhd['y']**2 )

Note, however, that if the user switches the data frame in a .outs file
to access data from a different epoch, these values will need to be
updated.

An exception to this is built-in ``calc_*`` methods, which perform common
MHD/fluid dynamic calculations (i.e., Alfven wave speed, vorticity, and
more.)  These values are updated when the data frame is switched (see the
``switch_frame`` method).

For example, to calculate number density and the fraction of the total
density for each species in a multi-fluid simulation,

>>> mhd = bats.Bats2d('tests/data/pybats_test/cut_multispecies.out')
>>> mhd.calc_ndens()

In an interactive session, users can tab-complete to explore the possible
built-in calculation methods.

Plotting
--------
While users can employ Matplotlib to plot values, a set of built-in
methods are available to expedite plotting.  These are the
``add_<plot type>`` methods.  These methods always have the following
keyword arguments that allow users to optionally build more complicated
plots: *target* and *loc*.  The *target* kwarg tells the plotting method
where to place the plot and can either be a Matplotlib figure or axes
object.  If it's an axes object, *loc* sets the subplot location using
the typical matplotlib syntax (e.g., ``loc=121``).  The default behavior is
to create a new figure and axes object.

This approach allows a user to over-plot contours, field lines, and
other plot artists as well as combine different subplots onto a single
figure.  Continuing with our example above, let's plot the grid layout
for our file as well as equatorial pressure and flow streamlines:

>>> import matplotlib.pyplot as plt
>>> fig = plt.Figure(figsize=(8,6))
>>> mhd.add_grid_plot(target=fig, loc=121)
>>> mhd.add_contour('x','y','p', target=fig, loc=122)
>>> mhd.add_stream_scatter('ux', 'uy', target=fig, loc=122)

Useful plotting methods include the following:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Plot Method
      - Description
    * - add_grid_plot
      - Create a quick-look diagram of the grid layout
    * - add_contour
      - Create a contour plot of a given variable
    * - add_pcolor
      - Add a p-color (no-interpolation contour) plot
    * - add_stream_scatter
      - Scatter stream traces (any vector field)
    * - add_b_magsphere
      - Add magnetic field lines for X-Z plane cuts
    * - add_planet
      - Add a simple black/white planet at the origin
    * - add_body
      - Add an inner boundary at the origin

Extracting and Stream Tracing
-----------------------------
Extracting values via interpolation to arbitrary points and creating
stream traces through any vector field (e.g., velocity or magnetic field)
are aided via the use of the following object methods:

.. list-table::
    :widths: 25 75
    :header-rows: 1

    * - Method
      - Description
    * - extract
      - Interpolate to arbitrary points and extract values
    * - get_stream
      - Integrate stream lines through vector fields


Be sure to read the docstring information of :class:`~spacepy.pybats.IdlFile` to
see how to handle multi-frame files (.outs) and for a list of critical
attributes.

.. versionchanged:: 0.5.0

    Unstructured data are now presented as in the files. See
    :class:`~spacepy.pybats.IdlFile` for details.