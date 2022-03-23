*************************************
SpacePy - A Quick Start Documentation
*************************************


The SpacePy Team
(Steve Morley, Josef Koller, Dan Welling, Brian Larsen, Jon Niehof,
Mike Henderson)


Installation
============

See :doc:`install`.

Toolbox - A Box Full of Tools
=============================

Contains tools that don't fit anywhere else but are, in general, quite
useful. The following functions are a selection of those implemented:

    * :func:`~spacepy.toolbox.windowMean`: windowing mean with variable window size and overlap
    * :func:`~spacepy.toolbox.dictree`: pretty prints the contents of dictionaries (recursively)
    * :func:`~spacepy.toolbox.loadpickle`: single line convenience routine for loading Python pickles
    * :func:`~spacepy.toolbox.savepickle`: same as loadpickle, but for saving
    * :func:`~spacepy.toolbox.update`: updates the OMNI database and the leap seconds database (internet connection required)
    * :func:`~spacepy.toolbox.tOverlap`: find interval of overlap between two time series
    * :func:`~spacepy.toolbox.tCommon`: find times common to two time series
    * :func:`~spacepy.toolbox.binHisto`: calculate number of bins for a histogram
    * :func:`~spacepy.toolbox.medAbsDev`: find the median absolute deviation of a data series
    * :func:`~spacepy.toolbox.normalize`: normalize a data series
    * :func:`~spacepy.toolbox.feq`: floating point equals

Import this module as::

>>> import spacepy.toolbox as tb

Examples:

>>> import spacepy.toolbox as tb
>>> a = {'entry1':'val1', 'entry2':2, 'recurse1':{'one':1, 'two':2}}
>>> tb.dictree(a)
+
|____entry1
|____entry2
|____recurse1
     |____one
     |____two
>>> import numpy as np
>>> dat = np.random.random_sample(100)
>>> tb.binHisto(dat)
(0.19151723370512266, 5.0)




Time and Coordinate Transformations
===================================

Import the modules as::

>>> import spacepy.time as spt
>>> import spacepy.coords as spc


Ticktock Class
--------------

The Ticktock class provides a number of time conversion routines and is
implemented as a container class built on the functionality of the Python
datetime module. The following time coordinates are provided

    * UTC: Coordinated Universal Time implemented as a :class:`datetime.datetime`
    * ISO: standard ISO 8601 format like ``2002-10-25T14:33:59``
    * TAI: International Atomic Time in units of seconds since Jan 1, 1958 (midnight) and includes leap seconds, i.e. every second has the same length
    * JD:  Julian Day
    * MJD: Modified Julian Day
    * UNX: UNIX time in seconds since Jan 1, 1970
    * RDT: Rata Die Time (Gregorian Ordinal Time) in days since Jan 1, 1 AD midnight
    * CDF: CDF Epoch time in milliseconds since Jan 1, year 0
    * DOY: Day of Year including fractions
    * leaps: Leap seconds according to ftp://maia.usno.navy.mil/ser7/tai-utc.dat

To access these time coordinates, you'll create an instance of a
Ticktock class, e.g.::

>>> t = spt.Ticktock('2002-10-25T12:30:00', 'ISO')

Instead of ISO you may use any of the formats listed above. You can also
use numpy arrays or lists of time points. ``t`` has now the class
attributes::

>>> t.dtype = 'ISO'
>>> t.data = '2002-10-25T12:30:00'

FYI ``t.UTC`` is added automatically.

If you want to convert/add a class attribute from the list above,
simply type e.g.::

>>> t.RTD

You can replace RTD with any from the list above.

You can find out how many leap seconds were used by issuing the command::

>>> t.getleapsecs()


Timedelta Class
---------------

You can add/subtract time from a Ticktock class instance by using an
instance of :class:`datetime.timedelta`::

>>> dt = datetime.timedelta(days=2.3)

Then you can add by e.g.::

>>> t+dt


Coords Class
------------

The spatial coordinate class includes the following coordinate systems in
Cartesian and spherical forms.

    * GZD:  (altitude, latitude, longitude) in km, deg, deg
    * GEO: cartesian, Re
    * GSM: cartesian, Re
    * GSE: cartesian, Re
    * SM: cartesian, Re
    * GEI: cartesian, Re
    * MAG: cartesian, Re
    * SPH: same as GEO but in spherical
    * RLL: radial distance, latitude, longitude, Re, deg, deg.

Create a Coords instance with spherical='sph' or cartesian='car'
coordinates::

>>> spaco = spc.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')

This will let you request, for example, all y-coordinates by ``spaco.y``
or if given in spherical coordinates by ``spaco.lati``. One can transform
the coordinates by ``newcoord = spaco.convert('GSM', 'sph')``.
This will return GSM coordinates in a spherical system. Since GSM
coordinates depend on time, you'll have to add first a Ticktock
vector with the name ``ticks`` like ``spaco.ticks = spt.Ticktock(['2002-02-02T12:00:00',
'2002-02-02T12:00:00'], 'ISO')``

Unit conversion will be implemented in the future.


The radbelt Module
==================

The radiation belt module currently includes a simple radial
diffusion code as a class. Import the module and instatiate a radbelt object::

>>> import spacepy.radbelt as sprb
>>> rb = sprb.RBmodel()

Add a time grid for a particular period that you are interested in::

>>> rb.setup_ticks('2002-02-01T00:00:00', '2002-02-10T00:00:00', 0.25)

This will automatically lookup required geomagnetic/solar wind conditions
for that period. Run the diffusion solver for that setup and plot the
results::

>>> rb.evolve()
>>> rb.plot()


The Data Assimilation Module
============================

This module includes data assimilation capabilities, through the
assimilation class. The class assimilates data for the radiation belt model
using the Ensemble Kalman Filter. The algorithm used is the SVD method
presented by Evensen in 2003 (Evensen, G., Ocean dynamics, 53, pp.343--367,
2003). To compensate for model errors, three inflation algorithms are
implemented. The inflation methodology is specified by the inflation
argument, where the options are the following:

   * inflation = 0: Add model error (perturbation for the ensemble) around
     model state values only where observations are available (DEFAULT).

   * inflation = 1: Add model error (perturbation for the ensemble) around
     observation values only where observations are available.

   * inflation = 2: Inflate around ensemble average for EnKF.

Prior to assimilation, a set of data values has to be specified by setting the
start and end dates, and time step, using the ``setup_ticks`` function of the
radiation belt model::

>>> import spacepy
>>> import datetime
>>> from spacepy import radbelt

>>> start = datetime.datetime(2002,10,23)
>>> end = datetime.datetime(2002,11,4)
>>> delta = datetime.timedelta(hours=0.5)
>>> rmod.setup_ticks(start, end, delta, dtype='UTC')

Once the dates and time step are specified, the data is added using the
``add_PSD`` function (NOTE: This requires a database available from the SpacePy team)::

>>> rmod.add_PSD()

The observations are averaged over the time windows, whose interval is give by
the time step. Once the dates and data are set, the assimilation is performed
using the ``assimilate`` function::

>>> rmod.assimilate(inflation=1)

This function will add the PSDa values, which are the analysis state of
the radiation belt using the observations within the dates. To plot the
analysis simply use the ``plot`` function::

>>> rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

Additionally, to create a summary plot of the observations use the ``plot_obs``
function within the radbelt module. For reference, the last closed drift shell,
Dst, and Kp are all included. These can be disabled individually using the
corresponding Boolean kwargs.

The clims kwarg can be used to manually set the color bar range.  To use, set
it equal to a two-element list containing minimum and maximum log :sub:`10` value to
plot.  Default action is to use [0,10] as the log :sub:`10` of the color range.  This
is good enough for most applications. The title of the top most plot defaults
to 'Summary Plot' but can be customized using the title kwarg.

The figure object and all three axis objects (PSD axis, Dst axis, and Kp axis)
are all returned to allow the user to further customize the plots as necessary.
If any of the plots are excluded, None is returned in their stead.

Example::

>>> rmod.plot_obs(clims=[-10,-6],Lmax=False,Kp=False,Dst=False,title='Observations Plot')

This command would create the summary plot with a color bar range of 10 :sup:`-10`
to 10 :sup:`-6`.  The Lmax line, Kp and Dst values would be excluded.  The title of
the topmost plot (phase space density) would be set to 'Observations Plot'.


OMNI Module
===========

The OMNI database is an hourly resolution, multi-source data set
with coverage from November 1963; higher temporal resolution versions of
the OMNI database exist, but with coverage from 1995. The primary data are
near-Earth solar wind, magnetic field and plasma parameters. However, a
number of modern magnetic field models require derived input parameters,
and Qin and Denton (2007) have used the publicly-available OMNI database to provide
a modified version of this database containing all parameters necessary
for these magnetic field models. These data are available through ViRBO  - the Virtual
Radiation Belt Observatory.

In SpacePy this data is made available, at 1-hourly resolution, on request 
on first import; if not downloaded when SpacePy is first used then any 
attempt to import the omni module will
ask the user whether they wish to download the data. Should the user
require the latest data, the toolbox.update function can
be used to fetch the latest files from ViRBO.

The following example fetches the OMNI data for the storms of
October and November, 2003.::

>>> import spacepy.time as spt
>>> import spacepy.omni as om
>>> import datetime as dt
>>> st = dt.datetime(2003,10,20)
>>> en = dt.datetime(2003,12,5)
>>> delta = dt.timedelta(days=1)
>>> ticks = spt.tickrange(st, en, delta, 'UTC')
>>> data = om.get_omni(ticks)

*data* is a dictionary containing all the OMNI data, by variable, for the timestamps
contained within the ``Ticktock`` object *ticks*. Now it is simple to plot Dst values
for instance::

>>> import pyplot as p
>>> p.plot(ticks.eDOY, data['Dst'])


The irbempy Module
==================

ONERA (Office National d'Etudes et Recherches Aerospatiales) initiated a
well-known FORTRAN library that provides routines to compute magnetic
coordinates for any location in the Earth's magnetic field, to perform
coordinate conversions, to compute magnetic field vectors in geospace for
a number of external field models, and to propagate satellite orbits in
time. Older versions of this library were called ONERA-DESP-LIB. Recently
the library has changed its name to IRBEM-LIB and is maintained by a number
of different institutions.

A number of key routines in IRBEM-LIB have been made available through the
module :mod:`~spacepy.irbempy`. Current functionality includes calls to calculate the local
magnetic field vectors at any point in geospace, calculation of the magnetic
mirror point for a particle of a given pitch angle (the angle between a
particle's velocity vector and the magnetic field line that it immediately
orbits such that a pitch angle of 90 degrees signifies gyration perpendicular
to the local field) anywhere in geospace, and calculation of electron drift
shells in the inner magnetosphere.::

>>> import spacepy.time as spt
>>> import spacepy.coordinates as spc
>>> import spacepy.irbempy as ib
>>> t = spt.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
>>> y = spc.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')
>>> ib.get_Bfield(t,y)
>>> # {'Blocal': array([  976.42565251,  3396.25991675]),
>>> #    'Bvec': array([[ -5.01738885e-01,  -1.65104338e+02,   9.62365503e+02], [  3.33497974e+02,  -5.42111173e+02,   3.33608693e+03]])}

One can also calculate the drift shell L* for a 90 degree pitch angle value by using::

>>> ib.get_Lstar(t,y, [90])
>>> # {'Bmin': array([  975.59122652,  3388.2476667 ]),
>>> #  'Bmirr': array([[  976.42565251], [ 3396.25991675]]),
>>> #  'Lm': array([[ 3.13508015], [ 2.07013638]]),
>>> #  'Lstar': array([[ 2.86958324], [ 1.95259007]]),
>>> #  'MLT': array([ 11.97222034,  12.13378624]),
>>> #  'Xj': array([[ 0.00081949], [ 0.00270321]])}

Other function wrapped with the IRBEM library include:

* :func:`~spacepy.irbempy.find_Bmirror`
* :func:`~spacepy.irbempy.find_magequator`
* :func:`~spacepy.irbempy.coord_trans`


pyCDF - Python Access to NASA CDF Library
=========================================

pycdf provides a "pythonic" interface to the NASA CDF library. It requires
that the NASA CDF C-library is properly installed.
The module can then be imported, e.g.::

>>> import spacepy.pycdf as cdf

To open and close a CDF file, we use the :class:`~pycdf.CDF` class::

>>> cdf_file = cdf.CDF('filename.cdf')
>>> cdf_file.close()

CDF files, like standard Python files, act as context managers::

>>> with cdf.CDF('filename.cdf') as cdf_file:
>>>     #do brilliant things with cdf_file
>>> #cdf_file is automatically closed here

CDF files act as Python dictionaries, holding CDF variables keyed
by the variable name::

>>> var_names = keys(cdf_file) #list of all variables
>>> for var_name in cdf_file:
>>>     print(len(cdf_file[var_name])) #number of records in each variable
>>> #list comprehensions work, too
>>> lengths = [len(cdf_file[var_name]) for var_name in cdf_file]

Each CDF variable acts like a numpy array, where the first dimension is the
record number. Multidimensional CDF variables can be subscripted using 
numpy's multidimensional slice notation. Many common list operations are also
implemented, where each record acts as one element of the list and can be
independently deleted, inserted, etc. Creating a Python :class:`~pycdf.Var`
object does not read the data from disc; data are only read as they are
accessed::

>>> epoch = cdf_file['Epoch'] #Python object created, nothing read from disc
>>> epoch[0] #time of first record in CDF (datetime object)
>>> a = epoch[...] #copy all times to list a
>>> a = epoch[-5:] #copy last five times to list a
>>> b_gse = cdf_file['B_GSE'] #B_GSE is a 1D, three-element array
>>> bz = b_gse[0,2] #Z component of first record
>>> bx = b_gse[:,0] #copy X component of all records to bx
>>> bx = cdf_file['B_GSE'][:,0] #same as above


The datamodel Module
====================

The SpacePy datamodel module implements classes that are designed to make implementing a standard
data model easy. The concepts are very similar to those used in standards like HDF5, netCDF and
NASA CDF.

The basic container type is analogous to a folder (on a filesystem; HDF5 calls this a
group): Here we implement this as a dictionary-like object, a :class:`datamodel.SpaceData` object, which
also carries attributes. These attributes can be considered to be global, i.e. relevant for the
entire folder. The next container type is for storing data and is based on a numpy array, this
class is :class:`datamodel.dmarray` and also carries attributes. The dmarray class is analogous to an
HDF5 dataset.


Guide for NASA CDF users
------------------------

By definition, a NASA CDF only has a single 'layer'. That is, a CDF contains a series of records
(stored variables of various types) and a set of attributes that are either global or local in
scope. Thus to use SpacePy's datamodel to capture the functionality of CDF the two basic data types
are all that is required, and the main constraint is that datamodel.SpaceData objects cannot be
nested (more on this later, if conversion from a nested datamodel to a flat datamodel is required).

This is best illustrated with an example. Imagine representing some satellite data within a CDF --
the global attributes might be the mission name and the instrument PI, the variables might be the
instrument counts [n-dimensional array], timestamps[1-dimensional array and an orbit number [scalar].
Each variable will have one attribute (for this example).

>>> import spacepy.datamodel as dm
>>> mydata = dm.SpaceData(attrs={'MissionName': 'BigSat1'})
>>> mydata['Counts'] = dm.dmarray([[42, 69, 77], [100, 200, 250]], attrs={'Units': 'cnts/s'})
>>> mydata['Epoch'] = dm.dmarray([1, 2, 3], attrs={'units': 'minutes'})
>>> mydata['OrbitNumber'] = dm.dmarray(16, attrs={'StartsFrom': 1})
>>> mydata.attrs['PI'] 'Prof. Big Shot'

This has now populated a structure that can map directly to a NASA CDF. To visualize our datamodel,
we can use the :meth:`~spacepy.datamodel.SpaceData.tree` method, which is equivalent to :func:`toolbox.dictree`
(which works for any dictionary-like object, including PyCDF file objects).

>>> mydata.tree(attrs=True)
+
:|____MissionName
:|____PI
|____Counts
    :|____Units
|____Epoch
    :|____units
|____OrbitNumber
    :|____StartsFrom
>>> import spacepy.toolbox as tb 
>>> tb.dictree(mydata, attrs=True)
+
:|____MissionName
:|____PI
|____Counts
    :|____Units
|____Epoch
    :|____units
|____OrbitNumber
    :|____StartsFrom


Attributes are denoted by a leading colon. The global attributes are those in the base level,
and the local attributes are attached to each variable.

If we have data that has nested 'folders', allowed by HDF5 but not by NASA CDF, then how can this be
represented such that the data structure can be mapped directly to a NASA CDF? The data will need to
be flattened so that it is single layered. Let us now store some ephemerides in our data structure:

>>> mydata['Ephemeris'] = dm.SpaceData()
>>> mydata['Ephemeris']['GSM'] = dm.dmarray([[1,3,3], [1.2,4,2.5], [1.4,5,1.9]])
>>> tb.dictree(mydata, attrs=True)
+
:|____MissionName
:|____PI
|____Counts
    :|____Units
|____Ephemeris
    |____GSM
|____Epoch
    :|____units
|____OrbitNumber
    :|____StartsFrom

Nested dictionary-like objects is not uncommon in Python (and can be exceptionally useful for representing
data, so to make this compatible with NASA CDF we call the :meth:`~spacepy.datamodel.SpaceData.flatten` method .

>>> mydata.flatten()
>>> tb.dictree(mydata, attrs=True)
+
:|____MissionName
:|____PI
|____Counts
    :|____Units
|____Ephemeris<--GSM
|____Epoch
    :|____units
|____OrbitNumber
    :|____StartsFrom

Note that the nested SpaceData has been moved to a variable with a new name reflecting its origin. The
data structure is now flat again and can be mapped directly to NASA CDF.


Converters to/from datamodel
----------------------------

Currently converters exist to read HDF5 and NASA CDF files directly to a SpacePy datamodel. This capability 
also exists for JSON-headed ASCII files (RBSP/AutoPlot-compatible). A converter from the datamodel to HDF5 
is now available and a converter to NASA CDF is under development. Also under development is the reverse of
the SpaceData.flatten method, so that flattened objects can be restored to their former glory.


Empiricals Module
=================

The empiricals module provides access to some useful empirical models.
As of SpacePy 0.1.2, the models available are:

    * :func:`~spacepy.empiricals.getLmax` An empirical parametrization of the L* of the last closed drift shell
      (Lmax)
    * :func:`~spacepy.empiricals.getPlasmaPause` The plasmapause location, following either Carpenter and Anderson
      (1992) or Moldwin et al. (2002)
    * :func:`~spacepy.empiricals.getMPstandoff` The magnetopause standoff location (i.e. the sub-solar point), using
      the Shue et al. (1997) model
    * :func:`~spacepy.empiricals.vampolaPA` A conversion of omnidirectional electron flux to pitch-angle dependent
      flux, using the sin :sup:`n` model of Vampola (1996)

Each of the first three models is called by passing it a Ticktock object (see above) which then
calculates the model output using the 1-hour Qin-Denton OMNI data (from the
OMNI module; see above). For example::

>>> import spacepy.time as spt
>>> import spacepy.empiricals as emp
>>> ticks = spt.tickrange('2002-01-01T12:00:00','2002-01-04T00:00:00',.25)

calls :func:`~spacepy.time.tickrange` and makes a Ticktock object
with times from midday on January 1st 2002 to midnight January 4th 2002,
incremented 6-hourly::

>>> Lpp = emp.getPlasmaPause(ticks)

then returns the model plasmapause location using the default setting of the
Moldwin et al. (2002) model. The Carpenter and Anderson model can be used by
setting the Lpp_model keyword to 'CA1992'.

The magnetopause standoff location can be called using this syntax, or can be
called for specific solar wind parameters (ram pressure, P, and IMF Bz) passed
through in a Python dictionary::

>>> data = {'P': [2,4], 'Bz': [-2.4, -2.4]}
>>> emp.getMPstandoff(data)
>>>   # array([ 10.29156018,   8.96790412])


SeaPy - Superposed Epoch Analysis in Python
===========================================

Superposed epoch analysis is a technique used to reveal consistent responses,
relative to some repeatable phenomenon, in noisy data . Time series of the variables
under investigation are extracted from a window around the epoch and all data
at a given time relative to epoch forms the sample of events at that lag. The
data at each time lag are then averaged so that fluctuations not
consistent about the epoch cancel. In many superposed epoch analyses the mean of
the data at each time *u* relative to epoch, is used to
represent the central tendency. In SeaPy we calculate both the mean and the median,
since the median is a more robust measure of central tendency and is less affected
by departures from normality. SeaPy also calculates a measure of spread at each time
relative to epoch when performing the superposed epoch analysis; the interquartile
range is the default, but the median absolute deviation and bootstrapped confidence
intervals of the median (or mean) are also available.

As an example we fetch OMNI data for 4 years and perform a superposed epoch analysis
of the solar wind radial velocity, with a set of epoch times read from a text file::

>>> import datetime as dt
>>> import spacepy.seapy as sea
>>> import spacepy.omni as om
>>> import spacepy.toolbox as tb
>>> import spacepy.time as spt
>>> # now read the epochs for the analysis (the path specified is the default 
>>> # install location on linux, different OS will have this elsewhere)
>>> epochs = sea.readepochs('~/.local/lib/python2.7/site-packages/spacepy/data/SEA_epochs_OMNI.txt')

The readepochs function can handle multiple formats by a user-specified format code.
ISO 8601 format is directly supported though it is not used here. The the readepochs docstring
for more information. As above, we use the get_omni function to retrieve the hourly data 
from the OMNI module::

>>> ticks = spt.tickrange(dt.datetime(2005,1,1), dt.datetime(2009,1,1), dt.timedelta(hours=1))
>>> omni1hr = om.get_omni(ticks)
>>> omni1hr.tree(levels=1, verbose=True)

::

    +
    |____ByIMF (spacepy.datamodel.dmarray (35065,))
    |____Bz1 (spacepy.datamodel.dmarray (35065,))
    |____Bz2 (spacepy.datamodel.dmarray (35065,))
    |____Bz3 (spacepy.datamodel.dmarray (35065,))
    |____Bz4 (spacepy.datamodel.dmarray (35065,))
    |____Bz5 (spacepy.datamodel.dmarray (35065,))
    |____Bz6 (spacepy.datamodel.dmarray (35065,))
    |____BzIMF (spacepy.datamodel.dmarray (35065,))
    |____DOY (spacepy.datamodel.dmarray (35065,))
    |____Dst (spacepy.datamodel.dmarray (35065,))
    |____G (spacepy.datamodel.dmarray (35065, 3))
    |____Hr (spacepy.datamodel.dmarray (35065,))
    |____Kp (spacepy.datamodel.dmarray (35065,))
    |____Pdyn (spacepy.datamodel.dmarray (35065,))
    |____Qbits (spacepy.datamodel.SpaceData [7])
    |____RDT (spacepy.datamodel.dmarray (35065,))
    |____UTC (spacepy.datamodel.dmarray (35065,))
    |____W (spacepy.datamodel.dmarray (35065, 6))
    |____Year (spacepy.datamodel.dmarray (35065,))
    |____akp3 (spacepy.datamodel.dmarray (35065,))
    |____dens (spacepy.datamodel.dmarray (35065,))

and these data are used for the superposed epoch analysis.
the temporal resolution is 1 hr and the window is +/- 3 days

>>> delta = dt.timedelta(hours=1)
>>> window= dt.timedelta(days=3)
>>> sevx = sea.Sea(omni1hr['velo'], omni1hr['UTC'], epochs, window, delta)
    #rather than quartiles, we calculate the 95% confidence interval on the median
>>> sevx.sea(ci=True)
>>> sevx.plot()
