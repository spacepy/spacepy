*************************************
SpacePy - A Quick Start Documentation
*************************************


Josef Koller and others


Installation
============

Option 1) to install it in the standard loction (depending on your system)::

	python setup.py install

Option 2) to install in custom location, e.g.::

	python setup.py install --home=/n/packages/lib/python


Toolbox - a box full of tools
=============================

Contains tools that don't fit anywhere else but are, in general, quite 
useful. The following function are implemented:

* item 1

* item 2

Time and Coordinate Transformations
===================================

Import the module as:: 

	import spacepy.spacetime as st
	

TickTock class
==============

The TickTock class provides as number of time conversion routines. 
The following time coordinates are provided

* UTC: Coordinated Universal Time implemented as a ``datetime.datetime`` class
* ISO: standard ISO 8601 format like ``2002-10-25T14:33:59``
* TAI: International Atomic Time in units of seconds since Jan 1, 1958 (midnight) and includes leap seconds, i.e. every second has the same length
* JD:  Julian Day
* MJD: Modified Julian Day
* UNX: UNIX time in seconds since Jan 1, 1970
* RDT: Rata Die Time (lat. fixed date) in days since Jan 1, 1 AD midnight
* CDF: CDF Epoch time in milliseconds since Jan 1, 0 
* DOY: Day of Year including fractions
* leaps: Leap seconds according to ftp://maia.usno.navy.mil/ser7/tai-utc.dat 

To access these time coordinates, you'll create an instance of a 
TickTock class, e.g.::

	t = st.TickTock('2002-10-25T12:30:00', 'ISO')

Instead of ISO you may use any of the formats listed above. You can also 
use numpy arrays or lists of time points. ``t`` has now the class 
attributes::

	t.dtype = 'ISO'
	t.data = '2002-10-25T12:30:00'

FYI ``t.UTC`` is added automatically.

If you want to convert/add a class attribute from the list above, 
simply type e.g.::

	t.RTD
	
You can replace RTD with any from the list above.

You can find out how many leap seconds were used by issuing the command::

	t.getleapsecs()


TickDelta class
===============

You can add/substract time from a TickTock class instance by creating a 
TickDelta instance first.::

	dt = st.TickDelta(days=2.3)

Then you can add by e.g.::

	t+dt 
	
	
SpaCo class
===========

The spatial coordinate class includes the following coordinate systems in 
cartesian and sphericals. 

* GZD:  (altitude, latitude, longitude in km, deg, deg
* GEO: cartesian, Re
* GSM: cartesian, Re
* GSE: cartesian, Re
* SM: cartesian, Re
* GEI: cartesian, Re
* MAG: cartesian, Re
* SPH: same as GEO but in spherical
* RLL: radial distance, latitude, longitude, Re, deg, deg.

Create a SpaCo instance with spherical='sph' or cartesian='car' 
coordinates::
 
 coord = st.SpaCo([[1,2,4],[1,2,2]], 'GEO', 'car')
 
This will let you request for example all y-coordinates by ``coord.y`` 
or if given in spherical coordinates by ``coord.lati``. One can transform 
the coordinates by ``newcoord = coord.convert('GSM', 'sph')``. 
This will return GSM coordinates in a spherical system. Since GSM 
coordinates depend on time, you'll have to add first a TickTock 
vector like ``coord.ticktock = st.TickTock(['2002-02-02T12:00:00', 
'2002-02-02T12:00:00'], 'ISO')``
 
 
RadBelt Module
==============

The radiation belt module currently include a simple radial 
diffusion code as a class. Import the module and create a class::

	import spacepy.radbelt as sprb
	rb = sprb.RBmodel()

Add a time grid for a particular period that you are interested in::

	rb.setup_ticks('2002-02-01T00:00:00', '2002-02-10T00:00:00', 0.25)

This will automatically lookup required geomagnetic/solar wind conditions 
for that period. Run the diffusion solver for that setup and plot the 
results.::

	rb.evolve()
	rb.plot()


OMNI Module
===========

bla bla


ONERA-DESP Module
=================

bla bla


pycdf - Python access to NASA CDF library
=========================================

pycdf provides a "pythonic" interface to the NASA CDF library (currently
read-only). It requires that the base C library be properly installed.
The module can then be imported, e.g.::

	import spacepy.pycdf as cdf

Extensive documentation is provided in epydoc format in docstrings.

To open and close a CDF file::

	cdf_file = cdf.CDF('filename.cdf')
	cdf_file.close()

CDF files, like standard Python files, act as context managers::

       with cdf.CDF('filename.cdf') as cdf_file:
	    #do brilliant things with cdf_file
       #cdf_file is automatically closed here

CDF files act as Python dictionaries, holding CDF variables keyed
by the variable name::

       var_names = keys(cdf_file) #list of all variables
       for var_name in cdf_file:
	   print(len(cdf_file[var_name])) #number of records in each variable
       #list comprehensions work, too
       lengths = [len(cdf_file[var_name]) for var_name in cdf_file]

Each CDF variable acts as a Python list, one element per record.
Multidimensional CDF variables are represented as nested lists and can be
subscripted using a multidimensional slice notation similar to numpy. Creating
a Python Var object does not read the data from disc; data are only read as
they are accessed::

       epoch = cdf_file['Epoch'] #Python object created, nothing read from disc
       epoch[0] #time of first record in CDF (datetime object)
       a = epoch[...] #copy all times to list a
       a = epoch[-5:] #copy last five times to list a
       b_gse = cdf_file['B_GSE'] #B_GSE is a 1D, three-element array
       bz = b_gse[0,2] #Z component of first record
       bx = b_gse[:,0] #copy X component of all records to bx
       bx = cdf_file['B_GSE'][:,0] #same as above


The testing.py module
=====================

Is supposed to test the implementation of spacepy modules.

