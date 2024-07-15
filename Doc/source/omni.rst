####
OMNI
####

Tools to read and process omni data (Qin-Denton, etc.)

See also the `full API documentation <spacepy.omni>`.

.. contents:: Table of Contents
    :depth: 2
    :local:

About omni
----------
The omni module primarily manages the hourly OMNI2 and Qin-Denton data, which
are sourced from the Virtual Radiation Belt Observatory (`ViRBO <http://virbo.org>`_), 
who maintain these data sources. The data can be kept up-to-date in SpacePy 
using the :func:`~spacepy.toolbox.update` function in the :mod:`spacepy.toolbox` module.

The OMNI2 data combines data from a variety of satellites that sample the solar
wind (notably ACE and Wind), and propagates the data to Earth's bow shock nose.
The `Qin-Denton <http://virbo.org/QinDenton>`_ data is derived from the OMNI2 
data and is designed for providing input to the Tsyganenko magnetic field 
models. The later Tsyganenko magnetic field models require subsidiary parameters
(G- and W-parameters) that are pre-calculated in the Qin-Denton data. Further,
the Qin-Denton data contains no data gaps -- all gaps are filled (for details on
the gap filling, see the paper by `Qin et al. <http://dx.doi.org/10.1029/2006SW000296>`_.)


Advanced features
-----------------

Higher resolution data, or custom data sources, can also be managed/accessed 
with this module, although this is considered an advanced use for this module.
This is achieved using custom names for the dbase keyword in get_omni, which
must be defined in the SpacePy configuration file (for a user-install on linux,
this is ~/.spacepy/spacepy.rc; see :doc:`/configuration`).
An example of the formatting required is

qd1min: /usr/somedir/QinDenton/YYYY/QinDenton_YYYYMMDD_1min.txt

In this example the custom data source name is qd1min. Wildcard substitutions 
can be made for the year (YYYY), month (MM) and day (DD). Future updates will
give more flexibility in data storage model, but currently we assume that all
custom data sources follow a convention in which the data files are daily, and
the files are organized into folders by year. The year, month and day must all
be specified in the filename.

Currently there are some restrictions on the data format for custom data 
sources. The stored data must currently be stored as JSON-headed ASCII.
If data conversions are required, then a valid dictionary of conversion
functions must be supplied via the convert keyword argument. See 
:func:`~spacepy.datamodel.readJSONheadedASCII` for details.
Additionally, by default this will interpolate the data to the requested time 
ticks. To return only the actual recorded data values for the specified time 
range set the keyword argument interp to False.