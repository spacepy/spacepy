====================
SpacePy Capabilities
====================

This page lists some capabilities of SpacePy and, in some cases, of
other packages that might be of interest to SpacePy users. It is
organized by topic; searching within this page is recommended. See the
:ref:`module reference <module_reference>` for every class/function
available in SpacePy, organized by module.

.. contents::
   :local:

Array manipulation
==================
Various :mod:`~spacepy.toolbox` functions are useful in manipulating
NumPy arrays. :mod:`~spacepy.datamanager` also contains many functions
for indexing arrays and manipulating them in ways that do not depend
on the interpretation of their contents.

Coordinate Transforms
=====================
:mod:`~spacepy.coordinates` provides a class for transforming among
most coordinate systems used in Earth magnetospheric and
ionospheric physics.

It also provides generalized coordinate transforms via quaternions.

File I/O
========
:mod:`~spacepy.pycdf` provides reading and writing of NASA CDF files,
with additional functionality for those with ISTP-compliant metadata.

:mod:`~spacepy.datamodel` provides easy reading and writing of HDF5
and most netCDF files. It also supports reading and writing ASCII-based
data files with rich JSON metadata, supported by tools such as
`Autoplot <http://autoplot.org>`_.

:func:`numpy.loadtxt` and related functions are helpful for reading
various "plain-text" files into numerical arrays.

:func:`scipy.io.readsav` reads IDL savesets.

:mod:`astropy.io.fits` supports FITS files.

Modeling
========
:mod:`~spacepy.ae9ap9` supports import and visualization of data from
the AE9/AP9 empirical radiation belt model.

:mod:`~spacepy.empiricals` implements several simple empirical and/or
analytic models for magnetospheric and solar wind phenomena, including
the plasmapause location, the Shue magnetopause model, and solar wind
temperature.

:mod:`~spacepy.pybats` supports output analysis and visualization of
many models compatible with the Space Weather Modeling Framework,
including the BATS-R-US global MHD model and the RAM-SCB ring current
model.

:mod:`~spacepy.omni` provides ready access to the `OMNI
<https://omniweb.gsfc.nasa.gov/>`_ near-Earth solar wind dataset,
useful for model inputs.

Statistics
==========
:mod:`~spacepy.poppy` supports determining confidence intervals on
 population metrics using the non-parametric bootstrap method.

Time conversions
================
:mod:`~spacepy.time` contains a class that easily allows time to be
represented in, and converted among, many representations, including
Python datetimes, ISO time strings, GPS time, TAI, etc.

Time series analysis and correlations
=====================================
:mod:`~spacepy.poppy` implements association analysis to determine the
relationship between point-in-time events.

:mod:`~spacepy.seapy` implements superposed epoch analysis, the
statistical evaluation of the time evolution of a system relative
to a set of starting epochs.

Visualization
=============
:mod:`~spacepy.plot` provides tools useful in making
publication-quality plots with the `matplotlib
<https://matplotlib.org/>`_ toolkit.
