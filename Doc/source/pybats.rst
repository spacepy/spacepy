
########################################
PyBats - SWMF & BATS-R-US Analysis Tools
########################################

PyBats!  An open source Python-based interface for reading, manipulating,
and visualizing BATS-R-US and SWMF output.
For more information on the SWMF, please visit the
`Center for Space Environment Modeling <http://csem.engin.umich.edu>`_.

See also the `full API documentation <spacepy.pybats>`.

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