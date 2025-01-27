#############################################################################################
Coordinates - Implementation of Coords class functions for coordinate transformations
#############################################################################################

.. contents:: Table of Contents
    :depth: 2
    :local:

.. currentmodule:: spacepy.coordinates


See also the `full API documentation <spacepy.coordinates>`.

Introduction
------------
The coordinate systems supported by this module cover the most commonly
used geophysical and magnetospheric systems. The naming conventions can
follow the names used by the popular IRBEM library, but for inertial
systems we use a more consistent, fine-grained naming convention that
clarifies the different systems.

Earth-centered Inertial Systems
-------------------------------
    * **ECI2000** Earth-centered Inertial, J2000 epoch
    * **ECIMOD** Earth-centered Inertial, mean-of-date
    * **ECITOD** Earth-centered Inertial, true-of-date
    * **GEI** Geocentric Equatorial Inertial (IRBEM approximation of TOD)

Magnetospheric Systems
----------------------
    * **GSM** Geocentric Solar Magnetospheric
    * **GSE** Geocentric Solar Ecliptic
    * **SM** Solar Magnetic
    * **MAG** Geomagnetic Coordinate System (aka CDMAG)

Earth-fixed Systems
-------------------
    * **GEO** Geocentric geographic, aka Earth-centered Earth-fixed
    * **GDZ** Geodetic coordinates

By convention *all* systems are treated as natively Cartesian except
geodetic (GDZ), which is defined in [altitude, latitude, longitude]
where altitude is relative to a reference ellipsoid. Similarly, distance
units are assumed to be Earth radii (Re) in all systems except GDZ, where
altitude is given in km. Conversions to GDZ will output altitude in km
regardless of the input distance units and conversions from GDZ will
output in Re regardless of input units. In all other cases, the distance
units will be preserved.

    .. versionchanged:: 0.3.0

    The new CTrans backend was added, which includes support for the names
    ``ECI2000``, ``ECIMOD``, ``ECITOD``, and ``CDMAG``. With the exception
    of ``ECIMOD``, these can be used with the existing IRBEM backend, and
    will be converted to their closest equivalents.

    .. versionchanged:: 0.4.0

    The default backend for coordinate transformations was changed from IRBEM
    to the CTrans-based SpacePy backend.

    .. versionchanged:: 0.8.0

    The IGRF model for the CTrans-based SpacePy backend was updated to IGRF14;
    the IRBEM backend was also updated to a version using IGRF14.

Differences between representations
--------------------------------------------
IRBEM's coordinate transformations are low-accuracy and were written for
a library with a driving philosophy of speed and robustness as priorities.
The coordinate transformations are therefore approximate. Further, most of
the geophysical systems (e.g., GSE, SM) are derived from an inertial
system. It is standard practice to use ECIMOD as this system. However,
IRBEM does not currently make ECIMOD available as one of its inertial
systems. IRBEM's default inertial system (called GEI) is consistent with
an approximation of ECITOD. Hence there will be small differences between
IRBEM's transformations and those using SpacePy's CTrans backend.
Further sources of difference include: IRBEM uses a low-order approximation
to the sidereal time and other parameters; the calculation of the Earth-Sun
vector differs between the representations; the definitions of an Earth
radius differ (SpacePy = 6378.137km; IRBEM = 6371.2 km). SpacePy's in-built
representation is higher accuracy and is comprehensively tested, including
tests for consistency with other high accuracy packages such as LANLGeoMag
and AstroPy. However, for use cases where the required precision is of order
1 percent the output can be considered equivalent.

Setting options for coordinate transformation
---------------------------------------------
The backend for coordinate transformations can be provided at
instantiation of a :class:`~spacepy.coordinates.Coords` object using a keyword
argument. However, for convenience and flexibility the options can be
set at the module level. Configurable options include the backend used
(:mod:`~spacepy.irbempy` or SpacePy's :mod:`~spacepy.ctrans`) and the
reference ellipsoid (only configurable for the SpacePy backend). A
warning will be raised if the backend is not set (either through the
defaults or the keyword argument). The final configurable option
(``itol``) is the maximum separation, in seconds, for which the
coordinate transformations will not be recalculated. To force all
transformations to use an exact transform for the time, set ``itol``
to zero. Values between 10s and 60s are recommended for speed while
also preserving accuracy, though different applications will require
different accuracies.  For example, assuming this module has been
imported as ``spc``, to set the SpacePy backend as the default and set
``itol`` to 5 seconds:

    >>> spc.DEFAULTS.set_values(use_irbem=False, itol=5)

Note also that the magnetospheric systems listed above require the
specification of a centered dipole axis. This is specified using the
first three coefficients of the IGRF, evaluated at the specified time
(specifically, the g[1][0], g[1][1], and h[1][1] coefficients).
For the SpacePy backend, this uses the bundled IGRF coefficients file
in the standard format. The coefficients from a different version of
the model can be used by, for example, downloading the desired model
version from `<https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/>`_,
and placing them in the ``data`` subdirectory of the
:doc:`.spacepy directory <configuration>`
with the name ``igrfcoeffs.txt``. This would be used, for example, if
the user wishes to reproduce a study using IGRF12. Note that the
IGRF coefficients are hard-coded in IRBEM thus changing the IGRF
version is non-trivial.
