=====================
SpacePy Configuration
=====================

SpacePy has a few tunable options that can be altered through the ``spacepy.rc``
config file. All options have defaults which will be used if not specified in
the config file. These defaults are usually fine for most people and may
change from version to version of SpacePy, so we do not recommend changing the
config file without substantial reason.

``spacepy.rc`` lives in the per-user SpacePy directory, called ``.spacepy``.
On Unix-like operating systems, it is in a user's home directory; on Windows, 
in the user's Document and Settings folder. If it doesn't exist, this directory
(and ``spacepy.rc``) is automatically created when SpacePy is imported.

``spacepy.rc`` has an INI-style format, parsed by :py:mod:`ConfigParser`. It
contains a single section, ``[spacepy]``.

    * `Available configuration options`_
    * `Developer documentation`_


Available configuration options
===============================
enable_deprecation_warning
  SpacePy raises :py:exc:`DeprecationWarning` when deprecated functions
  are called. Starting in Python 2.7, these are ignored. SpacePy adds a warnings
  filter to force display of deprecation warnings from SpacePy the first time a
  deprecated function is called. Set this option to False to retain the default
  Python behavior. (See :py:mod:`warnings` module for details on custom warning
  filters.)

leapsec_url
  URL of the leapsecond database used by time conversions.
  :py:func:`toolbox.update` will download from the URL.
  The default should almost always be acceptable.

ncpus
  Number of CPUs to use for computations that can be
  multithreaded/multiprocessed. By default, they will use the number of CPUs
  reported by :py:func:`multiprocessing.cpu_count`. You may wish to set this
  to a lower number if SpacePy makes other processes on your machine slow.

omni_url
  URL containing OMNI data.
  :py:func:`toolbox.update` will download from the URL.
  The default should almost always be acceptable.

psddata_url
  URL containing PSD data.
  :py:func:`toolbox.update` will download from the URL if requested.
  The default should almost always be acceptable.


Developer documentation
=======================
``spacepy.rc`` is loaded into a dictionary (``spacepy.config``) by SpacePy's
main ``__init__.py``. All options from the ``[spacepy]`` section are loaded,
with no developer intervention needed. Each key is the option's name; the
associated value is the option's value. To specify a default, add to the
``defaults`` dictionary at the top of ``_read_config``; each default, if not
overridden by the config file, will be included in the config dict. Values are
assumed to be strings. The ``caster`` dictionary is keyed by option name; the
value for each key is a function to be applied to the value with the same key
to produce a different type from a string.


--------------------------

:Release: |version|
:Doc generation date: |today|

For additions or fixes to this page, contact Jon Niehof at Los Alamos.
