=====================
SpacePy Configuration
=====================

SpacePy has a few tunable options that can be altered through the ``spacepy.rc``
configuration file. All options have defaults which will be used if not specified in
the configuration file. These defaults are usually fine for most people and may
change between SpacePy releases, so we do not recommend changing the
configuration file without substantial reason.

``spacepy.rc`` lives in the per-user SpacePy directory, called ``.spacepy``.
On Unix-like operating systems, it is in a user's home directory; on Windows, 
in the user's Documents and Settings folder. If it doesn't exist, this directory
(and ``spacepy.rc``) is automatically created when SpacePy is imported.

``spacepy.rc`` has an INI-style format, parsed by :py:mod:`ConfigParser`. It
contains a single section, ``[spacepy]``.

    * `The spacepy directory`_
    * `Available configuration options`_
    * `Developer documentation`_

The spacepy directory
=====================

When first imported, spacepy will create a ``.spacepy`` directory in
your ``$HOME`` folder. If you prefer a different location for this
directory, set the environment variable ``$SPACEPY`` to a location of
your choice. For example, with a ``csh``, or ``tcsh`` you would::

	setenv SPACEPY /a/different/dir

for the ``bash`` shell you would:

	export SPACEPY=/a/different/dir

If you change the default location, make sure you add the environment
variable ``$SPACEPY`` to your ``.cshrc, .tcshrc,`` or ``.bashrc``
script.

Available configuration options
===============================
enable_deprecation_warning
  SpacePy raises :py:exc:`~exceptions.DeprecationWarning` when deprecated functions
  are called. Starting in Python 2.7, these are ignored. SpacePy adds a warnings
  filter to force display of deprecation warnings from SpacePy the first time a
  deprecated function is called. Set this option to False to retain the default
  Python behavior. (See :py:mod:`warnings` module for details on custom warning
  filters.)

keepalive
  True to attempt to use HTTP keepalives when downloading data in
  :py:func:`~spacepy.toolbox.update` (default). This is faster when
  downloading many small files but may be fragile (e.g. if a proxy
  server is required). Set to False for a more robust and flexible,
  but slower, codepath.

leapsec_url
  URL of the leapsecond database used by time conversions.
  :py:func:`~spacepy.toolbox.update` will download from the URL.
  The default should almost always be acceptable.

ncpus
  Number of CPUs to use for computations that can be
  multithreaded/multiprocessed. By default, they will use the number of CPUs
  reported by :py:func:`multiprocessing.cpu_count`. You may wish to set this
  to a lower number if you need to reserve other processors on your machine.

notice
  True to display the SpacePy license and other information on import (default);
  False to omit.

omni2_url
  URL containing the OMNI2 data.
  :py:func:`~spacepy.toolbox.update` will download from the URL.
  The default should almost always be acceptable.

qindenton_url
  URL containing Qin-Denton packaging of OMNI data as as single file.
  :py:func:`~spacepy.toolbox.update` will download from the URL.
  The default should almost always be acceptable.

qd_daily_url
  URL containing Qin-Denton packaging of OMNI data in daily files,
  supplemental to ``qindenton_url``. :py:func:`~spacepy.toolbox.update`
  will download from the URL. The default should almost always be
  acceptable.

psddata_url
  URL containing PSD data.
  :py:func:`~spacepy.toolbox.update` will download from the URL if requested.
  The default should almost always be acceptable.

support_notice
  True to display a notice on import if not a release version of SpacePy
  (default); False to omit. Those regularly installing from git instead
  of a release may want to set this to False.

user_agent
  User Agent for network access. If this is set,
  :func:`~spacepy.toolbox.update` will use this User Agent string on all
  HTTP requests. Normally leaving this unset should be fine.


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

For additions or fixes to this page, contact the SpacePy Team at Los Alamos.
