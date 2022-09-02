######################################
pycdf - Python interface to CDF files
######################################

.. automodule:: spacepy.pycdf

Contents
--------

.. contents::
   :depth: 2
   :local:

Create a CDF
============
This example presents the entire sequence of creating a CDF and populating
it with some data; the parts are explained individually below.

>>> from spacepy import pycdf
>>> import datetime
>>> time = [datetime.datetime(2000, 10, 1, 1, val) for val in range(60)]
>>> import numpy as np
>>> data = np.random.random_sample(len(time))
>>> cdf = pycdf.CDF('MyCDF.cdf', '')
>>> cdf['Epoch'] = time
>>> cdf['data'] = data
>>> cdf.attrs['Author'] = 'John Doe'
>>> cdf.attrs['CreateDate'] = datetime.datetime.now()
>>> cdf['data'].attrs['units'] = 'MeV'
>>> cdf.close()

Import the pycdf module.

>>> from spacepy import pycdf

Make a data set of :class:`~datetime.datetime`. These will be converted into
CDF_TIME_TT2000 types.

>>> import datetime
>>> # make a dataset every minute for a hour
>>> time = [datetime.datetime(2000, 10, 1, 1, val) for val in range(60)]

.. warning::
    If you create a CDF in backwards compatibility mode (using
    :meth:`~spacepy.pycdf.Library.set_backward`),
    then :class:`~datetime.datetime` objects are degraded to CDF_EPOCH
    (millisecond resolution), not CDF_EPOCH16 (microsecond resolution).
    Use :meth:`~spacepy.pycdf.CDF.new` to specify a data type.

Create some random data.

>>> import numpy as np
>>> data = np.random.random_sample(len(time))

Create a new empty CDF.  The empty string, '', is the name of the CDF to use 
as a master; given an empty string, an empty CDF will be created, rather than 
copying from a master CDF.
If a master is used, data in the master will be copied to the new CDF.

>>> cdf = pycdf.CDF('MyCDF.cdf', '')

.. note::
    You cannot create a new CDF with a name that already exists on disk.
    It will throw a :exc:`~exceptions.NameError`


To put data into a CDF, assign it directly to an element of the CDF.
CDF objects behave like Python dictionaries.

>>> # put time into CDF variable Epoch
>>> cdf['Epoch'] = time
>>> # and the same with data (the smallest data type that fits the data is used by default)
>>> cdf['data'] = data

Adding attributes is done similarly. CDF attributes are also treated as dictionaries.

>>> # add some attributes to the CDF and the data
>>> cdf.attrs['Author'] = 'John Doe'
>>> cdf.attrs['CreateDate'] = datetime.datetime.now()
>>> cdf['data'].attrs['units'] = 'MeV'

Closing the CDF ensures the new data are written to disk:

>>> cdf.close()

CDF files, like standard Python files, act as context managers

>>> with cdf.CDF('filename.cdf', '') as cdf_file:
...     #do brilliant things with cdf_file
>>> #cdf_file is automatically closed here


Read a CDF
==========
Reading a CDF is very similar: the CDF object behaves like a dictionary.
The file is only accessed when data are requested. A full example using the above CDF:

>>> from spacepy import pycdf
>>> cdf = pycdf.CDF('MyCDF.cdf')
>>> print(cdf)
    Epoch: CDF_TIME_TT2000 [60]
    data: CDF_FLOAT [60]
>>> cdf['data'][4]
    0.8609974384307861
>>> data = cdf['data'][...] # don't forget the [...]
>>> cdf_dat = cdf.copy()
>>> cdf_dat.keys()
    ['Epoch', 'data']
>>> cdf.close()

Again import the pycdf module

>>> from spacepy import pycdf

Then open the CDF, this looks the same and creation, but without mention of a master CDF.

>>> cdf = pycdf.CDF('MyCDF.cdf')

The default ``__str__()`` and ``__repr__()`` behavior explains the contents, type, and size but not the data.

>>> print(cdf)
    Epoch: CDF_TIME_TT2000 [60]
    data: CDF_FLOAT [60]

To access the data one has to request specific elements of the variable, similar to a Python list.

>>> cdf['data'][4]
    0.8609974384307861
>>> data = cdf['data'][...] # don't forget the [...]

:func:`CDF.copy` will return the entire contents of a CDF, including
attributes, as a :class:`~spacepy.datamodel.SpaceData` object:

>>> cdf_dat = cdf.copy()

Since CDF objects behave like dictionaries they have a ``keys()`` method and iterations are over the names in ``keys()``

>>> cdf_dat.keys()
    ['Epoch', 'data']

Close the CDF when finished:

>>> cdf.close()


Modify a CDF
============
An example modifying the CDF created above:

>>> from spacepy import pycdf
>>> cdf = pycdf.CDF('MyCDF.cdf')
>>> cdf.readonly(False)
    False
>>> cdf['newVar'] = [1.0, 2.0]
>>> print(cdf)
    Epoch: CDF_TIME_TT2000 [60]
    data: CDF_FLOAT [60]
    newVar: CDF_FLOAT [2]
>>> cdf.close()

As before, each step in this example will now be individually explained.
Existing CDF files are opened in read-only mode and must be set to read-write
before modification:

>>> cdf.readonly(False)
    False

Then new variables can be added:

>>> cdf['newVar'] = [1.0, 2.0]

Or contents can be changed:

>>> cdf['data'][0] = 8675309

You can write all new data to an existing variable, leaving the
variable type, dimensionality, and attributes unchanged:

>>> cdf['Epoch'][...] = [datetime.datetime(2010, 10, 1, 1, val)
...     for val in range(60)]

This is the common usage when using a CDF file containing all the
variables and attributes but no data, sometimes called a "master
CDF". Although the ``[...]`` makes this explicit (writing new records
not a new variable), the same syntax as for a new variable can also be
used:

>>> # Either create a new variable or overwrite data in existing
>>> cdf['Epoch'] = [datetime.datetime(2010, 10, 1, 1, val)
...     for val in range(60)]

The new variables appear immediately:

>>> print(cdf)
    Epoch: CDF_TIME_TT2000 [60]
    data: CDF_FLOAT [60]
    newVar: CDF_FLOAT [2]

Closing the CDF ensures changes are written to disk:

>>> cdf.close()

Non record-varying
==================
Non record-varying (NRV) variables are usually used for data that does not vary
with time, such as the energy channels for an instrument.

NRV variables need to be created with :func:`CDF.new`, specifying the keyword 'recVary' as False.

>>> from spacepy import pycdf
>>> cdf = pycdf.CDF('MyCDF2.cdf', '')
>>> cdf.new('data2', [1], recVary=False)
    <Var:
    CDF_BYTE [1] NRV
    >
>>> cdf['data2'][...]
    [1]

Slicing and indexing
====================
Subsets of data in a variable can be easily referenced with Python's slicing
and indexing notation.

This example uses :py:mod:`bisect` to read a subset of the data from the
hourly data file created in earlier examples.

>>> from spacepy import pycdf
>>> cdf = pycdf.CDF('MyCDF.cdf')
>>> start = datetime.datetime(2000, 10, 1, 1, 9)
>>> stop = datetime.datetime(2000, 10, 1, 1, 35)
>>> import bisect
>>> start_ind = bisect.bisect_left(cdf['Epoch'], start)
>>> stop_ind = bisect.bisect_left(cdf['Epoch'], stop)
>>> # then grab the data we want
>>> time = cdf['Epoch'][start_ind:stop_ind]
>>> data = cdf['data'][start_ind:stop_ind]
>>> cdf.close()

The :class:`Var` documentation has several additional examples.

.. _pycdf_string_handling:

String handling
===============

   .. versionchanged:: 0.3.0

   Prior to SpacePy 0.3.0, pycdf treated all strings as ASCII-encoded, and
   would raise errors when writing or reading strings that were not valid
   ASCII.

Per the NASA CDF library, variable and attribute names must be in
ASCII. The contents of ``CDF_CHAR`` and ``CDF_UCHAR`` were redefined
to be UTF-8 as of CDF 3.8.1. As of SpacePy 0.3.0, pycdf treats all
``CHAR`` variables with a default encoding of UTF-8. This is true
regardless of the version of the underlying CDF library.

UTF-8 is a variable-length encoding, so the number of elements in the
variable may not correspond to the number of characters if data are
not restricted to the ASCII range.

A different encoding can be specified with the ``encoding`` argument
to :class:`~spacepy.pycdf.CDF.open` and this encoding will be used on
all reads and writes to that file. Opening a CDF read-write with
``encoding`` other than ``utf-8`` or ``ascii`` will issue a warning.

Writing strings which cannot be represented in the desired encoding
will raise an error. When reading from a CDF, characters which cannot
be decoded will be replaced with the Unicode "replacement character"
U+FFFD, which usually displays as a question mark.

It is always possible to write raw bytes data to a variable, if it is
desired to use a different encoding for one time. For arrays of data,
this will usually involve :func:`numpy.char.encode`:

   >>> cdf['Variable'] = data.encode('latin-1')
   >>> cdf['Variable'] = numpy.char.encode(data, encoding='latin-1')

All encoding and decoding can also be skipped using the
:meth:`~spacepy.pycdf.CDF.raw_var` method to access a variable;
however, without encoding, only :class:`bytes` can be written to
string variables.

Troubleshooting
===============
Cannot load CDF C library
^^^^^^^^^^^^^^^^^^^^^^^^^
pycdf requires the standard NASA CDF library; it can be installed
after SpacePy. See specific instructions for :ref:`Linux <linux_CDF>`,
:ref:`Mac <install_mac_cdf>`, and :ref:`Windows <windows_CDF>`.

The error ``Cannot load CDF C library`` indicates pycdf cannot find
this library. pycdf searches in locations where the library is
installed by default; if the library is not found, set the ``CDF_LIB``
environment variable to the directory containing the library file
(.dll, .dylib, or .so) before importing pycdf.

ZLIB_ERROR when opening a CDF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The error message ``ZLIB_ERROR: Error during ZLIB decompression`` most
commonly occurs when opening a CDF which has been compressed with
whole-file compression. In this case, it must be unzipped into a
temporary location (details are in the `CDF User's Guide
<https://cdf.gsfc.nasa.gov/html/cdf_docs.html>`_).

The temporary location is specified by environment variables, most
commonly ``CDF_TMP``. It appears that, particularly on Windows, some
installers of the library may set this to a location which is not
writeable. In that case, the solution is to change the environment
variable to a writeable location.

On Windows, environment variables are set in the System Properties
control panel. Click the "Environment Variables" button on the
Advanced tab. Usually a good value for ``CDF_TMP`` is
``C:\Users\USERNAME\AppData\Local\Temp``. If ``CDF_TMP`` is not
set, variables ``TMP`` and ``TEMP`` will be used, so those values are
worth checking. Values starting with ``C:\WINDOWS\system32\config``
are unlikely to work.

On Unix, including MacOS, ``CDF_TMP`` is used if set; otherwise
``TMPDIR``.

Access to CDF constants and the C library
=========================================
Constants defined in cdf.h and occasionally useful in accessing CDFs are
available in the :mod:`~spacepy.pycdf.const` module.

The underlying C library is represented by the :attr:`~spacepy.pycdf.lib`
variable.

Classes
=======

.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    CDF
    Var
    gAttrList
    zAttrList
    zAttr
    gAttr
    AttrList
    Attr
    Library
    CDFCopy
    VarCopy
    CDFError
    CDFException
    CDFWarning
    EpochError

Functions
=========

.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    concatCDF

Submodules
==========

.. autosummary::
    :toctree: autosummary  
    :template: clean_module.rst

    const
    istp

Data
====

.. attribute:: lib

    Module global :class:`Library` object.

    Initalized at :mod:`~spacepy.pycdf` load time so all classes have ready
    access to the CDF library and a common state. E.g:

    >>> from spacepy import pycdf
    >>> pycdf.lib.version
        (3, 3, 0, ' ')
