

######################################
pyCDF - Python interface to CDF files
######################################

.. automodule:: spacepy.pycdf

Contents
--------

- `Quickstart`_
    - `Create a CDF`_
    - `Read a CDF`_
    - `Modify a CDF`_
    - `Non record varying`_
    - `Slicing and indexing`_
- `Class reference`_

Quickstart
----------
Create a CDF
========
This quickstart guide should walk uses through the basics of CDF manipulation
using pyCDF

To create a new CDF from some data for your own use or to send to a colleague.  We will show the example then explain the parts.
    >>> from spacepy import pycdf
    >>> import numpy as np
    >>> import datetime
    >>> time = [datetime.datetime(2000, 10, 1, 1, val) for val in range(60)]
    >>> data = np.random.random_sample(len(time))
    >>> cdf = pycdf.CDF('MyCDF.cdf', '')
    >>> cdf['Epoch'] = time
    >>> cdf['data'] = data
    >>> cdf.attrs['Author'] = 'John Doe'
    >>> cdf.attrs['CreateDate'] = datetime.datetime.now()
    >>> cdf['data'].attrs['units'] = 'MeV'
    >>> cdf.close()

Import the pyCDF module.  This can be done however you like.
    >>> from spacepy import pycdf

Make a datetime data set, these are automatically converted into CDF_EPOCH types.
    >>> import datetime
    >>> # make a dataset every minute for a hour
    >>> time = [datetime.datetime(2000, 10, 1, 1, val) for val in range(60)]

.. warning::
    If you create a CDF in backwards compatibility mode (default) then :class:`datetime.datetime` objects are degraded to CDF_EPOCH, not CDF_EPOCH16 type.  This means millisecond resolution vs microsecond resolution.

Create some data of an arbitrary type
    >>> data = np.random.random_sample(len(time))

Create a new empty CDF.  The '' is the name of the CDF to use as a master.  Note that the data is copied form the master to the new CDF.
    >>> cdf = pycdf.CDF('MyCDF.cdf', '')
.. note::
    You cannot create a new CDF with a name that already exists on disk.  It will throw an :class:`exceptions.NameError`


To put that data into the CDF just do it, CDF objects behave like Python dictionaries.
    >>> # put time into the cdf as 'Epoch'
    >>> cdf['Epoch'] = time
    >>> # and the same with data (note that the smallest data type that fits the data is used by default)
    >>> cdf['data'] = data

Adding attributes is done the same way.  CDF variables are also treated as dictionaries.
    >>> # add some attributes to the variable data and to the global cdf
    >>> cdf.attrs['Author'] = 'John Doe'
    >>> cdf.attrs['CreateDate'] = datetime.datetime.now()
    >>> cdf['data'].attrs['units'] = 'MeV'

It is best to close the CDF manually to be sure you know that it will be written.
    >>> # and be sure to close the cdf to assure it is written
    >>> cdf.close()

CDF files, like standard Python files, act as context managers
    >>> with cdf.CDF('filename.cdf', '') as cdf_file:
    ...     #do brilliant things with cdf_file
    >>> #cdf_file is automatically closed here


Read a CDF
=======

Reading a CDF is done in much the same way, the CDF object behaves like a dictionary and only goes to disk when you request the data.  Shown here is a full example using the above CDF then, then explained (see also :py:func:`datamodel.fromCDF`).
    >>> from spacepy import pycdf
    >>> cdf = pycdf.CDF('MyCDF.cdf')
    >>> print(cdf)
        Epoch: CDF_EPOCH [60]
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
        Epoch: CDF_EPOCH [60]
        data: CDF_FLOAT [60]

To access the data one has to request specific elements of the variable that behaves like a list in this respect.
    >>> cdf['data'][4]
        0.8609974384307861
    >>> data = cdf['data'][...] # don't forget the [...]

One can also grab the entire contest of a CDF using the convenience routine :py:func:`pycdf.CDF.copy`
    >>> cdf_dat = cdf.copy()

Since CDF objects behave like dictionaries they have a ``keys()`` method and iterations are over the names in ``keys()``
    >>> cdf_dat.keys()
        ['Epoch', 'data']

and with writing it is Best to close the CDF (or use context managers)
    >>> cdf.close()

Modify a CDF
============
Again using the CDF created above a variable can be added or the contents of a variable changed.
    >>> from spacepy import pycdf
    >>> cdf = pycdf.CDF('MyCDF.cdf')
    >>> cdf.readonly(False)
        False
    >>> cdf['newVar'] = [1.0, 2.0]
    >>> print(cdf)
        Epoch: CDF_EPOCH [60]
        data: CDF_FLOAT [60]
        newVar: CDF_FLOAT [2]
    >>> cdf.close()

The parts of the example are straightforward.  A particular open CDF must be made write-able
    >>> cdf.readonly(False)
        False

Then new variables can be added
    >>> cdf['newVar'] = [1.0, 2.0]

Or contents changed
    >>> cdf['data'][0] = 8675309

And the new variable shows up
    >>> print(cdf)
        Epoch: CDF_EPOCH [60]
        data: CDF_FLOAT [60]
        newVar: CDF_FLOAT [2]

As with writing be sure to close the CDF
    >>> cdf.close()

Non record varying
==================
Creating a variable that is non record varying is really useful in the conversion of text files to CDF where whole columns do not change.


To create a variable that is non-record varying one has to manually create the variable using :py:func:`pycdf.CDF.new`.
    >>> from spacepy import pycdf
    >>> cdf = pycdf.CDF('MyCDF2.cdf', '')
    >>> # create a variable manually
    >>> cdf.new('data2', [1], recVary=False)
        <Var:
        CDF_BYTE [1] NRV
        >

Slicing and indexing
====================
This example is redundant to the above but worth a call out as it is a very common operation.

If one has the hourly data file created above and only wants to read in a portion of the data follow this recipe.  Using :py:mod:`bisect` can save a lot of disk I/O.
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


Class reference
===============



.. autosummary::
    :toctree:

    Attr
    zAttr
    gAttr
    AttrList
    gAttrList
    zAttrList
    const
    CDF
    CDFCopy
    Library
    Var
    VarCopy
    CDFError
    CDFException
    CDFWarning
    EpochError

