#########
Datamodel
#########

The datamodel classes constitute a data model implementation
meant to mirror the functionality of the data model output from pycdf, though
implemented slightly differently.

.. contents:: Table of Contents
    :depth: 2
    :local:

.. currentmodule:: spacepy.datamodel

See also the `full API documentation <spacepy.datamodel>`.

Documentation
=============

This contains the following classes:
 * :class:`~spacepy.datamodel.dmarray` - numpy arrays that support .attrs for information about the data
 * :class:`~spacepy.datamodel.SpaceData` - base class that extends dict, to be extended by others

Authors: Steve Morley and Brian Larsen

Additional Contributors: Charles Kiyanda and Miles Engel

Institution: Los Alamos National Laboratory

Contact: smorley@lanl.gov; balarsen@lanl.gov

Copyright 2010-2016 Los Alamos National Security, LLC.


About datamodel
---------------

The SpacePy datamodel module implements classes that are designed to make implementing a standard
data model easy. The concepts are very similar to those used in standards like HDF5, netCDF and
NASA CDF.

The basic container type is analogous to a folder (on a filesystem; HDF5 calls this a
group): Here we implement this as a dictionary-like object, a datamodel.SpaceData object, which
also carries attributes. These attributes can be considered to be global, i.e. relevant for the
entire folder. The next container type is for storing data and is based on a numpy array, this
class is datamodel.dmarray and also carries attributes. The dmarray class is analogous to an
HDF5 dataset.

In fact, HDF5 can be loaded directly into a SpacePy datamodel, carrying across all attributes,
using the function fromHDF5:

>>> import spacepy.datamodel as dm
>>> data = dm.fromHDF5('test.h5')

Functions are also available to directly load data and metadata into a
SpacePy datamodel from NASA CDF as well as JSON-headed ASCII. Writers also
exist to output a SpacePy datamodel directly to HDF5 or JSON-headed ASCII.
See :func:`~spacepy.datamodel.fromCDF`, :func:`~spacepy.datamodel.readJSONheadedASCII`,
:func:`~spacepy.datamodel.toHDF5`, and :func:`~spacepy.datamodel.toJSONheadedASCII` for more details.


Examples
--------

Imagine representing some satellite data within the global attributes might be
the mission name and the instrument PI, the variables might be the
instrument counts [n-dimensional array], timestamps[1-dimensional array and an orbit number [scalar].
Each variable will have one attribute (for this example).

>>> import spacepy.datamodel as dm
>>> mydata = dm.SpaceData(attrs={'MissionName': 'BigSat1'})
>>> mydata['Counts'] = dm.dmarray([[42, 69, 77], [100, 200, 250]], attrs={'Units': 'cnts/s'})
>>> mydata['Epoch'] = dm.dmarray([1, 2, 3], attrs={'units': 'minutes'})
>>> mydata['OrbitNumber'] = dm.dmarray(16, attrs={'StartsFrom': 1})
>>> mydata.attrs['PI'] 'Prof. Big Shot'

This has now populated a structure that can map directly to a NASA CDF, HDF5 or JSON-headed ASCII file.
To visualize our datamodel, we can use tree method (which can be applied to any dictionary-like object
using :func:`~spacepy.toolbox.dictree`).

>>> mydata.tree(attrs=True)

::

    +
    :|____MissionName
    :|____PI
    |____Counts
         :|____Units
    |____Epoch
         :|____units
    |____OrbitNumber
         :|____StartsFrom


Guide for NASA CDF users
------------------------
By definition, a NASA CDF only has a single 'layer'. That is, a CDF contains a series of records
(stored variables of various types) and a set of attributes that are either global or local in
scope. Thus to use SpacePy's datamodel to capture the functionality of CDF the two basic data types
are all that is required, and the main constraint is that datamodel.SpaceData objects cannot be
nested (more on this later, if conversion from a nested datamodel to a flat datamodel is required).


Opening a CDF and working directly with the contents can be easily done using the PyCDF module, however,
if you wish to load the entire contents of a CDF directly into a datamodel (complete with attributes)
the following will make life easier:

>>> import spacepy.datamodel as dm
>>> data = dm.fromCDF('inFile.cdf')


A quick guide to JSON-headed ASCII
----------------------------------
In many cases it is preferred to have a human-readable ASCII file, rather than a binary file like CDF
or HDF5. To make it easier to carry all the same metadata that is available in HDF5 or CDF we have
developed an ASCII data storage format that encodes the metadata using JSON (JavaScript Object Notation).
This notation supports two basic datatypes: key/value collections (like a SpaceData) and ordered lists
(which can represent arrays). JSON is human-readable, but if large arrays are stored in metadata is quickly
becomes difficult to read. For this reason we use JSON to encode the metadata (usually smaller datasets)
and store the data in a standard flat-ASCII format. The metadata is provided as a header that describes
the contents of the file.


To use JSON for storing only metadata associated with the data to be written to an ASCII file a minimal
metadata standard must be implemented. We use the following attribute names: DIMENSION and START_COLUMN.
We also recommend using the NASA ISTP metadata standard to assign attribute names. The biggest limitation
of flat ASCII is that sensibly formatting datasets of more than 2-dimensions (i.e. ranks greater than 2)
is not possible. For this reason if you have datasets of rank 3 or greater then we recommend using HDF5.
If text is absolutely required then it is possible to encode multi-dimensional arrays in the JSON metadata,
but this is not recommended.


This format is best understood by illustration. The following example builds a toy SpacePy datamodel and
writes it to a JSON-headed ASCII file. The contents of the file are then shown.

>>> import spacepy.datamodel as dm
>>> data = dm.SpaceData()
>>> data.attrs['Global'] = 'A global attribute'
>>> data['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
>>> data['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
>>> data['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
>>> dm.toJSONheadedASCII('outFile.txt', data, depend0='Var1', order=['Var1'])
#Note that not all field names are required, those not given will be listed
#alphabetically after those that are specified

The file looks like:

.. code-block:: none

    #{
    #    "MVar": {
    #        "Note": "Metadata",
    #        "VALUES": [7.8]
    #    },
    #    "Global": "A global attribute",
    #    "Var1": {
    #        "Local1": "A local attribute",
    #        "DIMENSION": [1],
    #        "START_COLUMN": 0
    #    },
    #    "Var2": {
    #        "DIMENSION": [2],
    #        "START_COLUMN": 2
    #    }
    #}
    1 8 9
    2 9 1
    3 3 4
    4 8 9
    5 7 8