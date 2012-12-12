#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The datamodel classes constitute a data model implementation
meant to mirror the functionality of the data model output from pycdf, though
implemented slightly differently.

This contains the following classes:
 * :py:class:`dmarray` - numpy arrays that support .attrs for information about the data
 * :py:class:`SpaceData` - base class that extends dict, to be extended by others

Currently used in GPScode and other projects

Authors: Steve Morley and Brian Larsen

Additional Contributors: Charles Kiyanda and Miles Engel

Institution: Los Alamos National Laboratory

Contact: smorley@lanl.gov; balarsen@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.


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

This has now populated a structure that can map directly to a NASA CDF or a HDF5. To visualize our datamodel,
we can use tree method (which can be applied to any dictionary-like object using the 
:py:func:`toolbox.dictree` function).

>>> mydata.tree(attrs=True)
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
By definition, a NASA CDF only has a single `layer'. That is, a CDF contains a series of records
(stored variables of various types) and a set of attributes that are either global or local in
scope. Thus to use SpacePy's datamodel to capture the functionality of CDF the two basic data types
are all that is required, and the main constraint is that datamodel.SpaceData objects cannot be
nested (more on this later, if conversion from a nested datamodel to a flat datamodel is required).


Opening a CDF and working directly with the contents can be easily done using the PyCDF module, however,
if you wish to load the entire contents of a CDF directly into a datamodel (complete with attributes)
the following will make life easier:

>>> from spacepy import pycdf
>>> with pycdf.CDF('test.cdf') as cdffile:
...     data = cdffile.copy()
"""

from __future__ import division
import copy, datetime, os, warnings
import re, json
try:
    import StringIO # can't use cStringIO as we might have unicode
except ImportError:
    import io as StringIO

from . import toolbox
import numpy


__contact__ = 'Steve Morley, smorley@lanl.gov'

class DMWarning(Warning):
    """
    Warnings class for datamodel, subclassed so it can be set to always
    """
    pass
warnings.simplefilter('always', DMWarning)

class dmarray(numpy.ndarray):
    """
    Container for data within a SpaceData object

    Raises
    ------
    NameError
        raised is the request name was not added to the allowed attributes list

    Examples
    --------
    >>> import spacepy.datamodel as datamodel
    >>> position = datamodel.dmarray([1,2,3], attrs={'coord_system':'GSM'})
    >>> position
    dmarray([1, 2, 3])
    >>> position.attrs
    {'coord_system': 'GSM'}a

    The dmarray, like a numpy ndarray, is versatile and can store
    any datatype; dmarrays are not just for arrays.

    >>> name = datamodel.dmarray('TestName')
    dmarray('TestName')

    To extract the string (or scalar quantity), use the tolist method

    >>> name.tolist()
    'TestName'

    .. currentmodule:: spacepy.datamodel
    .. autosummary::
        ~dmarray.addAttribute
    .. automethod:: addAttribute
    """
    Allowed_Attributes = ['attrs']

    def __new__(cls, input_array, attrs=None, dtype=None):
       # Input array is an already formed ndarray instance
       # We first cast to be our class type
       if not dtype:
           obj = numpy.asarray(input_array).view(cls)
       else:
           obj = numpy.asarray(input_array).view(cls).astype(dtype)
       # add the new attribute to the created instance
       if attrs != None:
           obj.attrs = attrs
       else:
           obj.attrs = {}
       # Finally, return the newly created object:
       return obj

    def __array_finalize__(self, obj):
       # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        for val in self.Allowed_Attributes:
            self.__setattr__(val, copy.deepcopy(getattr(obj, val, {})))

    def __reduce__(self):
        """This is called when pickling, see:
        http://www.mail-archive.com/numpy-discussion@scipy.org/msg02446.html
        for this particular example.
        Only the attributes in Allowed_Attributes can exist
        """
        object_state = list(numpy.ndarray.__reduce__(self))
        subclass_state = tuple([tuple([val, self.__getattribute__(val)]) for val in self.Allowed_Attributes])
        object_state[2] = (object_state[2],subclass_state)
        return tuple(object_state)

    def __setstate__(self, state):
        """Used for unpickling after __reduce__ the self.attrs is recovered from
        the way it was saved and reset.
        """
        nd_state, own_state = state
        numpy.ndarray.__setstate__(self,nd_state)
        for i, val in enumerate(own_state):
            if not val[0] in self.Allowed_Attributes: # this is attrs
                self.Allowed_Attributes.append(own_state[i][0])
            self.__setattr__(own_state[i][0], own_state[i][1])

    def __setattr__(self, name, value):
        """Make sure that .attrs is the only attribute that we are allowing
        dmarray_ne took 15.324803 s
        dmarray_eq took 15.665865 s
        dmarray_assert took 16.025478 s
        It looks like != is the fastest, but not by much over 10000000 __setattr__
        """
        if name == 'Allowed_Attributes':
            pass
        elif not name in self.Allowed_Attributes:
            raise(TypeError("Only attribute listed in Allowed_Attributes can be set"))
        super(dmarray, self).__setattr__(name, value)

    def addAttribute(self, name, value=None):
        """Method to add an attribute to a dmarray
        equivalent to
        a = datamodel.dmarray([1,2,3])
        a.Allowed_Attributes = a.Allowed_Attributes + ['blabla']
        """
        if name in self.Allowed_Attributes:
            raise(NameError('{0} is already an attribute cannot add again'.format(name)))
        self.Allowed_Attributes.append(name)
        self.__setattr__(name, value)

    def count(self, srchval):
        """
        Equivalent to count method on list

        """
        mask = self == srchval
        return int(mask.sum())


class SpaceData(dict):
    """
    Datamodel class extending dict by adding attributes.

    .. currentmodule:: spacepy.datamodel
    .. autosummary::
        ~SpaceData.flatten
        ~SpaceData.tree
    .. automethod:: flatten
    .. automethod:: tree
    """

    def __init__(self, *args, **kwargs):
        """
        Base class for "Data Model" representation data

        Abstract method, reimplement

        Attributes
        ----------
        attrs : dict
            dictionary of the attributes of the SpaceData object

        """
        #raise(ValueError("Abstract method called, reimplement __init__"))
        self.attrs = {}
        if 'attrs' in kwargs:
            if hasattr(kwargs['attrs'], '__getitem__'):
                self.attrs = kwargs['attrs']
            del kwargs['attrs']

        super(SpaceData, self).__init__(*args, **kwargs)

    def tree(self, **kwargs):
        '''Print the contents of the SpaceData object in a visual tree

        Examples
        --------
        >>> import spacepy.datamodel as dm
        >>> import spacepy.toolbox as tb
        >>> a = dm.SpaceData()
        >>> a['1'] = dm.SpaceData(dog = 5)
        >>> a['4'] = dm.SpaceData(cat = 'kitty')
        >>> a['5'] = 4
        >>> a.tree()
        +
        |____1
             |____dog
        |____4
             |____cat
        |____5
        '''
        toolbox.dictree(self, **kwargs)

    def flatten(self):
        '''
        Method to collapse datamodel to one level deep

        Examples
        --------
        >>> import spacepy.datamodel as dm
        >>> import spacepy.toolbox as tb
        >>> a = dm.SpaceData()
        >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        >>> a['4'] = dm.SpaceData(cat = 'kitty')
        >>> a['5'] = 4
        >>> tb.dictree(a)
        +
        |____1
             |____dog
             |____pig
                  |____fish
                       |____a
                       |____b
        |____4
             |____cat
        |____5

        >>> b = dm.flatten(a)
        >>> tb.dictree(b)
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5

        >>> a.flatten()
        >>> tb.dictree(a)
        +
        |____1<--dog
        |____1<--pig<--fish<--a
        |____1<--pig<--fish<--b
        |____4<--cat
        |____5

        '''
        flatobj = flatten(self)
        remkeys = [key for key in self]
        for key in remkeys:
            del self[key]
        for key in flatobj:
            self[key] = copy.copy(flatobj[key])


def convertKeysToStr(SDobject):
    if isinstance(SDobject, SpaceData):
        newSDobject = SpaceData()
        newSDobject.attrs = SDobject.attrs
    else:
        newSDobject = {}
    for key in SDobject:
        if not isinstance(key, str):
            if isinstance(SDobject[key], dict):
                newSDobject[str(key)] = convertKeysToStr(SDobject[key])
            else:
                newSDobject[str(key)] = SDobject[key]
        else:
            if isinstance(SDobject[key], dict):
                newSDobject[key] = convertKeysToStr(SDobject[key])
            else:
                newSDobject[key] = SDobject[key]

    return newSDobject


def flatten(dobj):
    '''Collapse datamodel to one level deep

    Examples
    --------

    >>> import spacepy.datamodel as dm
    >>> import spacepy.toolbox as tb
    >>> a = dm.SpaceData()
    >>> a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
    >>> a['4'] = dm.SpaceData(cat = 'kitty')
    >>> a['5'] = 4
    >>> tb.dictree(a)
    +
    |____1
         |____dog
         |____pig
              |____fish
                   |____a
                   |____b
    |____4
         |____cat
    |____5

    >>> b = dm.flatten(a)
    >>> tb.dictree(b)
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5

    >>> a.flatten()
    >>> tb.dictree(a)
    +
    |____1<--dog
    |____1<--pig<--fish<--a
    |____1<--pig<--fish<--b
    |____4<--cat
    |____5


    See Also
    --------
    SpaceData.flatten

    '''

    try:
        addme = dobj.__class__()
    except (TypeError):
        addme = SpaceData()
    remlist = []
    for key in dobj: #iterate over keys in SpaceData
        if isinstance(dobj[key], dict):
            remlist.append(key)
            newname = str(key) + '<--'
            for levkey in dobj[key]:
                if hasattr(dobj[key][levkey], 'keys'):
                    retdict = flatten(dobj[key][levkey])
                    for key2 in retdict:
                        addme[newname+levkey+'<--'+key2] = retdict[key2]
                else:
                    addme[newname+levkey] = copy.copy(dobj[key][levkey])
        else:
            addme[key] = copy.copy(dobj[key])
    return addme

def fromCDF(fname, **kwargs):
    '''
    Create a SpacePy datamodel representation of a NASA CDF file

    Parameters
    ----------
    file : string
        the name of the cdf file to be loaded into a datamodel

    Returns
    -------
    out : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> data = dm.fromCDF('test.cdf')

    See Also
    --------
    spacepy.pycdf.CDF.copy
    '''
    #TODO: add unflatten keyword and restore flattened variables
    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and SpacePy's pyCDF")

    with pycdf.CDF(fname) as cdfdata:
        return cdfdata.copy()

def fromHDF5(fname, **kwargs):
    '''
    Create a SpacePy datamodel representation of an HDF5 file

    Parameters
    ----------
    file : string
        the name of the HDF5 file to be loaded into a datamodel

    Returns
    -------
    out : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> data = dm.fromHDF5('test.hdf')

    Notes
    -----
    Known issues -- zero-sized datasets will break in h5py
    This is kluged by returning a dmarray containing a None
    '''
    def hdfcarryattrs(SDobject, hfile, path):
        if hasattr(hfile[path],'attrs'):
            for key, value in hfile[path].attrs.iteritems():
                try:
                    SDobject.attrs[key] = value
                except:
                    warnings.warn('The following key:value pair is not permitted\n' +
                                    'key = {0} ({1})\n'.format(key, type(key)) +
                                    'value = {0} ({1})'.format(value, type(value)), DMWarning)

    try:
        import h5py as hdf
    except ImportError:
        raise ImportError('HDF5 converter requires h5py')

    if type(fname) == str:
        hfile = hdf.File(fname, mode='r')
    else:
        hfile = fname
        #should test here for HDF file object

    if 'path' not in kwargs:
        path = '/'
    else:
        path = kwargs['path']

    SDobject = SpaceData()
    allowed_elems = [hdf.Group, hdf.Dataset]
    ##carry over the attributes
    hdfcarryattrs(SDobject, hfile, path)
    ##carry over the groups and datasets
    for key, value in hfile[path].iteritems():
        try:
            if type(value) is allowed_elems[0]: #if a group
                SDobject[key] = SpaceData()
                SDobject[key] = fromHDF5(hfile, path=path+'/'+key)
            elif type(value) is allowed_elems[1]: #if a dataset
                try:
                    SDobject[key] = dmarray(value)
                except (TypeError, ZeroDivisionError): #ZeroDivisionError catches zero-sized DataSets
                    SDobject[key] = dmarray(None)
                hdfcarryattrs(SDobject[key], hfile, path+'/'+key)
        except:
            raise ValueError('HDF5 file contains type other than Group or Dataset')
    if path=='/': hfile.close()
    return SDobject

def toHDF5(fname, SDobject, **kwargs):
    '''
    Create an HDF5 file from a SpacePy datamodel representation

    Parameters
    ----------
    fname : str
        Filename to write to

    SDobject : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    overwrite : bool (optional)
        allow overwrite of an existing target file (default True)
    mode : str (optional)
        HDF5 file open mode (a, w, r) (default 'a')

    Returns
    -------
    None
    '''
    def SDcarryattrs(SDobject, hfile, path, allowed_attrs):
        if hasattr(SDobject, 'attrs'):
            for key, value in SDobject.attrs.iteritems():
                dumval, dumkey = copy.copy(value), copy.copy(key)
                if type(value) in allowed_attrs:
                    #test for datetimes in iterables
                    if hasattr(value, '__iter__'):
                        dumval = [b.isoformat() if isinstance(b, datetime.datetime) else b for b in value]
                    truth = False
                    try:
                        if value.nbytes: truth = True #empty arrays of any dimension are nbytes=0
                    except AttributeError: #not an array
                        if value or value is 0: truth = True

                    if truth:
                        if type(key) is unicode:
                            dumkey = str(key)
                        if type(value) is unicode:
                            dumval = str(value)
                        try:
                            hfile[path].attrs[dumkey] = dumval
                        except:
                            hfile[path].attrs[dumkey] = str(dumval)
                            warnings.warn('The following value is not permitted\n' +
                                    'key, value = {0} ({1})\n'.format(key, value) +
                                    'value has been converted to a string for output', DMWarning)
                    else:
                        hfile[path].attrs[dumkey] = ''
                elif isinstance(value, datetime.datetime):
                    dumval = value.isoformat()
                    if type(key) is unicode:
                        dumkey = str(key)
                    hfile[path].attrs[dumkey] = dumval
                else:
                    #TODO: add support for arrays(?) in attrs (convert to isoformat)
                    warnings.warn('The following key:value pair is not permitted\n' +
                                    'key = {0} ({1})\n'.format(key, type(key)) +
                                    'value type {0} is not in the allowed attribute list'.format(type(value)),
                                        DMWarning)

    try:
        import h5py as hdf
    except ImportError:
        raise ImportError('h5py is required to use HDF5 files')

    #mash these into a defaults dict...
    if 'mode' not in kwargs:
        wr_mo = 'a'
    else:
        wr_mo = kwargs['mode']

    if 'overwrite' not in kwargs: kwargs['overwrite'] = True
    if type(fname) == str:
        if os.path.isfile(fname) and not kwargs['overwrite']:
            raise(IOError('Cannot write HDF5, file exists (see overwrite) "{0!s}"'.format(fname)))
        if os.path.isfile(fname) and kwargs['overwrite']:
            os.remove(fname)
        hfile = hdf.File(fname, mode=wr_mo)
    else:
        hfile = fname
        #should test here for HDF file object

    if 'path' in kwargs:
        path = kwargs['path']
    else:
        path = '/'

    allowed_attrs = [int, long, float, str, unicode, numpy.ndarray, list, tuple, numpy.string_]
    for v in numpy.typecodes['AllInteger']:
        allowed_attrs.append(numpy.typeDict[v])
    for v in numpy.typecodes['AllFloat']:
        allowed_attrs.append(numpy.typeDict[v])

    allowed_elems = [SpaceData, dmarray]

    #first convert non-string keys to str
    SDobject = convertKeysToStr(SDobject)

    SDcarryattrs(SDobject,hfile,path,allowed_attrs)
    for key, value in SDobject.iteritems():
        if isinstance(value, allowed_elems[0]):
            hfile[path].create_group(key)
            toHDF5(hfile, SDobject[key], path=path+'/'+key)
        elif isinstance(value, allowed_elems[1]):
            try:
                hfile[path].create_dataset(key, data=value)
            except:
                dumval = value.copy()
                if isinstance(value[0], datetime.datetime):
                    for i, val in enumerate(value): dumval[i] = val.isoformat()
                hfile[path].create_dataset(key, data=dumval.astype('|S35'))
                #else:
                #    hfile[path].create_dataset(key, data=value.astype(float))
            SDcarryattrs(SDobject[key], hfile, path+'/'+key, allowed_attrs)
        else:
            warnings.warn('The following data is not being written as is not of an allowed type\n' +
                           'key = {0} ({1})\n'.format(key, type(key)) +
                              'value type {0} is not in the allowed data type list'.format(type(value)),
                                  DMWarning)
    if path=='/': hfile.close()


def toHTML(fname, SDobject, attrs=(),
           varLinks=False, linkFormat=None, echo=False, tableTag='<table border="1">'):
    """
    Create an HTML dump of the structure of a spacedata

    Parameters
    ----------
    fname : str
        Filename to write to

    SDobject : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    overwrite : bool (optional)
        allow overwrite of an existing target file (default True)
    mode : str (optional)
        HDF5 file open mode (a, w, r) (default 'a')
    echo : bool
        echo the html to the screen
    varLinks : bool
        make the variable name a link to a stub page

    """
    output = StringIO.StringIO() # put the output into a StringIO
    keys = SDobject.keys()
    keys.sort()

    output.write(tableTag)
    output.write('\n')
    output.write('<tr><th>{0}</th>'.format('Variable'))
    for attr in attrs:
        output.write('<th>{0}</th>'.format(attr))
    output.write('</tr>')

    for ii, key in enumerate(keys):
        if ii % 2 == 0:
            output.write('<tr>')
        else:
            output.write('<tr class="alt">')
        output.write('<td>')
        if varLinks:
            output.write('<a href="{0}.html">'.format(key))
        output.write('{0}'.format(key))
        if varLinks:
            output.write('</a>')
        output.write('</td>')

        for attr in attrs:
            try:
                if not isinstance(SDobject[key].attrs[attr], (str, unicode)):
                    tmp = str(SDobject[key].attrs[attr])
                    output.write('<td>{0}</td>'.format(_idl2html(tmp)))
                else:
                    output.write('<td>{0}</td>'.format(_idl2html(SDobject[key].attrs[attr])))
            except KeyError:
                output.write('<td></td>')
        output.write('</tr>\n')
    output.write('</table>\n')
    with open(fname, 'w') as fp:
        fp.write(output.getvalue())
    if echo:
        print(output.getvalue())
    output.close()

def _idl2html(idl):
    """
    given an idl format string for text change it to html

    Parameters
    ==========
    idl : str
        idl formated string

    Returns
    =======
    out : str
        html formatted string
    """
    html = idl
    conv = {'!!': '!',
            '!E': '<sup>',
            '!I': '<sub>'}
    while True: # hate it but for now
        ind = html.find('!')
        if ind == -1:
            break
        code = html[ind:ind+2]
        html = html.replace(code, conv[code])
        if code == '!I':
            if '!N' in html:
                html = html.replace('!N', '</sub>', 1 ) # just replace 1
            else:
                html = html + '</sub>'
        elif code == '!E':
            if '!N' in html:
                html = html.replace('!N', '</sup>', 1 ) # just replace 1
            else:
                html = html + '</sup>'
    return html

def readJSONMetadata(fname, **kwargs):
    '''Read JSON metadata from an ASCII data file

    Parameters
    ----------
    fname : str
        Filename to read metadata from

    Other Parameters
    ----------------
    verbose : bool (optional)
        set verbose output so metadata tree prints on read (default False)

    Returns
    -------
    mdata: spacepy.datamodel.SpaceData
        SpaceData with the metadata from the file
    '''
    if hasattr(fname, 'read'):
        lines = fname.read()
    else:
        with open(fname, 'r') as f:
            lines = f.read()

    # isolate header
    p_srch = re.compile(r"^#(.*)$", re.M)
    hreg = re.findall(p_srch, lines)
    header = "".join(hreg)

    # isolate JSON field
    srch = re.search( r'\{\s*(.*)\s*\}', header )
    if isinstance(srch, type(None)):
        raise IOError('The input file has no valid JSON header. Must be valid JSON bounded by braces "{ }".')
    js = srch.group(1)
    inx = js.rfind('end JSON')

    if inx == -1:
        js = ' '.join(('{', js, '}'))
        mdatadict = json.loads(js)
    else:
        js = ' '.join(('{', js[:inx]))
        mdatadict = json.loads(js)

    mdata = SpaceData()
    for key in mdatadict:
       if 'START_COLUMN' in mdatadict[key]:
           mdata[key] = SpaceData(attrs=mdatadict[key])
       elif 'VALUES' in mdatadict[key]:
           dum = mdatadict[key].pop('VALUES')
           mdata[key] = dmarray(dum, attrs=mdatadict[key])
       else:
           mdata.attrs[key] = mdatadict[key]

    if 'verbose' in kwargs:
        if kwargs['verbose']:
            mdata.tree(verbose=True, attrs=True)
    return mdata

def readJSONheadedASCII(fname, mdata=None, comment='#', convert=False):
    '''read JSON-headed ASCII data files into a SpacePy datamodel

    Parameters
    ----------
    fname : str or list
        Filename(s) to read data from

    Other Parameters
    ----------------
    mdata : spacepy.datamodel.SpaceData (optional)
        supply metadata object, otherwise is read from fname (default None)
    comment: str (optional)
        comment string in file to be read; lines starting with comment are
        ignored (default '#')
    convert: bool or dict-like (optional)
        If True, uses common names to try conversion from string. If a dict-
        like then uses the functions specified as the dict values to convert
        each element of 'key' to a non-string

    Returns
    -------
    mdata: spacepy.datamodel.SpaceData
        SpaceData with the data and metadata from the file
    '''
    import dateutil.parser as dup
    filelike = False
    try:
        if isinstance(fname, (str, unicode)):
            fname=[fname]
        elif hasattr(fname, 'readlines'):
            fname = [fname]
            filelike = True
    except NameError: # for Py3
        if isinstance(fname, str):
            fname=[fname]
    if not mdata:
        mdata = readJSONMetadata(fname[0])
    mdata_copy = dmcopy(mdata)
    def innerloop(fh, mdata, mdata_copy):
        line = fh.readline()
        while (line and line[0]==comment):
            line = fh.readline()
        fh.seek(-len(line), os.SEEK_CUR) # fixes the missing first data bug
        alldata = fh.readlines()
        if not alldata:
            return mdata
        ncols = len(alldata[0].rstrip().split())
        # fixes None in the data from empty lines at the end
        for row in xrange(len(alldata)): # reverse order
            if not alldata[-1].rstrip(): # blank line (or al white space)
                alldata.pop(-1)
            else:
                break
        nrows = len(alldata)
        data = numpy.empty((nrows, ncols), dtype=object)
        for ridx, line in enumerate(alldata):
            for cidx, el in enumerate(line.rstrip().split()):
                data[ridx, cidx] = el
        keys = mdata.keys()
        for key in keys:
            if 'START_COLUMN' in mdata_copy[key].attrs:
                st = mdata_copy[key].attrs['START_COLUMN']
                if 'DIMENSION' in mdata_copy[key].attrs:
                    varDims = numpy.array(mdata_copy[key].attrs['DIMENSION'])
                    singleDim = True
                    if len(varDims)>1 or varDims[0]>1:
                        singleDim = False
                if ('DIMENSION' in mdata_copy[key].attrs) and not singleDim:
                    en = int(mdata_copy[key].attrs['DIMENSION'][0]) + int(st)
                    try:
                        assert mdata[key]=={}
                        mdata[key] = data[:,int(st):int(en)]
                    except AssertionError:
                        mdata[key] = numpy.vstack((mdata[key], data[:,int(st):int(en)]))
                else:
                    try:
                        assert mdata[key]=={}
                        mdata[key] = data[:,int(st)]
                    except AssertionError:
                        mdata[key] = numpy.hstack((mdata[key], data[:,int(st)]))
        return mdata
    for fn in fname:
        if not filelike:
            with open(fn, 'rb') as fh: # fixes windows bug with seek()
                mdata = innerloop(fh, mdata, mdata_copy)
        else:
            mdata = innerloop(fh, mdata, mdata_copy)
    #now add the attributres to the variables
    keys = mdata.keys()
    for key in keys:
        if isinstance(mdata[key], SpaceData):
            mdata[key] = dmarray(None, attrs=mdata_copy[key].attrs)
        else:
            mdata[key] = dmarray(mdata[key], attrs=mdata_copy[key].attrs)

    if convert:
        if isinstance(convert, dict):
            conversions=convert
        else:
            conversions = {'DateTime': lambda x: dup.parse(x, ignoretz=True),
                           'ExtModel': lambda x: str(x)}
        for conkey in conversions:
            try:
                name = keys.pop(keys.index(conkey)) #remove from keylist
                for i,element in numpy.ndenumerate(mdata[name]):
                    mdata[name][i] = conversions[name](element)
            except:
                print('Key {0} for conversion not found in file'.format(conkey))
                #this should be a warning, not a print
        for remkey in keys:
            try:
                mdata[remkey] = numpy.asanyarray(mdata[remkey], dtype=float)
            except ValueError:
                pass #this will skip any unspecified string fields
    return mdata

def writeJSONMetadata(fname, insd, depend0=None, order=None, verbose=False, returnString=False):
    '''Scrape metadata from SpaceData object and make a JSON header

    Parameters
    ----------
    fname : str
        Filename to write to (can also use a file-like object)
        None can be given in conjunction with the returnString keyword to skip writing output

    insd : spacepy.datamodel.SpaceData
        SpaceData with associated attributes and variables in dmarrays

    Other Parameters
    ----------------
    depend0 : str (optional)
        variable name to use to indicate parameter on which other data depend (e.g. Time)
    order : list (optional)
        list of key names in order of start column in output JSON file
    verbose: bool (optional)
        verbose output
    returnString: bool (optional)
        return JSON header as string instead of returning None

    Returns
    -------
    None  (unless returnString keyword is True)
    '''
    js_out = {}

    def stripNL(text):
        out = text.group().replace('\n','').replace('  ','')
        return out

    #if required, identify depend0 for deciding what's data/metadata
    if depend0 is None:
        #search for DEPEND_0 in metadata
        for key in insd:
            if insd[key].attrs.has_key('DEPEND_0'):
                depend0 = insd[key].attrs['DEPEND_0']
                if not isinstance(depend0, str):
                    #assume it's a singleton list
                    depend0 = depend0[0]
                if not isinstance(depend0, str):
                    depend0 = None #Failed to get a depend0
                else:
                    break #we're done here
        if depend0 is None:
            #fall back to most common var length
            tmp, keylist = [], insd.keys()
            for key in keylist:
                tmp.append(len(insd[key]))
            depend0 = keylist[tmp.index(numpy.bincount(tmp).argmax())]
            #TODO Set using Time, or Epoch, or similar...
    else:
        if not insd.has_key(depend0): raise KeyError('Invalid key supplied for ordering metadata on write')
    datalen = len(insd[depend0])

    #start with global attrs
    #TODO: check for datetime objs in attributes
    if insd.attrs:
        glattr = _dateToISO(insd.attrs)
        for key in glattr:
            js_out[key] = dmcopy(glattr[key])
            #TODO Mark these as global somehow (by omission of some metadata?)
            try:
                js_out[key] = js_out[key].tolist()
            except:
                pass
    #collect keys and put in order for output
    #TODO first check for extant START_COLUMN
    #then check dimensionality so that start column and dims can be added, if not present
    if hasattr(order, '__iter__'):
        keylist = order
        #now make sure that all missing keys are added to end
        for key in sorted(insd.keys()):
            if key not in order: keylist.append(key)
    else:
        ##TODO do we want to have DEPEND0 first in order by default?
        keylist = sorted(insd.keys())

    idx = 0
    for key in keylist:
        js_out[key] = dmcopy(_dateToISO(insd[key].attrs))
        if len(insd[key]) == datalen: #is data
            if verbose: print('data: {0}'.format(key))
            try:
                js_out[key]['DIMENSION'] = list(insd[key].shape[1:])
            except AttributeError: #AttrErr if just metadata
                js_out[key]['DIMENSION'] = insd[key].attrs['DIMENSION']
            if not js_out[key]['DIMENSION']: js_out[key]['DIMENSION'] = [1]
            js_out[key]['START_COLUMN'] = idx
            dims = js_out[key]['DIMENSION']
            idx += int(dims[0])
            if len(dims)>1:
                l1 = 'The data cannot be properly represented in JSON-headed ASCII as it has too high a rank\n'
                l2 = 'key = {0} ({1})\n'.format(key, insd[key].shape)
                l3 = 'Maximum allowed number of dimensions is 2\n'
                warnings.warn(''.join([l1, l2, l3]), DMWarning)
        else: #is metadata
            if verbose: print('metadata: {0}'.format(key))
            js_out[key]['VALUES'] = dmcopy(_dateToISO(insd[key]))
        for kk in js_out[key]:
            try:
                js_out[key][kk] = js_out[key][kk].tolist()
            except:
                pass
    json_str = json.dumps(js_out, indent=4)
    reob = re.compile('\[.*?\]', re.DOTALL)
    json_str = re.sub(reob, stripNL, json_str) #put lists back onto one line
    #add comment field for header
    json_str = ''.join(['#', json_str])
    json_str = '\n#'.join(json_str.split('\n'))
    json_str = ''.join([json_str,'\n'])

    if isinstance(fname, str):
        with open(fname,'w') as fh:
            fh.writelines(json_str)
    elif hasattr(fname, 'writelines'):
        fname.writelines(json_str)
    elif (fname is None) and (returnString):
        return json_str

    if returnString: return json_str


def _dateToISO(indict):
    retdict = dmcopy(indict)
    if isinstance(indict, dict):
        for key in indict:
            if isinstance(indict[key], datetime.datetime):
                retdict[key] = retdict[key].isoformat()
            elif hasattr(indict[key], '__iter__'):
                for idx, el in enumerate(indict[key]):
                    if isinstance(el, datetime.datetime):
                        retdict[idx] = el.isoformat()
    else:
        if isinstance(indict, datetime.datetime):
            retdict = retdict.isoformat()
        elif hasattr(indict, '__iter__'):
            dum = numpy.asanyarray(indict)
            for idx, el in numpy.ndenumerate(indict):
                if isinstance(el, datetime.datetime):
                    retdict[idx] = el.isoformat()
    return retdict
    

def toJSONheadedASCII(fname, insd, metadata=None, depend0=None, order=None, **kwargs):
    ''' '''
    kwarg_dict = {'sep': ' '}
    for key in kwarg_dict.keys():
        if key in kwargs:
            kwarg_dict[key] = kwargs[key]
    import StringIO
    if not metadata:
        metadata = StringIO.StringIO()
        writeJSONMetadata(metadata, insd, depend0=depend0, order=order)
        metadata.seek(0) #rewind StringIO object to start
    hdr = readJSONMetadata(metadata)

    datlist = []
    for key in hdr:
        if hdr[key].attrs.has_key('START_COLUMN'):
            #add to list of (start_col, keyname) pairs
            datlist.append((hdr[key].attrs['START_COLUMN'], key, hdr[key].attrs['DIMENSION'][0]))
            #also use for data length
            datlen = len(insd[key])
            if datlen==0:
                raise ValueError('No data present to write: Use writeJSONmetadata')
                #TODO: Set this to just default to writing the header out and raise a warning
    datlist.sort()
    ncols = datlist[-1][0]+datlist[-1][-1]

    #now open file (file-like) and for each line in len(data)
    #write the line using start_column, name, dimension
    data = numpy.zeros([datlen,ncols], dtype=object)
    for stcol, name, dim in datlist:
        if dim==1:
            data[:, stcol] = _dateToISO(insd[name])
        else:
            data[:, stcol:stcol+dim] = _dateToISO(insd[name])    
    hdstr = writeJSONMetadata(None, hdr, depend0=depend0, order=order, returnString=True)
    with open(fname, 'w') as fh:
        fh.writelines(hdstr)
        for line in data:
            prline = kwarg_dict['sep'].join([str(el) for el in line])
            fh.write(''.join([prline,'\n']))



def dmcopy(dobj):
    '''Generic copy utility to return a copy of a (datamodel) object

    Parameters
    ----------
    dobj : object
        object to return a copy of

    Returns
    -------
    copy_obj: object (same type as input)
        copy of input oibject

    Examples
    --------
    >>> import spacepy.datamodel as dm
    >>> dat = dm.dmarray([2,3], attrs={'units': 'T'})
    >>> dat1 = dm.dmcopy(dat)
    >>> dat1.attrs['copy': True]
    >>> dat is dat1
    False
    >>> dat1.attrs
    {'copy': True, 'units': 'T'}
    >>> dat.attrs
    {'units': 'T'}
    '''
    if isinstance(dobj, (SpaceData, dmarray)):
        return copy.deepcopy(dobj)
    elif isinstance(dobj, numpy.ndarray):
        return numpy.copy(dobj)
    else:
        return copy.copy(dobj)
