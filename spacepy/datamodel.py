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
we can use the :py:func:`toolbox.dictree` function (which works for any dictionary-like object).

>>> import spacepy.toolbox as tb
>>> tb.dictree(mydata, attrs=True)
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

from .toolbox import dictree
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


class SpaceData(dict):
    """
    Datamodel class extending dict

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
        dictree(self, **kwargs)

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

    .. deprecated:: 0.1.3
        See :meth:`~spacepy.pycdf.CDF.copy` in :class:`~spacepy.pycdf.CDF`.

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
    '''
    #TODO: add unflatten keyword and restore flattened variables
    try:
        from spacepy import pycdf
    except ImportError:
        raise ImportError("CDF converter requires NASA CDF library and SpacePy's pyCDF")
    warnings.warn('fromCDF is deprecated, see pycdf.CDF.copy',
                  DeprecationWarning)
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
                except TypeError:
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
                    if value or value is 0:
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
            raise(IOError('Cannot write HDF5, file exists (see overwrite)'))
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

    allowed_attrs = [int, long, float, str, unicode, numpy.ndarray, list, tuple,
                     numpy.float32, numpy.float64, numpy.string_]
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
            #pretty-print config_dict
            config_dict.tree(verbose=True, attrs=True)

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
    if isinstance(fname, str):
        fname=[fname]
    if not mdata:
        mdata = readJSONMetadata(fname[0])
    mdata_copy = copy.copy(mdata)
    for fn in fname:
        with open(fn, 'rb') as fh: # fixes windows bug with seek()
            line = fh.readline()
            while line[0]==comment:
                line = fh.readline()
            fh.seek(-len(line), os.SEEK_CUR) # fixes the missing first data bug
            alldata = fh.readlines()
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
                        en = int(mdata_copy[key].attrs['DIMENSION'][0]) + int(st)
                        try:
                            mdata[key] = numpy.vstack((mdata[key], data[:,int(st):int(en)]))
                        except ValueError:
                            mdata[key] = data[:,int(st):int(en)]
                    else:
                        try:
                            assert mdata[key]=={}
                            mdata[key] = data[:,int(st)]
                        except AssertionError:
                            mdata[key] = numpy.hstack((mdata[key], data[:,int(st)]))

    #now add the attributres to the variables
    for key in keys:
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
