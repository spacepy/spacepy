#!/usr/bin/env python

"""pyCDF implementation

This module contains the implementation details of pyCDF, the
python interface to U{NASA's CDF library<http://cdf.gsfc.nasa.gov/>}.
"""

__version__ = '0.11'
__author__ = 'Jonathan Niehof <jniehof@lanl.gov>'

#Main implementation file for pycdf, automatically imported

import collections
import ctypes
import ctypes.util
import datetime
import os
import os.path
import shutil
import sys
import warnings

from . import const


try:
    str_classes = (str, bytes, unicode)
except NameError:
    str_classes = (str, bytes)


class Library(object):
    """Abstraction of the base CDF C library and its state.

    Not normally intended for end-user use. An instance of this class
    is created at package load time as the L{lib} variable, providing
    access to the underlying C library if necessary.

    Calling the C library directly requires knowledge of the
    U{ctypes<http://docs.python.org/library/ctypes.html>} package.

    L{__init__} searches for and loads the C library; details on the
    search are documented there.

    @ivar _del_middle_rec_bug: does this version of the library have a bug
                               when deleting a record from the middle of a
                               variable?
    @type _del_middle_rec_bug: boolean
    @ivar _library: C{ctypes} connection to the library
    @type _library: ctypes.WinDLL or ctypes.CDLL
    @ivar version: version of the CDF library, in order version, release,
                   increment, subincrement
    @type version: tuple
    """
    def __init__(self):
        """Load the CDF C library.

        Searches for the library in the order:
            1. Appropriately-named file in CDF_LIB
            2. Appropriately-named file in CDF_BASE
            3. Standard library search path
        @raise CDFError: BAD_DATA_TYPE if can't map types properly
        """

        if 'CDF_LIB' in os.environ:
            libdir = os.environ['CDF_LIB']
        elif 'CDF_BASE' in os.environ:
            libdir = os.path.join(os.environ['CDF_BASE'], 'lib')
        else:
            libdir = None
            libpath = None
        if libdir:
            if sys.platform == 'win32':
                libpath = os.path.join(libdir, 'cdf.dll')
            elif sys.platform == 'darwin':
                libpath = os.path.join(libdir, 'libcdf.dylib')
                if not os.path.exists(libpath):
                    libpath = os.path.join(libdir, 'cdf.dylib')
            else:
                libpath = os.path.join(libdir, 'libcdf.so')
            if not os.path.exists(libpath):
                libpath = None

        if libpath == None:
            if sys.platform == 'win32':
                libpath = ctypes.util.find_library('cdf.dll')
            else:
                libpath = ctypes.util.find_library('cdf')
            if not libpath:
                raise Exception('Cannot find CDF C library. ' + \
                                'Try os.putenv("CDF_LIB", library_directory) ' + \
                                'before import.')

        if sys.platform == 'win32':
            self._library = ctypes.WinDLL(libpath)
        else:
            self._library = ctypes.CDLL(libpath)
        self._library.CDFlib.restype = ctypes.c_long #commonly used, so set it up here
        self._library.EPOCHbreakdown.restype = ctypes.c_long
        self._library.computeEPOCH.restype = ctypes.c_double
        self._library.computeEPOCH.argtypes = [ctypes.c_long] * 7
        self._library.computeEPOCH16.restype = ctypes.c_double
        self._library.computeEPOCH16.argtypes = [ctypes.c_long] * 10 + \
            [ctypes.POINTER(ctypes.c_double * 2)]

        #Set up the dictionary for CDF type - ctypes lookup
        c_types = {}
        c_types['unsigned'] = [ctypes.c_ubyte, ctypes.c_ushort, ctypes.c_uint,
                      ctypes.c_ulong, ctypes.c_ulonglong]
        c_types['signed'] = [ctypes.c_byte, ctypes.c_short, ctypes.c_int,
                     ctypes.c_long, ctypes.c_longlong]
        c_types['float'] = [ctypes.c_float, ctypes.c_double, ctypes.c_longdouble]
        c_sizes = {}
        for i in c_types:
            c_sizes[i] = [ctypes.sizeof(j) for j in c_types[i]]
        types_wanted = {const.CDF_BYTE.value: 'signed',
                        const.CDF_CHAR.value: 'signed',
                        const.CDF_INT1.value: 'signed',
                        const.CDF_UCHAR.value: 'unsigned',
                        const.CDF_UINT1.value: 'unsigned',
                        const.CDF_INT2.value: 'signed',
                        const.CDF_UINT2.value: 'unsigned',
                        const.CDF_INT4.value: 'signed',
                        const.CDF_UINT4.value: 'unsigned',
                        const.CDF_FLOAT.value: 'float',
                        const.CDF_REAL4.value: 'float',
                        const.CDF_DOUBLE.value: 'float',
                        const.CDF_REAL8.value: 'float',
                        const.CDF_EPOCH.value: 'float',
                        }
        self.sizedict = {const.CDF_BYTE.value: 1,
                         const.CDF_CHAR.value: 1,
                         const.CDF_INT1.value: 1,
                         const.CDF_UCHAR.value: 1,
                         const.CDF_UINT1.value: 1,
                         const.CDF_INT2.value: 2,
                         const.CDF_UINT2.value: 2,
                         const.CDF_INT4.value: 4,
                         const.CDF_UINT4.value: 4,
                         const.CDF_FLOAT.value: 4,
                         const.CDF_REAL4.value: 4,
                         const.CDF_DOUBLE.value: 8,
                         const.CDF_REAL8.value: 8,
                         const.CDF_EPOCH.value: 8,
                         }
        self.ctypedict = {}
        for i in types_wanted:
            type_wanted = types_wanted[i]
            size_wanted = self.sizedict[i]
            try:
                self.ctypedict[i] = c_types[type_wanted][
                    c_sizes[type_wanted].index(size_wanted)
                    ]
            except ValueError:
                raise CDFError(const.BAD_DATA_TYPE)
        #EPOCH16 is a double[2] (NOT doubledouble) and needs special handling
        self.sizedict[const.CDF_EPOCH16.value] = 16
        self.ctypedict[const.CDF_EPOCH16.value] = self.ctypedict[
            const.CDF_EPOCH.value] * 2
        self.cdftypenames = {const.CDF_BYTE.value: 'CDF_BYTE',
                             const.CDF_CHAR.value: 'CDF_CHAR',
                             const.CDF_INT1.value: 'CDF_INT1',
                             const.CDF_UCHAR.value: 'CDF_UCHAR',
                             const.CDF_UINT1.value: 'CDF_UINT1',
                             const.CDF_INT2.value: 'CDF_INT2',
                             const.CDF_UINT2.value: 'CDF_UINT2',
                             const.CDF_INT4.value: 'CDF_INT4',
                             const.CDF_UINT4.value: 'CDF_UINT4', 
                             const.CDF_FLOAT.value: 'CDF_FLOAT',
                             const.CDF_REAL4.value: 'CDF_REAL4',
                             const.CDF_DOUBLE.value: 'CDF_DOUBLE',
                             const.CDF_REAL8.value: 'CDF_REAL8',
                             const.CDF_EPOCH.value: 'CDF_EPOCH',
                             const.CDF_EPOCH16.value: 'CDF_EPOCH16',
                             }

        #Get CDF version information
        ver = ctypes.c_long(0)
        rel = ctypes.c_long(0)
        inc = ctypes.c_long(0)
        sub = ctypes.c_char(' ')
        self.call(const.GET_, const.LIB_VERSION_, ctypes.byref(ver),
                  const.GET_, const.LIB_RELEASE_, ctypes.byref(rel),
                  const.GET_, const.LIB_INCREMENT_, ctypes.byref(inc),
                  const.GET_, const.LIB_subINCREMENT_, ctypes.byref(sub))
        ver = ver.value
        rel = rel.value
        inc = inc.value
        sub = sub.value
        self.version = (ver, rel, inc, sub)
        self._del_middle_rec_bug = ver < 3 or (ver == 3 and
                                               (rel < 3 or
                                                (rel == 3 and inc < 1)))

    def check_status(self, status, ignore=()):
        """Raise exception or warning based on return status of CDF call

        @param status: status returned by the C library, equivalent to C{CDFStatus}
        @type status: int
        @param ignore: CDF statuses to ignore. If any of these
                       is returned by CDF library, any related warnings or
                       exceptions will I{not} be raised. (Default none).
        @type ignore: sequence of ctypes.c_long
        @raise CDFError: if status < CDF_WARN, indicating an error
        @raise CDFWarning: if CDF_WARN <= status < CDF_OK, indicating a warning,
                           I{and} interpreter is set to error on warnings.
        @return: L{status} (unchanged)
        @rtype: int
        """
        if status == const.CDF_OK or status in ignore:
            return status
        if status < const.CDF_WARN:
            raise CDFError(status)
        else:
            warning = CDFWarning(status)
            warning.warn()
            return status

    def call(self, *args, **kwargs):
        """Call the CDF internal interface

        Passes all parameters directly through to the CDFlib routine of the
        CDF library's C internal interface. Checks the return value with
        L{check_status}.

        @param args: Passed directly to the CDF library interface. Useful
                     constants are defined in the L{const} module of this package.
        @type args: various, see C{ctypes}.
        @keyword ignore: sequence of CDF statuses to ignore. If any of these
                         is returned by CDF library, any related warnings or
                         exceptions will I{not} be raised.
        @return: CDF status from the library
        @rtype: int
        @note: Terminal NULL_ is automatically added to L{args}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        if 'ignore' in kwargs:
            return self.check_status(self._library.CDFlib(
                *(args + (const.NULL_, ))
                ), kwargs['ignore'])
        else:
            return self.check_status(self._library.CDFlib(
                *(args + (const.NULL_, ))
                ))

    def epoch_to_datetime(self, epoch):
        """Converts a CDF epoch value to a datetime

        @param epoch: epoch value from CDF
        @type epoch: float
        """
        yyyy = ctypes.c_long(0)
        mm = ctypes.c_long(0)
        dd = ctypes.c_long(0)
        hh = ctypes.c_long(0)
        min = ctypes.c_long(0)
        sec = ctypes.c_long(0)
        msec = ctypes.c_long(0)
        self._library.EPOCHbreakdown(ctypes.c_double(epoch),
                                     ctypes.byref(yyyy), ctypes.byref(mm),
                                     ctypes.byref(dd),
                                     ctypes.byref(hh), ctypes.byref(min),
                                     ctypes.byref(sec), ctypes.byref(msec))
        return datetime.datetime(yyyy.value, mm.value, dd.value,
                                 hh.value, min.value, sec.value,
                                 msec.value * 1000)

    def datetime_to_epoch(self, dt):
        """Converts a Python datetime to a CDF Epoch value

        @param dt: date and time to convert
        @type dt: datetime.datetime
        @return: epoch corresponding to L{dt}
        @rtype: float
        """
        if dt.tzinfo != None and dt.utcoffset() != None:
            dt = dt - dt.utcoffset()
        dt.replace(tzinfo=None)
        micro = dt.microsecond % 1000
        if micro >= 500:
            dt += datetime.timedelta(0, 0, 1000)
        return self._library.computeEPOCH(dt.year, dt.month, dt.day, dt.hour,
                                          dt.minute, dt.second,
                                          int(dt.microsecond / 1000))

    def epoch16_to_datetime(self, epoch):
        """Converts a CDF epoch16 value to a datetime

        @param epoch: epoch16 value from CDF
        @type epoch: list of two floats
        @raise EpochError: if input invalid
        """
        try:
            if len(epoch) != 2:
                raise EpochError('EPOCH16 values must be a pair of doubles.')
        except TypeError:
            raise EpochError('EPOCH16 values must be a pair of doubles.')
        yyyy = ctypes.c_long(0)
        mm = ctypes.c_long(0)
        dd = ctypes.c_long(0)
        hh = ctypes.c_long(0)
        min = ctypes.c_long(0)
        sec = ctypes.c_long(0)
        msec = ctypes.c_long(0)
        usec = ctypes.c_long(0)
        nsec = ctypes.c_long(0)
        psec = ctypes.c_long(0)
        self._library.EPOCH16breakdown((ctypes.c_double * 2)(*epoch),
                                     ctypes.byref(yyyy), ctypes.byref(mm),
                                     ctypes.byref(dd),
                                     ctypes.byref(hh), ctypes.byref(min),
                                     ctypes.byref(sec), ctypes.byref(msec),
                                     ctypes.byref(usec), ctypes.byref(nsec),
                                     ctypes.byref(psec))
        micro = int(float(msec.value) * 1000 + float(usec.value) +
                    float(nsec.value) / 1000 + float(psec.value) / 1e6 + 0.5)
        if micro < 1000000:
            return datetime.datetime(yyyy.value, mm.value, dd.value,
                                     hh.value, min.value, sec.value,
                                     micro)
        else:
            add_sec = int(micro / 1000000)
            try:
                return datetime.datetime(yyyy.value, mm.value, dd.value,
                                         hh.value, min.value, sec.value,
                                         micro - add_sec * 1000000) + \
                                         datetime.timedelta(seconds=add_sec)
            except OverflowError:
                return datetime.datetime(datetime.MAXYEAR, 12, 31,
                                         23, 59, 59,
                                         999999)

    def datetime_to_epoch16(self, dt):
        """Converts a Python datetime to a CDF Epoch16 value

        @param dt: date and time to convert
        @type dt: datetime.datetime
        @return: epoch16 corresponding to L{dt}
        @rtype: list of float
        """
        if dt.tzinfo != None and dt.utcoffset() != None:
            dt = dt - dt.utcoffset()
        dt.replace(tzinfo=None)
        epoch16 = (ctypes.c_double * 2)(0.0, 0.0)
        self._library.computeEPOCH16(dt.year, dt.month, dt.day, dt.hour,
                                     dt.minute, dt.second,
                                     int(dt.microsecond / 1000),
                                     dt.microsecond % 1000, 0, 0,
                                     epoch16)
        return [epoch16[0], epoch16[1]]


lib = Library()
"""Module global library object.

Initalized at module load time so all classes have ready
access to the CDF library and a common state.
"""


class CDFException(Exception):
    """Base class for errors or warnings in the CDF library.

    Not normally used directly (see subclasses L{CDFError} and L{CDFWarning}).

    Error messages provided by this class are looked up from the underlying
    C library.

    @ivar status: CDF library status code
    @type status: ctypes.c_long
    @ivar string: CDF library error message for L{status}
    @type string: string
    """

    def __init__(self, status):
        """Create a CDF Exception

        Uses CDF C library to look up an appropriate error messsage.

        @param status: CDF status
        @type status: ctypes.c_long
        """
        self.status = status
        self.string = 'CDF error ' + repr(status) + ', unable to get details.'
        message = ctypes.create_string_buffer(const.CDF_STATUSTEXT_LEN + 1)
        try:
            retval = lib._library.CDFlib(const.SELECT_, const.CDF_STATUS_,
                                         ctypes.c_long(status),
                                         const.GET_, const.STATUS_TEXT_, message,
                                         const.NULL_)
            if retval == const.CDF_OK:
                if isinstance(message.value, str):
                    self.string = message.value
                elif isinstance(message.value, bytes):
                    self.string = message.value.decode()
        except:
            pass

    def __str__(self):
        """Error string associated with the library error.

        @return: Error message from the CDF library.
        @rtype: string
        """
        return self.string


class CDFError(CDFException):
    """Raised for an error in the CDF library."""
    pass


class CDFWarning(CDFException, UserWarning):
    """Used for a warning in the CDF library."""

    def warn(self, level=4):
        """Issues a warning based on the information stored in my exception

        Intended for use in check_status or similar wrapper function.

        @param level: optional (default 3), how far up the stack the warning
                      should be reported. Passed directly to C{warnings.warn}.
        @type level: int
        """
        warnings.warn(self, self.__class__, level)


class EpochError(Exception):
    """Used for errors in epoch routines"""
    pass


def _compress(obj, comptype=None, param=None):
    """Set or check the compression of a L{CDF} or L{Var}

    @param obj: object on which to set or check compression
    @type obj: L{CDF} or L{Var}
    @param comptype: type of compression to change to, see CDF C reference
                     manual section 4.10. Constants for this parameter
                     are in L{const}. If not specified, will not change
                     compression.
    @type comptype: ctypes.c_long
    @param param: Compression parameter, see CDF CRM 4.10 and L{const}.
                  If not specified, will choose reasonable default (5 for
                  gzip; other types have only one possible parameter.)
    @type param: ctypes.c_long
    @return: (comptype, param) currently in effect
    @rtype: tuple
    """
    if isinstance(obj, CDF):
        COMPRESSION_ = const.CDF_COMPRESSION_
    elif isinstance(obj, Var):
        COMPRESSION_ = const.zVAR_COMPRESSION_
    else:
        raise ValueError('Must specify a CDF or Var type.')

    validparams = {const.NO_COMPRESSION.value: [ctypes.c_long(0)],
                   const.RLE_COMPRESSION.value: [const.RLE_OF_ZEROs],
                   const.HUFF_COMPRESSION.value:
                       [const.OPTIMAL_ENCODING_TREES],
                   const.AHUFF_COMPRESSION.value:
                       [const.OPTIMAL_ENCODING_TREES],
                   const.GZIP_COMPRESSION.value: [ctypes.c_long(5),
                                                  ctypes.c_long(1),
                                                  ctypes.c_long(2),
                                                  ctypes.c_long(3),
                                                  ctypes.c_long(4),
                                                  ctypes.c_long(6),
                                                  ctypes.c_long(7),
                                                  ctypes.c_long(8),
                                                  ctypes.c_long(9),
                                                  ],
                   }
    comptypes = [const.NO_COMPRESSION, const.RLE_COMPRESSION,
                 const.HUFF_COMPRESSION, const.AHUFF_COMPRESSION,
                 const.GZIP_COMPRESSION]
    comptypevalues = [i.value for i in comptypes]

    if comptype != None:
        if not hasattr(comptype, 'value'):
            comptype = ctypes.c_long(comptype)
        if param == None:
            if not comptype.value in validparams:
                raise CDFError(const.BAD_COMPRESSION)
            param = validparams[comptype.value][0]
        paramlist = (ctypes.c_long * 1)(param)
        obj._call(const.PUT_, COMPRESSION_,
                   comptype, paramlist)
    params = (ctypes.c_long *
              const.CDF_MAX_PARMS)(*([0] * const.CDF_MAX_PARMS))
    comptype = ctypes.c_long(0)
    percent = ctypes.c_long(0)
    obj._call(const.GET_, COMPRESSION_,
               ctypes.byref(comptype), ctypes.byref(params),
               ctypes.byref(percent))
    param = params[0]
    if not comptype.value in comptypevalues:
        raise CDFError(const.BAD_COMPRESSION)
    validparamvalues = [i.value for i in validparams[comptype.value]]
    if not param in validparamvalues:
        raise CDFError(const.BAD_COMPRESSION_PARM)
    comptype = comptypes[comptypevalues.index(comptype.value)]
    if comptype in (const.RLE_COMPRESSION, const.HUFF_COMPRESSION,
                    const.AHUFF_COMPRESSION):
        param = validparams[comptype.value][validparamvalues.index(param)]
    return (comptype, param)


class CDF(collections.Mapping):
    """Python object representing a CDF file.

    Opening existing
    ================
    Open a CDF by creating a CDF object, e.g.:
    C{cdffile = pycdf.CDF('cdf_filename.cdf')}
    Be sure to L{close} or L{save} when done.

    CDF supports the U{with
    <http://docs.python.org/tutorial/inputoutput.html#methods-of-file-objects>}
    keyword, like other file objects, so::
        with pycdf.CDF('cdf_filename.cdf') as cdffile:
            do brilliant things with the CDF
    will open the CDF, execute the indented statements, and close the CDF when
    finished or when an error occurs. The U{python docs
    <http://docs.python.org/reference/compound_stmts.html#with>} include more
    detail on this 'context manager' ability.

    CDF objects behave like a python
    U{dictionary
    <http://docs.python.org/tutorial/datastructures.html#dictionaries>},
    where the keys are names of variables in the CDF, and the values,
    L{Var} objects. As a dictionary, they are also U{iterable
    <http://docs.python.org/tutorial/classes.html#iterators>} and it is easy
    to loop over all of the variables in a file. Some examples:
      1. List the names of all variables in the open CDF cdffile::
         cdffile.keys()
      2. Get a L{Var} object corresponding to the variable named Epoch::
         epoch = cdffile['Epoch']
      3. Determine if a CDF contains a variable named B_GSE::
         if 'B_GSE' in cdffile:
           print 'B_GSE is in the file'
         else:
           print 'B_GSE is not in the file'
      4. Find how many variables are in the file::
         print len(cdffile)
      5. Open the CDF named C{cdf_filename.cdf}, read I{all} the data from all
         variables into it, and close it when done or if an error occurs::
           with pycdf.CDF('cdf_filename.cdf') as cdffile:
                data = cdffile.copy()
    This last example can be very inefficient as it reads the entire CDF.
    Normally it's better to treat the CDF as a dictionary and access only
    the data needed, which will be pulled transparently from disc. See
    L{Var} for more subtle examples.

    Potentially useful dictionary methods and related functions:
      - U{in<http://docs.python.org/reference/expressions.html#in>}
      - U{keys<http://docs.python.org/tutorial/datastructures.html#dictionaries>}
      - U{len<http://docs.python.org/library/functions.html#len>}
      - U{list comprehensions
        <http://docs.python.org/tutorial/datastructures.html#list-comprehensions>}
      - U{sorted<http://docs.python.org/library/functions.html#sorted>}
    See also toolbox.dictree in SpacePy.

    The L{attrs} Python attribute acts as a dictionary referencing CDF
    attributes (do not confuse the two); all the dictionary methods above
    also work on the attribute dictionary. See L{gAttrList} for more on the
    dictionary of global attributes.

    Creating New
    ============
    Creating a new CDF from a master (skeleton) CDF has similar syntax to
    opening one::
      cdffile = pycdf.CDF('cdf_filename.cdf', 'master_cdf_filename.cdf')
    This creates and opens C{cdf_filename.cdf} as a copy of
    C{master_cdf_filename.cdf}

    Using a skeleton CDF is recommended over making a CDF entirely from
    scratch, but this is possible by specifying a blank master::
      cdffile = pycdf.CDF('cdf_filename.cdf', '')

    When CDFs are created in this way, they are opened read-write, see
    L{readonly} to change.

    Add variables by direct assignment, which will automatically set type
    and dimension based on the data provided::
      cdffile['new_variable_name'] = [1, 2, 3, 4]
    or, if more control is needed over the type and dimensions, use L{new}.

    @ivar _handle: file handle returned from CDF library open functions.
    @type _handle: ctypes.c_void_p
    @ivar pathname: filename of the CDF file
    @type pathname: string
    @ivar attrs: All global attributes for this CDF
    @type attrs: L{gAttrList}
    @note: CDF is opened read-only by default, see L{readonly} to change.
    """

    def __init__(self, pathname, masterpath=None):
        """Open or create a CDF file.

        @param pathname: name of the file to open or create
        @type pathname: string
        @param masterpath: name of the master CDF file to use in creating
                           a new file. If not provided, an existing file is
                           opened; if provided but evaluates to C{False}
                           (e.g. C{''}), an empty new CDF is created.
        @type masterpath: string
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        self.pathname = pathname.encode()
        self._handle = ctypes.c_void_p(None)
        if masterpath == None:
            self._open()
        elif masterpath:
            self._from_master(masterpath.encode())
        else:
            self._create()
        lib.call(const.SELECT_, const.CDF_zMODE_, ctypes.c_long(2))
        self.attrs = gAttrList(self)

    def __del__(self):
        """Destructor

        Close CDF file if there is still a valid handle.
        @note: To avoid data loss, explicitly call L{close} or L{save}.
        """
        if self._handle:
            self.close()

    def __delitem__(self, name):
        """Delete a zVariable in this CDF, by name or number

        @param name: Name or number of the CDF variable
        @type name: string or int
        @note: variable numbers may change if variables are added or removed.
        """
        self[name]._delete()

    def __enter__(self):
        """Context manager entrance function."""
        return self

    def __exit__(self, type, value, traceback):
        """Context manager exit function.

        Close CDF file.
        """
        self.close()

    def __getitem__(self, name):
        """Gets a zVariable in this CDF, by name or number

        The CDF acts like a dict

        @param name: Name or number of the CDF variable
        @type name: string or int
        @return: CDF variable named or numbered L{name}
        @rtype: L{Var}
        @raise KeyError: for pretty much any problem in lookup
        @note: variable numbers may change if variables are added or removed.
        """
        try:
            return Var(self, name)
        except CDFException:
            (type, value, traceback) = sys.exc_info()
            raise KeyError(name + ': ' + str(value))

    def __setitem__(self, name, data):
        """Writes data to a zVariable in this CDF

        If the zVariable does not exist, will create one matching
        L{data}. If it does exist, will attempt to write L{data}
        to it without changing the type or dimensions.

        @param name: name or number of the variable to write
        @type name: str or int
        @param data: data to write, or a L{Var} to copy
        """
        if isinstance(data, Var):
            self.clone(data, name)
        elif name in self:
            self[name][...] = data
        else:
            self.new(name, data)

    def __iter__(self, current = 0):
        """Iterates over zVars in CDF

        Iterators for dicts return keys
        @note: Returned in variable-number order
        """
        while current < self.__len__():
            name = self[current].name()
            value = (yield name)
            if value == None:
                current += 1
            else:
                current = self[value]._num()
                current += 1

    def __len__(self):
        """Implements 'length' of CDF (number of zVars)

        @return: number of zVars in the CDF
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.CDF_NUMzVARS_, ctypes.byref(count))
        return count.value

    def __contains__(self, key):
        """Determines whether a particular variable name is in the CDF

        @note: Essentially an efficiency function; L{__iter__} is called
               if this isn't defined
        @param key: key/variable name to check
        @type key: string
        @return: True if L{key} is the name of a variable in CDF, else False
        @rtype: boolean
        """
        try:
            foo = self[key]
            return True
        except KeyError:
            (type, value, traceback) = sys.exc_info()
            expected = "'" + str(key) + \
               ": NO_SUCH_VAR: Named variable not found in this CDF.'"
            if str(value) == expected:
                return False
            raise

    def __repr__(self):
        """Returns representation of CDF

        Cannot return anything that can be eval'd to create a copy of the
        CDF, so just wrap the informal representation in angle brackets.
        @return: all the data in this list of attributes
        @rtype: str
        """
        return '<CDF:\n' + str(self) + '\n>'

    def __str__(self):
        """Returnss a string representation of the CDF

        This is an 'informal' representation in that it cannot be evaluated
        directly to create a L{CDF}, just the names, types, and sizes of all
        variables. (Attributes are not listed.)

        @return: description of the variables in the CDF
        @rtype: str
        """
        return '\n'.join([key + ': ' + str(value)
                          for (key, value) in self.items()])

    def _open(self):
        """Opens the CDF file (called on init)

        Will open an existing CDF file read/write.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        @note: Not intended for direct call; pass parameters to L{CDF}
               constructor.
        """

        lib.call(const.OPEN_, const.CDF_, self.pathname, ctypes.byref(self._handle))
        self.readonly(True)

    def _create(self):
        """Creates (and opens) a new CDF file

        Created at L{pathname}.
        Assumes zero-dimension r variables
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        @note: Not intended for direct call; pass parameters to L{CDF}
               constructor.
        """

        lib.call(const.CREATE_, const.CDF_, self.pathname, ctypes.c_long(0),
                              (ctypes.c_long * 1)(0), ctypes.byref(self._handle))

    def _from_master(self, master_path):
        """Creates a new CDF from a master CDF file

        L{master_path} is copied to L{pathname} and opened.

        @param master_path: location of the master CDF file
        @type master_path: string
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        @note: Not intended for direct call; pass parameters to L{CDF}
               constructor.
        """

        if os.path.exists(self.pathname):
            raise CDFError(const.CDF_EXISTS)
        shutil.copy2(master_path, self.pathname)
        self._open()

    def _call(self, *args, **kwargs):
        """Select this CDF as current and call the CDF internal interface

        Adds call to select this CDF to L{args} and passes all parameters
        directly through to the CDFlib routine of the CDF library's C internal
        interface. Checks the return value with L{Library.check_status}.

        @param args: Passed directly to the CDF library interface. Useful
                     constants are defined in the L{const} module of this package.
        @type args: various, see C{ctypes}.
        @return: CDF status from the library
        @rtype: ctypes.c_long
        @note: Terminal NULL_ is automatically added to L{args}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        return lib.call(const.SELECT_, const.CDF_, self._handle,
                        *args, **kwargs)

    def clone(self, zVar, name=None, data=True):
        """Clone a zVariable (from another CDF or this) into this CDF

        @param zVar: variable to clone
        @type zVar: L{Var}
        @param name: Name of the new variable (default: name of the original)
        @type name: str
        @param data: Copy data, or only type, dimensions, variance, attributes?
                     (default: copy data as well)
        @type data: bool
        """
        if name == None:
            name = zVar.name()
        if name in self:
            del self[name]
        self.new(name, cdftype=zVar.type(), recVary=zVar.rv(),
                 dimVarys=zVar.dv(), dims=zVar._dim_sizes(),
                 n_elements=zVar._nelems())
        self[name].compress(*zVar.compress())
        self[name].attrs.clone(zVar.attrs)
        if data:
            self[name][...] = zVar[...]

    def col_major(self, new_col=None):
        """Finds the majority of this CDF file

        @param new_col: Specify True to change to column-major, False to
                        change to row major, or do not specify
                        to leave majority alone (check only)
        @type new_col: bool
        @returns: True if column-major, false if row-major
        @rtype: bool
        """
        if new_col != None:
            new_maj = const.COLUMN_MAJOR if new_col else const.ROW_MAJOR
            self._call(const.PUT_, const.CDF_MAJORITY_, new_maj)
        maj = ctypes.c_long(0)
        self._call(const.GET_, const.CDF_MAJORITY_, ctypes.byref(maj))
        if not maj.value in (const.ROW_MAJOR.value, const.COLUMN_MAJOR.value):
            raise CDFError(const.BAD_MAJORITY)
        return maj.value == const.COLUMN_MAJOR.value

    def readonly(self, ro=None):
        """Sets or check the readonly status of this CDF

        @param ro: True to set the CDF readonly,
                   False to set it read/write,
                   otherwise to check
        @type ro: boolean
        @return: True if CDF is read-only, else False
        @rtype: boolean
        @note: If the CDF has been changed since opening,
               setting readonly mode will have no effect.
        """
        if ro == True:
            self._call(const.SELECT_, const.CDF_READONLY_MODE_,
                       const.READONLYon)
        elif ro == False:
            self._call(const.SELECT_, const.CDF_READONLY_MODE_,
                       const.READONLYoff)
        mode = ctypes.c_long(0)
        self._call(const.CONFIRM_, const.CDF_READONLY_MODE_,
                   ctypes.byref(mode))
        if mode.value == const.READONLYon.value:
            return True
        elif mode.value == const.READONLYoff.value:
            return False
        else:
            raise CDFError(const.BAD_READONLY_MODE.value)

    def checksum(self, new_val=None):
        """Set or check the checksum status of this CDF

        @param new_val: True to enable checksum, False to disable, or leave out
                        to simply check.
        @type new_val: bool
        @return: True if the checksum is enabled or False if disabled
        @rtype: bool
        """
        if new_val != None:
            self._call(const.PUT_, const.CDF_CHECKSUM_,
                       const.MD5_CHECKSUM if new_val else const.NO_CHECKSUM)
        chk = ctypes.c_long(0)
        self._call(const.GET_, const.CDF_CHECKSUM_, ctypes.byref(chk))
        if not chk.value in (const.MD5_CHECKSUM.value,
                             const.NO_CHECKSUM.value):
            raise CDFError(const.BAD_CHECKSUM)
        return chk.value == const.MD5_CHECKSUM.value

    def close(self):
        """Closes the CDF file

        @note: Although called on object destruction (from L{__del__}),
               to ensure all data are saved the user should explicitly call
               L{close} or L{save}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        self._call(const.CLOSE_, const.CDF_)
        self._handle = ctypes.c_void_p(None)

    def compress(self, comptype=None, param=None):
        """Set or check the compression of this CDF

        Sets compression on entire I{file}, not per-variable
        (see L{Var.compress}).

        @param comptype: type of compression to change to, see CDF C reference
                         manual section 4.10. Constants for this parameter
                         are in L{const}. If not specified, will not change
                         compression.
        @type comptype: ctypes.c_long
        @param param: Compression parameter, see CDF CRM 4.10 and L{const}.
                      If not specified, will choose reasonable default (5 for
                      gzip; other types have only one possible parameter.)
        @type param: ctypes.c_long
        @return: (comptype, param) currently in effect
        @rtype: tuple
        """
        return _compress(self, comptype, param)

    def new(self, name, data=None, cdftype=None, recVary=True, dimVarys=None,
            dims=None, n_elements=None):
        """Create a new zVariable in this CDF

        @param name: name of the new variable
        @type name: str
        @param data: data to store in the new variable
        @param cdftype: CDF type of the variable, from L{const}
        @type cdftype: ctypes.c_long
        @param recVary: record variance of the variable
        @type recVary: bool
        @param dimVarys: dimension variance of each dimension
        @type dimVarys: list of bool
        @param dims: size of each dimension of this variable
        @type dims: list of int
        @param n_elements: number of elements, should be 1 except for
                           CDF_CHAR, for which it's the length of the string.
        @type n_elements: int
        @return: the newly-created zVariable
        @rtype: L{Var}
        @raise ValueError: if neither data nor sufficient typing information
                           is provided.
        """
        if data == None:
            if cdftype == None:
                raise ValueError('Must provide either data or a CDF type.')
            if dims == None:
                raise ValueError('Must provide either data or dimension list.')
            if n_elements == None:
                n_elements = 1
        else:
            (guess_dims, guess_types, guess_elements) = _Hyperslice.types(data)
            if dims == None:
                if recVary:
                    dims = guess_dims[1:]
                else:
                    dims = guess_dims
            if cdftype == None:
                cdftype = guess_types[0]
            if n_elements == None:
                n_elements = guess_elements
        if dimVarys == None:
            dimVarys = [True for i in dims]
        recVary = const.VARY if recVary else const.NOVARY
        dimVarys = [const.VARY if dimVary else const.NOVARY
                    for dimVary in dimVarys]
        if not hasattr(cdftype, 'value'):
            cdftype = ctypes.c_long(cdftype)
        new_var = Var(self, name, cdftype, n_elements, dims, recVary, dimVarys)
        if data != None:
            new_var[...] = data
        return new_var

    def save(self):
        """Saves the CDF file but leaves it open.

        If closing the CDF, L{close} is sufficient, there is no need to call
        C{save} before C{close}.

        @note: Relies on an undocumented call of the CDF C library, which
               is also used in the Java interface.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        self._call(const.SAVE_, const.CDF_)

    def copy(self):
        """Make a copy of all data and attributes in this CDF

        @return: dict of all data
        @rtype: L{CDFCopy}
        """
        return CDFCopy(self)


class CDFCopy(dict):
    """A copy of all data and attributes in a L{CDF}

    Data are L{VarCopy} objects, keyed by variable name (i.e.
    data are accessed much like from a L{CDF}).

    @ivar attrs: attributes for the CDF
    @type attrs: dict
    """

    def __init__(self, cdf):
        """Copies all data and attributes from a CDF

        @param cdf: CDF to take data from
        @type cdf: L{CDF}
        """
        self.attrs = cdf.attrs.copy()
        super(CDFCopy, self).__init__((key, var.copy())
                                      for (key, var) in cdf.items())


class Var(collections.Sequence):
    """A CDF variable.

    Reading
    =======
    This object does not directly store the data from the CDF; rather,
    it provides access to the data in a format that looks like a Python
    list. General list information is available in the python docs:
    U{1<http://docs.python.org/tutorial/introduction.html#lists>},
    U{2<http://docs.python.org/tutorial/datastructures.html#more-on-lists>},
    U{3<http://docs.python.org/library/stdtypes.html#typesseq>}.

    A record-varying variable's data are viewed as a hypercube of dimensions
    n_dims+1 and are indexed in row-major fashion, i.e. the last
    index changes most frequently / is contiguous in memory. If
    the CDF is column-major, the data are transformed to row-major
    before return. The extra dimension is the record number.

    Non record-varying variables are similar, but do not have the extra
    dimension of record number.

    Variables can be subscripted by a multidimensional index to return the
    data. Indices are in row-major order with the first dimension
    representing the record number. If the CDF is column major,
    the data are reordered to row major. Each dimension is specified
    by standard Python
    U{slice<http://docs.python.org/tutorial/introduction.html#strings>}
    notation, with dimensions separated by commas. The ellipsis fills in
    any missing dimensions with full slices. The returned data are
    lists; Python represents multidimensional arrays as nested lists.
    The innermost set of lists represents contiguous data.

    Degenerate dimensions are 'collapsed', i.e. no list of only one
    element will be returned if a single subscript is specified
    instead of a range. (To avoid this, specify a slice like 1:2,
    which starts with 1 and ends before 2).

    Two special cases:
      1. requesting a single-dimension slice for a
         record-varying variable will return all data for that
         record number (or those record numbers) for that variable.
      2. Requests for multi-dimensional variables may skip the record-number
         dimension and simply specify the slice on the array itself. In that
         case, the slice of the array will be returned for all records.
    In the event of ambiguity (i.e. single-dimension slice on a one-dimensional
    variable), case 1 takes priority.
    Otherwise, mismatch between the number of dimensions specified in
    the slice and the number of dimensions in the variable will cause
    an IndexError to be thrown.

    This all sounds very complicated but it's essentially attempting
    to do the 'right thing' for a range of slices.

    As a list type, variables are also U{iterable
    <http://docs.python.org/tutorial/classes.html#iterators>}; iterating
    over a variable returns a single complete record at a time.

    This is all clearer with examples. Consider a variable B_GSM, with
    three elements per record (x, y, z components) and fifty records in
    the CDF. Then:
      1. C{B_GSM[0, 1]} is the y component of the first record.
      2. C{B_GSM[10, :]} is a three-element list, containing x, y, and z
         components of the 11th record. As a shortcut, if only one dimension
         is specified, it is assumed to be the record number, so this
         could also be written C{B_GSM[10]}.
      3. C{B_GSM[...]} reads all data for B_GSM and returns it as a
         fifty-element list, each element itself being a three-element
         list of x, y, z components.

    Multidimensional example: consider fluxes stored as a function of
    pitch angle and energy. Such a variable may be called Flux and
    stored as a two-dimensional array, with the first dimension
    representing (say) ten energy steps and the second, eighteen
    pitch angle bins (ten degrees wide, centered from 5 to 175 degrees).
    Assume 100 records stored in the CDF (i.e. 100 different times).
      1. C{Flux[4]} is a list of ten elements, one per energy step,
         each element being a list of 18 fluxes, one per pitch bin.
         All are taken from the fifth record in the CDF.
      2. C{Flux[4, :, 0:4]} is the same record, all energies, but
         only the first four pitch bins (roughly, field-aligned).
      3. C{Flux[..., 0:4]} is a 100-element list (one per record),
         each element being a ten-element list (one per energy step),
         each containing fluxes for the first four pitch bins.
    This slicing notation is very flexible and allows reading
    specifically the desired data from the CDF.

    All data are, on read, converted to appropriate Python data
    types; EPOCH and EPOCH16 types are converted to
    U{datetime<http://docs.python.org/library/datetime.html>}.

    Potentially useful list methods and related functions:
      - U{count<http://docs.python.org/tutorial/datastructures.html#more-on-lists>}
      - U{in<http://docs.python.org/reference/expressions.html#in>}
      - U{index<http://docs.python.org/tutorial/datastructures.html#more-on-lists>}
      - U{len<http://docs.python.org/library/functions.html#len>}
      - U{list comprehensions
        <http://docs.python.org/tutorial/datastructures.html#list-comprehensions>}
      - U{sorted<http://docs.python.org/library/functions.html#sorted>}

    The topic of array majority can be very confusing; good background material
    is available at U{IDL Array Storage and Indexing
    <http://www.dfanning.com/misc_tips/colrow_major.html>}. In brief,
    I{regardless of the majority stored in the CDF}, pycdf will always present
    the data in the native Python majority, row-major order, also known as
    C order. This is the default order in U{numPy
    <http://docs.scipy.org/doc/numpy/reference/arrays.ndarray.html
    #internal-memory-layout-of-an-ndarray>}.
    However, packages that render image data may expect it in column-major
    order. If the axes seem 'swapped' this is likely the reason.

    The L{attrs} Python attribute acts as a dictionary referencing zAttributes
    attributes (do not confuse the two); all the dictionary methods above
    also work on the attribute dictionary. See L{zAttrList} for more on the
    dictionary of global attributes.

    Writing
    =======
    As with reading, every attempt has been made to match the behavior of
    Python lists. You can write one record, many records, or even certain
    elements of all records. There is one restriction: only the record
    dimension (i.e. dimension 0) can be resized by write, as all records
    in a variable must have the same dimensions. Similarly, only whole
    records can be deleted.

    For these examples, assume Flux has 100 records and dimensions [2, 3].
      1. C{Flux[0] = [[1, 2, 3], [4, 5, 6]]} rewrites the first record
         without changing the rest.
      2. C{Flux[...] = [[1, 2, 3], [4, 5, 6]]} writes a new first record
         and deletes all the rest.
      3. C{Flux[99:] = [[[1, 2, 3], [4, 5, 6]],  [[11, 12, 13], [14, 15, 16]]]}
         writes a new record in the last position and adds a new record after.
      4. C{Flux[5:6] = [[[1, 2, 3], [4, 5, 6]],  [[11, 12, 13], [14, 15, 16]]]}
         inserts two new records between the current number 5 and 6. This
         operation can be quite slow, as it requires reading and rewriting the
         entire variable. (CDF does not directly support record insertion.)
      5. C{Flux[0:2, 0, 0] = [1, 2]} changes the first element of the first
         two records but leaves other elements alone.
      6. C{del Flux[0]} removes the first record.
      7. C{del Flux[5]} removes record 5 (the sixth). Due to the need to work
         around a bug in the CDF library, this operation can be quite slow.
      8. C{del Flux[...]} deletes I{all data} from C{Flux}, but leaves the
         variable definition intact.

    @ivar cdf_file: the CDF file containing this variable
    @type cdf_file: L{CDF}
    @ivar attrs: All variable attributes for this variable
    @type attrs: L{zAttrList}
    @ivar _name: name of this variable
    @type _name: string
    @raise CDFError: if CDF library reports an error
    @raise CDFWarning: if CDF library reports a warning and interpreter
                       is set to error on warnings.
    @note: Not intended to be created directly; use methods of L{CDF}
           to gain access to a variable.
    @note: Although this interface only directly supports zVariables, zMode is
           set on opening the CDF so rVars appear as zVars. See p.24 of the
           CDF user's guide; pyCDF uses zMode 2.
    """

    def __init__(self, cdf_file, var_name, *args):
        """Create or locate a variable

        @param cdf_file: CDF file containing this variable
        @type cdf_file: L{CDF}
        @param var_name: name of this variable
        @type var_name: string
        @param args: additional arguments passed to L{_create}. If none,
                     opens an existing variable. If provided, creates a
                     new one.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        self.cdf_file = cdf_file
        self._name = None
        if len(args) == 0:
            self._get(var_name)
        else:
            self._create(var_name, *args)
        self.attrs = zAttrList(self)

    def __getitem__(self, key):
        """Returns a slice from the data array. Details under L{Var}.

        @return: The data from this variable
        @rtype: list-of-lists of appropriate type.
        @raise IndexError: if L{key} is out of range, mismatches dimensions,
                           or simply unparseable.
        @raise CDFError: for errors from the CDF library
        """
        hslice = _Hyperslice(self, key)
        if hslice.counts[0] == 0:
            return []
        buffer = hslice.create_buffer()
        hslice.select()
        lib.call(const.GET_, const.zVAR_HYPERDATA_, ctypes.byref(buffer))
        result = hslice.unpack_buffer(buffer)
        
        if self.type() == const.CDF_EPOCH.value:
            if isinstance(result, float): #single value
                result = lib.epoch_to_datetime(result)
            else:
                hslice.transform_each(result, lib.epoch_to_datetime)
        elif self.type() == const.CDF_EPOCH16.value:
            if isinstance(result[0], float): #single value
                result = lib.epoch16_to_datetime(result)
            else:
                hslice.transform_each(result, lib.epoch16_to_datetime)
            
        return hslice.convert_array(result)

    def __delitem__(self, key):
        """Removes a record (or set of records) from the CDF

        Only whole records can be deleted, so the del call must either specify
        only one dimension or it must specify all elements of the non-record
        dimensions. This is I{not} a way to resize a variable!

        Deleting records from the middle of a variable may be very slow in
        some circumstances. To work around a bug in CDF library versions
        before 3.3.1, all the data must be read in, the requested deletions
        done, and then all written back out.

        @param key: index or slice to delete
        @type key: int or slice
        @raise TypeError: if an attempt is made to delete from a non
                          record-varying variable, or to delete below
                          the record level
        """
        if not self.rv():
            raise TypeError('Cannot delete records from non-record-varying '
                            'variable.')
        hslice = _Hyperslice(self, key)
        if hslice.dims > 1 and hslice.counts[1:] != hslice.dimsizes[1:]:
            raise TypeError('Can only delete entire records.')
        if hslice.counts[0] == 0:
            return
        start = hslice.starts[0]
        count = hslice.counts[0]
        interval = hslice.intervals[0]
        dimsize = hslice.dimsizes[0]

        self._call()
        dangerous_delete = False
        if lib._del_middle_rec_bug and \
               (interval != 1 or (start != 0 and start + count < dimsize)):
            #delete from middle is dangerous if only have one index entry
            entries = ctypes.c_long(0)
            lib.call(const.GET_, const.zVAR_nINDEXENTRIES_,
                     ctypes.byref(entries))
            dangerous_delete = (entries.value == 1)
            
        if dangerous_delete:
            data = self[...]
            del data[start:start + count * interval:interval]
            self[0:dimsize - count] = data
            first_rec = dimsize - count
            last_rec = dimsize - 1
            lib.call(const.DELETE_, const.zVAR_RECORDS_,
                     ctypes.c_long(first_rec), ctypes.c_long(last_rec))
        elif interval == 1:
            first_rec = ctypes.c_long(start)
            last_rec = ctypes.c_long(start + count - 1)
            lib.call(const.DELETE_, const.zVAR_RECORDS_,
                     first_rec, last_rec)
        else:
            self._call()
            #delete from end to avoid renumbering of records
            for recno in range(start + (count - 1) * interval,
                               start - 1, -1 * interval):
                lib.call(const.DELETE_, const.zVAR_RECORDS_,
                         ctypes.c_long(recno), ctypes.c_long(recno))

    def __setitem__(self, key, data):
        """Puts a slice into the data array. Details under L{Var}.

        @param key: index or slice to store
        @type key: int or slice
        @param data: data to store
        @type data: list
        @raise IndexError: if L{key} is out of range, mismatches dimensions,
                           or simply unparseable. IndexError will 
        @raise CDFError: for errors from the CDF library
        """
        hslice = _Hyperslice(self, key)
        n_recs = hslice.counts[0]
        hslice.expand(data)
        buff = hslice.create_buffer()
        try:
            hslice.pack_buffer(buff, data)
        except (IndexError, TypeError):
            #Can be thrown if data isn't same size as slice
            #This is an expensive check, so only do it if something went wrong
            data_dims = _Hyperslice.dimensions(data)
            expected = hslice.expected_dims()
            if data_dims != expected:
                raise ValueError('attempt to assign data of dimensions ' +
                                 str(data_dims) + ' to slice of dimensions ' +
                                 str(expected))
            else: #Something else went wrong
                raise
        if hslice.counts[0] > n_recs and \
               hslice.starts[0] + n_recs < hslice.dimsizes[0]:
            #Specified slice ends before last record, so insert in middle
            saved_data = self[hslice.starts[0] + n_recs:]
        if hslice.counts[0] > 0:
            hslice.select()
            lib.call(const.PUT_, const.zVAR_HYPERDATA_, ctypes.byref(buff))
        if hslice.counts[0] < n_recs:
            first_rec = hslice.starts[0] + hslice.counts[0]
            last_rec = hslice.dimsizes[0] - 1
            lib.call(const.DELETE_, const.zVAR_RECORDS_,
                     ctypes.c_long(first_rec), ctypes.c_long(last_rec))
        elif hslice.counts[0] > n_recs and \
               hslice.starts[0] + n_recs < hslice.dimsizes[0]:
            #Put saved data in after inserted data
            self[hslice.starts[0] + hslice.counts[0]:] = saved_data

    def _create(self, var_name, datatype, n_elements = 1, dims = (),
               recVary = const.VARY, dimVarys = None):
        """Creates a new zVariable

        @param var_name: name of this variable
        @type var_name: string
        @param datatype: CDF data type
        @type datatype: ctypes.c_long
        @param n_elements: number of elements (should be 1 except for
                           CDF_CHAR variables).
        @type n_elements: long
        @param dims: size of each dimension for multi-dimensional variable,
                     or empty for a zero-dimensional
        @type dims: sequence of long
        @param recVary: record variance for this variable (VARY/NOVARY)
        @type recVary: long
        @param dimVarys: array of VARY or NOVARY, variance for each dimension
        @type dimVarys: sequence of long
        @return: new variable with this name
        @rtype: L{Var}
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        @note: Not intended to be used directly; use L{CDF.new}.
        """

        dim_array = (ctypes.c_long * len(dims))(*dims)
        enc_name = var_name.encode('ascii')
        if dimVarys == None:
            dim_vary_array = (ctypes.c_long * (len(dims) if len(dims) > 0 else 1))(const.VARY)
        else:
            dim_vary_array = (ctypes.c_long * len(dims))(*dimVarys)
        varNum = ctypes.c_long(0)
        self.cdf_file._call(const.CREATE_, const.zVAR_, enc_name, datatype,
                 ctypes.c_long(n_elements), ctypes.c_long(len(dims)), dim_array,
                 recVary, dim_vary_array, ctypes.byref(varNum))
        self._name = enc_name

    def _delete(self):
        """Removes this zVariable from the CDF

        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        self.cdf_file._call(const.SELECT_, const.zVAR_NAME_, self._name,
                 const.DELETE_, const.zVAR_)
        self._name = None

    def _get(self, var_name):
        """Gets an existing zVariable

        @param var_name: name of this variable
        @type var_name: string
        @return: variable with this name
        @rtype: L{Var}
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        @note: Not intended to be used directly; use L{CDF.__getitem__}.
        """

        if isinstance(var_name, str_classes):
            try:
                enc_name = var_name.encode('ascii')
            except AttributeError:
                enc_name = var_name #already in ASCII
            #This call simply 'touches' the CDF to cause an error if the name isn't there
            varNum = ctypes.c_long(0)
            self.cdf_file._call(const.GET_, const.zVAR_NUMBER_, enc_name, ctypes.byref(varNum))
            self._name = enc_name
        else:
            name = ctypes.create_string_buffer(const.CDF_VAR_NAME_LEN256+1)
            self.cdf_file._call(const.SELECT_, const.zVAR_, ctypes.c_long(var_name),
                     const.GET_, const.zVAR_NAME_, name)
            self._name = name.value

    def _num(self):
        """Returns the zVar number for this variable

        @return: number of this zVar
        @rtype: int
        """
        varNum = ctypes.c_long(0)
        self.cdf_file._call(const.GET_, const.zVAR_NUMBER_, self._name, ctypes.byref(varNum))
        return varNum.value

    def __len__(self):
        """Get number of records for this variable in this file

        @return: Number of records
        @rtype: long
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_MAXREC_, ctypes.byref(count))
        return (count.value + 1)

    def __repr__(self):
        """Returns representation of the variable

        Cannot return anything that can be eval'd to create a copy,
        so just wrap the informal representation in angle brackets.
        @return: info on this zVar
        @rtype: str
        """
        return '<Var:\n' + str(self) + '\n>'
        
    def __str__(self):
        """Returns a string representation of the variable

        This is an 'informal' representation in that it cannot be evaluated
        directly to create a L{Var}.

        @return: info on this zVar, CDFTYPE [dimensions] NRV
                 (if not record-varying)
        @rtype: str
        """
        cdftype = self.type()
        chartypes = (const.CDF_CHAR.value, const.CDF_UCHAR.value)
        rv = self.rv()
        typestr = lib.cdftypenames[cdftype] + \
                  ('*' + str(self._nelems()) if cdftype in chartypes else '' )
        if rv:
            sizestr = str([len(self)] + self._dim_sizes())
        else:
            sizestr = str(self._dim_sizes())
        return typestr + ' ' + sizestr + ('' if rv else ' NRV')

    def _n_dims(self):
        """Get number of dimensions for this variable

        @return: the number of dimensions
        @rtype: long
        """
        n_dims = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_NUMDIMS_, ctypes.byref(n_dims))
        return n_dims.value

    def _dim_sizes(self):
        """Get the dimension sizes for this variable

        @return: sequence of sizes
        @rtype: sequence of long
        @note: This will always be in Python order (i.e. row major,
               last index iterates most quickly), I{regardless} of the
               majority of the CDF.
        """
        sizes = (ctypes.c_long * const.CDF_MAX_DIMS)(0)
        self._call(const.GET_, const.zVAR_DIMSIZES_, sizes)
        sizes = sizes[0:self._n_dims()]
        return sizes

    def rv(self, new_rv=None):
        """Gets or sets whether this variable has record variance

        @param new_rv: True to change to record variance, False to change
                       to NRV (unspecified to simply check variance.)
        @type new_rv: boolean
        @return: True if record variance, False if NRV
        @rtype: boolean
        @note: If the variance is unknown, True is assumed
               (this replicates the apparent behaviour of the
               CDF library on variable creation).
        """
        if new_rv != None:
            self._call(const.PUT_, const.zVAR_RECVARY_,
                       const.VARY if new_rv else const.NOVARY)
        vary = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_RECVARY_, ctypes.byref(vary))
        return vary.value != const.NOVARY.value

    def dv(self, new_dv=None):
        """Gets or sets dimension variance of each dimension of variable.

        @param new_dv: Each element True to change that dimension to dimension
                       variance, False to change to not dimension variance.
                       (Unspecified to simply check variance.)
        @type new_dv: list of bool
        @return: True if that dimension has variance, else false.
        @rtype: list of bool
        @note: If the variance is unknown, True is assumed
               (this replicates the apparent behaviour of the
               CDF library on variable creation).
        """
        ndims = self._n_dims()
        if new_dv != None:
            if len(new_dv) != ndims:
                raise ValueError('Must specify variance for ' +
                                 str(ndims) + 'dimensions.')
            varies = (ctypes.c_long * ndims)(
                *[const.VARY if dv else const.NOVARY for dv in new_dv])
            self._call(const.PUT_, const.zVAR_DIMVARYS_,
                       varies)
        if ndims == 0:
            return []
        varies = (ctypes.c_long * const.CDF_MAX_DIMS)()
        self._call(const.GET_, const.zVAR_DIMVARYS_, varies)
        return [dv != const.NOVARY.value for dv in varies[0:ndims]]

    def _call(self, *args, **kwargs):
        """Select this CDF and variable and call the CDF internal interface

        Adds call to select this CDF to L{args} and passes all parameters
        directly through to the CDFlib routine of the CDF library's C internal
        interface. Checks the return value with L{Library.check_status}.

        @param args: Passed directly to the CDF library interface. Useful
                     constants are defined in the L{const} module of this package.
        @type args: various, see C{ctypes}.
        @return: CDF status from the library
        @rtype: ctypes.c_long
        @note: Terminal NULL_ is automatically added to L{args}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        return self.cdf_file._call(const.SELECT_, const.zVAR_NAME_, self._name,
                                   *args, **kwargs)

    def _c_type(self):
        """Returns the C type of this variable

        @return: ctypes type that will hold value from this variable
        @rtype: type
        @raise CDFError: for library-reported error or failure to find C type
        """
        cdftype = self.type()
        if cdftype == const.CDF_CHAR.value or cdftype == const.CDF_UCHAR.value:
            return ctypes.c_char * self._nelems()
        try:
            return lib.ctypedict[cdftype]
        except KeyError:
            raise CDFError(const.BAD_DATA_TYPE)

    def type(self, new_type=None):
        """Returns or sets the CDF type of this variable

        @param new_type: the new type, see L{const}
        @type new_type: ctypes.c_long
        @return: CDF type
        @rtype: int
        """
        if new_type != None:
            if not hasattr(new_type, 'value'):
                new_type = ctypes.c_long(new_type)
            n_elements = ctypes.c_long(self._nelems())
            self._call(const.PUT_, const.zVAR_DATASPEC_,
                       new_type, n_elements)
        cdftype = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_DATATYPE_,
                   ctypes.byref(cdftype))
        return cdftype.value

    def _nelems(self):
        """Number of elements for each value in this variable

        This is the length of strings for CHAR and UCHAR,
        should be 1 otherwise.
        @return: length of strings
        @rtype: int
        """
        nelems = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_NUMELEMS_, ctypes.byref(nelems))
        return nelems.value

    def name(self):
        """Returns the name of this variable

        @return: variable's name
        @rtype: string
        """
        if isinstance(self._name, str):
            return self._name
        elif isinstance(self._name, bytes):
            return self._name.decode()

    def compress(self, comptype=None, param=None):
        """Set or check the compression of this variable

        @param comptype: type of compression to change to, see CDF C reference
                         manual section 4.10. Constants for this parameter
                         are in L{const}. If not specified, will not change
                         compression.
        @type comptype: ctypes.c_long
        @param param: Compression parameter, see CDF CRM 4.10 and L{const}.
                      If not specified, will choose reasonable default (5 for
                      gzip; other types have only one possible parameter.)
        @type param: ctypes.c_long
        @return: (comptype, param) currently in effect
        @rtype: tuple
        @note: Compression may not be changeable on variables with data already
               written; even deleting the data may not permit the change.
        """
        return _compress(self, comptype, param)

    def copy(self):
        """Copies all data and attributes from this variable

        @return: list of all data in record order
        @rtype: L{VarCopy}
        """
        return VarCopy(self)

    def rename(self, new_name):
        """Renames this variable

        @param new_name: the new name for this variable
        @type new_name: str
        """
        try:
            enc_name = new_name.encode('ascii')
        except AttributeError:
            enc_name = new_name
        if len(enc_name) > const.CDF_VAR_NAME_LEN256:
            raise CDFError(const.BAD_VAR_NAME)
        self._call(const.PUT_, const.zVAR_NAME_, enc_name)
        self._name = enc_name


class VarCopy(list):
    """A copy of the data and attributes in a L{Var}

    Data are in the list elements, accessed much like L{Var}

    @ivar attrs: attributes for the variable
    @type attrs: dict
    @ivar _dims: dimensions of this variable
    @type _dims: list
    """

    def __init__(self, zVar):
        """Copies all data and attributes from a zVariable

        @param zVar: variable to take data from
        @type zVar: L{Var}
        """
        self.attrs = zVar.attrs.copy()
        super(VarCopy, self).__init__(zVar[...])
        self._dims = [len(zVar)] + zVar._dim_sizes()

    def __getitem__(self, key):
        """Returns a subset of the data in this copy"""
        if not hasattr(key, '__len__') and key != Ellipsis:
            return super(VarCopy, self).__getitem__(key)
        key = _Hyperslice.expand_ellipsis(key, len(self._dims))
        result = super(VarCopy, self).__getitem__(key[0])
        for subkey in key[1:]:
            result = result[subkey]
        return result


class _Hyperslice(object):
    """Represents a CDF 'slice' used for the hyper CDF functions

    For internal module use only.

    @ivar dims: number of dimensions to this slice, usually
                number of dimensions to the variable plus one
                for the record, which represents the 0th
                (least rapidly varying) dimension.
    @type dims: int
    @ivar dimsizes: size of each dimension (0th is number of records)
    @type dimsizes: list of int
    @ivar starts: index of the start value for each dimension
                  ('dimension indices' in CDF speak)
    @type starts: list of int
    @ivar counts: number of values to get from each dimension.
                  Final result will be the product of everything
                  in counts.
                  ('dimension counts' in CDF speak)
    @type counts: list of int
    @ivar intervals: interval between successive indices
                     to use for each dimension.
                     ('dimension invervals' in CDF speak)
    @type intervals: list of int
    @ivar degen: is this dimension degenerate, i.e. should be
                 removed in the returned dataset. A 3D array
                 with one dimension degenerate will be returned
                 as a 2D array (i.e. list-of-lists.)
    @type degen: list of boolean
    @ivar rev: should this dimension be returned in reverse order?
    @type rev: list of boolean
    @ivar column: is this slice in column-major mode (if false, row-major)
    @type column: boolean
    @ivar zvar: what CDF variable this object slices on
    @type zvar: L{Var}
    @ivar expanded_key: fully-expanded version of the key passed to the
                        constructor (all dimensions filled in)
    @type expanded_key: tuple
    @note: All dimension-related variables are stored row-major
           (Python order)
    """

    def __init__(self, zvar, key):
        """Create a Hyperslice

        @param zvar: zVariable that this slices
        @type zvar: L{Var}
        @param key: Python multi-dimensional slice as passed to
                    __getitem__
        @type key: tuple of slice and/or int
        @raise IndexError: if slice is out of range, mismatches dimensions, or
                           otherwise unparsable.
        @raise ValueError: if slice has invalid values
        """

        self.zvar = zvar
        self.rv = self.zvar.rv()
        self.dims = zvar._n_dims() + 1
        self.dimsizes = [len(zvar)] + \
                        zvar._dim_sizes()
        self.starts = [0] * self.dims
        self.counts = [1] * self.dims
        self.intervals = [1] * self.dims
        self.degen = [False] * self.dims
        self.rev = [False] * self.dims
        #key is:
        #1. a single value (integer or slice object) if called 1D
        #2. a tuple (of integers and/or slice objects) if called nD
        #3. Each item is either a single value (degenerate dim)
        #   or a slice object.

        if not hasattr(key, '__len__'): #Not a container object, pack in tuple
            key = (key, )
        key = self.expand_ellipsis(key,
                                   (self.dims + 1) if self.rv else self.dims)
        if len(key) == 1 and self.rv: #get all data for this record(s)
            key = self.expand_ellipsis(key + (Ellipsis, ), self.dims + 1)
        elif len(key) == self.dims - 1 or not self.rv: #NRV is always rec 0
            if self.rv: #get all records
                key = (slice(None, None, None), ) +key
            else: #NRV, so get 0th record (degenerate)
                key = (0, ) + key
        self.expanded_key = key
        if len(key) == self.dims:
            for i in range(self.dims):
                idx = key[i]
                if hasattr(idx, 'start'): #slice
                    (self.starts[i], self.counts[i],
                     self.intervals[i], self.rev[i]) = \
                     self.convert_range(idx.start, idx.stop,
                                              idx.step, self.dimsizes[i])
                else: #Single degenerate value
                    if idx < 0:
                        idx += self.dimsizes[i]
                    if idx != 0 and (idx >= self.dimsizes[i] or idx < 0):
                        raise IndexError('list index out of range')
                    self.starts[i] = idx
                    self.degen[i] = True
        else:
            raise IndexError('Slice does not match dimensions for zVar ' +
                             str(zvar._name))

        self.column = zvar.cdf_file.col_major()

    def expected_dims(self, data=None):
        """Calculate size of non-degenerate dimensions

        Figures out size, in each dimension, of expected input data

        @return: size of each dimension for this slice, excluding degnerate
        @rtype: list of int
        """
        return [self.counts[i] for i in range(self.dims) if not self.degen[i]]

    def expand(self, data):
        """Expands the record dimension of this slice to hold a set of data

        If the length of data (outermost dimension) is larger than the record
        count (counts[0]) for this slice, expand the slice to hold all the data.
        This requires that the record dimension of the slice not be degenerate,
        and also that it not have been completely specified when the hyperslice
        was created (i.e. record dimension either ellipsis or no specified
        stop.)

        Does I{not} expand any other dimension, since that's Very Hard in CDF.

        @param data: the data which are intended to be stored in this slice
        @type data: list
        """
        rec_slice = self.expanded_key[0]
        if isinstance(data, str_classes) or self.degen[0] or \
               not hasattr(rec_slice, 'stop'):
            return
        if len(data) < self.counts[0]: #Truncate to fit data
            if rec_slice.stop == None and rec_slice.step in (None, 1):
                self.counts[0] = len(data)
        elif len(data) > self.counts[0]: #Expand to fit data
            if rec_slice.step in (None, 1):
                self.counts[0] = len(data)            

    def expand_ellipsis(self, slices):
        """Expands ellipses into the correct number of full-size slices

        @param slices: tuple of slices, integers, or ellipse objects
        @type slices: tuple
        @return: key with ellipses replaced by appropriate number of blank slices
        @rtype: tuple
        @raise IndexError: if ellipses specified when already have enough dimensions
        """
        if not Ellipsis in slices:
            return slices

        size = len(slices)
        extra = self.dims - size #how many dims to replace ellipsis
        if self.rv:
            extra += 1
        if extra < 0:
            raise IndexError('Too many dimensions specified to use ellipsis.')
        idx = slices.index(Ellipsis)
        result = slices[0:idx] + \
                 tuple([slice(None, None, None) for i in range(extra)]) + \
                 slices[idx+1:]
        if Ellipsis in result:
            raise IndexError('Ellipses can only be used once per slice.')
        return result

    def create_buffer(self):
        """Creates a ctypes array to hold the data from this slice

        @return: array sized, typed, and dimensioned to hold data from
                 this slice
        @rtype: instance of subclass of ctypes.Array
        """
        counts = self.counts
        degens = self.degen
        if self.column:
            counts = self.reorder(counts)
            degens = self.reorder(degens)
        cdftype = self.zvar.type()
        constructor = self.zvar._c_type()
        #Build array from innermost out
        for count, degen in zip(counts[-1::-1], degens[-1::-1]):
            if not degen:
                constructor *= count
        return constructor()

    def convert_array(self, array):
        """Converts a nested list-of-lists to format of this slice

        Takes a list-of-lists returned in CDF order
        Converts to the order specified in this slice:
          1. Row-major (if CDF was column-major)
          2. Reversed indices if specified in slice

        @param array: data to reorder
        @type array: list (of lists)
        @return: array, reordered as necessary
        @rtype: list (of lists)
        """
        if self.column:
            if self.degen[0]:
                array = self.flip_majority(array)
            else:
                #Record-number dimension is not degenerate, so keep it first
                array = [self.flip_majority(i) for i in array]

        if True in self.rev:
            #Remove degenerate dimensions
            rev = [self.rev[i] for i in range(self.dims) if not self.degen[i]]
            for i in range(len(rev)):
                if not rev[i]: #no need to reverse
                    continue
                if i == 0:
                    array.reverse()
                else:
                    #make a flattened representation that goes to one level
                    #above what we want to flip
                    flatter = array
                    for j in range(i - 1):
                        flatter = [item for sublist in flatter for item in sublist]
                    for list in flatter:
                        list.reverse()
        return array

    def unpack_buffer(self, buffer):
        """Unpacks a buffer of data from this slice into pythonic form

        @param buffer: ctypes object created by L{create_buffer}
        @type buffer: subclass of ctypes.Array
        @return: list-of-lists of the data
        """
        cdftype = self.zvar.type()
        if cdftype == const.CDF_CHAR.value or cdftype == const.CDF_UCHAR.value:
            finaltype = self.zvar._c_type()
            for count, degen in zip(self.counts[-1::-1], self.degen[-1::-1]):
                if not degen:
                    finaltype *= count
            buffer = ctypes.cast(buffer, ctypes.POINTER(finaltype)).contents
            buffer = self.c_array_to_list(buffer)
            if str != bytes: #Need to decode output
                if isinstance(buffer, bytes):
                    buffer = buffer.decode()
                else:
                    dec = lambda x:x.decode()
                    self.transform_each(buffer, dec)
            return buffer
        else:
            return self.c_array_to_list(buffer)

    def pack_buffer(self, buff, data):
        """Packs data in the form of this slice into a buffer


        @param buff: ctypes object created by L{create_buffer}
        @type buff: subclass of ctypes.Array
        @param data: data to pack into L{buff}
        @type data: list-of-lists
        """
        counts = self.counts
        degen = self.degen
        rev = self.rev
        if self.column:
            counts = self.reorder(counts)
            degen = self.reorder(degen)
            rev = self.reorder(rev)
        cdf_type = self.zvar.type()

        #Handle completely degenerate case, i.e. scalar
        if min(degen):
            if cdf_type == const.CDF_EPOCH16.value:
                buff[:] = lib.datetime_to_epoch16(data)
            elif cdf_type == const.CDF_EPOCH.value:
                buff.value = lib.datetime_to_epoch(data)
            else:
                buff.value = data
            return

        cdf_counts = [] #non-degenerate no. of elements in each dim., CDF order
        cdf_rev = [] #is this dimension reversed? CDF order
        for i in range(len(counts)):
            if not degen[i]:
                cdf_counts.append(counts[i])
                cdf_rev.append(rev[i])
        lastdim = len(cdf_counts) - 1
        record_degen = self.degen[0]
        reordered = max(cdf_rev) or self.column #python and CDF order differ?

        cdf_idx = [0 for i in range(len(cdf_counts))]
        while cdf_idx[0] < cdf_counts[0]:
            #Calculate the Python index
            py_idx = cdf_idx[:]
            if reordered:
                for i in range(lastdim + 1):
                    if self.column and (record_degen or i != 0):
                        if cdf_rev[i]:
                            py_idx[lastdim - i] = \
                                           cdf_counts[i] - cdf_idx[i] - 1
                        else:
                            py_idx[lastdim - i] = cdf_idx[i]
                    else:
                        if cdf_rev[i]:
                            py_idx[i] = cdf_counts[i] - cdf_idx[i] - 1
                        else:
                            py_idx[i] = cdf_idx[i]

            py_obj = data
            cdf_obj = buff
            for i in range(lastdim):
                py_obj = py_obj[py_idx[i]]
                cdf_obj = cdf_obj[cdf_idx[i]]
            if cdf_type == const.CDF_EPOCH16.value:
                cdf_obj[cdf_idx[lastdim]][:] = lib.datetime_to_epoch16(
                    py_obj[py_idx[lastdim]])
            elif cdf_type == const.CDF_EPOCH.value:
                cdf_obj[cdf_idx[lastdim]] = lib.datetime_to_epoch(
                    py_obj[py_idx[lastdim]])
            elif cdf_type == const.CDF_CHAR.value or \
                     cdf_type == const.CDF_UCHAR.value:
                cdf_obj[cdf_idx[lastdim]].value = py_obj[py_idx[lastdim]]
            else:
                cdf_obj[cdf_idx[lastdim]] = py_obj[py_idx[lastdim]]

            #Switch to the next index
            currdim = lastdim
            while currdim >= 0:
                cdf_idx[currdim] += 1
                if cdf_idx[currdim] >= cdf_counts[currdim] and currdim > 0:
                    cdf_idx[currdim] = 0
                    currdim -= 1
                else:
                    break

    def c_array_to_list(self, array):
        """Converts a ctypes array type to a python nested list

        @param array: the array to convert
        @type array: instance of ctypes.Array subclass
        @return: contents of array
        @rtype: list (of lists)
        """
        if hasattr(array, 'value'):
            return array.value

        dimsizes = []
        counts = self.counts
        degen = self.degen
        if self.column:
            counts = self.reorder(counts)
            degen = self.reorder(degen)
        for count, degen in zip(counts, degen):
            if not degen:
                dimsizes.append(count)

        if self.zvar.type() == const.CDF_EPOCH16.value:
            basetype = lib.ctypedict[const.CDF_EPOCH.value]
            dimsizes.append(2)
        else:
            basetype = self.zvar._c_type()

        n_elements = 1
        for j in dimsizes:
            n_elements *= j
        flat_type = basetype * n_elements
        flat = ctypes.cast(array,
                           ctypes.POINTER(flat_type)
                           ).contents

        dims = [i for i in reversed(dimsizes)]
        if len(dims) == 1:
            if hasattr(flat[0], 'value'):
                return [i.value for i in flat]
            else:
                return [i for i in flat]

        for i in range(len(dims) - 1):
            size = dims[i]
            n_lists = 1
            for j in dims[i + 1:]:
                n_lists *= j
            if i == 0:
                if hasattr(flat[0], 'value'):
                    result = [[k.value for k in flat[j * size:(j + 1) * size]]
                              for j in range(n_lists)]
                else:
                    result = [[k for k in flat[j * size:(j + 1) * size]]
                              for j in range(n_lists)]
            else:
                result = [result[j * size:(j + 1) * size] for j in range(n_lists)]
        return result

    def transform_each(self, data, func):
        """Transforms each element of an n-D array

        L{data} must be dimensioned the same as this slice. Iterates
        over each element in L{data}, applying L{func} to the element
        and assigning the result back to the element.

        L{data} must I{not} be a scalar, as then it cannot be changed.

        @param data: the data to transform
        @type data: list
        @param func: the function to call for each element in L{data}
        @type func: callable
        """
        counts = self.counts
        degen = self.degen
        if self.column:
            counts = self.reorder(counts)
            degen = self.reorder(degen)
        cdf_counts = [] #non-degenerate no. of elements in each dim., CDF order
        for i in range(len(counts)):
            if not degen[i]:
                cdf_counts.append(counts[i])

        lastdim = len(cdf_counts) - 1
        idx = [0 for i in range(len(cdf_counts))]
        while idx[0] < cdf_counts[0]:
            parent = data
            for i in idx[:-1]:
                parent = parent[i]
            parent[idx[-1]] = func(parent[idx[-1]])

            currdim = lastdim
            while currdim >= 0:
                idx[currdim] += 1
                if idx[currdim] >=  cdf_counts[currdim] and currdim > 0:
                    idx[currdim] = 0
                    currdim -= 1
                else:
                    break

    def select(self):
        """Selects this hyperslice in the CDF

        Calls the CDF library to select the CDF, variable, records, and
        array elements corresponding to this slice.
        """
        args = (const.SELECT_, const.zVAR_RECNUMBER_, ctypes.c_long(self.starts[0]),
                const.SELECT_, const.zVAR_RECCOUNT_, ctypes.c_long(self.counts[0]),
                const.SELECT_, const.zVAR_RECINTERVAL_,
                ctypes.c_long(self.intervals[0]))
        if self.dims > 1:
            dims = self.dims - 1
            args += (const.SELECT_, const.zVAR_DIMINDICES_,
                     (ctypes.c_long * dims)(*self.starts[1:]),
                     const.SELECT_, const.zVAR_DIMCOUNTS_,
                     (ctypes.c_long * dims)(*self.counts[1:]),
                     const.SELECT_, const.zVAR_DIMINTERVALS_,
                     (ctypes.c_long * dims)(*self.intervals[1:]))
        self.zvar._call(*args)

    @staticmethod
    def expand_ellipsis(slices, n_dims):
        """Expands any ellipses into correct number of full-size slices

        @param slices: tuple of slices, integers, or ellipse objects
        @type slices: tuple
        @param n_dims: number of dimensions this slice is over
        @type n_dims: int
        @return: L{slices} with ellipses replaced by appropriate number of
                 full-dimension slices
        @rtype: tuple
        @raise IndexError: if ellipses specified when already have enough
                           dimensions
        """
        if slices == Ellipsis:
            return tuple([slice(None, None, None)
                          for i in range(n_dims)])
        if not Ellipsis in slices:
            return slices

        extra = n_dims - len(slices) #how many dims to replace ellipsis
        if extra < 0:
            raise IndexError('Too many dimensions specified to use ellipsis.')
        idx = slices.index(Ellipsis)
        result = slices[0:idx] + \
                 tuple([slice(None, None, None) for i in range(extra)]) + \
                 slices[idx+1:]
        if Ellipsis in result:
            raise IndexError('Ellipses can only be used once per slice.')
        return result

    @staticmethod
    def flip_majority(array):
        """Inverts the array majority of an existing array

        @param array: data to change majority of
        @type array: list (of lists)
        @return: array, with majority changed (return[j][i] == array[i][j])
        @rtype: list (of lists)
        @note: Assumes the fundamental type is not iterable--
               this will not work for e.g. an array of tuples
        @note: array must be 'square' / 'regular' --i.e.
               all lists at a particular level of dimensionality
               must have the same size
        @raise RuntimeError: if dimensionality of result not what it
               should be (i.e. error in I{this} function).
        """
        try:
            dims = [len(array)] #dimensions of array
        except TypeError:
            return array #scalar

        flat = array #this is progressively flattened (i.e. dims stripped off)
        while True:
            try:
                if isinstance(flat[0], str_classes):
                    break
                lengths = [len(i) for i in flat]
            except TypeError: #Now completely flat
                break
            if min(lengths) != max(lengths):
                raise TypeError('Array dimensions not regular')
            dims.append(lengths[0])
            flat = [item for sublist in flat for item in sublist]

        result = flat
        for i in range(len(dims) - 1):
            stride = 1
            for j in dims[i+1:]:
                stride *= j
            result = [result[j::stride] for j in range(stride)]
        if len(result) != dims[-1]:
            raise RuntimeError('Dimensionality mismatch: ' +
                               len(result) + dims)
        return result

    @staticmethod
    def dimensions(data):
        """Finds the dimensions of a nested list-of-lists

        @param data: data of which dimensions are desired
        @type data: list (of lists)
        @return: dimensions of L{data}, in order outside-in
        @rtype: list of int
        @raise ValueError: if L{data} has irregular dimensions
        """
        return _Hyperslice.types(data)[0]

    @staticmethod
    def types(data):
        """Find dimensions and valid types of a nested list-of-lists

        Any given data may be representable by a range of CDF types; infer
        the CDF types which can represent this data. This breaks down to:
          1. Proper kind (numerical, string, time)
          2. Proper range (stores highest and lowest number)
          3. Sufficient resolution (EPOCH16 required if datetime has
             microseconds or below.)

        If more than one value satisfies the requirements, types are returned
        in preferred order:
          1. Type that matches precision of data first, then
          2. integer type before float type, then
          3. Smallest type first, then
          4. signed type first, then
          5. specifically-named (CDF_BYTE) vs. generically named (CDF_INT1)
        So for example, EPOCH_16 is preferred over EPOCH if L{data} specifies
        below the millisecond level (rule 1), but otherwise EPOCH is preferred
        (rule 2).

        For floats, four-byte is preferred unless eight-byte is required:
          1. absolute values between 0 and 3e-39
          2. absolute values greater than 1.7e38
        This will switch to an eight-byte double in some cases where four bytes
        would be sufficient for IEEE 754 encoding, but where DEC formats would
        require eight.

        @param data: data for which dimensions and CDF types are desired
        @type data: list (of lists)
        @return: dimensions of L{data}, in order outside-in;
                 CDF types which can represent this data;
                 number of elements required (i.e. length of longest string)
        @rtype: 3-tuple of lists ([int], [ctypes.c_long], [int])
        @raise ValueError: if L{data} has irregular dimensions
        """
        if isinstance(data, str_classes):
            dims = []
        else:
            try:
                dims = [len(data)]
            except TypeError:
                dims = []
        if dims:
            flat = data
        else:
            flat = [data]

        while True:
            try:
                if isinstance(flat[0], str_classes):
                    break
                lengths = [len(i) for i in flat]
            except TypeError: #Now completely flat
                break
            if min(lengths) != max(lengths):
                raise ValueError('Data irregular in dimension ' +
                                str(len(dims)))
            dims.append(lengths[0])
            flat = [item for sublist in flat for item in sublist]

        elements = 1
        if isinstance(flat[0], str_classes):
            types = [const.CDF_CHAR, const.CDF_UCHAR]
            elements = max([len(s) for s in flat])
        elif isinstance(flat[0], datetime.datetime):
            if max([dt.microsecond % 1000 for dt in flat]) > 0:
                types = [const.CDF_EPOCH16, const.CDF_EPOCH]
            else:
                types = [const.CDF_EPOCH, const.CDF_EPOCH16]
        elif max([hasattr(i, 'is_integer') for i in flat]):
            absolutes = [abs(i) for i in flat if i != 0]
            if len(absolutes) > 0 and \
                   (max(absolutes) > 1.7e38 or min(absolutes) < 3e-39):
                types = [const.CDF_DOUBLE, const.CDF_REAL8]
            else:
                types = [const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8]
        else:
            minval = min(flat)
            maxval = max([abs(i) for i in flat])
            if minval < 0:
                types = [const.CDF_BYTE, const.CDF_INT1,
                         const.CDF_INT2, const.CDF_INT4,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8]
                cutoffs = [2 ** 7, 2 ** 7, 2 ** 15, 2 ** 31,
                           1.7e38, 1.7e38, 8e307, 8e307]
            else:
                types = [const.CDF_BYTE, const.CDF_INT1, const.CDF_UINT1,
                         const.CDF_INT2, const.CDF_UINT2,
                         const.CDF_INT4, const.CDF_UINT4,
                         const.CDF_FLOAT, const.CDF_REAL4,
                         const.CDF_DOUBLE, const.CDF_REAL8]
                cutoffs = [2 ** 7, 2 ** 7, 2 ** 8,
                           2 ** 15, 2 ** 16, 2 ** 31, 2 ** 32,
                           1.7e38, 1.7e38, 8e307, 8e307]
            types = [t for (t, c) in zip(types, cutoffs) if c > maxval]
        types = [t.value for t in types]
        return (dims, types, elements)

    @staticmethod
    def reorder(seq):
        """Reorders seq to switch array majority

        Used to take an array of subscripts between row
        and column majority. First element is not touched,
        being the record number.

        @param seq: a sequence of I{subscripts}
        @type seq: sequence of integers
        @return: seq with all but element 0 reversed in order
        @rtype: sequence of integers
        """
        if hasattr(seq, '__getitem__'):
            return seq[0:1] + seq[-1:0:-1]
        else:
            return seq

    @staticmethod
    def convert_range(start, stop, step, size):
        """Converts a start/stop/step range to start/count/interval

        (i.e. changes from Python-style slice to CDF-style)
        @param start: index to start a slice at, may be none or negative
        @type start: int
        @param stop: index at end of slice (one-past, standard Python),
                     may be none or negative
        @type stop: int
        @param step: interval for stepping through stlice
        @type step: int
        @param size: size of list to slice
        @type size: int
        @return: (start, count, interval, rev) where:
                   1. start is the start index, normalized to be within
                      the size of the list and negatives handled
                   2. count is the number of records in the slice,
                      guaranteed to stop before the end
                   3. interval is the skip between records
                   4. rev indicates whether the sequence should be reversed
        @rtype: (int, int, int, boolean)
        """
        (start, stop, step) = slice(start, stop, step).indices(size)
        if step < 0:
            step *= -1
            count = int((start - stop + step - 1) / step)
            start = start - (count - 1) * step
            rev = True
        else:
            count = int((stop - start + step - 1) / step)
            rev = False
        if count < 0:
            count = 0
            start = 0
        return (start, count, step, rev)


class AttrList(collections.Mapping):
    """Object representing a list of attributes.

    Only used in its subclasses, L{gAttrList} and L{zAttrList}

    @ivar _cdf_file: CDF these attributes are in
    @type _cdf_file: L{CDF}
    @ivar special_entry: callable which returns a "special"
                         entry number, used to limit results
                         for zAttrs to those which match the zVar
    @type special_entry: callable
    @ivar AttrType: type of attribute in this list, L{zAttr} or L{gAttr}
    @type AttrType: type
    @ivar attr_name: name of attribute type, zAttribute or gAttribute
    @type attr_name: str
    @ivar global_scope: is this list scoped global (True) or variable (False)
    @type global_scope: bool
    """

    def __init__(self, cdf_file, special_entry=None):
        """Initialize the attribute collection

        @param cdf_file: CDF these attributes are in
        @type cdf_file: L{CDF}
        @param special_entry: callable which returns a "special"
                              entry number, used to limit results
                              for zAttrs to those which match the zVar
        @type special_entry: callable
        """
        self._cdf_file = cdf_file
        self.special_entry = special_entry

    def __getitem__(self, name):
        """Find an Attribute by name

        @param name: name of the Attribute to return
        @type name: str
        @return: attribute named L{name}
        @rtype: L{Attr}
        @raise KeyError: if there is no attribute named L{name}
        @raise CDFError: other errors in CDF library
        """
        try:
            attrib = self.AttrType(self._cdf_file, name)
        except CDFError:
            (t, v, tb) = sys.exc_info()
            if v.status == const.NO_SUCH_ATTR:
                raise KeyError(name + ': ' + str(v))
            else:
                raise
        if attrib.global_scope() != self.global_scope:
            raise KeyError(name + ': no ' + self.attr_name + ' by that name.')
        return attrib

    def __setitem__(self, name, data):
        """Create an Attribute or change its entries

        @param name: name of Attribute to change
        @type name: str
        @param data: Entries to populate this Attribute with.
                     Any existing Entries will be deleted!
                     Another C{Attr} may be specified, in which
                     case all its entries are copied.
        @type data: scalar, list, or L{Attr}
        """
        if isinstance(data, AttrList):
            if name in self:
                del self[name]
            attr = self._get_or_create(name)
            for entryno in range(data.max_idx()):
                if data.has_entry(entryno):
                    attr.new(data[entryno], data.type(entryno), entryno)
        else:
            attr = self._get_or_create(name)
            if isinstance(data, str_classes):
                data = [data]
            else:
                try:
                    junk = len(data)
                except TypeError:
                    data = [data]
            attr[:] = data
            del attr[len(data):]

    def __delitem__(self, name):
        """Delete an Attribute (and all its entries)

        @param name: name of Attribute to delete
        @type name: str
        """
        try:
            attr = self.AttrType(self._cdf_file, name)
        except CDFError:
            (t, v, tb) = sys.exc_info()
            if v.status == const.NO_SUCH_ATTR:
                raise KeyError(name + ': ' + str(v))
            else:
                raise
            if attr.global_scope() != self.global_scope:
                raise KeyError(name + ': not ' + self.attr_name)
        attr._delete()

    def __iter__(self, current=0):
        """Iterates over all Attr in this CDF or variable

        Returns name of one L{Attr} at a time until reaches the end.
        @note: Returned in number order.
        """
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMATTRS_,
                             ctypes.byref(count))
        while current < count.value:
            candidate = self.AttrType(self._cdf_file, current)
            if candidate.global_scope() == self.global_scope:
                if self.special_entry == None or \
                        candidate.has_entry(self.special_entry()):
                    if str == bytes:
                        value = yield(candidate._name)
                    else:
                        value = yield(candidate._name.decode())
                    if value != None:
                        current = self[value].number()
            current += 1

    def __repr__(self):
        """Returns representation of attribute list

        Cannot return anything that can be eval'd to create a copy of the
        list, so just wrap the informal representation in angle brackets.
        @return: all the data in this list of attributes
        @rtype: str
        """
        return '<' + self.__class__.__name__ + ':\n' + str(self) + '\n>'

    def __str__(self):
        """Returns a string representation of the attribute list

        This is an 'informal' representation in that it cannot be evaluated
        directly to create an L{AttrList}.

        @return: all the data in this list of attributes
        @rtype: str
        """
        return '\n'.join([key + ': ' + (
            ('\n' + ' ' * (len(key) + 2)).join(
            [str(value[i]) + ' [' + lib.cdftypenames[value.type(i)] + ']'
             for i in range(len(value))])
            if isinstance(value, Attr)
            else str(value) +
            ' [' + lib.cdftypenames[self.type(key)] + ']'
            )
            for (key, value) in self.items()])

    def clone(self, master, name=None, new_name=None):
        """Clones this attribute list, or one attribute in it, from another

        @param master: the attribute list to copy from
        @type master: L{AttrList}
        @param name: name of attribute to clone (default: clone entire list)
        @type name: str
        @param new_name: name of the new attribute, default L{name}
        @type new_name: str
        """
        if name == None:
            self._clone_list(master)
        else:
            self._clone_attr(master, name, new_name)

    def copy(self):
        """Create a copy of this attribute list

        @return: copy of the entries for all attributes in this list
        @rtype: dict
        """
        return dict((key, value[:] if isinstance(value, Attr) else value)
                    for (key, value) in self.items())

    def new(self, name, data=None, type=None):
        """Create a new Attr in this AttrList

        @param name: name of the new Attribute
        @type name: str
        @param data: data to put into the first entry in the new Attribute
        @param type: CDF type of the first entry from L{const}. Only used
                     if L{data} are specified.
        @raise KeyError: if L{name} already exists in this list
        """
        if name in self:
            raise KeyError(name + ' already exists.')
        attr = self._get_or_create(name)
        if data != None:
            if self.special_entry == None:
                attr.new(data, type)
            else:
                attr.new(data, type, self.special_entry())

    def rename(self, old_name, new_name):
        """Rename an attribute in this list

        Renaming a zAttribute renames it for I{all} zVariables in this CDF!

        @param old_name: the current name of the attribute
        @type old_name: str
        @param new_name: the new name of the attribute
        @type new_name: str
        """
        AttrList.__getitem__(self, old_name).rename(new_name)

    def _clone_attr(self, master, name, new_name=None):
        """Clones a single attribute from one in this list or another

        Copies data and types from the master attribute to the new one

        @param master: attribute list to copy attribute from
        @type master: L{AttrList}
        @param name: name of attribute to copy
        @type name: str
        @param new_name: name of the new attribute, default L{name}
        @type new_name: str
        """
        if new_name == None:
            new_name = name
        self[new_name] = master[name]

    def _clone_list(self, master):
        """Clones this attribute list from another

        @param master: the attribute list to copy from
        @type master: L{AttrList}
        """
        for name in master:
            self._clone_attr(master, name)
        for name in list(self): #Can't iterate over a list we're changing
            if not name in master:
                del self[name]    

    def _get_or_create(self, name):
        """Retrieve L{Attr} or create it if it doesn't exist

        @param name: name of the attribute to look up or create
        @type name: str
        @return: attribute with this name
        @rtype: L{Attr}
        """
        attr = None
        try:
            attr = self.AttrType(self._cdf_file, name)
        except CDFError:
            (t, v, tb) = sys.exc_info()
            if v.status != const.NO_SUCH_ATTR:
                raise
        if attr == None:
            attr = self.AttrType(self._cdf_file, name, True)
        elif attr.global_scope() != self.global_scope:
                raise KeyError(name + ': not ' + self.attr_name)
        return attr


class gAttrList(AttrList):
    """Object representing I{all} the gAttributes in a CDF.

    Normally accessed as an attribute of an open L{CDF}::
        global_attribs = cdffile.attrs

    Appears as a dictionary: keys are attribute names; each value is an
    attribute represented by a L{gAttr} object. To access the global
    attribute TEXT::
        text_attr = cdffile.attrs['TEXT']
    """

    def __init__(self, cdf_file):
        """Initialize the attribute collection

        @param cdf_file: CDF these attributes are in
        @type cdf_file: L{CDF}
        """
        self.AttrType = gAttr
        self.attr_name = 'gAttribute'
        self.global_scope = True
        super(gAttrList, self).__init__(cdf_file)

    def __len__(self):
        """Number of gAttributes in this CDF

        @return: number of gAttributes in the CDF
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMgATTRS_,
                             ctypes.byref(count))
        return count.value


class zAttrList(AttrList):
    """Object representing I{all} the zAttributes in a zVariable.

    Normally access as an attribute of a L{Var} in an open CDF::
        epoch_attribs = cdffile['Epoch'].attrs

    Appears as a dictionary: keys are attribute names, values are
    the value of the zEntry associated with the appropriate zVariable.
    Each vAttribute in a CDF may only have a I{single} entry associated
    with each variable. The entry may be a string, a single numerical value,
    or a series of numerical values. Entries with multiple values are returned
    as an entire list; direct access to the individual elements is not
    possible.

    Example: finding the first dependency of (ISTP-compliant) variable
    Flux::
        print cdffile['Flux'].attrs['DEPEND_0']

    zAttributes are shared among zVariables, one zEntry allowed per zVariable.
    (pyCDF hides this detail.) Deleting the last zEntry for a zAttribute will
    delete the underlying zAttribute.

    zEntries are created and destroyed by the usual dict methods on the
    zAttrlist::
        epoch_attribs['new_entry'] = [1, 2, 4] #assign a list to new zEntry
        del epoch_attribs['new_entry'] #delete the zEntry
    L{__setitem__} describes how the type of an zEntry is determined.

    @ivar _zvar: zVariable these attributes are in
    @type _zvar: L{Var}
    @ivar _cdf_file: CDF these attributes are in
    @type _cdf_file: L{CDF}
    """

    def __init__(self, zvar):
        """Initialize the attribute collection

        @param zvar: zVariable these attributes are in
        @param zvar: L{Var}
        """
        self.AttrType = zAttr
        self.attr_name = 'zAttribute'
        self.global_scope = False
        super(zAttrList, self).__init__(zvar.cdf_file, zvar._num)
        self._zvar = zvar

    def __getitem__(self, name):
        """Find an zEntry by name

        @param name: name of the zAttribute to return
        @type name: str
        @return: attribute named L{name}
        @rtype: L{zAttr}
        @raise KeyError: if there is no attribute named L{name} associated
                         with this zVariable
        @raise CDFError: other errors in CDF library
        """
        attrib = super(zAttrList, self).__getitem__(name)
        zvar_num = self._zvar._num()
        if attrib.has_entry(zvar_num):
            return attrib[zvar_num]
        else:
            raise KeyError(name + ': no such attribute for variable ' +
                           self._zvar.name())

    def __delitem__(self, name):
        """Delete an zEntry by name

        @param name: name of the zEntry to delete
        @type name: str
        @raise KeyError: if there is no attribute named L{name} associated
                         with this zVariable
        @raise CDFError: other errors in CDF library
        @note: If this is the only remaining entry, the Attribute will be
               deleted.
        """
        attrib = super(zAttrList, self).__getitem__(name)
        zvar_num = self._zvar._num()
        if not attrib.has_entry(zvar_num):
            raise KeyError(str(name) + ': no such attribute for variable ' +
                           str(self._zvar._name))
        del attrib[zvar_num]
        if len(attrib) == 0:
            attrib._delete()

    def __setitem__(self, name, data):
        """Sets a zEntry by name

        The type of the zEntry is guessed from L{data}. The type is chosen to
        match the data; subject to that constraint, it will try to match
        (in order):
          1. existing zEntry corresponding to this zVar
          2. other zEntries in this zAttribute
          3. the type of this zVar
          4. data-matching constraints described in L{_Hyperslice.types}

        @param name: name of zAttribute; zEntry for this zVariable will be set
                     in zAttribute by this name
        @type name: str
        @raise CDFError: errors in CDF library
        @raise ValueError: if unable to find a valid CDF type matching L{data},
                           or if L{data} is the wrong dimensions.
        """
        try:
            attr = super(zAttrList, self).__getitem__(name)
        except KeyError:
            attr = zAttr(self._cdf_file, name, True)
        zvar_num = self._zvar._num()
        attr[zvar_num] = data

    def __len__(self):
        """Number of zAttributes in this variable

        @return: number of zAttributes in the CDF
                 which have entries for this variable.
        @rtype: int
        """
        length = 0
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMATTRS_,
                             ctypes.byref(count))
        current = 0
        while current < count.value:
            candidate = zAttr(self._cdf_file, current)
            if not candidate.global_scope():
                if candidate.has_entry(self._zvar._num()):
                    length += 1
            current += 1
        return length

    def type(self, name, new_type=None):
        """Find or change the CDF type of a zEntry in this zVar

        @param name: name of the zAttr to check or change
        @type name: str
        @param new_type: type to change it to, see L{const}
        @type new_type: ctypes.c_long
        @return: CDF variable type, see L{const}
        @rtype: int
        @note: If changing types, old and new must be equivalent, see CDF
               User's Guide section 2.5.5 pg. 57
        """
        attrib = super(zAttrList, self).__getitem__(name)
        zvar_num = self._zvar._num()
        if not attrib.has_entry(zvar_num):
            raise KeyError(name + ': no such attribute for variable ' +
                           self._zvar.name())
        return attrib.type(zvar_num, new_type)

    def _clone_attr(self, master, name, new_name=None):
        """Clones a single attribute from one in this list or another

        Copies data and types from the master attribute to the new one

        @param master: attribute list to copy attribute from
        @type master: L{zAttrList}
        @param name: name of attribute to copy
        @type name: str
        @param new_name: name of the new attribute, default L{name}
        @type new_name: str
        """
        if new_name == None:
            new_name = name
        if new_name in self:
            del self[new_name]
        self.new(new_name, master[name], master.type(name))


class Attr(collections.Sequence):
    """An attribute, z or g, for a CDF

    This class should not be used directly.

    Represents a CDF attribute, providing access to the Entries in a format
    that looks like a Python
    list. General list information is available in the python docs:
    U{1<http://docs.python.org/tutorial/introduction.html#lists>},
    U{2<http://docs.python.org/tutorial/datastructures.html#more-on-lists>},
    U{3<http://docs.python.org/library/stdtypes.html#typesseq>}.

    Each element of the list is a single Entry of the appropriate type.
    The index to the elements is the Entry number.

    Multi-dimensional slicing is I{not} supported; an Entry with multiple
    elements will have all elements returned (and can thus be sliced itself).
    Example::
        first_three = attribute[5, 0:3] #will fail
        first_three = attribute[5][0:3] #first three elements of 5th Entry

    @ivar _cdf_file: CDF file containing this attribute
    @type _cdf_file: L{CDF}
    @ivar _name: Name of the attribute
    @type _name: bytes
    """

    def __init__(self, cdf_file, attr_name, create=False):
        """Initialize this attribute

        @param cdf_file: CDF file containing this attribute
        @type cdf_file: L{CDF}
        @param attr_name: Name of this attribute
        @type attr_name: str
        @param create: True to create attribute,
                       False to look up existing.
        @type create: bool
        """
        self._cdf_file = cdf_file
        if isinstance(attr_name, str_classes):
            try:
                self._name = attr_name.encode('ascii')
            except AttributeError:
                self._name = attr_name
            if create:
                attrno = ctypes.c_long(0)
                self._cdf_file._call(const.CREATE_, const.ATTR_,
                                     self._name, self.SCOPE,
                                     ctypes.byref(attrno))
            else:
                self._cdf_file._call(const.CONFIRM_, const.ATTR_EXISTENCE_,
                                     self._name)
        else:
            name = ctypes.create_string_buffer(const.CDF_ATTR_NAME_LEN256 + 1)
            self._cdf_file._call(const.SELECT_, const.ATTR_,
                                 ctypes.c_long(attr_name))
            self._cdf_file._call(const.GET_, const.ATTR_NAME_, name)
            self._name = name.value

    def __getitem__(self, key):
        """Return a slice of Entries.

        Because Attributes may be sparse, a multi-element slice will return
        None for those elements which do not have associated Entries.

        @param key: index or range of Entry number to return
        @type key: slice or int
        @return: a list of entries, appropriate type.
        @raise IndexError: if L{key} is an int and that Entry number does not
                           exist.
        """
        if hasattr(key, 'indices'):
            idx = range(*key.indices(self.max_idx() + 1))
            return [self._get_entry(i) if self.has_entry(i) else None
                    for i in idx]
        else:
            if self.has_entry(key):
                return self._get_entry(key)
            else:
                raise IndexError('list index ' + str(key) + ' out of range.')

    def __setitem__(self, key, data):
        """Set a slice of Entries.

        @param key: index or range of Entry numbers to set
        @type key: slice or int
        @param data: the data to set these entries to.
                     Normally each entry should be a sequence; if
                     a scalar is provided, it is treated as a single-element
                     list.
        @type data: scalar or list
        @raise ValueError: if size of {data} does not match size of L{key}
        @note: Attributes do not 'grow' or 'shrink' as entries are added
               or removed. Indexes of entries never change and there is no
               way to 'insert'.
        """
        if key == Ellipsis:
            key = slice(None, None, None)
        if not hasattr(key, 'indices'):
            #Single value, promote everything a dimension
            idx = (key, key + 1, 1)
            data = [data]
        else:
            idx = key.indices(self.max_idx() + 1)
            if key.step == None or key.step > 0:
                #Iterating forward, extend slice to match data
                if len(data) > len(range(*idx)):
                    idx = (idx[0], idx[0] + idx[2] * len(data), idx[2])

        #get, and check, types and sizes for all data
        #checks first so don't have error after changing half the Entries
        data_idx = -1
        typelist = []
        for i in range(*idx):
            data_idx += 1
            if data_idx >= len(data):
                continue
            datum = data[data_idx]
            if datum == None:
                typelist[i] = (None, None, None)
                continue
            (dims, types, elements) = _Hyperslice.types(datum)
            if len(types) <= 0:
                raise ValueError('Cannot find a matching CDF type.')
            if len(dims) > 1:
                raise ValueError('Entries must be scalar or 1D.')
            elif len(dims) == 1 and isinstance(datum[0], str_classes):
                raise ValueError('Entry strings must be scalar.')
            entrytypes = []
            if self.has_entry(i):
                entrytype = self.type(i)
                if entrytype in types:
                    entrytypes = [entrytype]
            if not entrytypes:
                for num in range(self.max_idx()):
                    if self.has_entry(num):
                        entrytype = self.type(num)
                        if entrytype in types:
                            entrytypes.append(entrytype)
            if entrytypes:
                #Of those types in entrytypes, find the one which is earliest
                # in types, i.e. the preferred type
                entry_type = types[
                    min([types.index(entrytype) for entrytype in entrytypes])
                    ]
            elif self.ENTRY_ == const.zENTRY_:
                vartype = self._cdf_file[i].type()
                if vartype in types:
                    entry_type = vartype
                else:
                    entry_type = types[0]
            else:
                entry_type = types[0]
            if not entry_type in lib.ctypedict:
                raise ValueError('Cannot find a matching ctypes type.')
            typelist.append((dims, entry_type, elements))

        data_idx = -1
        for i in range(*idx):
            data_idx += 1
            if data_idx >= len(data) or data[data_idx] == None:
                if self.has_entry(i):
                    del self[i]
                continue
            datum = data[data_idx]
            (dims, entry_type, elements) = typelist[data_idx]
            self._write_entry(i, datum, entry_type, dims, elements)

    def __delitem__(self, key):
        """Delete a slice of Entries.

        @param key: index or range of Entry numbers to delete
        @type key: slice or int
        @note: Attributes do not 'grow' or 'shrink' as entries are added
               or removed. Indexes of entries never change and there is no
               way to 'insert'.
        """
        if key == Ellipsis:
            key = slice(None, None, None)
        if not hasattr(key, 'indices'):
            idx = (key, key + 1, 1)
        else:
            idx = key.indices(self.max_idx() + 1)
        for i in range(*idx):
            self._call(const.SELECT_, self.ENTRY_, ctypes.c_long(i),
                       const.DELETE_, self.ENTRY_)

    def __iter__(self, current=0):
        """Iterates over all entries in this Attribute

        Returns data from one entry at a time until reaches the end.
        @note: Returned in entry-number order.
        """
        while current <= self.max_idx():
            if self.has_entry(current):
                value = yield(self._get_entry(current))
                if value != None:
                    current = value
            current += 1

    def __reversed__(self, current=None):
        """Iterates over all entries in this Attribute

        Returns data from one entry at a time, starting at end and going
        to beginning.
        @note: Returned in entry-number order.
        """
        if current == None:
            current = self.max_idx()
        while current >= 0:
            if self.has_entry(current):
                value = yield(self._get_entry(current))
                if value != None:
                    current = value
            current -= 1

    def __len__(self):
        """Number of Entries for this Attr. NOT same as max Entry number.

        @return: Number of Entries
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, self.ATTR_NUMENTRIES_, ctypes.byref(count))
        return count.value

    def __repr__(self):
        """Returns representation of an attribute

        Cannot return anything that can be eval'd to create a copy of the
        attribtute, so just wrap the informal representation in angle brackets.
        @return: all the data in this attribute
        @rtype: str
        """
        return '<\n' + str(self) + '\n>'

    def __str__(self):
        """Returns a string representation of the attribute

        This is an 'informal' representation in that it cannot be evaluated
        directly to create an L{Attr}.

        @return: all the data in this attribute
        @rtype: str
        """
        return '\n'.join([str(item) for item in self])

    def _call(self, *args, **kwargs):
        """Select this CDF and Attr and call the CDF internal interface

        @param args: Passed directly to the CDF library interface.
        @type args: various, see C{ctypes}.
        @return: CDF status from the library
        @rtype: ctypes.c_long
        @note: Terminal NULL_ is automatically added to L{args}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        return self._cdf_file._call(
            const.SELECT_, const.ATTR_NAME_, self._name,
            *args, **kwargs)

    def _entry_len(self, number):
        """Number of elements in an Entry

        @param number: number of Entry
        @type number: int
        @return: number of elements
        @rtype: int
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        count = ctypes.c_long(0)
        self._call(
            const.SELECT_, self.ENTRY_, number,
            const.GET_, self.ENTRY_NUMELEMS_, ctypes.byref(count))
        return count.value

    def type(self, number, new_type=None):
        """Find or change the CDF type of a particular Entry number
        
        @param number: number of Entry to check or change
        @type number: int
        @param new_type: type to change it to, see L{const}
        @type new_type: ctypes.c_long
        @return: CDF variable type, see L{const}
        @rtype: int
        @note: If changing types, old and new must be equivalent, see CDF
               User's Guide section 2.5.5 pg. 57
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        if new_type != None:
            if not hasattr(new_type, 'value'):
                new_type = ctypes.c_long(new_type)
            size = ctypes.c_long(self._entry_len(number))
            self._call(const.SELECT_, self.ENTRY_, number,
                       const.PUT_, self.ENTRY_DATASPEC_, new_type, size)
        cdftype = ctypes.c_long(0)
        self._call(const.SELECT_, self.ENTRY_, number,
                   const.GET_, self.ENTRY_DATATYPE_, ctypes.byref(cdftype))
        return cdftype.value

    def has_entry(self, number):
        """Check if this attribute has a particular Entry number

        @param number: number of Entry to check
        @type number: int
        @return: True if L{number} is a valid entry number; False if not
        @rtype: bool
        """
        status = self._call(const.CONFIRM_, self.ENTRY_EXISTENCE_,
                            ctypes.c_long(number),
                            ignore=(const.NO_SUCH_ENTRY, ))
        return not status == const.NO_SUCH_ENTRY

    def max_idx(self):
        """Maximum index of Entries for this Attr

        @return: maximum Entry number
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, self.ATTR_MAXENTRY_, ctypes.byref(count))
        return count.value

    def new(self, data, cdftype=None, number=None):
        """Create a new Entry in this Attribute

        @param data: data to put in the Entry
        @param cdftype: type of the new Entry (otherwise guessed from L{data})
        @param number: Entry number to write, default is lowest available number.
        @note: Will overwrite an existing Entry.
        """
        if number == None:
            number = 0
            while self.has_entry(number):
                number += 1
        (dims, types, elements) = _Hyperslice.types(data)
        if cdftype == None:
            cdftype = types[0]
        elif hasattr(cdftype, 'value'):
            cdftype = cdftype.value
        self._write_entry(number, data, cdftype, dims, elements)
                
    def number(self):
        """Find the attribute number for this attribute

        @return: attribute number
        @rtype: int
        """
        no = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.ATTR_NUMBER_,
                             self._name, ctypes.byref(no))
        return no.value

    def global_scope(self):
        """Determine scope of this attribute.

        @return: True if global (i.e. gAttr)
                 False if zAttr
        @rtype: bool
        """
        scope = ctypes.c_long(0)
        self._call(const.GET_, const.ATTR_SCOPE_, ctypes.byref(scope))
        if scope.value == const.GLOBAL_SCOPE.value:
            return True
        elif scope.value == const.VARIABLE_SCOPE.value:
            return False
        else:
            raise CDFError(const.BAD_SCOPE)

    def rename(self, new_name):
        """Rename this attribute

        Renaming a zAttribute renames it for I{all} zVariables in this CDF!

        @param new_name: the new name of the attribute
        @type new_name: str
        """
        try:
            enc_name = new_name.encode('ascii')
        except AttributeError:
            enc_name = new_name
        if len(enc_name) > const.CDF_ATTR_NAME_LEN256:
            raise CDFError(const.BAD_ATTR_NAME)
        self._call(const.PUT_, const.ATTR_NAME_, enc_name)
        self._name = enc_name

    def _get_entry(self, number):
        """Read an Entry associated with this L{Attr}

        @param number: number of Entry to return
        @type number: int
        @return: data from entry numbered L{number}
        @rtype: list or str
        """
        #Make a big enough buffer
        length = self._entry_len(number)
        cdftype = self.type(number)
        if cdftype in (const.CDF_CHAR.value, const.CDF_UCHAR.value):
            buff = ctypes.create_string_buffer(length)
        else:
            try:
                buff = (lib.ctypedict[cdftype] * length)()
            except KeyError:
                raise CDFError(const.BAD_DATA_TYPE)

        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        self._call(const.SELECT_, self.ENTRY_, number,
                   const.GET_, self.ENTRY_DATA_, ctypes.byref(buff))

        #decode
        if cdftype in (const.CDF_CHAR.value, const.CDF_UCHAR.value):
            result = buff.value
            if not isinstance(result, str):
                result = result.decode()
        else:
            if cdftype == const.CDF_EPOCH.value:
                result = [lib.epoch_to_datetime(item) for item in buff]
            elif cdftype == const.CDF_EPOCH16.value:
                result = [lib.epoch16_to_datetime(item[:]) for item in buff]
            else:
                #subscripting c_array usually returns Python type, not ctype
                result = [item for item in buff]
            if length == 1:
                result = result[0]

        return result

    def _write_entry(self, number, data, cdftype, dims, elements):
        """Write an Entry to this Attr.

        @param number: number of Entry to write
        @type number: int
        @param data: data to write
        @param cdftype: the CDF type to write, from L{const}
        @param dims: dimensions of L{data}
        @type dims: list
        @param elements: number of elements in L{data}, 1 unless it is a string
        @type elements: int
        """
        if isinstance(data, str_classes):
            buff = ctypes.create_string_buffer(elements)
            buff.value = data
        else:
            if len(dims) == 0:
                elements = 1
                data = [data]
            else:
                elements = dims[0]
            buff = (lib.ctypedict[cdftype] * elements)()
            if cdftype == const.CDF_EPOCH16.value:
                for j in range(elements):
                    buff[j][:] = lib.datetime_to_epoch16(data[j])
            elif cdftype == const.CDF_EPOCH.value:
                for j in range(elements):
                    buff[j] = lib.datetime_to_epoch(data[j])
            else:
                for j in range(elements):
                    buff[j] = data[j]
        self._call(const.SELECT_, self.ENTRY_, ctypes.c_long(number),
                   const.PUT_, self.ENTRY_DATA_, ctypes.c_long(cdftype),
                   ctypes.c_long(elements), ctypes.byref(buff))

    def _delete(self):
        """Delete this Attribute

        Also deletes all Entries associated with it.
        """
        self._call(const.DELETE_, const.ATTR_)


class zAttr(Attr):
    """zAttribute for zVariables within a CDF.

    Do not use directly.
    """

    def __init__(self, *args, **kwargs):
        """Initialize this attribute"""
        self.ENTRY_ = const.zENTRY_
        self.ENTRY_DATA_ = const.zENTRY_DATA_
        self.SCOPE = const.VARIABLE_SCOPE
        self.ENTRY_EXISTENCE_ = const.zENTRY_EXISTENCE_
        self.ATTR_NUMENTRIES_ = const.ATTR_NUMzENTRIES_
        self.ATTR_MAXENTRY_ = const.ATTR_MAXzENTRY_
        self.ENTRY_NUMELEMS_ = const.zENTRY_NUMELEMS_
        self.ENTRY_DATATYPE_ = const.zENTRY_DATATYPE_
        self.ENTRY_DATASPEC_ = const.zENTRY_DATASPEC_
        super(zAttr, self).__init__(*args, **kwargs)


class gAttr(Attr):
    """Global Attribute for a CDF

    Represents a CDF attribute, providing access to the gEntries in a format
    that looks like a Python
    list. General list information is available in the python docs:
    U{1<http://docs.python.org/tutorial/introduction.html#lists>},
    U{2<http://docs.python.org/tutorial/datastructures.html#more-on-lists>},
    U{3<http://docs.python.org/library/stdtypes.html#typesseq>}.

    Normally accessed by providing a key to a L{gAttrList}, e.g.::
        attribute = cdffile.attrs['attribute_name']
        first_gentry = attribute[0]

    Each element of the list is a single gEntry of the appropriate type.
    The index to the elements is the gEntry number.

    A gEntry may be either a single string or a 1D array of numerical type.
    Entries of numerical type (everything but CDF_CHAR and CDF_UCHAR)
    with a single element are returned as scalars; multiple-element entries
    are returned as a list. No provision is made for accessing below
    the entry level; the whole list is returned at once (but Python's
    slicing syntax can be used to extract individual items from that list.)

    Multi-dimensional slicing is I{not} supported; an entry with multiple
    elements will have all elements returned (and can thus be sliced itself).
    Example::
        first_three = attribute[5, 0:3] #will fail
        first_three = attribute[5][0:3] #first three elements of 5th Entry

    gEntries are I{not} necessarily contiguous; a gAttribute may have an
    entry 0 and entry 2 without an entry 1. C{len} will return the I{number}
    of gEntries; use L{max_idx} to find the highest defined gEntry number and
    L{has_entry} to determine if a particular gEntry number exists. Iterating
    over all entries is also supported::
        entrylist = [entry for entry in attribute]

    Deleting gEntries will leave a "hole"::
        attribute[0:3] = [1, 2, 3]
        del attribute[1]
        attribute.has_entry(1) #False
        attribute.has_entry(2) #True
        print attribute[0:3] #[1, None, 3]

    Multi-element slices over nonexistent gEntries will return None where
    no entry exists. Single-element indices for nonexistent gEntries will
    raise IndexError. Assigning None to a gEntry will delete it.

    When assigning to a gEntry, the type is chosen to match the data; subject
    to that constraint, it will try to match (in order):
      1. existing gEntry of the same number in this gAttribute
      2. other gEntries in this gAttribute
      3. data-matching constraints described in L{_Hyperslice.types}
    """

    def __init__(self, *args, **kwargs):
        """Initialize this attribute"""
        self.ENTRY_ = const.gENTRY_
        self.ENTRY_DATA_ = const.gENTRY_DATA_
        self.SCOPE = const.GLOBAL_SCOPE
        self.ENTRY_EXISTENCE_ = const.gENTRY_EXISTENCE_
        self.ATTR_NUMENTRIES_ = const.ATTR_NUMgENTRIES_
        self.ATTR_MAXENTRY_ = const.ATTR_MAXgENTRY_
        self.ENTRY_NUMELEMS_ = const.gENTRY_NUMELEMS_
        self.ENTRY_DATATYPE_ = const.gENTRY_DATATYPE_
        self.ENTRY_DATASPEC_ = const.gENTRY_DATASPEC_
        super(gAttr, self).__init__(*args, **kwargs)
