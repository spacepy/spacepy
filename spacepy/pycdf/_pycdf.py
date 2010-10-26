#!/usr/bin/env python

"""pyCDF implementation

This module contains the implementation details of pyCDF, the
python interface to U{NASA's CDF library<http://cdf.gsfc.nasa.gov/>}.
"""

__version__ = '0.5'
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

    @ivar _library: C{ctypes} connection to the library
    @type _library: ctypes.WinDLL or ctypes.CDLL
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
        if not hasattr(ignore, '__iter__'):
            ignore = (ignore, )
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
        @keyword check: If True, check returned CDF status and raise exceptions
                        or warnings as necessary (default).
                        If False, do not check.
        @keyword ignore: sequence of CDF statuses to ignore. If any of these
                         is returned by CDF library, any related warnings or
                         exceptions will I{not} be raised.
        @type check: bool
        @return: CDF status from the library
        @rtype: int
        @note: Terminal NULL_ is automatically added to L{args}.
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        if 'check' in kwargs and not kwargs['check']:
            return self._library.CDFlib(*(args + (const.NULL_, )))
        else:
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

    def epoch16_to_datetime(self, epoch):
        """Converts a CDF epoch16 value to a datetime

        @param epoch: epoch16 value from CDF
        @type epoch: list of two floats
        #raise EpochError: if input invalid
        """
        if len(epoch) != 2:
            raise EpochError
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

    def warn(self, level = 3):
        """Issues a warning based on the information stored in my exception

        Intended for use in check_status or similar wrapper function.

        @param level: optional (default 3), how far up the stack the warning
                      should be reported. Passed directly to C{warnings.warn}.
        @type level: int
        """
        warnings.warn(self.__str__(), self.__class__, level)


class EpochError(Exception):
    """Used for errors in epoch routines"""
    pass


class CDF(collections.Mapping):
    """Python object representing a CDF file.

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

    The L{attrs} Python attribute acts as a dictionary referencing CDF
    attributes (do not confuse the two); all the dictionary methods above
    also work on the attribute dictionary. See L{gAttrList} for more on the
    dictionary of global attributes.

    @ivar _handle: file handle returned from CDF library open functions.
    @type _handle: ctypes.c_void_p
    @ivar pathname: filename of the CDF file
    @type pathname: string
    @ivar attrs: All global attributes for this CDF
    @type attrs: L{gAttrList}
    @note: Write support has not been implemented, so CDF is opened
           read-only by default.
    """

    def __init__(self, pathname, masterpath = None):
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
        self._readonly(True)
        self.attrs = gAttrList(self)

    def __del__(self):
        """Destructor

        Close CDF file if there is still a valid handle.
        @note: To avoid data loss, explicitly call L{close} or L{save}.
        """
        if self._handle:
            self.close()

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

    def _majority(self):
        """Finds the majority of this CDF file

        @returns: const.COLUMN_MAJOR or const.ROW_MAJOR
        @rtype: ctypes.clong
        """
        maj = ctypes.c_long(0)
        self._call(const.GET_, const.CDF_MAJORITY_, ctypes.byref(maj))
        if maj.value != const.ROW_MAJOR.value and \
               maj.value != const.COLUMN_MAJOR.value:
            raise CDFError(const.BAD_MAJORITY)
        return maj

    def _new_var(self, var_name, *args):
        """Add a new zVariable to this CDF
        
        @param var_name: Name of the variable.
        @type var_name: string
        @param args: Additional arguments passed to L{Var._create}
        @type args: various
        @return: a new variable named L{var_name}.
        @rtype: L{Var}
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """

        return Var(self, var_name, *args)

    def _readonly(self, ro = None):
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
        @raise IndexError: TODO detail this
        @raise CDFError: for errors from the CDF library
        """
        hslice = _Hyperslice(self, key)
        self._call(const.SELECT_, const.zVAR_RECNUMBER_, ctypes.c_long(hslice.starts[0]),
                  const.SELECT_, const.zVAR_RECCOUNT_, ctypes.c_long(hslice.counts[0]),
                  const.SELECT_, const.zVAR_RECINTERVAL_,
                  ctypes.c_long(hslice.intervals[0]))
        if hslice.dims > 1:
            dims = hslice.dims - 1
            lib.call(const.SELECT_, const.zVAR_DIMINDICES_,
                     (ctypes.c_long * dims)(*hslice.starts[1:]),
                     const.SELECT_, const.zVAR_DIMCOUNTS_,
                     (ctypes.c_long * dims)(*hslice.counts[1:]),
                     const.SELECT_, const.zVAR_DIMINTERVALS_,
                     (ctypes.c_long * dims)(*hslice.intervals[1:]))
        buffer = hslice.create_buffer()
        lib.call(const.GET_, const.zVAR_HYPERDATA_, ctypes.byref(buffer))
        result = hslice.unpack_buffer(buffer)
        
        if self._cdf_type() == const.CDF_EPOCH.value:
            flat = result
            while True:
                try:
                    flat = [item for sublist in flat for item in sublist]
                except TypeError:
                    break
            try:
                for i in range(len(flat)):
                    #NOT a list comprehension, changing values so they'll be
                    #different in result as well
                    flat[i] = lib.epoch_to_datetime(flat[i])
            except TypeError:
                result = lib.epoch_to_datetime(result)
        elif self._cdf_type() == const.CDF_EPOCH16.value:
            flat = result
            if isinstance(flat[0], float): #single value
                result = lib.epoch16_to_datetime(result)
            else:
                while not isinstance(flat[0][0], float):
                    #Flatten list until just two-element floats
                    flat = [item for sublist in flat for item in sublist]
                for i in range(len(flat)):
                    flat[i] = lib.epoch16_to_datetime(flat[i])
            
        return hslice.convert_array(result)

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
        @note: Not intended to be used directly; use L{CDF._new_var}.
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
        majority = self.cdf_file._majority()
        if majority == const.COLUMN_MAJOR:
            sizes.reverse()
        return sizes

    def _rec_vary(self):
        """Gets whether this variable has record variance

        @return: True if record variance, False if NRV
        @rtype: boolean
        @note: If the variance is unknown, True is assumed
               (this replicates the apparent behaviour of the
               CDF library on variable creation).
        """
        vary = ctypes.c_long(0)
        self._call(const.GET_, const.zVAR_RECVARY_, ctypes.byref(vary))
        if vary.value == const.NOVARY.value:
            return False
        else:
            return True

    def _del_recs(self, start, count):
        """Delete a series of records from this variable

        @param start: zero-base index of first record to remove
        @type start: long
        @param count: number of records to remove
        @type count: long
        @raise CDFError: if CDF library reports an error
        @raise CDFWarning: if CDF library reports a warning and interpreter
                           is set to error on warnings.
        """
        #TODO: sanity-checking...
        first = ctypes.c_long(start)
        last = ctypes.c_long(start + count - 1)
        self._call(const.DELETE_, const.zVAR_RECORDS_, first, last)

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
        cdftype = self._cdf_type()
        if cdftype == const.CDF_CHAR.value or cdftype == const.CDF_UCHAR.value:
            return ctypes.c_char * self._nelems()
        try:
            return lib.ctypedict[cdftype]
        except KeyError:
            raise CDFError(const.BAD_DATA_TYPE)

    def _cdf_type(self):
        """Returns the CDF type of this variable

        @return: CDF type
        @rtype: int
        """
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

    def copy(self):
        """Copies all data and attributes from this variable

        @return: list of all data in record order
        @rtype: L{VarCopy}
        """
        return VarCopy(self)


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
        key = self._expand_ellipsis(key)
        result = super(VarCopy, self).__getitem__(key[0])
        for subkey in key[1:]:
            result = result[subkey]
        return result

    def _expand_ellipsis(self, slices):
        """Expands ellipses into the correct number of full-sized slices

        @param slices: tuple of slices, integers, or ellipse objects
        @type slices: tuple
        @return: key with ellipses replaced by appropriate number of blank
                 slices
        @rtype: tuple
        @raise IndexError: if ellipses specified when already have enough
                           dimensions
        """
        if not Ellipsis in slices:
            return slices

        size = len(slices)
        extra = len(self._dims) - size
        if extra < 0:
            raise IndexError('Too many dimensions specified to use ellipsis.')
        idx = slices.index(Ellipsis)
        result = slices[0:idx] + \
                 tuple([slice(None, None, None) for i in range(extra)]) + \
                 slices[idx+1:]
        if Ellipsis in result:
            raise IndexError('Ellipses can only be used once per slice.')
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
    @ivar counts: number of values to get from each dimension
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
    @note: All variables are stored in the majority of the CDF
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
        self.rv = self.zvar._rec_vary()
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
        key = self.expand_ellipsis(key)
        if len(key) == 1 and self.rv: #get all data for this record(s)
            key = self.expand_ellipsis(key + (Ellipsis, ))
        elif len(key) == self.dims - 1 or not self.rv: #NRV is always rec 0
            if zvar._rec_vary(): #get all records
                key = (slice(None, None, None), ) +key
            else: #NRV, so get 0th record (degenerate)
                key = (0, ) + key
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
                    if idx >= self.dimsizes[i] or idx < 0:
                        raise IndexError('list index out of range')
                    self.starts[i] = idx
                    self.degen[i] = True
        else:
            raise IndexError('Slice does not match dimensions for zVar ' +
                             zvar._name)
                
        self.column = (zvar.cdf_file._majority().value ==
                       const.COLUMN_MAJOR.value)

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
        cdftype = self.zvar._cdf_type()
        if cdftype == const.CDF_CHAR.value or cdftype == const.CDF_UCHAR.value:
            buffsize = self.zvar._nelems()
            for count, degen in zip(self.counts[-1::-1], self.degen[-1::-1]):
                if not degen:
                    buffsize *= count            
            return ctypes.create_string_buffer(buffsize)
        else:
            constructor = self.zvar._c_type()
            #Build array from innermost out
            for count, degen in zip(self.counts[-1::-1], self.degen[-1::-1]):
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
        cdftype = self.zvar._cdf_type()
        if cdftype == const.CDF_CHAR.value or cdftype == const.CDF_UCHAR.value:
            finaltype = self.zvar._c_type()
            for count, degen in zip(self.counts[-1::-1], self.degen[-1::-1]):
                if not degen:
                    finaltype *= count
            buffer = ctypes.cast(buffer, ctypes.POINTER(finaltype)).contents
            buffer = self.c_array_to_list(buffer)
            #Flatten the list and convert each one to Unicode if necessary
            flat = buffer
            if not isinstance(flat, str) and not isinstance(flat, bytes):
                while not isinstance(flat[0], str) and not isinstance(flat[0], bytes):
                    flat = [item for sublist in flat for item in sublist]
                for i in range(len(flat)):
                    #NOT a list comprehension, changing values so they'll be
                    #different in result as well
                    if not isinstance(flat[i], str):
                        flat[i] = flat[i].decode()
            else:
                if not isinstance(buffer, str):
                    buffer = buffer.decode()
            return buffer
        else:
            return self.c_array_to_list(buffer)

    def c_array_to_list(self, array):
        """Converts a ctypes array type to a python nested list

        @param array: the array to convert
        @type array: instance of ctypes.Array subclass
        @return: contents of array
        @rtype: list (of lists)
        """
        if hasattr(array, 'value'):
            return array.value
        if not isinstance(array, ctypes.Array):
            return array

        dimsizes = []
        counts = self.counts
        degen = self.degen
        if self.column:
            counts = self.reorder(counts)
            degen = self.reorder(degen)
        for count, degen in zip(counts, degen):
            if not degen:
                dimsizes.append(count)

        if self.zvar._cdf_type() == const.CDF_EPOCH16.value:
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
            return [i.value if hasattr(i, 'value') else i for i in flat]
        
        for i in range(len(dims) - 1):
            size = dims[i]
            n_lists = 1
            for j in dims[i + 1:]:
                n_lists *= j
            if i == 0:
                result = [[k.value if hasattr(k, 'value') else k
                           for k in flat[j * size:(j + 1) * size]]
                          for j in range(n_lists)]
            else:
                result = [result[j * size:(j + 1) * size] for j in range(n_lists)]
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

        try:
            string_classes = (str, bytes, unicode)
        except NameError:
            string_classes = (str, bytes)            

        flat = array #this is progressively flattened (i.e. dims stripped off)
        while True:
            try:
                if isinstance(flat[0], string_classes):
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
        """Convertss a start/stop/step range to start/count/interval

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
        #Normalize step
        if step == 0:
            raise ValueError('slice step cannot be zero')
        if step == None:
            step = 1

        #Everything is different based on direction
        if step < 0:
            rev = True
            step *= -1
            (start, stop) = (stop, start) #swap "beginning" and "end"
            if start == None:
                #didn't specify an end point, start from beginning of string
                start = 0
            elif start >= 0:
                #end points *past* end, so shift start past specified end
                start += 1
                if start > size:
                    start = size
            else:
                #specified a point from END of string
                #translate to point from START, and do the off-by-one shift
                start += size
                start += 1
                if start < 0:
                    start = 0
            if stop == None:
                #didn't specify a start point, go to end of string
                stop = size
            elif stop >= 0:
                #Shift stop by one to get it past the start pointed to
                stop += 1
                if stop > size:
                    stop = size
            else:
                #specified point from END of string
                #translate to point from START, and do off-by-one
                stop += size
                stop += 1
                if stop < 0:
                    stop = 0
            #Convert from start/stop to counts
            count = int((stop - start + step - 1) / step)
            if count < 0:
                count = 0
                start = 0
            #"Nail" the sequence to the STOP, not the start
            start = (stop - 1) - step * (count - 1)
        else:
            rev = False
            if start == None:
                #didn't specify a start point, start from beginning of string
                start = 0
            elif start >= 0:
                if start > size:
                    start = size
            else:
                #specified a point from END of string
                #translate to point from START
                start += size
                if start < 0:
                    start = 0
            if stop == None:
                #didn't specify an end point, go to end of string
                stop = size
            elif stop >= 0:
                if stop > size:
                    stop = size
            else:
                #specified point from END of string
                #translate to point from START
                stop += size
                if stop < 0:
                    stop = 0
            #Convert from start/stop to counts
            #This is basically a poor man's floor
            count = int((stop - start + step - 1) / step)
            if count < 0:
                count = 0
                start = 0
        return (start, count, step, rev)


class gAttrList(collections.Mapping):
    """Object representing I{all} the gAttributes in a CDF.

    Normally access as an attribute of an open L{CDF}::
        global_attribs = cdffile.attrs

    Appears as a dictionary: keys are attribute names, values are
    L{gAttr} objects representing that attribute. Accessing the global
    attribute TEXT::
        text_attr = cdffile.attrs['TEXT']

    @ivar _cdf_file: CDF these attributes are in
    @type _cdf_file: L{CDF}
    """

    def __init__(self, cdf_file):
        """Initialize the attribute collection

        @param cdf_file: CDF these attributes are in
        @type cdf_file: L{CDF}
        """
        self._cdf_file = cdf_file

    def __getitem__(self, name):
        """Find an gAttribute by name

        @param name: name of the gAttribute to return
        @type name: str
        @return: attribute named L{name}
        @rtype: L{gAttr}
        @raise KeyError: if there is no attribute named L{name}
        @raise CDFError: other errors in CDF library
        """
        try:
            attrib = gAttr(self._cdf_file, name)
        except CDFError:
            (t, v, tb) = sys.exc_info()
            if v.status == const.NO_SUCH_ATTR:
                raise KeyError(name + ': ' + str(v))
            else:
                raise
        if not attrib.global_scope():
            raise KeyError(name + ': no gAttribute by that name.')
        return attrib

    def __len__(self):
        """Number of gAttributes in this CDF

        @return: number of gAttributes in the CDF
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMgATTRS_,
                             ctypes.byref(count))
        return count.value

    def __iter__(self, current=0):
        """Iterates over all gAttr in this CDF

        Returns name of one L{gAttr} at a time until reaches the end.
        @note: Returned in number order.
        """
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMATTRS_,
                             ctypes.byref(count))
        while current < count.value:
            candidate = gAttr(self._cdf_file, current)
            if candidate.global_scope():
                value = yield(candidate._name)
                if value != None:
                    current = self[value].number()
            current += 1

    def copy(self):
        """Create a copy of this attribute list

        @return: copy of all entries for all attributes in this list
        @rtype: dict of lists
        """
        return dict((key, value[:]) for (key, value) in self.items())


class zAttrList(collections.Mapping):
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
        self._cdf_file = zvar.cdf_file
        self._zvar = zvar

    def __getitem__(self, name):
        """Find an zAttribute by name

        @param name: name of the zAttribute to return
        @type name: str
        @return: attribute named L{name}
        @rtype: L{zAttr}
        @raise KeyError: if there is no attribute named L{name} associated
                         with this zVariable
        @raise CDFError: other errors in CDF library
        """
        try:
            attrib = zAttr(self._cdf_file, name)
        except CDFError:
            (t, v, tb) = sys.exc_info()
            if v.status == const.NO_SUCH_ATTR:
                raise KeyError(name + ': ' + str(v))
            else:
                raise
        if attrib.global_scope():
            raise KeyError(str(name) + ': no zAttribute by that name.')
        zvar_num = self._zvar._num()
        if attrib.has_entry(zvar_num):
            return attrib[zvar_num]
        else:
            raise KeyError(str(name) + ': no such attribute for variable ' +
                           self._zvar._name)

    def __len__(self):
        """Number of zAttributes in this CDF

        @return: number of zAttributes in the CDF
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

    def __iter__(self, current=0):
        """Iterates over all zAttr in this CDF

        Returns name of one L{zAttr} at a time until reaches the end.
        @note: Returned in number order.
        """
        count = ctypes.c_long(0)
        self._cdf_file._call(const.GET_, const.CDF_NUMATTRS_,
                             ctypes.byref(count))
        while current < count.value:
            candidate = zAttr(self._cdf_file, current)
            if not candidate.global_scope():
                if candidate.has_entry(self._zvar._num()):
                    value = yield(candidate._name)
                    if value != None:
                        current = self[value].number()
            current += 1

    def copy(self):
        """Create a copy of this attribute list

        @return: copy of the entries for all attributes in this list
        @rtype: dict
        """
        return dict((key, value) for (key, value) in self.items())


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

    def __init__(self, cdf_file, attr_name):
        """Initialize this attribute

        @param cdf_file: CDF file containing this attribute
        @type cdf_file: L{CDF}
        @param attr_name: Name of this attribute
        @type attr_name: str
        """
        self._cdf_file = cdf_file
        if isinstance(attr_name, str_classes):
            try:
                self._name = attr_name.encode('ascii')
            except AttributeError:
                self._name = attr_name
            self._cdf_file._call(const.CONFIRM_, const. ATTR_EXISTENCE_,
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

    def _get_entry(self, number):
        """Read an Entry associated with this L{Attr}

        @param number: number of Entry to return
        @type number: int
        @return: data from entry numbered L{number}
        @rtype: list or str
        """
        #Make a big enough buffer
        length = self._entry_len(number)
        cdftype = self.entry_type(number)
        if cdftype in (const.CDF_CHAR.value, const.CDF_UCHAR.value):
            buff = ctypes.create_string_buffer(length)
        else:
            try:
                buff = (lib.ctypedict[cdftype] * length)()
            except KeyError:
                raise CDFError(const.BAD_DATA_TYPE)

        self._entry_data(number, buff)

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


class zAttr(Attr):
    """zAttribute for zVariables within a CDF.
    
    Do not use directly.
    """

    def __len__(self):
        """Number of zEntries for this zAttr. NOT same as max zEntry number.

        @return: Number of zEntries
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.ATTR_NUMzENTRIES_, ctypes.byref(count))
        return count.value

    def max_idx(self):
        """Maximum index of zEntries for this zAttr

        @return: maximum zEntry number
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.ATTR_MAXzENTRY_, ctypes.byref(count))
        return count.value

    def _entry_data(self, number, buff):
        """Reads all data from a zEntry in the CDF

        @param number: number of zEntry to read
        @type number: int
        @param buff: buffer large enough to hold the data
        @type buff: ctypes.Array
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        self._call(const.SELECT_, const.zENTRY_, number,
                   const.GET_, const.zENTRY_DATA_, ctypes.byref(buff))

    def _entry_len(self, number):
        """Number of elements in a zEntry

        @param number: number of zEntry
        @type number: int
        @return: number of elements
        @rtype: int
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        count = ctypes.c_long(0)
        self._call(
            const.SELECT_, const.zENTRY_, number,
            const.GET_, const.zENTRY_NUMELEMS_, ctypes.byref(count))
        return count.value

    def has_entry(self, number):
        """Check if this attribute has a particular Entry number

        @param number: number of zEntry to check
        @type number: int
        @return: True if L{number} is a valid entry number; False if not
        @rtype: bool
        """
        status = self._call(const.CONFIRM_, const.zENTRY_EXISTENCE_,
                            ctypes.c_long(number),
                            ignore=(const.NO_SUCH_ENTRY))
        return not status == const.NO_SUCH_ENTRY

    def entry_type(self, number):
        """Find the CDF type of a particular Entry number
        @param number: number of zEntry to check
        @type number: int
        @return: CDF variable type, see L{const}
        @rtype: int
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        cdftype = ctypes.c_long(0)
        self._call(const.SELECT_, const.zENTRY_, number,
                   const.GET_, const.zENTRY_DATATYPE_, ctypes.byref(cdftype))
        return cdftype.value


class gAttr(Attr):
    """Global Attribute for a CDF
    
    Represents a CDF attribute, providing access to the Entries in a format
    that looks like a Python
    list. General list information is available in the python docs:
    U{1<http://docs.python.org/tutorial/introduction.html#lists>},
    U{2<http://docs.python.org/tutorial/datastructures.html#more-on-lists>},
    U{3<http://docs.python.org/library/stdtypes.html#typesseq>}.

    Normally accessed by providing a key to a L{gAttrList}, e.g.::
        attribute = cdffile.attrs['attribute_name']

    Each element of the list is a single entry of the appropriate type.
    The index to the elements is the entry number.

    Entries of numerical type (i.e. everything but CDF_CHAR and CDF_UCHAR)
    with a single element are returned as scalars; multiple element entries
    are returned as a list. No provision is made for accessing below
    the entry level; the whole list is returned at once (but Python's
    slicing syntax can be used to extract individual items from that list.)

    Multi-dimensional slicing is I{not} supported; an entry with multiple
    elements will have all elements returned (and can thus be sliced itself).
    Example::
        first_three = attribute[5, 0:3] #will fail
        first_three = attribute[5][0:3] #first three elements of 5th Entry
    """

    def __len__(self):
        """Number of gEntries for this zAttr. NOT same as max gEntry number.

        @return: Number of gEntries
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.ATTR_NUMgENTRIES_, ctypes.byref(count))
        return count.value

    def max_idx(self):
        """Maximum index of gEntries for this gAttr

        @return: maximum gEntry number
        @rtype: int
        """
        count = ctypes.c_long(0)
        self._call(const.GET_, const.ATTR_MAXgENTRY_, ctypes.byref(count))
        return count.value

    def _entry_data(self, number, buff):
        """Reads all data from a gEntry in the CDF

        @param number: number of gEntry to read
        @type number: int
        @param buff: buffer large enough to hold the data
        @type buff: ctypes.Array
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        self._call(const.SELECT_, const.gENTRY_, number,
                   const.GET_, const.gENTRY_DATA_, ctypes.byref(buff))


    def _entry_len(self, number):
        """Number of elements in a gEntry

        @param number: number of gEntry
        @type number: int
        @return: number of elements
        @rtype: int
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        count = ctypes.c_long(0)
        self._call(
            const.SELECT_, const.gENTRY_, number,
            const.GET_, const.gENTRY_NUMELEMS_, ctypes.byref(count))
        return count.value

    def has_entry(self, number):
        """Check if this attribute has a particular Entry number

        @param number: number of gEntry to check
        @type number: int
        @return: True if L{number} is a valid entry number; False if not
        @rtype: bool
        """
        status = self._call(const.CONFIRM_, const.gENTRY_EXISTENCE_,
                            ctypes.c_long(number),
                            ignore=(const.NO_SUCH_ENTRY))
        return not status == const.NO_SUCH_ENTRY

    def entry_type(self, number):
        """Find the CDF type of a particular Entry number
        @param number: number of gEntry to check
        @type number: int
        @return: CDF variable type, see L{const}
        @rtype: int
        """
        if not self.has_entry(number):
            raise IndexError('list index ' + str(number) + ' out of range.')
        cdftype = ctypes.c_long(0)
        self._call(const.SELECT_, const.gENTRY_, number,
                   const.GET_, const.gENTRY_DATATYPE_, ctypes.byref(cdftype))
        return cdftype.value
