# -*- coding: utf-8 -*-
'''
PyBats!  An open source Python-based interface for reading, manipulating,
and visualizing BATS-R-US and SWMF output.
For more information on the SWMF, please visit the
`Center for Space Environment Modeling <http://csem.engin.umich.edu>`_. 

Introduction
------------

At its most fundamental level, PyBats provides access to output files written
by the Space Weather Modeling Framework and the codes contained within.  
The key task performed by PyBats is loading simulation data into a Spacepy
data model object so that the user can move on to the important tasks of 
analyzing and visualizing the values.  The secondary goal of PyBats is to make
common tasks performed with these data as easy as possible.  The result is that
most SWMF output can be opened and visualized using only a few lines of code.
Many complicated tasks, such as field line integration, is included as well.

Organization
------------

Many output files in the SWMF share a common format.  Objects to handle broad
formats like these are found in the base module.  The base module has classes
to handle SWMF *input* files, as well.

The rest of the classes are organized by code, i.e. classes and functions 
specifically relevant to BATS-R-US can be found in 
:mod:`spacepy.pybats.bats`.  Whenever a certain code included in the SWMF
requires a independent class or a subclass from the PyBats base module, it
will receive its own submodule.

Conventions and Prefixes
------------------------

Nearly every class in PyBats inherits from :class:`spacepy.datamodel.SpaceData`,
so it is important for users to understand how to employ and explore SpaceData
objects.  There are a few exceptions, so always pay close attention to the
docstrings and examples.  Legacy code that does not adhere to this pattern is
slowly being brought up-to-date with each release.

Visualization methods have two prefixes: *plot_* and *add_*.  Whenever a method
begins with *plot_*, a quick-look product will be created that is not highly-
configurable.  These methods are meant to yeild either simple
diagnostic plots or static, often-used products.  There are few methods that use
this prefix.  The prefix *add_* is always followed by *plot_type*; it indicates
a plotting method that is highly configurable and meant to be combined with
other *add_*-like methods and matplotlib commands.

Common calculations, such as calculating Alfven wave speeds of MHD results,
are strewn about PyBats' classes.  They are always given the method prefix
*calc_*, i.e. *calc_alfven*.  Methods called *calc_all* will search for all
class methods with the *calc_* prefix and call them.

Copyright ©2010 Los Alamos National Security, LLC.


Submodules
----------
		
There are submodules for most models included within the SWMF.  The classes
and methods contained within are code-specific, yielding power and
convenience at the cost of flexibility.  A few of the submodules are helper
modules- they are not code specific, but rather provide functionality not
related to an SWMF-included code.

.. autosummary::
   :template: clean_module.rst
   :toctree: autosummary

   bats
   dgcpm
   dipole
   gitm
   kyoto
   pwom
   ram
   rim
   trace2d
   
Top-Level Classes & Functions
-----------------------------

Top-level PyBats classes handle common-format input and output from the SWMF
and are very flexible.  However, they do little beyond open files for the user.

There are several functions found in the top-level module.  These are mostly
convenience functions for customizing plots.

.. rubric:: Classes	    
.. autosummary::
   :template: clean_class.rst
   :toctree: autosummary

   IdlFile
   ImfInput
   LogFile
   NgdcIndex
   PbData
   SatOrbit

.. rubric:: Functions
.. autosummary::
   :template: clean_function.rst
   :toctree: autosummary

   add_body
   add_planet
   parse_tecvars

'''

__contact__ = 'Dan Welling, dwelling@umich.edu'

# Global imports (used ubiquitously throughout this module.
from spacepy.datamodel import dmarray, SpaceData
import spacepy.plot.apionly
import spacepy.plot as spu
import numpy as np

# Some common, global functions.
def parse_filename_time(filename):
    '''
    Given an SWMF file whose name follows the usual standards (see below), 
    attempt to parse the file name to extract information about the iteration,
    runtime and date/time at which the file was written.  All three are
    returned to the caller in that order.  If any cannot be found in the 
    file name, "None" is returned in its place.

    Information on the format of time strings contained within SWMF output
    file names can be found in the `SWMF User Manual 
    <http://herot.engin.umich.edu/~gtoth/SWMF/doc/HTML/SWMF/node225.html>`_.

    
    Parameters
    ==========
    filename : string
        An SWMF output file name.

    Other Parameters
    ================
    None

    Returns
    =======
    i_iter : integer or None
        The iteration at which the file was written, if found.
    runtime : float or None
        The run time (in seconds) at which the file was written, if found.
    time : datetime.datetime or None
        Either a datetime object or *None*, depending on the success of the
        file name parsing.

    Examples
    ========
    >>> parse_filename_time('z=0_mhd_2_t00050000_n00249620.out')
    (249620, 50000.0, None)

    >>> parse_filename_time('mag_grid_e20130924-232600.out')
    (None, None, datetime.datetime(2013, 9, 24, 23, 26))

    '''

    from dateutil.parser import parse
    from datetime import timedelta
    import re

    filename = filename.split('/')[-1]
    
    # Look for date/time:
    if '_e' in filename:
        t_string = re.search('_e(\d{8}\-\d{6})', filename).groups()[-1]
        time = parse(t_string)
        #time += timedelta()
    else:
        time = None

    # Look for run time:
    if '_t' in filename:
        raw     = re.search('_t(\d+)', filename).groups()[-1]
        runtime = 3600*float(raw[:-4]) + 60*float(raw[-4:-2]) + float(raw[-2:])
    else:
        runtime = None

    # Look for file iteration:
    if '_n' in filename:
        i_iter = int(re.search('_n(\d+)', filename).groups()[-1])
    else:
        i_iter = None

    return i_iter, runtime, time

        
def mhdname_to_tex(varname):
    '''
    Convert common MHD variable and unit names into LaTeX-formated strings.
    '''

    import re
    
    match_u = re.search('(.*)u([xyz])', varname)
    match_b = re.match('([bj])([xyz])', varname)
    
    if 'rho' in varname.lower():
        out = r'$\rho_{'+varname.replace('rho', '')+'}$'
    elif varname.lower() == 'nt':
        out = '$nT$'
    elif varname[:2] == 'dB':
        subscript=varname[2]+', '*(len(varname)>3)+varname[3:]
        out = r'$\Delta B _{'+subscript+'}$'
    elif varname.lower()[-1] == 'p':
        out = '$P_{'+varname[:-1]+'}$'
    elif match_u:
        out = '$U_{'+match_u.group(2)+', '*bool(match_u.group(1)) \
              +match_u.group(1)+'}$'
    elif match_b:
        out = '$'+match_b.group(1).upper()+'_{'+match_b.group(2)+'}$'
    elif 'j' == varname[0].lower():
        out = '$J_{'+varname[1:]+'}$'
    else:
        out=varname

    return out

def parse_tecvars(line):
    '''
    Parse the VARIABLES line from a TecPlot-formatted ascii data file.  Create
    a list of name-unit tuples for each variable.
    '''
    
    import re

    # Create return list.
    ret = []

    # Ensure that the input line is a TecPlot variable line.
    if "VARIABLES" not in line:
        raise ValueError('Input line is not a TecPlot VARIABLES line.')

    # Strip out "VARIABLES = "
    line=re.sub('(^\s*VARIABLES\s*\=\s*)|(")', '', line)

    # break into individual vars using commas.
    for s in line.split(','):
        m=re.match('\s*(\w+)\s*(\[(.*)\])*\s*', s)
        ret.append( (m.group(1).lower(), m.group(3)) )

    return ret


def add_planet(ax, rad=1.0, ang=0.0, add_night=True, zorder=1000,
               **extra_kwargs):
    '''
    Creates a circle of ``radius=self.para['rbody']`` and returns the
    MatPlotLib Ellipse patch object for plotting.  If an axis is specified
    using the *ax* keyword, the patch is added to the plot.
    
    Unlike the add_body method, the circle is colored half white (dayside)
    and half black (nightside) to coincide with the direction of the 
    sun. Additionally, because the size of the planet is not intrinsically
    known to the MHD file, the kwarg "rad", defaulting to 1.0, sets the
    size of the planet.  *add_night* can turn off this behavior.
    
    Extra keywords are handed to the Ellipse generator function.

    
    Parameters
    ==========
    ax : Matplotlib Axes object
       Set the axes on which to place planet.
        
    Other Parameters
    ================
    rad : float
       Set radius of planet.  Defaults to 1.
    ang : float
       Set the rotation of the day-night terminator from the y-axis, in degrees.
       Defaults to zero (terminator is aligned with Y-axis.)
    add_night : boolean
       Add night hemisphere.  Defaults to **True**
    zorder : int
       Set the matplotlib zorder of the patch to set how other plot
       elements order with the inner boundary patch. Defaults to 1000,
       nightside patch is given zorder of *zorder+5*.
    
    '''

    from matplotlib.patches import Circle, Wedge

    body = Circle((0,0), rad, fc='w', zorder=zorder, **extra_kwargs)
    arch = Wedge((0,0), rad, 90+ang, -90+ang, fc='k', 
                 zorder=zorder+5, **extra_kwargs)
        
    ax.add_artist(body)
    if add_night: ax.add_artist(arch)

    return body, arch

def add_body(ax, rad=2.5, facecolor='lightgrey', show_planet=True, 
             ang=0.0, add_night=True, zorder=1000, **extra_kwargs):
    '''
    Creates a circle of radius=self.attrs['rbody'] and returns the
    MatPlotLib Ellipse patch object for plotting.  If an axis is specified
    using the "ax" keyword, the patch is added to the plot.
    Default color is light grey; extra keywords are handed to the Ellipse
    generator function.
    
    Because the body is rarely the size of the planet at the center of 
    the modeling domain, add_planet is automatically called.  This can
    be negated by using the show_planet kwarg.
    
    Parameters
    ==========
    ax : Matplotlib Axes object
       Set the axes on which to place planet.
        
    Other Parameters
    ================
    rad : float
       Set radius of the inner boundary.  Defaults to 2.5.
    facecolor : string
       Set color of face of inner boundary circle via Matplotlib color 
       selectors (name, hex, etc.)  Defaults to 'lightgrey'.
    show_planet : boolean
       Turns on/off planet indicator inside inner boundary.  
       Defaults to **True**
    ang : float
       Set the rotation of the day-night terminator from the y-axis, in degrees.
       Defaults to zero (terminator is aligned with Y-axis.)
    add_night : boolean
       Add night hemisphere.  Defaults to **True**
    zorder : int
       Set the matplotlib zorder of the patch to set how other plot
       elements order with the inner boundary patch. Defaults to 1000.
       If a planet is added, it is given a zorder of *zorder*+5.

    '''
    from matplotlib.patches import Ellipse
    
    dbody = 2.0 * rad
    body = Ellipse((0,0),dbody,dbody,facecolor=facecolor, zorder=zorder,
                   **extra_kwargs)

    if show_planet:
        add_planet(ax, ang=ang, add_night=add_night, zorder=zorder+5)

    ax.add_artist(body)
        

def _read_idl_ascii(pbdat, header='units', keep_case=True):
    '''
    Load a SWMF IDL ascii output file and load into a pre-existing PbData
    object.  This should only be called by :class:`IdlFile`.

    The input object must have the name of the file to be opened and read 
    stored in its attributes list and named 'file'.

    The kwarg *header* dictates how the header line will be handled.  In some
    output files, the header is simply the list of units.  In others, it
    contains other data.  If *header* is set to 'units', the header is
    assumed to be a list of units for each variable contained within the
    file.  If set to **None** or not recognized, the header will be saved
    in the object's attribute list under 'header'.

    Parameters
    ==========
    pbdat : PbData object
        The object into which the data will be loaded.

    Other Parameters
    ================
    header : string or **None**
        A string indicating how the header line will be handled; see above.
    keep_case : boolean
        If set to True, the case of variable names will be preserved.  If
        set to False, variable names will be set to all lower case.
    
    Returns
    =======
    True : Boolean
        Returns True on success.

    '''
    
    # Open the file:
    infile = open(pbdat.attrs['file'], 'r')

    # Read the top header line:
    headline = infile.readline().strip()

    # Read & convert iters, runtime, etc. from next line:
    parts = infile.readline().split()
    pbdat.attrs['iter']   = int(parts[0])
    pbdat.attrs['runtime']= float(parts[1])
    pbdat.attrs['ndim']   = int(parts[2])
    pbdat.attrs['nparam'] = int(parts[3])
    pbdat.attrs['nvar']   = int(parts[4])

    # Read & convert grid dimensions.
    grid = [int(x) for x in infile.readline().split()]
    pbdat['grid'] = dmarray(grid)

    # Data from generalized (structured but irregular) grids can be 
    # detected by a negative ndim value.  Unstructured grids (e.g.
    # BATS, AMRVAC) are signified by negative ndim values AND
    # the grid size is always [x, 1(, 1)]
    # Here, we set the grid type attribute to either Regular, 
    # Generalized, or Unstructured.  Let's set that here.
    pbdat['grid'].attrs['gtype'] = 'Regular'
    pbdat['grid'].attrs['npoints']  = abs(pbdat['grid'].prod())
    if pbdat.attrs['ndim'] < 0: 
        if any(pbdat['grid'][1:] > 1): 
            pbdat['grid'].attrs['gtype'] = 'Generalized'
        else:
            pbdat['grid'].attrs['gtype']   = 'Unstructured'
            pbdat['grid'].attrs['npoints'] = pbdat['grid'][0]
    pbdat.attrs['ndim'] = abs(pbdat.attrs['ndim'])

        # Quick ref vars:
    time=pbdat.attrs['runtime']
    gtyp=pbdat['grid'].attrs['gtype']
    npts=pbdat['grid'].attrs['npoints']
    ndim=pbdat['grid'].size
    nvar=pbdat.attrs['nvar']
    npar=pbdat.attrs['nparam']

     # Read parameters stored in file.
    para = np.zeros(npar)
    if npar>0:
        para[:] = infile.readline().split()

    # Read variable names.  Preserve or destroy case
    # via keyword argument selection:
    if keep_case:
        names = infile.readline().split()
    else:
        names = infile.readline().lower().split()

    # Now that we know the number of variables, we can properly handle
    # the headline and units based on the kwarg *header*:
    pbdat.attrs['header']=headline
    if header == 'units':
        # If headline is just units:
        units = headline.split()
    else:
        # If headline is NOT just units, create blank units:
        units = [''] * (len(names)-npar)

    # For some reason, there are often more units than variables
    # in these files.  It looks as if there are more grid units
    # than grid vectors (e.g. 'R R R' implies X, Y, and Z data
    # in file but only X and Y are present.)  Let's try to work
    # around this rather egregious error.
    nSkip=len(units)+npar-len(names)
    if nSkip<0: nSkip=0
          
    # Save grid names (e.g. 'x' or 'r') and save associated params.
    pbdat['grid'].attrs['dims']=tuple(names[0:ndim])
    for name, para in zip(names[(nvar+ndim):], para):
        pbdat.attrs[name]=para
        
    # Create string representation of time.
    pbdat.attrs['strtime']='%4.4ih%2.2im%06.3fs'%\
        (np.floor(time/3600.), np.floor(time%3600. / 60.0),
         time%60.0)
        
    # Create containers for the rest of the data:
    for v, u in zip(names, units[nSkip:]):
        pbdat[v] = dmarray(np.zeros(npts), {'units':u})

    # Load grid points and data:
    for i, line in enumerate(infile.readlines()):
        parts=line.split()
        for j, p in enumerate(parts):
            pbdat[names[j]][i] = p

    # Close the file:
    infile.close()

    # Arrange data into multidimentional arrays if necessary.
    gridnames = names[:ndim]
    if gtyp == 'Irregular':
        for v in names:
            if v not in pbdat: continue
            pbdat[v] = dmarray(np.reshape(pbdat[v], pbdat['grid'], order='F'),
                               attrs=pbdat[v].attrs)
    elif gtyp == 'Regular':
        # Put coords into vectors:
        prod = [1]+pbdat['grid'].cumprod().tolist()
        for i,x in enumerate(pbdat['grid'].attrs['dims']):
            pbdat[x] = dmarray(pbdat[x][0:prod[i+1]-prod[i]+1:prod[i]],
                attrs=pbdat[x].attrs)
        for v in names:
            if v not in pbdat.keys(): continue
            if v not in pbdat['grid'].attrs['dims']:
                pbdat[v] = dmarray(np.reshape(pbdat[v], pbdat['grid'],
                                              order='F'), attrs=pbdat[v].attrs)

def readarray(f, dtype=np.float32, inttype=np.int32):
    '''
    Read an array from an unformatted binary file written out by a 
    Fortran program.

    Parameters
    ==========
    f : Binary file object
        The file from which to read the array of values.

    Other Parameters
    ================
    dtype : data type
        The data type used for data conversion.
    inttype : Numpy integer type
        Set the precision for the integers that store the size of each
        entry.   Defaults to numpy.int32
    
    Returns
    =======
    numpy.array : Numpy array
        Return the data read from file as a 1-dimensional array.

    '''

    if dtype is str:
        dtype_size_bytes=1
    else:
        dtype_size_bytes=dtype.itemsize

    # Get the record length
    rec_len=np.fromfile(f,dtype=inttype,count=1)

    # Check that the record length is consistent with the data type
    if rec_len%dtype_size_bytes!=0:
        raise ValueError('Read error: Data type inconsistent with record' +
                         ' length (data type size is {0:d} bytes, record' +
                         ' length is {1:d} bytes'.format(
                             int(dtype_size_bytes),int(n)))

    if len(rec_len)==0:
        # Zero-length record...file may be truncated
        raise EOFError('Zero-length read at start marker')

    try:
        startpos = f.tell()
    except IOError:
        seekable = False
    else:
        seekable = True
    if seekable:
        f.seek(0, 2)
        endpos = f.tell()
        if endpos - startpos < (int(rec_len[0]) + np.dtype(inttype).itemsize):
            raise EOFError('File is shorter than expected data')
        f.seek(startpos + int(rec_len[0]), 0)
        rec_len_end = np.fromfile(f,dtype=inttype,count=1)
        f.seek(startpos, 0)
        if rec_len_end != rec_len:
            raise ValueError((
                'Read error: End marker length ({0:d}) does not match start '
                'marker length ({1:d}).').format(rec_len[0], rec_len_end[0]) + 
                'This indicates incorrect endiannes, wrong file type, '
                'or file is corrupt.')
    # Read the data
    if dtype is str:
        A=f.read(rec_len[0])
    else:
        A=np.fromfile(f,dtype=dtype,count=int(rec_len[0]/dtype_size_bytes))

    # Check the record length marker at the end
    rec_len_end=np.fromfile(f,dtype=inttype,count=1)
    if len(rec_len_end)==0:
        # Couldn't read, file may be truncated
        raise EOFError('Zero-length read at end marker')

    if rec_len_end!=rec_len:
        # End marker is inconsistent with start marker. Something is wrong.
        raise ValueError(
            'Read error: End marker does not match start marker '+
            '(start marker says record length is {0:d} bytes, end marker says {1:d} bytes).'.format(int(rec_len),int(rec_len_end)) +
            'This indicates incorrect endiannes, wrong file type, '+
            'or file is corrupt'.format(int(rec_len),int(rec_len_end)))

    return A

def _read_idl_bin(pbdat, header='units', keep_case=True, headeronly=False):
    '''
    Load a SWMF IDL binary output file and load into a pre-existing PbData
    object.  This should only be called by :class:`IdlFile`.

    The input object must have the name of the file to be opened and read 
    stored in its attributes list and named 'file'.

    The kwarg *header* dictates how the header line will be handled.  In some
    output files, the header is simply the list of units.  In others, it
    contains other data.  If *header* is set to 'units', the header is
    assumed to be a list of units for each variable contained within the
    file.  If set to **None** or not recognized, the header will be saved
    in the object's attribute list under 'header'.

    Parameters
    ==========
    pbdat : PbData object
        The object into which the data will be loaded.

    Other Parameters
    ================
    header : string or **None**
        A string indicating how the header line will be handled; see above.
    keep_case : boolean
        If set to True, the case of variable names will be preserved.  If
        set to False, variable names will be set to all lower case.

    Returns
    =======
    True : Boolean
        Returns True on success.

    '''

    # Open, read, and parse the file into numpy arrays.
    # Note that Fortran writes integer buffers around records, so
    # we must parse those as well.
    with open(pbdat.attrs['file'], 'rb') as infile:
        # On the first try, we may fail because of wrong-endianess.
        # If that is the case, swap that endian and try again.
        endian='little'

        inttype=np.dtype(np.int32)
        EndChar='<'
        inttype.newbyteorder(EndChar)

        try:
            headline=readarray(infile,str,np.int32)
        except (ValueError,EOFError):
            endian='big'
            EndChar='>'
            inttype.newbyteorder(EndChar)
            infile.seek(0)
            headline=readarray(infile,str,)
        headline=headline.decode('utf-8')
        
        pbdat.attrs['endian']=endian

        # detect double-precision file.
        pos=infile.tell()

        RecLen=np.fromfile(infile,dtype=inttype,count=1)
        infile.seek(pos)

        # Set data types
        if RecLen > 20:
            floattype=np.dtype(np.float64)
        else:
            floattype=np.dtype(np.float32)
        floattype.newbyteorder(EndChar)

        # Parse rest of header
        header_fields_dtype=np.dtype([
            ('it',np.int32),('t',floattype),('ndim',np.int32),
            ('npar',np.int32),('nvar',np.int32)])
        header_fields_dtype.newbyteorder(EndChar)

        (pbdat.attrs['iter'], pbdat.attrs['runtime'],
         pbdat.attrs['ndim'], pbdat.attrs['nparam'], pbdat.attrs['nvar']) = \
            readarray(infile,
                      dtype=header_fields_dtype,
                      inttype=inttype)[0]

        # Get gridsize
        pbdat['grid']=dmarray(readarray(infile,inttype,inttype))

        # Data from generalized (structured but irregular) grids can be
        # detected by a negative ndim value.  Unstructured grids (e.g.
        # BATS, AMRVAC) are signified by negative ndim values AND
        # the grid size is always [x, 1(, 1)]
        # Here, we set the grid type attribute to either Regular,
        # Generalized, or Unstructured.  Let's set that here.
        pbdat['grid'].attrs['gtype'] = 'Regular'
        pbdat['grid'].attrs['npoints']  = abs(pbdat['grid'].prod())
        if pbdat.attrs['ndim'] < 0:
            if any(pbdat['grid'][1:] > 1):
                pbdat['grid'].attrs['gtype'] = 'Generalized'
            else:
                pbdat['grid'].attrs['gtype']   = 'Unstructured'
                pbdat['grid'].attrs['npoints'] = pbdat['grid'][0]
        pbdat.attrs['ndim'] = abs(pbdat.attrs['ndim'])

        # Quick ref vars:
        time=pbdat.attrs['runtime']
        gtyp=pbdat['grid'].attrs['gtype']
        npts=pbdat['grid'].attrs['npoints']
        ndim=pbdat['grid'].size
        nvar=pbdat.attrs['nvar']
        npar=pbdat.attrs['nparam']

        # Read parameters stored in file.
        para  = np.zeros(npar)
        if npar>0:
            para[:] = readarray(infile,floattype,inttype)

        names = readarray(infile,str,inttype).decode('utf-8')

        # Preserve or destroy original case of variable names:
        if not keep_case: names = names.lower()

        names.strip()
        names = names.split()

        # Now that we know the number of variables, we can properly handle
        # the headline and units based on the kwarg *header*:
        pbdat.attrs['header']=headline
        if header == 'units':
            # If headline is just units:
            units = headline.split()
        else:
            # If headline is NOT just units, create blank units:
            units = [''] * (len(names)-npar)
        
        # For some reason, there are often more units than variables
        # in these files.  It looks as if there are more grid units
        # than grid vectors (e.g. 'R R R' implies X, Y, and Z data
        # in file but only X and Y are present.)  Let's try to work
        # around this rather egregious error.
        nSkip = len(units)+npar-len(names)
        if nSkip<0: nSkip = 0
        
        # Save grid names (e.g. 'x' or 'r') and save associated params.
        pbdat['grid'].attrs['dims']=tuple(names[0:ndim])
        for name, para in zip(names[(nvar+ndim):], para):
            pbdat.attrs[name]=para
                
                
        # Create string representation of time.
        pbdat.attrs['strtime'] = '{0:04d}h{1:02d}m{2:06.3f}s'.format(
            int(time//3600), int(time%3600//60), time%60)

        # Get the grid points...
        prod = [1] + pbdat['grid'].cumprod().tolist()

        # Read the data into a temporary array
        griddata = readarray(infile,floattype,inttype)
        for i in range(0,ndim):
            # Get the grid coordinates for this dimension
            tempgrid = griddata[npts*i:npts*(i+1)]

            # Unstructred grids get loaded as vectors.
            if gtyp == 'Unstructured':
                pbdat[names[i]] = dmarray(tempgrid)
            # Irregularly gridded items need multidimensional grid arrays:
            elif gtyp == 'Irregular':
                pbdat[names[i]] = dmarray(np.reshape(tempgrid, pbdat['grid'],
                                                     order='F'))
            # Regularly gridded ones need vector grid arrays:
            elif gtyp == 'Regular':
                pbdat[names[i]] = dmarray(np.zeros(pbdat['grid'][i]))
                for j in range(int(pbdat['grid'][i])):
                    pbdat[names[i]][j] = tempgrid[j*int(prod[i])]
            else:
                raise ValueError('Unknown grid type: {0}'.format(pbdat.gridtype))
            # Add units to grid.
            if units: pbdat[names[i]].attrs['units'] = units.pop(nSkip)

        # Get the actual data and sort.
        for i in range(ndim,nvar+ndim):
            pbdat[names[i]] = dmarray(readarray(infile,floattype,inttype))
            if units: pbdat[names[i]].attrs['units']=units.pop(nSkip)
            if gtyp != 'Unstructured':
                # Put data into multidimensional arrays.
                pbdat[names[i]] = pbdat[names[i]].reshape(pbdat['grid'], order='F')

        # Unstructured data can be in any order, so let's sort it.
        if gtyp == 'Unstructured':
            gridtotal = np.zeros(npts)
            offset = 0.0  # The offset ensures no repeating vals while sorting.
            for key in pbdat['grid'].attrs['dims']:
                gridtotal = gridtotal + offset + pbdat[key]
                offset = offset + np.pi/2.0
                SortIndex = np.argsort(gridtotal)
            for key in list(pbdat.keys()):
                if key=='grid': continue
                pbdat[key] = pbdat[key][SortIndex]

class PbData(SpaceData):
    '''
    The base class for all PyBats data container classes.  Inherits from
    :class:`spacepy.datamodel.SpaceData` but has additional methods for quickly 
    exploring an SWMF dataset.

    Just like :class:`spacepy.datamodel.SpaceData` objects, *PbData* objects
    work just like dictionaries except they have special **attr** dictionary
    attributes for both the top-level object and most values.  This means that
    the following syntax can be used to explore a generic *PbData* object:
    
    >>>print obj.keys()
    >>>print obj.attrs
    >>>value = obj[key]

    Printing *PbData* objects will produce a tree of contents and attributes;
    calling ``self.listunits()`` will print all values that have the 'units'
    attribute and the associated units.  Hence, it is often most instructive
    to use the following two lines to quickly learn a *PbData*'s contents:
    
    >>>print obj
    >>>obj.listunits()

    *PbData* is the main organizational tool for Pybats datasets, so the
    information here is applicable to nearly all Pybats classes.
    '''

    def __init__(self, *args, **kwargs):
        super(PbData, self).__init__(*args, **kwargs)  # Init as SpaceData.

    def __repr__(self):
        return 'PyBats data object'

    def __str__(self):
        '''
        Display contents of container.
        '''
        print(type(self))
        self.tree(attrs=True, verbose=True)
        return ''

    def listunits(self):
        '''
        List all variables and associated units.
        '''
        keys=list(self.keys())
        keys.sort()
        length=0
        for key in keys:
            if len(key)>length: length=len(key)
        form="%%%is:%%s"%length
        for key in keys:
            if 'units' in self[key].attrs:
                print(form%(key, self[key].attrs['units']))

    def timeseries_append(self, obj):
        '''
        If *obj* is of the same type as self, and both have a 'time'
        entry, append all arrays within *obj* that have the same size as
        'time' to the corresponding arrays within self.  This is useful
        for combining time series data of consecutive data.
        '''

        # check to make sure that the append appears feasible:
        if type(obj) != type(self):
            raise TypeError('Cannot append type {} to type {}.'.format(
                type(obj), type(self)))
        if 'time' not in obj:
            raise ValueError('No time vector in {} object'.format(type(obj)))
        if 'time' not in self:
            raise ValueError('No time fector in {} object'.format(type(self)))

        npts = self['time'].size

        for v in self:
            # Only combine vectors that are the same size as time and are
            # in both objects.
            if v not in obj: continue
            if self[v].size != npts: continue
            # Append away.
            self[v] = dmarray(np.append(self[v], obj[v]), self[v].attrs)
        
class IdlFile(PbData):
 
    '''
    An object class that reads/parses an IDL-formatted output file from the 
    SWMF and places it into a :class:`spacepy.pybats.PbData` object.

    Usage:
    >>>data = spacepy.pybats.IdlFile('binary_file.out')

    See :class:`spacepy.pybats.PbData` for information on how to explore
    data contained within the returned object.

    This class serves as a parent class to SWMF component-specific derivative
    classes that do more preprocessing of the data before returning the 
    object.  Hence, using this class to read binary files is typically not
    the most efficient way to proceed.  Look for a PyBats sub module that suits
    your specific needs, or use this base object to write your own.

    A note on byte-swapping: PyBats assumes little endian byte ordering because
    this is what most machines use.  However, there is an autodetect feature
    such that, if PyBats doesn't make sense of the first read (a record length
    entry, or RecLen), it will proceed using big endian ordering.  If this
    doesn't work, the error will manifest itself through the "struct" package
    as an "unpack requires a string of argument length 'X'".
    '''

    def __init__(self, filename,format=None,header='units',
                 keep_case=True, *args,**kwargs):
        super(IdlFile, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file']   = filename   # Save file name.
        self.attrs['format'] = format     # Save file format.
        self.read(header, keep_case)   # Read file.

    def __repr__(self):
        return 'SWMF IDL-Binary file "%s"' % (self.attrs['file'])
    
    def read(self, header, keep_case):
        '''
        This method reads an IDL-formatted BATS-R-US output file and places
        the data into the object.  The file read is self.filename which is
        set when the object is instantiation.
        '''

        if self.attrs['format'] is None:
            try:
                _read_idl_bin(self, header=header, keep_case=keep_case)
            except (ValueError, EOFError, MemoryError):
                _read_idl_ascii(self, header=header, keep_case=keep_case)
        elif self.attrs['format'][:3] == 'bin':
            _read_idl_bin(self, header=header, keep_case=keep_case)
        elif self.attrs['format'][:3] == 'asc':
            _read_idl_ascii(self, header=header, keep_case=keep_case)
        else:
            raise ValueError('Unrecognized file format: {}'.format(
                self.attrs['format']))
        
class LogFile(PbData):
    ''' An object to read and handle SWMF-type logfiles.

    *LogFile* objects read and hold all information in an SWMF
    ascii time-varying logfile.  The file is read upon instantiation.  
    Many SWMF codes
    produce flat ascii files that can be read by this class; the most
    frequently used ones have their own classes that inherit from this.

    See :class:`spacepy.pybats.PbData` for information on how to explore
    data contained within the returned object.

    Usage:
    >>>data = spacepy.pybats.LogFile('filename.log')

    =========  ===========================================
    kwarg      Description
    ---------  -------------------------------------------
    starttime  Manually set the start time of the data.
    =========  ===========================================

    Time is handled by Python's datetime package.  Given that time
    may or may not be given in the logfile, there are three options
    for how time is returned:
        
        1. if the full date and time are listed in the file, ``self['time']``
        is an array of datetime objects corresponding to the entries.  The
        starttime kwarg is ignored.
        
        2. If only the runtime (seconds from start of simulation) is given,
        self['time'] is an array of datetime objects that starts from
        the given *starttime* kwarg which defaults to  1/1/1 00:00UT.
        
        3. If neither of the above are given, the time is assumed to 
        advance one second per line starting from either the starttime kwarg
        or from 2000/1/1 00:00UT + the first iteration (if given in file.)
        As you can imagine, this is sketchy at best.

    This time issue is output dependent: some SWMF models place the full date
    and time into the log by default while others will almost never include the
    full date and time.  The variable ``self['runtime']`` contains the more
    generic seconds from simulation start values.
    
    Example usage:

    >>> import spacepy.pybats as pb
    >>> import pylab as plt
    >>> import datetime as dt
    >>> time1 = dt.datetime(2009,11,30,9,0)
    >>> file1 = pb.logfile('satfile_n000000.dat', starttime=time1)
    >>> plt.plot(file1['time'], file1['dst'])

    '''

    import datetime as dt

    def __init__(self, filename, starttime=(2000,1,1,0,0,0), keep_case=True,
                 *args, **kwargs):
        super(LogFile, self).__init__(*args, **kwargs)
        self.attrs['file'] = filename
        self.read(starttime, keep_case)

    def read(self, starttime, keep_case=True):
        '''
        Load the ascii logfile located at self.filename.
        This method is automatically called upon instantiation.
        '''
        import numpy as np
        import datetime as dt

        # Convert starttime from tuple to datetime object.
        if type(starttime) != dt.datetime:
            if len(starttime) != 6:
                raise ValueError('starttime must be a length 6 Tuple ' +
                                 'or a datetime.datetime object')
            starttime=dt.datetime(starttime[0], starttime[1], starttime[2],
                                  starttime[3], starttime[4], starttime[5])

        # Slurp in entire file.
        infile = open(self.attrs['file'], 'r')
        raw = infile.readlines()
        infile.close()

        # Parse the header.
        self.attrs['descrip'] = raw.pop(0)
        raw_names = raw.pop(0)
        if not keep_case: raw_names = raw_names.lower()
        names = raw_names.split()
        loc={}
        # Keep track of in which column each data vector lies.
        for i, name in enumerate(names):
            loc[name] = i

        # Use the header info to set up arrays 
        # for time, iteration, and the data.
        npts = len(raw)
        # If opening an incomplete file, we must skip the last line.
        if len(raw[-1].split()) < len(names):
            npts=npts-1
        self.attrs['npts']=npts

        # Pop time/date/iteration names off of Namevar.
        for key in ['year','mo','dy','hr','mn','sc','msc',
                    't','it','yy','mm','dd','hh','ss','ms']:
            if key in loc: names.pop(names.index(key))

        # Create containers for data:
        time=dmarray(np.zeros(npts, dtype=object))
        runtime=dmarray(np.zeros(npts), attrs={'units':'s'})
        self['iter']=dmarray(np.zeros(npts))
        for name in names:
            self[name] = dmarray(np.zeros(npts))

        for i in range(npts):
            vals = raw[i].split()
            # Set time:
            if 'year' in loc or 'yy' in loc:
                # If "year" or "yy" is listed, we have the full datetime.
                if 'year' in loc:
                    # BATS date format
                    time[i]=(dt.datetime(
                            int(vals[loc['year']]), # Year
                            int(vals[loc['mo']  ]), # Month
                            int(vals[loc['dy']]), # Day
                            int(vals[loc['hr']]), # Hour
                            int(vals[loc['mn']]), # Minute
                            int(vals[loc['sc']]), # Second
                            int(vals[loc['msc']]) * 1000 #microsec
                            ))
                elif 'yy' in loc:
                    # RIM date format
                    time[i]=(dt.datetime(
                            int(vals[1]), # Year
                            int(vals[2]), # Month
                            int(vals[3]), # Day
                            int(vals[4]), # Hour
                            int(vals[5]), # Minute
                            int(vals[6]), # Second
                            int(vals[7]) * 1000 #microsec
                            ))
                diffT = time[i] - time[0]
                runtime[i]=diffT.days*24.0*3600.0 + \
                    diffT.seconds + \
                    diffT.microseconds*1E-6
            elif 't' in loc:
                # If "t" but no "year" is listed, we only have runtime.
                # Check to ensure number of seconds doesn't
                # exceed 24 hours.
                nowsecs = float(vals[loc['t']])
                nowdays = int(nowsecs / (24.0 * 3600.0))
                nowsecs = nowsecs - (24.0 * 3600.0 * nowdays)
                delta = dt.timedelta(\
                    days    = nowdays,
                    seconds = nowsecs
                    )
                newtime = starttime + delta
                time[i]=newtime
                runtime[i]=nowdays*24.0*3600.0 + nowsecs
            elif 'it' in loc:
                # Check to ensure number of seconds doesn't
                # exceed 24 hours.
                nowsecs = float(vals[loc['it']])
                nowdays = int(nowsecs / (24.0 * 3600.0))
                nowsecs = nowsecs - (24.0 * 3600.0 * nowdays)
                delta = dt.timedelta(\
                    days    = nowdays,
                    seconds = nowsecs
                    )
                newtime = starttime + delta
                time[i]=newtime
                runtime[i]=nowdays*24.*3600.+nowsecs
            else:
                time[i]=starttime + dt.timedelta(float(i))
            # Set iteration:
            if 'it' in loc:
                if '*' in vals[loc['it']]:
                    self['iter'][i] = -1
                else:
                    self['iter'][i] = int(vals[loc['it']]) 
            else: self['iter'][i] = i
            # Collect data
            for j, name in enumerate(names):
                self[name][i] = float(vals[loc[name]])

        # Convert time and runtime to dmarrays.
        self['time'] = time
        self['runtime'] = runtime


class NgdcIndex(PbData):
    '''
    Many models incorporated into the SWMF rely on National Geophysical Data
    Center (NGDC) provided index files (especially F10.7 and Kp).  These
    files, albeit simple ascii, have a unique format and expansive header
    that can be cumbersome to handle.  Data files can be obtained from
    http://spidr.ngdc.noaa.gov .

    NgdcIndex objects aid in reading, creating, and visualizing these files.

    Creating an :class:`NgdcIndex` object is simple:
    
    >>> from spacepy import pybats
    >>> obj=pybats.NgdcIndex(filename='ngdc_example.dat')


    Upon instantiation, if *filename* is a valid file AND kwarg *load* is set
    to boolean True, the contents of *filename* are loaded into the object
    and no other work needs to be done.

    If *filename* is False or *load* is False, a blank :class:`NgdcIndex`
    is created for the user to manipulate.  
    The user can set the time array and the ssociated data values to any values 
    desired and use the method *obj.write()* to dump the contents to an NGDC
    formatted input file.  See the documentation for the write method for 
    more details.

    This class is a work-in-progress.  It is especially tuned for SWMF-needs
    and cannot be considered a general function for the handling of generic
    NGDC files.
    
    =========== ============================================================
    Kwarg       Description
    ----------- ------------------------------------------------------------
    filename    Set the input/output file name.
    load        Read file upon instantiation?  Defaults to **True**
    =========== ============================================================

    '''

    def __init__(self, filename=None, load=True, *args, **kwargs):
        
        super(NgdcIndex, self).__init__(*args, **kwargs)
        # Set key information.
        self.attrs['file'] = filename
        self.attrs['comments'] = []
        
        if load:
            self._read()

    def _read(self):
        '''
        Read an existing NGDC file into memory.  Should only be called upon
        instantiation.
        '''
        from dateutil.parser import parse

        # Start by reading the header.
        infile = open(self.attrs['file'], 'r')
        temp = infile.readline()
        while temp != '#'+50*'-'+'\n': #Detect start of file.
            # Skip blank comments, save substantive ones.
            if temp == '#\n': 
                temp = infile.readline()
                continue
            self.attrs['comments'].append(temp)
            temp = infile.readline()
            
        # Skip "data start" marker.
        infile.readline()

        # Continue to read and sort all data in file.
        def read_one_var():
            attrs = {}
            # Get current variable.  Strip out all but name.
            varname = infile.readline().split(':')[-1].strip()
            temp = infile.readline()
            while temp != '#>\n':
                if temp != '#\n':
                    parts = temp[1:].split(':')
                    attrs[parts[0].strip()] = parts[-1].strip()
                temp = infile.readline()

            # Create some dummy arrays.
            t = []
            d = []

            # Skip the sub-header.
            temp = infile.readline()
            temp = infile.readline()

            # Parse data until next variable or EOF.
            while (temp != '#>\n') and (temp != ''):
                parts = temp.split()
                t.append(parse(' '.join(parts[0:2])))
                d.append(float(parts[2]))
                temp = infile.readline()

            # Put the variable together.
            self[varname] = dmarray([t, d], attrs=attrs)

            # Make judgement as to whether there is more data or not.
            return (temp == '#>\n')

        IsData = True
        while IsData:
            IsData = read_one_var()


    def write(self, outfile=False):
        '''
        Write the :class:`NgdcIndex` object to file.  Kwarg *outfile* can be
        used to specify the path of the output file; if it is not set, 
        *self.attrs['file']* is used.  If this is not set, default to 
        "ngdc_index.dat".
        '''
        
        import datetime as dt
    
        if not outfile:
            if self.attrs['file']!=None:
                outfile=self.attrs['file']
            else:
                outfile='ngdc_index.dat'

        out = open(outfile, 'w')
        
        # Write header.
        for c in self.attrs['comments']:
            out.write(c)
        out.write(2*'#\n'+'#'+50*'-'+'\n')
        
        # Write variables.
        for k in self:
            out.write('#>\n')
            out.write('#%s: %s\n' % ('Element', k))
            for a in self[k].attrs:
                out.write('#%s: %s\n' % (a, self[k].attrs[a]))
            out.write('#>\n#yyyy-MM-dd HH:mm value qualifier description\n')
            for i in range(len(self[k][0,:])):
                t = self[k][0,i]; d = self[k][1,i]
                out.write('%04i-%02i-%02i %02i:%02i%7.1f\t""\t""\n' % 
                          (t.year,t.month,t.day,t.hour,t.minute, d))

class ImfInput(PbData):
    '''
    A class to read, write, manipulate, and visualize solar wind upstream
    input files for SWMF simulations.  More about such files can be found
    in the SWMF/BATS-R-US documentation for the \#SOLARWINDFILE command.

    Creating an :class:`ImfInput` object is simple:
    
    >>> from spacepy import pybats
    >>> obj=pybats.ImfInput(filename='test.dat', load=True)

    Upon instantiation, if *filename* is a valid file AND kwarg *load* is set
    to boolean True, the contents of *filename* are loaded into the object
    and no other work needs to be done.

    If *filename* is False or *load* is False, a blank :class:`ImfInput file`
    is created for the user to manipulate.  
    The user can set the time array and the 
    associated data values (see *obj.attrs['var']* for a list) to any values 
    desired and use the method *obj.write()* to dump the contents to an SWMF
    formatted input file.  See the documentation for the write method for 
    more details.

    Like most :mod:`~spacepy.pybats` objects, you may interact with
    :class:`ImfInput` objects as if they were specialized dictionaries.
    Access data like so:
    
    >>> obj.keys()
    ['bx', 'by', 'bz', 'vx', 'vy', 'vz', 'rho', 'temp']
    >>> density=obj['rho']

    Adding new data entries is equally simple so long as you have the values
    and the name for the values::

    >>> import numpy as np
    >>> v = np.sqrt(obj['vx']**2 + obj['vy']**2 + obj['vz']**2)
    >>> obj['v']=v
       
    =========== ============================================================
    Kwarg       Description
    ----------- ------------------------------------------------------------
    filename    Set the input/output file name.
    load        Read file upon instantiation?  Defaults to **True**
    npoints     For empty data sets, sets number of points (default is 0)
    =========== ============================================================
    '''
    
    def __init__(self, filename=False, load=True, npoints=0, *args, **kwargs):
        from numpy import zeros

        # Initialize data object and required attributes.
        super(ImfInput, self).__init__(*args, **kwargs)
        self.attrs['var']= ['bx', 'by', 'bz', 'ux', 'uy', 'uz', 'rho', 'temp']
        self.attrs['std_var']=True
        self.attrs['coor']='GSM'
        self.attrs['satxyz']=[None, None, None]
        self.attrs['zerobx']=False
        self.attrs['reread']=False
        self.attrs['delay']=None
        self.attrs['plane']=[None, None]
        self.attrs['header']=[]
        self['time']=dmarray(zeros(npoints, dtype=object))
            
        # Set Filename.
        if filename:
            self.attrs['file'] = filename
        else:
            self.attrs['file'] = None

        # Load/create data vectors.
        if filename and load:  # Load contents from existing file.
            self.read(filename)
        else:
            units   = ['nT', 'nT', 'nT', 'km/s', 'km/s', 'km/s', 'cm^-3', 'K']
            for i, key in enumerate(self.attrs['var']):
                self[key]=dmarray(zeros(npoints), attrs={'units':units[i]})

        # Determine the density variable, which can either be "n" or "rho".
        if "n" in self.attrs['var']:
            self._denvar="n"
        elif "rho" in self.attrs['var']:
            self._denvar="rho"
        else:
            raise ValueError('Could not find density variable in file.')
        
            
        self.calc_pram()


    def calc_pram(self):
        '''
        Calculate ram pressure from SW conditions.  Output units in nPa.
        If object was instantiated via an existing imf file, this value
        is calculated automatically.
        '''
        n = self._denvar

        self['pram']=dmarray(self['ux']**2.*self[n]*1.67621E-6,{'units':'nPa'})


    def calc_u(self):
        '''
        Calculate the magnitude of the total solar wind bulk velocity.  Store
        internally as self['u'].
        '''

        self['u'] = dmarray( np.sqrt(self['ux']**2+self['uy']**2+self['uz']**2),
                             {'units':'km/s'} )

        return True

    def calc_b(self):
        '''
        Calculate the magnitude of the IMF in nT.  Store as self['b'].
        '''
        
        self['b'] = dmarray( np.sqrt(self['bx']**2+self['by']**2+self['bz']**2),
                             {'units':'nT'} )
        
        return True

    def calc_alf(self):
        '''
        Calculate the solar wind Alfven speed in $km/s$.  Result is stored
        internally as self['vAlf']
        '''

        if 'b'    not in self: self.calc_b()

        # Const: nT->T, m->km, mu_0, proton mass, cm-3->m-3.
        const = 1E-12/np.sqrt(4.*np.pi*10**-7*1.67E-27*100**3)
        
        self['vAlf']=dmarray(const*self['b']/np.sqrt(self['rho']),
                            {'units':'km/s'})

        return True
        
    def calc_alfmach(self):
        '''
        Calculate the Alvenic Mach number and save as self['machA'].
        Units default to $km/s$.
        '''

        if 'vAlf' not in self: self.calc_alf()
        if 'u'    not in self: self.calc_u()
        self['machA']=dmarray(self['u']/self['vAlf'], {'units':None})

        return True
        
    def varcheck(self):
        '''
        Ensure that the variable list, which gives the order of the
        variables as they are written to file, correctly corresponds
        to the variables stored in the ImfInput object.

        Returns True on success, returns warnings and False on fail.
        '''
        # Convenience:
        var=self.attrs['var']
        key=list(self.keys())
        key.remove('time')

        # Number of variables check:
        if len(var) > len(key):
            print('Not enough variables in IMF object:')
            print('\t%i listed, %i actual.\n' % (len(var),len(key)))
            return False
        # Each variable corresponds to only one in the dict 
        # and occurs only once:
        for v in var:
            if v not in key:
                print('Variable %s listed but does not exist.\n' % v)
                return False
            if var.count(v) != 1:
                print('Variable %s listed multiple times.\n' % v)
                return False
        # Success!
        return True

    def read(self, infile):
        '''
        Read an SWMF IMF/solar wind input file into a newly
        instantiated imfinput object.
        '''
        import datetime as dt
        
        in_header = True
        with open(infile, 'r') as f:
            while True:
                line = f.readline()
                # All non-blank lines before first Param are header
                if in_header:
                    if line.strip() and line[0] != '#':
                        self.attrs['header'].append(line)
                        continue
                    else:
                        in_header = False
                # Parse all params.
                # Grab line, continue if it's not a Param.
                param = line.strip()
                if not param or param[0] != '#':
                    continue
                # For all possible Params, set object attributes/info.
                if param[:5] == '#COOR':
                    self.attrs['coor'] = f.readline()[0:3]
                elif param[:7] == '#REREAD':
                    self.attrs['reread']=True
                elif param[:7] == '#ZEROBX':
                    setting = f.readline()[0]
                    self.attrs['zerobx'] = (setting == 'T')
                elif param[:4] == '#VAR':
                    self.attrs['var'] = f.readline().split()
                    self.attrs['std_var']=False
                elif param[:6] == '#PLANE':
                    xp = float(f.readline().split()[0])
                    yp = float(f.readline().split()[0])
                    self.attrs['plane']=[xp, yp]
                elif param[:9] == '#POSITION':
                    yp = float(f.readline().split()[0])
                    zp = float(f.readline().split()[0])
                    self.attrs['satxyz'][1:]=(yp, zp)
                elif param[:13] == '#SATELLITEXYZ':
                    xp = float(f.readline().split()[0])
                    yp = float(f.readline().split()[0])
                    zp = float(f.readline().split()[0])
                    self.attrs['satxyz']=[xp, yp, zp]
                elif param[:10] == '#TIMEDELAY':
                    self.attrs['delay']=float(f.readline().split()[0])
                elif param[:6] == '#START':
                    break
                else:
                    raise Exception('Unknown file parameter: ' + param)

            # Read data
            indata = np.fromfile(f, sep=' ').reshape(
                -1, 7 + len(self.attrs['var']))
        npoints = indata.shape[0]
        # Create containers for data.
        self['time'] = dmarray(np.empty(npoints, dtype=object))
        for key in self.attrs['var']:
            self[key] = dmarray(np.empty(npoints, dtype=np.float64))
        indata[:, 6] *= 1000 # to microseconds
        self['time'][:] = np.frompyfunc(dt.datetime, 7, 1)(
            *np.require(indata[:, 0:7], dtype=np.int).transpose())
        for i, name in enumerate(self.attrs['var']):
            self[name][:] = indata[:, i + 7]

    def write(self, outfile=False):
        '''
        Write the :class:`ImfInput` object to file.  Kwarg *outfile* can be
        used to specify the path of the output file; if it is not set, 
        *self.attrs['file']* is used.  If this is not set, default to 
        "imfinput.dat".
        '''
        
        import datetime as dt
    
        # Check that the var attribute agrees with the data dictionary.
        if not self.varcheck():
            raise Exception('Number of variables does not match variable order.')

        if not outfile:
            if self.attrs['file']!=None:
                outfile=self.attrs['file']
            else:
                outfile='imfinput.dat'

        with open(outfile, 'wb') as out:
        
            # Convenience variable:
            var=self.attrs['var']

            # Write the header:
            out.write('File created on {}\n'.format(dt.datetime.now().isoformat())
                      .encode())
            for head in self.attrs['header']:
                out.write(head.encode())

            # Handle Params:
            if self.attrs['coor']:
                out.write('#COOR\n{}\n\n'.format(self.attrs['coor']).encode())
            if self.attrs['zerobx']:
                out.write(b'#ZEROBX\nT\n\n')
            if self.attrs['reread']:
                out.write(b'#REREAD')
            if not self.attrs['std_var']:
                out.write('#VAR\n{}\n\n'.format(' '.join(var)).encode())
            if self.attrs['satxyz'].count(None)<3:
                xyz = self.attrs['satxyz']
                if (xyz[0]==None) and (None not in xyz[1:]):
                    out.write('#POSITION\n{0[1]:-6.2f}\n{0[2]:-6.2f}\n\n'
                              .format(xyz).encode())
                elif None not in xyz:
                    out.write('#SATELLITEXYZ\n{}\n'.format(
                        ''.join("{:-6.2f}\n".format(n) for n in xyz)).encode())
            if self.attrs['delay']:
                out.write('#DELAY\n{:-9.2f}\n\n'.format(self.attrs['delay'])
                          .encode())
            if None not in self.attrs['plane']:
                out.write('#PLANE\n{}\n'.format(
                    ''.join('{:-6.2f}\n'.format(n)
                            for n in self.attrs['plane'])).encode())

            # Write the data:
            out.write(b'\n#START\n')
            # Round time to millisecond and format it
            timestr = np.vectorize(
                lambda t: (t.replace(microsecond=0)
                + dt.timedelta(microseconds=int(round(t.microsecond, -3))))
                .strftime('%Y %m %d %H %M %S %f')[:-3] + ' ',
                otypes=[bytes])(self['time'])
            outarray = np.column_stack([timestr] + [
                    np.char.mod('%10.2f', self[key]) for key in var])
            np.savetxt(out, outarray, delimiter=' ', fmt='%s')

    def add_pram_bz(self, target=None, loc=111, pcol='#CC3300', bcol='#3333CC',
                    xlim=None, plim=None, blim=None, epoch=None):
        '''
        Plot, on a single set of axes, both ram pressure and IMF Bz for a 
        quick-look summary of an event. 

        The figure and main axes object are returned.

        Parameters
        ==========

        Other Parameters
        ================
        target : Figure or Axes
            If None (default), a new figure is generated from scratch.
            If a matplotlib Figure object, a new axis is created
            to fill that figure.
            If a matplotlib Axes object, the plot is placed
            into that axis.
        loc : int
            Use to specify the subplot placement of the axis
            (e.g. loc=212, etc.) Used if target is a Figure or None.
            Default 111 (single plot).
        pcol : matplotlib color specifier
            Set the line and label color for ram pressure.
        bcol : matplotlib color specifier
            Set the line and label color for IMF Bz.
        xlim : 2-element list of datetimes
            Set the time range for the plot.  Defaults to full range.
        plim : 2-element list
            Set the y-axis range for dynamic pressure.
        blim : 2-element list
            Set the y-axis range for IMF Bz.
        epoch : datetime.datetime
            Create an epoch marker (vertical line) at time=epoch.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object
        '''

        import matplotlib.pyplot as plt
        from spacepy.plot import set_target, applySmartTimeTicks
        
        # Set ax and fig based on given target.
        fig, a1 = set_target(target, figsize=(8,4), loc=loc)
        a2 = a1.twinx()
        
        self.calc_pram()

        # Plot lines:
        a1.plot(self['time'], self['bz'],   color=bcol)
        a2.plot(self['time'], self['pram'], color=pcol, lw=1.5)

        # Restrict x-range:
        if not xlim: xlim=self['time']
        applySmartTimeTicks(a1, xlim, dolabel=True)

        # Zero line for IMF Bz:
        a1.hlines(0, xlim[0], xlim[-1], color=bcol, linestyles='dashed')

        # Label axes:
        a1.set_ylabel('IMF B$_{Z}$ ($nT$)', color=bcol, size=16)
        a2.set_ylabel('P$_{dyn}$ ($nPa$)',  color=pcol, size=16)
        plt.setp(a1.get_yticklabels(), color=bcol)
        plt.setp(a2.get_yticklabels(), color=pcol)

        # Add epoch marker.
        if epoch:
            ymin, ymax = a1.get_ylim()
            a1.vlines(epoch, ymin, ymax, linestyles='solid', color='g', 
                      linewidths=2.)
            a1.set_ylim([ymin, ymax])

        return fig, a1
        
    def quicklook(self, timerange=None):
        '''
        Generate a quick-look plot of solar wind conditions driving the
        SWMF.  Default values show IMF, number density, and Earthward velocity.
        Returns a figure object containing the plots.
        '''

        import matplotlib.pyplot as plt
        from spacepy.plot import applySmartTimeTicks
        
        if not timerange:
            timerange = [self['time'][0], self['time'][-1]]

        def adjust_plots(ax, ylab, xlab=False, Zero=True):
            ax.grid(True)
            applySmartTimeTicks(ax,timerange)
            ax.set_ylabel(ylab)
            labels =ax.get_yticklabels()
            labels[-1].set_visible(False)
            labels[0].set_visible(False)
            if Zero:
                ax.plot(timerange, [0,0], 'k--')
            if xlab:
                ax.set_xlabel('Universal Time from %s' % 
                              timerange[0].isoformat())
            else:
                ax.xaxis.set_ticklabels([])

        fig = plt.figure(figsize=(8,10))
        fig.subplots_adjust(hspace=0.025, top=0.95, bottom=0.05, right=0.95)
        
        a1 = fig.add_subplot(511)
        a1.plot(self['time'], self['bx'], lw=1.25, c='#003366')
        adjust_plots(a1, 'IMF $B_{X}$ ($nT$)')
        a1.set_title('Solar Wind Drivers (%s Coordinates)' 
                     % (self.attrs['coor']))

        a2 = fig.add_subplot(512)
        a2.plot(self['time'], self['by'], lw=1.25, c='#333399')
        adjust_plots(a2, 'IMF $B_{Y}$ ($nT$)')

        a3 = fig.add_subplot(513)
        a3.plot(self['time'], self['bz'], lw=1.25, c='#0033CC')
        adjust_plots(a3, 'IMF $B_{Z}$ ($nT$)')

        a4 = fig.add_subplot(514)
        a4.plot(self['time'], self[self._denvar], lw=1.25, c='red')
        adjust_plots(a4, 'Density ($cm^{-3}$)', Zero=False)

        a5 = fig.add_subplot(515)
        a5.plot(self['time'], -1.0*self['ux'], lw=1.25, c='green')
        adjust_plots(a5, '$|V_{X}|$ ($km/s$)', Zero=False, xlab=True)

        return fig

class SatOrbit(PbData):
    '''
    An class to load, read, write, and handle BATS-R-US satellite orbit
    input files.  These files are used to fly virtual satellites through
    the MHD domain.  Note that the output files should be handled by 
    the :class:`LogFile` and not this satorbit object.  The object's 
    required and always present attributes are:

    ============ ==============================================================
    Attribute    Description
    ============ ==============================================================
    head         A list of header lines for the file that contain comments.
    coor         The three-letter code (see SWMF doc) of the coord system.
    file         Location of the file to read/write.
    ============ ==============================================================


    The object should always have the following two data keys:

    ============ ==============================================================
    Key          Description
    ============ ==============================================================
    time         A list or numpy vector of datetime objects
    xyz          A 3 x len(time) numpy array of x,y,z coordinates associated
                 with the time vector.
    ============ ==============================================================

    A "blank" instantiation will create an empty object for the user to fill.
    This is desired if the user wants to create a new orbit, either from
    data or from scratch, and save it in a properly formatted file.  Here's
    an example with a "stationary probe" type orbit where the attributes are
    filled and then dropped to file:

    >>> from spacepy.pybats import SatOrbit
    >>> import datetime as dt
    >>> import numpy as np
    >>> sat = SatOrbit()
    >>> sat['time'] = [ dt.datetime(2000,1,1), dt.datetime(2000,1,2) ]
    >>> pos = np.zeros( (3,2) )
    >>> pos[:,0]=[6.6, 0, 0]
    >>> pos[:,1]=[6.6, 0, 0]
    >>> sat['xyz'] = pos
    >>> sat.attrs['coor'] = 'SMG'
    >>> sat.attrs['file'] = 'noon_probe.dat'
    >>> sat.write()

    If instantiated with a file name, the name is loaded into the object.
    For example,

    >>> sat=SatOrbit('a_sat_orbit_file.dat')

    ...will populate all fields with the data in *a_sat_orbit_file.dat*.
    '''

    def __init__(self, filename=None, *args, **kwargs):
        '''
        Create a satorbit object.  If I{file} is not given, buid an
        empty satorbit object is instantiated for the user to fill
        with values of their own creation.
        '''
        import numpy as np

        super(SatOrbit,self).__init__(*args, **kwargs)

        self.attrs['file'] = filename
        self.attrs['head']=[]
        self.attrs['coor'] = 'GSM'
        self['time'] = dmarray(np.zeros(0, dtype=object))

        if filename:
            try:
                self.read()
            except IOError as reason:
                self['xyz'] = dmarray(np.zeros( (3,1) ))
                raise IOError(reason)
        else:
            # fill with empty stuff.
            self['xyz'] = dmarray(np.zeros( (3,1) ))

    def read(self):
        '''
        Read and parse a satellite orbit input file into the satorbit object.
        '''
        import numpy as np
        import datetime as dt
        # Try to open the file.  If we fail, pass the exception
        # to the caller.
        infile = open(self.attrs['file'], 'r')

        # Slurp contents
        raw = infile.readlines()
        infile.close()

        # Parse header.  Save coordinate system and any comments.
        while 1:
            line = raw.pop(0).strip()
            if line:
                if line == '#START':
                    break
                elif line == '#COOR':
                    self.attrs['coor'] = raw.pop(0).strip()
                else:
                    self.attrs['head'].append(line)

        # Read and store all values.
        npts = len(raw)
        self['xyz'] = dmarray(np.zeros( (3,npts) ))
        self['time']=dmarray(np.zeros(npts, dtype=object))
        for i, line in enumerate(raw):
            parts = line.split()
            self['time'][i]=dt.datetime(
                    int(parts[0]), #year
                    int(parts[1]), #month
                    int(parts[2]), #day
                    int(parts[3]), #hour
                    int(parts[4]), #min
                    int(parts[5]), #sec
                    int(parts[6]) * 1000 #micro seconds
                    )
            self['xyz'][:,i] = parts[7:]

    def write(self):
        '''Write a L{satorbit object<pybats.satorbit>} to file using the 
        correct SWMF-input format.  The file will be written to the 
        location specified by self.filename.  An error is thrown if this
        attribute is not set.
        '''

        try:
            outfile = open(self.attrs['file'], 'w')
        except:
            raise Exception('Could not open self.filename!')
        
        # Start by writing header, coordinate system, and then #START.
        for line in self.attrs['head']:
            outfile.write(line+'\n')
        outfile.write('\n')
        outfile.write('#COOR\n{}\n\n'.format(self.attrs['coor']))
        outfile.write('#START\n')

        # Write the rest of the orbit.
        npts=len(self['time'])
        for i in range(npts):
            #Time:
            outfile.write('{:%Y %m %d %H %M %S} {:03.0f} '.format(
                self['time'][i], self['time'][i].microsecond/1000.0) )
            #Position:
            outfile.write('{:13.7E} {:13.7E} {:13.7E}\n'.format(
                self['xyz'][0,i], self['xyz'][1,i], self['xyz'][2,i]) )

        outfile.close()

