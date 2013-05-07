# -*- coding: utf-8 -*-
'''
PyBats!  An open source Python-based interface for reading, manipulating,
and visualizing BATS-R-US and SWMF output.
For more information on the SWMF, please visit the
`Center for Space Environment Modeling<http://csem.engin.umich.edu>`_. 

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
:module:`spacepy.pybats.bats`.  Whenever a certain code included in the SWMF
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
'''

__contact__ = 'Dan Welling, dwelling@umich.edu'

# Global imports (used ubiquitously throughout this module.
from spacepy.datamodel import dmarray, SpaceData
import numpy as np

# Some common, global functions.
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

def smart_timeticks(time):
    '''
    Given a list or array of time values, create intelligent timeticks based
    on the max and min of the input.
    Return three items: major tick locator, minor tick locator, and a 
    format string.

    Example:
    >>>Mtick, mtick, fmt = smart_timeticks([time1, time2])

    One could then apply these to an axis object as follows:
    >>>ax.xaxis.set_major_locator(Mtick)
    >>>ax.xaxis.set_minor_locator(mtick)
    >>>ax.xaxis.set_major_formatter(fmt)

    Often, it is much easier to just use 
    :func:`spacepy.pybats.apply_smart_timeticks`.
    '''
    
    import matplotlib.dates as mdt

    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if nHours < 1:
        Mtick=mdt.MinuteLocator(byminute=[0,15,30,45])
        mtick=mdt.MinuteLocator(byminute=range(60), interval=5)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick=mdt.MinuteLocator(byminute=[0,30])
        mtick=mdt.MinuteLocator(byminute=range(60), interval=10)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick=mdt.HourLocator(byhour=range(24), interval=2)
        mtick=mdt.MinuteLocator(byminute=[0,15,30,45])
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 24:
        Mtick=mdt.HourLocator(byhour=[0,3,6,9,12,15,18,21])
        mtick=mdt.HourLocator(byhour=range(24))
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick=mdt.HourLocator(byhour=[0,6,12,18])
        mtick=mdt.HourLocator(byhour=range(24))
        fmt = mdt.DateFormatter('%H:%M UT')
    elif deltaT.days < 15:
        Mtick=mdt.DayLocator(bymonthday=range(1,32))
        mtick=mdt.HourLocator(byhour=[0,6,12,18])
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 32:
        Mtick=mdt.DayLocator(bymonthday=range(5,35,5))
        mtick=mdt.HourLocator(byhour=[0,6,12,18])
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 60:
        Mtick=mdt.MonthLocator()
        mtick=mdt.DayLocator(bymonthday=range(5,35,5))
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 731:
        Mtick=mdt.MonthLocator()
        mtick=mdt.DayLocator(bymonthday=15)
        fmt = mdt.DateFormatter('%b %Y')
    else:
        Mtick=mdt.YearLocator()
        mtick=mdt.MonthLocator(bymonth=7)
        fmt = mdt.DateFormatter('%Y')
    return(Mtick, mtick, fmt)

def apply_smart_timeticks(ax, time, dolimit=True, dolabel=False):
    '''
    Given an axis *ax* and a list/array of datetime objects, *time*, 
    use the smart_timeticks function to built smart time ticks and
    then immediately apply them to the give axis.

    The range of the *time* input value will be used to set the limits
    of the x-axis as well.  Set kwarg *dolimit* to False to override 
    this behavior.

    If you wish to auto-label the X axis with ISO-formatted date and time
    taken from ``time[0]``, use the kwarg *dolabel* (defaults to **False**).
    '''

    Mtick, mtick, fmt = smart_timeticks(time)
    ax.xaxis.set_major_locator(Mtick)
    ax.xaxis.set_minor_locator(mtick)
    ax.xaxis.set_major_formatter(fmt)
    if dolimit:
        ax.set_xlim([time[0], time[-1]])
    if dolabel:
        ax.set_xlabel('Time from %s' % time[0].isoformat())
    return True

def add_planet(ax, rad=1.0, ang=0.0, **extra_kwargs):
    '''
    Creates a circle of ``radius=self.para['rbody']`` and returns the
    MatPlotLib Ellipse patch object for plotting.  If an axis is specified
    using the *ax* keyword, the patch is added to the plot.
    
    Unlike the add_body method, the circle is colored half white (dayside)
    and half black (nightside) to coincide with the direction of the 
    sun. Additionally, because the size of the planet is not intrinsically
    known to the MHD file, the kwarg "rad", defaulting to 1.0, sets the
    size of the planet.
    
    Extra keywords are handed to the Ellipse generator function.
    '''

    from matplotlib.patches import Circle, Wedge

    body = Circle((0,0), rad, fc='w', zorder=1000, **extra_kwargs)
    arch = Wedge((0,0), rad, 90+ang, -90+ang, fc='k', 
                 zorder=1001, **extra_kwargs)
        
    ax.add_artist(body)
    ax.add_artist(arch)

    return body, arch

def add_body(ax, rad=2.5, facecolor='lightgrey', DoPlanet=True, 
             ang=0.0, **extra_kwargs):
    '''
    Creates a circle of radius=self.para['rbody'] and returns the
    MatPlotLib Ellipse patch object for plotting.  If an axis is specified
    using the "ax" keyword, the patch is added to the plot.
    Default color is light grey; extra keywords are handed to the Ellipse
    generator function.
    
    Because the body is rarely the size of the planet at the center of 
    the modeling domain, add_planet is automatically called.  This can
    be negated by using the DoPlanet kwarg.
    '''
    from matplotlib.patches import Ellipse
    
    dbody = 2.0 * rad
    body = Ellipse((0,0),dbody,dbody,facecolor=facecolor, zorder=999,
                   **extra_kwargs)

    if DoPlanet:
        add_planet(ax, ang=ang)

    ax.add_artist(body)
        
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
        print type(self)
        self.tree(attrs=True, verbose=True)
        return ''

    def listunits(self):
        '''
        List all variables and associated units.
        '''
        keys=self.keys(); keys.sort()
        length=0
        for key in keys:
            if len(key)>length: length=len(key)
        form="%%%is:%%s"%length
        for key in keys:
            if self[key].attrs.has_key('units'):
                print(form%(key, self[key].attrs['units']))
        
class IdlBin(PbData):
 
    '''
    An object class that reads/parses a binary output file from the SWMF
    and places it into a :class:spacepy.pybats.PbData object.

    Usage:
    >>>data = spacepy.pybats.IdlBin('binary_file.out')

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

    def __init__(self, filename, *args, **kwargs):
        super(IdlBin, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file'] = filename   # Save file name.
        self.read()   # Read binary file.

    def __repr__(self):
        return 'SWMF IDL-Binary file "%s"' % (self.attrs['file'])
    
    def read(self):
        '''
        This method reads an IDL-formatted BATS-R-US output file and places
        the data into the object.  The file read is self.filename which is
        set when the object is instantiation.
        '''
        import numpy as np
        import struct

        # Open, read, and parse the file into numpy arrays.
        # Note that Fortran writes integer buffers around records, so
        # we must parse those as well.
        infile = open(self.attrs['file'], 'rb')

        # On the first try, we may fail because of wrong-endianess.
        # If that is the case, swap that endian and try again.
        EndChar = '<' # Endian marker (default: little.)
        self.attrs['endian']='little'
        RecLenRaw = infile.read(4)

        RecLen = ( struct.unpack(EndChar+'l', RecLenRaw) )[0]
        if (RecLen > 10000) or (RecLen < 0):
            EndChar = '>'
            self.attrs['endian']='big'
            RecLen = ( struct.unpack(EndChar+'l', RecLenRaw) )[0]

        header = ( struct.unpack(EndChar+'%is'%RecLen,
                                 infile.read(RecLen)) )[0]    
        header.strip()
        units = header.split()

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        format = 'f'
        # parse header; detect double-precision file.
        if RecLen > 20: format = 'd'
        (self.attrs['iter'], self.attrs['time'],
         self.attrs['ndim'], self.attrs['nparam'], self.attrs['nvar']) = \
            struct.unpack(EndChar+'l%s3l' % format, infile.read(RecLen))
        # Get gridsize
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        self['grid']=dmarray(struct.unpack(EndChar+'%il' % 
                                           abs(self.attrs['ndim']), 
                                           infile.read(RecLen)))
        # Data from generalized (structured by irregular) grids can be 
        # detected by a negative ndim value.  Unstructured grids (e.g.
        # BATS, AMRVAC) are signified by negative ndim values AND
        # the grid size is always [x, 1(, 1)]
        # Here, we set the grid type attribute to either Regular, 
        # Generalized, or Unstructured.  Let's set that here.
        self['grid'].attrs['gtype'] = 'Regular'
        self['grid'].attrs['npoints']  = abs(self['grid'].prod())
        if self.attrs['ndim'] < 0: 
            if any(self['grid'][1:] > 1): 
                self['grid'].attrs['gtype'] = 'Generalized'
            else:
                self['grid'].attrs['gtype']   = 'Unstructured'
                self['grid'].attrs['npoints'] = self['grid'][0]
        self.attrs['ndim'] = abs(self.attrs['ndim'])

        # Quick ref vars:
        time=self.attrs['time']
        gtyp=self['grid'].attrs['gtype']
        npts=self['grid'].attrs['npoints']
        ndim=self['grid'].size
        nvar=self.attrs['nvar']
        npar=self.attrs['nparam']

        # Read parameters stored in file.
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        para = np.zeros(npar)
        para[:] = struct.unpack(EndChar+'%i%s' % (npar,format), 
                                      infile.read(RecLen))

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        names = ( struct.unpack(EndChar+'%is' % RecLen, 
                                infile.read(RecLen)) )[0].lower()
        names.strip()
        names = names.split()       

        # For some reason, there are often more units than variables
        # in these files.  It looks as if there are more grid units
        # than grid vectors (e.g. 'R R R' implies X, Y, and Z data
        # in file but only X and Y are present.)  Let's try to work
        # around this rather egregious error.
        nSkip=len(units)+npar-len(names)
        if nSkip<0: nSkip=0
        # Some debug crap:
        #print "nSkip=", nSkip
        #for n, u in zip(names, units[nSkip:]):
        #    print n, u

        # Save grid names (e.g. 'x' or 'r') and save associated params.
        self['grid'].attrs['dims']=names[0:ndim]
        for name, para in zip(names[(nvar+ndim):], para):
            self.attrs[name]=para

        # Create string representation of time.
        self.attrs['strtime']='%4.4ih%2.2im%06.3fs'%\
            (np.floor(time/3600.), np.floor(time%3600. / 60.0),
             time%60.0)

        # Get the grid points...
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        prod = [1] + self['grid'].cumprod().tolist()
        for i in range(0,ndim):
            tempgrid = np.array(struct.unpack(
                    EndChar+'%i%s' % (npts, format), 
                    infile.read(RecLen/ndim) ) )
            # Unstructred grids get loaded as vectors.
            if gtyp == 'Unstructured':
                self[names[i]] = dmarray(tempgrid)
            # Irregularly gridded items need multidimensional grid arrays:
            elif gtyp == 'Irregular':
                self[names[i]] = dmarray(
                    np.reshape(tempgrid, self['grid']))
            # Regularly gridded ones need vector grid arrays:
            elif gtyp == 'Regular':
                self[names[i]] = dmarray(np.zeros(self['grid'][i]))
                for j in range(int(self['grid'][i])):
                    self[names[i]][j] = tempgrid[j*int(prod[i])]
            else:
                raise ValueError('Unknown grid type: %s'%self.gridtype)
            # Add units to grid.
            self[names[i]].attrs['units']=units.pop(nSkip)

                    
        # Get the actual data and sort.
        for i in range(ndim,nvar+ndim):
            (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
            self[names[i]] = dmarray(\
                np.array(struct.unpack(EndChar+'%i%s' % (npts, format), 
                                       infile.read(RecLen))) )
            self[names[i]].attrs['units']=units.pop(nSkip)
            if gtyp != 'Unstructured':
                # Put data into multidimensional arrays.
                self[names[i]] = self[names[i]].reshape(self['grid'])

        # Unstructured data can be in any order, so let's sort it.
        if gtyp == 'Unstructured':
            gridtotal = np.zeros(npts)
            offset = 0.0  # The offset ensures no repeating vals while sorting.
            for key in self['grid'].attrs['dims']:
                gridtotal = gridtotal + offset + self[key]
                offset = offset + np.pi/2.0
                SortIndex = np.argsort(gridtotal)
            for key in self.keys():
                if key=='grid': continue
                self[key] = self[key][SortIndex]

        infile.close()

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

    def __init__(self, filename, starttime=(2000,1,1,0,0,0), *args, **kwargs):
        super(LogFile, self).__init__(*args, **kwargs)
        self.attrs['file'] = filename
        self.read(starttime)

    def read(self, starttime):
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
        names = (raw.pop(0)).split()
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
        if 'year'in loc: names.pop(names.index('year'))
        if 'mo'  in loc: names.pop(names.index('mo'))
        if 'dy'  in loc: names.pop(names.index('dy'))
        if 'hr'  in loc: names.pop(names.index('hr'))
        if 'mn'  in loc: names.pop(names.index('mn'))
        if 'sc'  in loc: names.pop(names.index('sc'))
        if 'msc' in loc: names.pop(names.index('msc'))
        if 't'   in loc: names.pop(names.index('t'))
        if 'it'  in loc: names.pop(names.index('it'))

        # Create containers for data:
        time=dmarray(np.zeros(npts, dtype=object))
        runtime=dmarray(np.zeros(npts), attrs={'units':'s'})
        self['iter']=dmarray(np.zeros(npts))
        for name in names:
            self[name] = dmarray(np.zeros(npts))

        for i in range(npts):
            vals = raw[i].split()
            # Set time:
            if 'year' in loc:
                # If "year" is listed, we have the full datetime.
                time[i]=(dt.datetime(
                        int(vals[loc['year']]), # Year
                        int(vals[loc['mo']  ]), # Month
                        int(vals[loc['dy']]), # Day
                        int(vals[loc['hr']]), # Hour
                        int(vals[loc['mn']]), # Minute
                        int(vals[loc['sc']]), # Second
                        int(vals[loc['msc']]) * 1000 #microsec
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
                self['iter'][i] = int(vals[loc['it']]) 
            else: self['iter'][i] = i
            # Collect data
            for j, name in enumerate(names):
                self[name][i] = float(vals[loc[name]])

            # Convert time and runtime to dmarrays.
            self['time']   =time
            self['runtime']=runtime

    def add_dst_quicklook(self, target=None, loc=111, showObs=True):
        '''
        Create a quick-look plot of Dst (if variable present in file) 
        and compare against observations.

        Like all *add_\* * methods in Pybats, the *target* kwarg determines
        where to place the plot.
        If kwarg *target* is **None** (default), a new figure is 
        generated from scratch.  If *target* is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot location
        *loc* (defaults to 111).  If target is a matplotlib Axes object, 
        the plot is placed into that axis at subplot location *loc*.

        Observed Dst is automatically fetched from the Kyoto World Data Center
        via the :mod:`spacepy.pybats.kyoto` module.  The associated 
        :class:`spacepy.pybats.kyoto.KyotoDst` object, which holds the observed
        Dst, is stored as *self.obs_dst* for future use.

        The figure and axes objects are returned to the user.
        '''
        
        import matplotlib.pyplot as plt
        
        if 'dst' not in self.keys():
            return None, None

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(loc)

        ax.plot(self['time'], self['dst'], 
                label='BATS-R-US $D_{ST}$ (Biot-Savart)')
        ax.hlines(0.0, self['time'][0], self['time'][-1], 
                  'k', ':', label='_nolegend_')
        apply_smart_timeticks(ax, self['time'])
        ax.set_ylabel('Dst ($nT$)')
        ax.set_xlabel('Time from '+ self['time'][0].isoformat()+' UTC')

        if(showObs):
            try:
                import spacepy.pybats.kyoto as kt
            except ImportError:
                return fig, ax
        
            try:
                stime = self['time'][0]; etime = self['time'][-1]
                if not hasattr(self, 'obs_dst'):
                    self.obs_dst = kt.fetch('dst', stime, etime)

            except BaseException as args:
                print('WARNING! Failed to fetch Kyoto Dst: '+ args)
            else:
                ax.plot(self.obs_dst['time'], self.obs_dst['dst'], 
                        'k--', label='Obs. Dst')
                ax.legend(loc='best')
                apply_smart_timeticks(ax, self['time'])
        else:
            ax.legend(loc='best')
            

        return fig, ax

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
                outfile='imfinput.dat'

        out = open(outfile, 'w')
        
        # Write header.
        for c in self.attrs['comments']:
            out.write(c)
        out.write(2*'#\n'+'#'+50*'-'+'\n')
        
        # Write variables.
        for k in self.keys():
            out.write('#>\n')
            out.write('#%s: %s\n' % ('Element', k))
            for a in self[k].attrs.keys():
            #for a in ['Table','Description','Measure units','Origin']:
                out.write('#%s: %s\n' % (a, self[k].attrs[a]))
            #out.write(2*'#\n')
            #for a in ['Sampling', 'Missing value']:
            #    out.write('#%s: %s\n' % (a, self[k].attrs[a]))
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

    Like most :module:`pybats` objects, you may interact with :class:`ImfInput`
    objects as if they were specialized dictionaries.  Access data like so:
    
    >>> obj.keys()
    ['bx', 'by', 'bz', 'vx', 'vy', 'vz', 'dens', 'temp']
    >>> density=obj['dens']

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
        self.attrs['var']= ['bx', 'by', 'bz', 'vx', 'vy', 'vz', 'dens', 'temp']
        self.attrs['std_var']=True
        self.attrs['coor']='GSM'
        self.attrs['satxyz']=[None, None, None]
        self.attrs['zerobx']=False
        self.attrs['reread']=False
        self.attrs['delay']=None
        self.attrs['plane']=[None, None]
        self.attrs['header']=[]
        self['time']=dmarray(zeros(npoints, dtype=object))
        units   = ['nT', 'nT', 'nT', 'km/s', 'km/s', 'km/s', 'cm^-3', 'K']
        for i, key in enumerate(self.attrs['var']):
            self[key]=dmarray(zeros(npoints), attrs={'units':units[i]})
            
        # Set Filename.
        if filename:
            self.attrs['file'] = filename
        else:
            self.attrs['file'] = None

        # Load/create data vectors.
        if filename and load:  # Load contents from existing file.
            self.read(filename)
            self.calc_pram()

    def calc_pram(self):
        '''
        Calculate ram pressure from SW conditions.  Output units in nPa.
        If object was instantiated via an existing imf file, this value
        is calculated automatically.
        '''
        self['pram']=dmarray(self['vx']**2.*self['dens']*1.67621E-6, 
                             {'units':'nPa'})

    def varcheck(self):
        '''
        Ensure that the variable list, which gives the order of the
        variables as they are written to file, correctly corresponds
        to the variables stored in the ImfInput object.

        Returns True on success, returns warnings and False on fail.
        '''
        # Convenience:
        var=self.attrs['var']
        key=self.keys()
        key.remove('time')

        # Number of variables check:
        if len(var) != len(key):
            print('Number of listed variables is incorrect:')
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
        from numpy import zeros
        import datetime as dt

        # Slurp lines into memory.
        f = open(infile, 'r')
        lines = f.readlines()
        f.close()
        
        # Read header.  All non-blank lines before first Param are header.
        self.attrs['header']=[]
        while 1:
            if (lines[0].strip() != '') and lines[0][0] != '#':
                self.attrs['header'].append(lines.pop(0)) 
            else:
                break

        # Parse all Params.
        while len(lines)>0:
            # Grab line, continue if it's not a Param.
            param=lines.pop(0).strip()
            if param=='': continue
            if param[0] != '#': continue
            # For all possible Params, set object attributes/info.
            if param == '#COOR':
                self.attrs['coor']=lines.pop(0)[0:3]
            elif param == '#REREAD':
                self.attrs['reread']=True
            elif param == '#ZEROBX':
                setting=lines.pop(0)[0]
                if setting=='T':
                    self.attrs['zerobx']=True
                else:
                    self.attrs['zerobx']=False
            elif param == '#VAR':
                self.attrs['var']=lines.pop(0).split()
                self.attrs['std_var']=False
            elif param == '#PLANE':
                xp = float(lines.pop(0).split()[0])
                yp = float(lines.pop(0).split()[0])
                self.attrs['plane']=[xp, yp]
            elif param == '#POSITION':
                yp = float(lines.pop(0).split()[0])
                zp = float(lines.pop(0).split()[0])
                self.attrs['satxyz'][1:]=(yp, zp)
            elif param == '#SATELLITEXYZ':
                Xp = float(lines.pop(0).split()[0])
                yp = float(lines.pop(0).split()[0])
                zp = float(lines.pop(0).split()[0])
                self.attrs['satxyz']=[xp, yp, zp]
            elif param == '#TIMEDELAY':
                self.attrs['delay']=float(lines.pop(0).split()[0])
            elif param == '#START':
                break
            else:
                raise Exception('Unknown file parameter: ' + param)

        # Create containers for data.
        npoints = len(lines)
        self['time']=dmarray(zeros(npoints, dtype=object))
        for key in self.attrs['var']:
            self[key]=dmarray(zeros(npoints))

        # Parse data.
        for i, line in enumerate(lines):
            parts = line.split()
            self['time'][i]=(dt.datetime(
                    int(parts[0]), #year
                    int(parts[1]), #month
                    int(parts[2]), #day
                    int(parts[3]), #hour
                    int(parts[4]), #min
                    int(parts[5]), #sec
                    int(parts[6]) * 1000 #micro seconds
                    )) 
            for j, name in enumerate(self.attrs['var']):
                self[name][i] = float(parts[7+j])

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

        out = open(outfile, 'w')
        
        # Convenience variable:
        var=self.attrs['var']

        # Write the header:
        out.write('File created on %s\n' % (dt.datetime.now().isoformat()))
        for head in self.attrs['header']:
            out.write(head)

        # Handle Params:
        if self.attrs['zerobx']:
            out.write('#ZEROBX\nT\n\n')
        if self.attrs['reread']:
            out.write('#REREAD')
        if not self.attrs['std_var']:
            out.write('#VAR\n%s\n\n'%(' '.join(var)))
        if self.attrs['satxyz'].count(None)<3:
            if (self.attrs['satxyz'][0]==None) and (
                None not in self.attrs['satxyz'][1]):
                out.write('#POSITION\n%-6.2f\n%-6.2f\n\n'%
                          (self.attrs['satxyz'][2],self.attrs['satxyz'][3]))
            elif None not in self.attrs['satxyz']:
                out.write('#SATELLITEXYZ\n%s\n' % 
                          (''.join("%-6.2f\n"%n for n in self.attrs['satxyz'])))
        if self.attrs['delay']:
            out.write('#DELAY\n%-9.2f\n\n' % (self.attrs['delay']))
        if None not in self.attrs['plane']:
            out.write('#PLANE\n%s\n'%
                      (''.join('%-6.2f\n'%n for n in self.attrs['plane'])))

        # Write the data:
        out.write('\n#START\n')
        for i in range(len(self['time'])):
            out.write('%04d %02d %02d %02d %02d %02d %03d' % 
                          (self['time'][i].year, 
                           self['time'][i].month,
                           self['time'][i].day,
                           self['time'][i].hour,
                           self['time'][i].minute,
                           self['time'][i].second,
                           self['time'][i].microsecond/1000.0 ) )
            out.write(' %s\n' % ' '.join('%10.2f'%self[key][i] for key in var))
            #for key in var:
            #    out.write(' %9.2f' % (self[key][i]))
            #out.write('\n')
                #out.write('  %11.4E\n' % (self['temp'][i]))
        out.close()

    def quicklook(self, timerange=None):
        '''
        Generate a quick-look plot of solar wind conditions driving the
        SWMF.  Default values show IMF, number density, and Earthward velocity.
        Returns a figure object containing the plots.
        '''

        import matplotlib.pyplot as plt

        if not timerange:
            timerange = [self['time'][0], self['time'][-1]]

        def adjust_plots(ax, ylab, xlab=False, Zero=True):
            ax.grid(True)
            apply_smart_timeticks(ax,timerange)
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
        a4.plot(self['time'], self['dens'], lw=1.25, c='red')
        adjust_plots(a4, 'Density ($cm^{-3}$)', Zero=False)

        a5 = fig.add_subplot(515)
        a5.plot(self['time'], -1.0*self['vx'], lw=1.25, c='green')
        adjust_plots(a5, '$V_{X}$ ($km/s$)', Zero=False, xlab=True)

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
            outfile.write(line)
        outfile.write('\n')
        outfile.write('#COOR\n%s\n\n' % self.attrs['coor'])
        outfile.write('#START\n')

        # Write the rest of the orbit.
        npts=len(self['time'])
        for i in range(npts):
            #Time:
            outfile.write('%04d %02d %02d %02d %02d %02d %03d ' % 
                          (self['time'][i].year, 
                           self['time'][i].month,
                           self['time'][i].day,
                           self['time'][i].hour,
                           self['time'][i].minute,
                           self['time'][i].second,
                           self['time'][i].microsecond/1000.0 ) )
            #Position:
            outfile.write('%13.7E %13.7E %13.7E\n' % 
                          (self['xyz'][0,i], 
                           self['xyz'][1,i],
                           self['xyz'][2,i]) )

        outfile.close()

