# -*- coding: utf-8 -*-
'''
PyBats!  An open source Python-based interface for reading, manipulating,
and visualizing BATS-R-US and SWMF output.

Copyright ©2010 Los Alamos National Security, LLC.
'''

def load(pbfile):
    '''
    Load and return a pickled pybats object.

    Example:

    >>> pyobj = pybats.load('pybats_file.pb')

    '''
    try:
        from cPickle import load
    except:
        from pickle import load
    infile = open(pbfile, 'r')
    pbobj = load(infile)
    infile.close
    return pbobj

def save(pbobj, filename=None):
    '''
    Use cPickle to save a pybats object as a binary file.
    If kwarg filename is set, save to that file name.  
    Otherwise, use pbobj.filename to create a path with a new extension, '.pb'.
    Example:
    
    >>> pbobj = pybats.IdlBin('examplefile.out')
    >>> pybats.save(pbobj, 'newfilename.pb')
    '''
    try:
        from cPickle import dump
    except:
        from pickle import dump

    if filename == None:
        i = pbobj.filename.rfind('.')
        filename = pbobj.filename[0:i]+'.pb'

    outfile = open(filename, 'w')
    dump(pbobj, outfile)
    outfile.close()

class IdlBin(object):
 
    '''
    An object class to hold information from an IDL-format SWMF binary
    output file.  Upon instantiation, the file name given to __init__() will
    be opened and read into the class instance.

    This class serves as a parent class to SWMF component-specific derivative
    classes that do more preprocessing of the data before returning the 
    object.  Hence, using this class to read binary files is typically not
    the most efficient way to proceed.  Look for a PyBats sub module that suits
    your specific needs, or use this base object to write your own!

    A note on byte-swapping: PyBats assumes little endian byte ordering because
    this is what most machines use.  However, there is an autodetect feature
    such that, if PyBats doesn't make sense of the first read (a record length
    entry, or RecLen), it will proceed using big endian ordering.  If this
    doesn't work, the error will manifest itself through the "struct" package
    as an "unpack requires a string of argument length 'X'".

    Class Attibutes:
    filename: File opened and read into object.
    iter: Simulation iteration when file was created.
    time: Time elapsed, in seconds, from beginning of simulation.
    ndim, nvar, npara: number of dimensions, variables, and parameters.
    namevar: list of string names of variables saved in file.
    namepara: list of string names of parameters saved in file.
    gridtype: Either Regular, Generalized, or Unstructured.
    gridsize: An array listing the number of elements in each dimension.
    grid: A dictionary of dimension name and corresponding gridpoints.
    data: A dictionary of variable names and corresponding data arrays.
    para: A dictionary of parameter names and corresponding values.
    '''

    def __init__(self, file):
        self.filename = file
        self.read()   # Read binary file.

    def __repr__(self):
        return 'SWMF IDL-Binary file %s' % (self.filename)

    def __str__(self):
        return 'SWMF IDL-Binary file "%s"' % (self.filename)

    def read(self):
        '''
        This method reads an IDL-formatted BATS-R-US output file and places
        the data into the object.  The file read is self.filename which is
        set when the object is instantiation.
        '''
        from math import pi
        import numpy as np
        import struct

        # Open, read, and parse the file into numpy arrays.
        # Note that Fortran writes integer buffers around records, so
        # we must parse those as well.
        infile = open(self.filename, 'rb')

        # On the first try, we may fail because of wrong-endianess.
        # If that is the case, swap that endian and try again.
        EndChar = '<' # Endian marker (default: little.)
        RecLenRaw = infile.read(4)

        RecLen = ( struct.unpack(EndChar+'l', RecLenRaw) )[0]
        if (RecLen > 10000) or (RecLen < 0):
            print('Confusing data (RecLen=%i), Switching Endian...' % RecLen)
            EndChar = '>'
            RecLen = ( struct.unpack(EndChar+'l', RecLenRaw) )[0]

        header = ( struct.unpack(EndChar+'%is'%RecLen,
                                 infile.read(RecLen)) )[0]    
        header.strip()
        self.units = header.split()

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        format = 'f'
        # parse header; detect double-precision file.
        if RecLen > 20: format = 'd'
        (self.iter, self.time, self.ndim, self.npara, self.nvar) = \
            struct.unpack(EndChar+'l%s3l' % format, infile.read(RecLen))
        # Get gridsize
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        self.gridsize = np.zeros(abs(self.ndim))
        self.gridsize[:] = struct.unpack(EndChar+'%il' % abs(self.ndim), 
                                         infile.read(RecLen))
        # Data from generalized (structured by irregular) grids can be 
        # detected by a negative ndim value.  Unstructured grids (e.g.
        # BATS, AMRVAC) are signified by negative ndim values AND
        # the grid size is always [x, 1(, 1)]
        # Here, we set the grid type attribute to either Regular, 
        # Generalized, or Unstructured.  Let's set that here.
        self.gridtype = 'Regular'
        self.npoints  = abs(self.gridsize.prod())
        if self.ndim < 0: 
            if any(self.gridsize[1:] > 1): 
                self.gridtype = 'Generalized'
            else:
                self.gridtype = 'Unstructured'
                self.npoints = self.gridsize[0]
        self.ndim = abs(self.ndim)

        # Read parameters stored in file.
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        para = np.zeros(self.npara)
        para[:] = struct.unpack(EndChar+'%i%s' % (self.npara,format), 
                                      infile.read(RecLen))

        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        names = ( struct.unpack(EndChar+'%is' % RecLen, 
                                infile.read(RecLen)) )[0].lower()
        names.strip()
        names = names.split()
        
        gridnames = names[0:self.ndim]
        self.namevar  = names[self.ndim:self.nvar+self.ndim]
        self.namepara = names[self.nvar+self.ndim:]

        # Initialize dictionaries to hold data.
        self.grid = {}
        self.data = {}
        self.para = {}
        # Get the grid points...
        (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
        prod = [1] + self.gridsize.cumprod().tolist()
        for i in range(0,self.ndim):
            tempgrid = np.array(struct.unpack(
                    EndChar+'%i%s' % (self.npoints, format), 
                    infile.read(RecLen/self.ndim) ) )
            # Unstructred grids get loaded as vectors.
            if self.gridtype == 'Unstructured':
                self.grid[gridnames[i]] = tempgrid
            # Irregularly gridded items need multidimensional grid arrays:
            elif self.gridtype == 'Irregular':
                self.grid[gridnames[i]] = \
                    np.reshape(self.grid[gridnames[i]], self.gridsize)
            # Regularly gridded ones need vector grid arrays:
            elif self.gridtype == 'Regular':
                self.grid[gridnames[i]] = np.zeros(self.gridsize[i])
                for j in range(int(self.gridsize[i])):
                    self.grid[gridnames[i]][j] = tempgrid[j*int(prod[i])]
            else:
                raise ValueError('Unknown grid type: %s'%self.gridtype)


                    
        # Get the actual data and sort.
        for i in range(0,self.nvar):
            (OldLen, RecLen) = struct.unpack(EndChar+'2l', infile.read(8))
            self.data[self.namevar[i]] = \
                np.array(struct.unpack(EndChar+'%i%s' % (self.npoints, format), 
                                       infile.read(RecLen) ) )
            if self.gridtype != 'Unstructured':
                # Put data into multidimensional arrays.
                self.data[self.namevar[i]] = \
                    np.reshape(self.data[self.namevar[i]], self.gridsize)

        # Sort the parameters into a dictionary, too.
        for i in range(0,self.npara):
            self.para[self.namepara[i]]= para[i]

        # Unstructured data can be in any order, so let's sort it!
        if self.gridtype == 'Unstructured':
            gridtotal = np.zeros(self.npoints)
            offset = 0.0  # The offset ensures no repeating vals while sorting.
            for keys in self.grid:
                gridtotal = gridtotal + offset + self.grid[keys]
                offset = offset + pi/2.0
                SortIndex = np.argsort(gridtotal)
            for keys in self.grid:
                self.grid[keys] = self.grid[keys][SortIndex]
            for keys in self.data:
                self.data[keys] = self.data[keys][SortIndex]

        infile.close()

class LogFile(object):
    ''' An object to read and handle SWMF-type logfiles.

    logfile objects read and hold all information in an SWMF
    ascii logfile.  The file is read upon instantiation.

    Time is handled by Python's datetime package.  Given that time
    may or may not be given in the logfile, there are three options
    for how time is returned:
        1. if the full date and time are listed in the file, self.time
        is a list of datetime objects corresponding to the entries.  The
        starttime kwarg is ignored.
        2. If only the runtime (seconds from start of simulation) is given,
        self.time is a list of datetime objects that either start from
        the given starttime kwarg or starting from 1/1/1 00:00UT.
        3. If neither of the above are given, the time is assumed to 
        advance one second per line starting from either the starttime kwarg
        or from 1/1/1 00:00UT + the first iteration (if given in file.)
        As you can imagine, this is sketchy at best.
    
    Example usage:

    >>> import pybats
    >>> file1 = pybats.logfile('satfile_n000000.dat')
    >>> file1.namevar
    ['it', 't', 'rho', 'dst','bx', 'by', 'bz']
    
    This file does not have the time clearly set, only time from the 
    beginning of the run.  Let's set the start time and plot Dst:

    >>> import pybats
    >>> import pylab as plt
    >>> import datetime as dt
    >>> time1 = dt.datetime(2009,11,30,9,0)
    >>> file1 = pybats.logfile('satfile_n000000.dat', starttime=time1)
    >>> plt.plot(file1.time, file1.data['dst'])
    
    Just like other PyBats objects, data is referenced by its name from
    the file.  

    '''

    import datetime as dt

    def __init__(self, file, starttime=dt.datetime(1,1,1,0,0,0)):
        self.filename = file
        self.read(starttime)

    def __repr__(self):
        print('Yeah, working on that.')

    def read(self, starttime):
        '''
        Load the ascii logfile located at self.filename.
        This method is automatically called upon instantiation.
        '''
        import numpy as np
        import datetime as dt

        # Slurp in entire file.
        infile = open(self.filename, 'r')
        raw = infile.readlines()
        infile.close()

        # Parse the header.
        self.description = raw.pop(0)
        self.namevar = (raw.pop(0)).split()
        loc={}
        for i, name in enumerate(self.namevar):
            loc[name] = i

        # Use the header info to set up arrays 
        # for time, iteration, and the data.
        self.npoints = len(raw) 

        # Pop time/date/iteration names off of Namevar.
        if 'year'in loc: self.namevar.pop(self.namevar.index('year'))
        if 'mo'  in loc: self.namevar.pop(self.namevar.index('mo'))
        if 'dy'  in loc: self.namevar.pop(self.namevar.index('dy'))
        if 'hr'  in loc: self.namevar.pop(self.namevar.index('hr'))
        if 'mn'  in loc: self.namevar.pop(self.namevar.index('mn'))
        if 'sc'  in loc: self.namevar.pop(self.namevar.index('sc'))
        if 'msc' in loc: self.namevar.pop(self.namevar.index('msc'))
        if 't'   in loc: self.namevar.pop(self.namevar.index('t'))
        if 'it'  in loc: self.namevar.pop(self.namevar.index('it'))

        # Create containers for data:
        self.data = {}
        self.time = []
        self.iter = np.zeros(self.npoints)
        for name in self.namevar:
            self.data[name] = np.zeros(self.npoints)

        for i, line in enumerate(raw):
            vals = line.split()
            # Set time:
            if 'year' in loc:
                self.time.append(dt.datetime(
                        int(vals[loc['year']]), # Year
                        int(vals[loc['mo']  ]), # Month
                        int(vals[loc['dy']]), # Day
                        int(vals[loc['hr']]), # Hour
                        int(vals[loc['mn']]), # Minute
                        int(vals[loc['sc']]), # Second
                        float(vals[loc['msc']]) * 1000 #microsec
                        ))
            elif 't' in loc:
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
                self.time.append(newtime)
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
                self.time.append(newtime)
            else:
                self.time.append(starttime + dt.timedelta(float(i)))
            # Set iteration:
            if 'it' in loc:
                self.iter[i] = int(vals[loc['it']]) 
            else: self.iter[i] = i
            # Collect data
            for j, name in enumerate(self.namevar):
                self.data[name][i] = float(vals[loc[name]])


class ImfInput(object):
    '''
    An object to hold the contents for an SWMF IMF input file.

    Currently, this object simply reads in a file upon instantiation.
    It is a candidate for strong revision in the future, so be careful
    using this class in your programs.  Should be a subclass of logfile.
    '''

    def __init__(self, infile):
        # It may be a good idea to make file name a kwarg, 
        # therefore we could create this object w/o reading
        # a file.
        self.namevar = ['bx', 'by', 'bz', 'vx', 'vy', 'vz', 'dens', 'temp']
        self.units   = ['nT', 'nT', 'nT', 'km/s', 'km/s', 'km/s', 'cm^-3', 'K']
        self.filename = infile
        self.read(infile)

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
        
        # Read header.
        self.header=[]
        while 1:
            self.header.append(lines.pop(0))
            if self.header[-1].strip() == '#START':break
            
        # Create containers for data.
        npts = len(lines)
        self.time = []
        self.data = {}
        for name in self.namevar:
            self.data[name] = zeros(npts)
         
        # Parse data.
        for i, line in enumerate(lines):
            parts = line.split()
            self.time.append(dt.datetime(
                    int(parts[0]), #year
                    int(parts[1]), #month
                    int(parts[2]), #day
                    int(parts[3]), #hour
                    int(parts[4]), #min
                    int(parts[5]), #sec
                    int(parts[6]) * 1000 #micro seconds
                    ) )
            for j, name in enumerate(self.namevar):
                self.data[name][i] = float(parts[7+j])

    def quicklook(self, timerange=None):
        '''
        Generate a quick-look plot of solar wind conditions driving the
        SWMF.  Default values show IMF, number density, and Earthward velocity.
        Returns a figure object containing the plots.
        '''

        from .rampy import apply_smart_timeticks
        import matplotlib.pyplot as plt
        
        if not timerange:
            timerange = [self.time[0], self.time[-1]]

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
                              timerange[-1].isoformat())
            else:
                ax.xaxis.set_ticklabels([])

        fig = plt.figure(figsize=(8,10))
        fig.subplots_adjust(hspace=0.025)
        
        a1 = fig.add_subplot(511)
        a1.plot(self.time, self.data['bx'], 'b-')
        adjust_plots(a1, 'IMF $B_{X}$ ($nT$)')
        a1.set_title('Solar Wind Drivers')

        a2 = fig.add_subplot(512)
        a2.plot(self.time, self.data['by'], 'b-')
        adjust_plots(a2, 'IMF $B_{Y}$ ($nT$)')

        a3 = fig.add_subplot(513)
        a3.plot(self.time, self.data['bz'], 'b-')
        adjust_plots(a3, 'IMF $B_{Z}$ ($nT$)')

        a4 = fig.add_subplot(514)
        a4.plot(self.time, self.data['dens'], 'r-')
        adjust_plots(a4, 'Density ($cm^{-3}$)', Zero=False)

        a5 = fig.add_subplot(515)
        a5.plot(self.time, -1.0*self.data['vx'], 'g-')
        adjust_plots(a5, '$V_{X}$ ($km/s$)', Zero=False, xlab=True)

        return fig

class SatOrbit(object):
    '''
    An class to load, read, write, and handle BATS-R-US satellite orbit
    input files.  These files are used to fly virtual satellites through
    the MHD domain.  Note that the output files should be handled by 
    the L{logfile object<pybats.logfile>} and not this satorbit object.
    '''

    def __init__(self, filename=None):
        '''
        Create a satorbit object.  If I{file} is not given, buid an
        empty satorbit object is instantiated for the user to fill
        with values of their own creation.
        '''
        import numpy as np

        self.filename = filename
        self.header=[]
        self.npoints=0
        self.coord_system = None
        self.time = []

        if filename:
            try:
                self.read()
            except IOError:
                self.position = np.zeros( (3,1) )
                raise

        else:
            # fill with empty stuff.
            self.position = np.zeros( (3, 1) )


    def read(self):
        '''
        Read and parse a satellite orbit input file into the satorbit object.
        '''
        import numpy as np
        import datetime as dt
        # Try to open the file.  If we fail, pass the exception
        # to the caller.
        try:
            infile = open(self.filename, 'r')
        except IOError as reason:
            raise

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
                    self.coord_system = raw.pop(0).strip()
                else:
                    self.header.append(line)

        # Read and store all values.
        self.npoints = len(raw)
        self.position = np.zeros( (3,self.npoints) )
        for i, line in enumerate(raw):
            parts = line.split()
            self.time.append(dt.datetime(
                    int(parts[0]), #year
                    int(parts[1]), #month
                    int(parts[2]), #day
                    int(parts[3]), #hour
                    int(parts[4]), #min
                    int(parts[5]), #sec
                    int(parts[6]) * 1000 #micro seconds
                    ) )
            self.position[:,i] = parts[7:]

    def write(self):
        '''Write a L{satorbit object<pybats.satorbit>} to file using the 
        correct SWMF-input format.  The file will be written to the 
        location specified by self.filename.  An error is thrown if this
        attribute is not set.
        '''

        try:
            outfile = open(self.filename, 'w')
        except:
            print('Could not open self.filename!')
            raise
        
        # Start by writing header, coordinate system, and then #START.
        for line in self.header:
            outfile.write(line)
        outfile.write('\n')
        outfile.write('#COOR\n%s\n\n' % self.coord_system)
        outfile.write('#START\n')

        # Write the rest of the orbit.
        for i in range(self.npoints):
            #Time:
            outfile.write('%04d %02d %02d %02d %02d %02d %03d ' % 
                          (self.time[i].year, 
                           self.time[i].month,
                           self.time[i].day,
                           self.time[i].hour,
                           self.time[i].minute,
                           self.time[i].second,
                           self.time[i].microsecond/1000.0 ) )
            #Position:
            outfile.write('%13.7E %13.7E %13.7E\n' % 
                          (self.position[0,i], 
                           self.position[1,i],
                           self.position[2,i]) )

        outfile.close()
