'''
A PyBats module for handling input, output, and visualization of 
binary SWMF output files taylored to BATS-R-US-type data.
'''

import numpy as np
import spacepy.pybats

class Stream(object):
    '''
    A class for streamlines.  Contains all of the information about
    the streamline, including extracted variables.

    Upon instantiation, the object will trace through the vector
    field determined by the "[x/y]field" values and the Bats object
    "bats".

    Line colors are set using the style keyword.  For more information,
    see the set_style method.

    The integration method is set by the method kwarg.  The default is
    Runge-Kutta 4 (rk4), but others are available.  RK4 gives a good
    blend of speed and accuracy; see the test functions in pybats.trace2d
    for more info.
    '''
    
    def __init__(self, bats, xstart, ystart, xfield, yfield,
                 style = 'mag', type='streamline', method='rk4',
                 extract=False):
        # Key values:
        self.xstart = xstart #X and Y starting
        self.ystart = ystart #points in the field.
        self.xvar = xfield #Strings that list the variables
        self.yvar = yfield #that will be used for tracing.
        # Descriptors:
        self.type   = type
        self.open   = True
        self.method = method

        # Do tracing:
        self.trace(bats, extract)

        # Set style
        self.set_style(style)


    def __repr__(self):
        pass

    def __str__(self):
        pass

    def set_style(self, style):
        '''
        Set the line style either using a simple matplotlib-type style
        string or using a preset style type.  Current types include:
        
        'mag' : treat line as a magnetic field line.  Closed lines are
                white, other lines are black.
        '''
        
        if style == 'mag':
            if self.open:
                self.style = 'k-'
            else:
                self.style = 'w-'
        else:
            self.style = style

    def trace(self, bats, extract=False):
        '''
        Trace through the vector field.
        '''
        from numpy import array, sqrt
        if self.method == 'euler':
            from spacepy.pybats.trace2d import trace2d_eul as trc
        elif self.method == 'rk4':
            from spacepy.pybats.trace2d import trace2d_rk4 as trc
        
        # Get name of dimensions in order.
        grid = list(bats.grid.keys())

        # Trace forward
        x1, y1 = trc(bats.data[self.xvar], bats.data[self.yvar], 
                     self.xstart, self.ystart,
                     bats.grid[grid[0]], bats.grid[grid[1]])
        # Trace backward
        x2, y2 = trc(bats.data[self.xvar], bats.data[self.yvar], 
                     self.xstart, self.ystart,
                     bats.grid[grid[0]], bats.grid[grid[1]], ds=-0.1)
        # Join lines together such that point 0 is beginning of line
        # and point -1 is the end (when moving parallel to line.)
        self.x = array(x2[::-1].tolist() + x1[1:].tolist())
        self.y = array(y2[::-1].tolist() + y1[1:].tolist())

        # Check if line is closed to body.
        if 'rbody' in bats.para:
            r1 = sqrt(self.x[0]**2.0  + self.y[0]**2.0)
            r2 = sqrt(self.x[-1]**2.0 + self.y[-1]**2.0)
            if (r1 < bats.para['rbody']) and (r2 < bats.para['rbody']):
                self.open = False

    def plot(self, ax):
        '''
        Add streamline to axes object "ax". 
        '''
        ax.plot(self.x, self.y, self.style)

class Bats2d(spacepy.pybats.IdlBin):
    '''
    An object class of parent pybats.idlbin taylored to BATS-R-US output.
    This function requires the Matplotlib griddata function.
    
    '''

    ####################
    # CALCULATIONS
    ####################
    def calc_temp(self, units='eV'):
        '''
        If pressure and number density are available, calculate temperature
        for each species using P = nkT.  Current version works for single
        species only.

        Use the units kwarg to set output units.  Current choices are
        KeV, eV, and K.  Default is eV.
        '''
        units = units.lower()

        # Check for necessary variables:
        try:
            self.namevar.index('rho')
        except ValueError:
            raise AttributeError("Number density not found.")

        try:
            self.namevar.index('p')
        except ValueError:
            raise AttributeError("Pressure not found.")
        
        # Calculate value then adjust units.
        self.data['temp'] = self.data['p'] / self.data['rho']
        if units == 'ev':
            self.data['temp'] = self.data['temp'] * 6250.0
        elif units == 'kev':
            self.data['temp'] = self.data['temp'] * 6.250
        elif units == 'k':
            self.data['temp'] = self.data['temp'] * 72432275.82

    def calc_all(self):
        '''
        Perform all variable calculations (e.g. calculations that
        begin with 'calc_').  Any exceptions raised by functions that
        could not be peformed (typicaly from missing variables) are
        discarded.
        '''
        for command in dir(self):
            if (command[0:5] == 'calc_') and (command != 'calc_all'):
                try:
                    eval('self.'+command+'()')
                except AttributeError as Error:
                    print('WARNING: Did not perform %s: %s' % (command, Error))

    def regrid(self, cellsize=1.0, dim1range=-1, dim2range=-1, debug=False):
        '''
        Re-bin data to regular grid of spacing cellsize.  Action is 
        performed on all data entries in the bats2d object.

        '''
        from matplotlib.mlab import griddata

        if self.gridtype == 'Regular': return
        
        # Order our dimensions alphabetically.
        dims = []
        for key in sorted(self.grid.keys()):
            dims.append(key)
        if debug: print("Ordered dimensions: ", dims)

        # Check to see if dimranges are 2-element lists.
        # If not, either set defaults or raise exceptions.
        if dim1range == -1:
            dim1range = [min(self.grid[dims[0]]), max(self.grid[dims[0]])]
        else:
            if isinstance(dim1range, ( type(()), type([]) ) ):
                if len(dim1range) != 2:
                    raise ValueError('dim1range must have two elements!')
            else: raise TypeError('dim1range must be a tuple or list!')
        if dim2range == -1:
            dim2range = [min(self.grid[dims[1]]), max(self.grid[dims[1]])]
        else:
            if isinstance(dim2range, ( type(()), type([]) ) ):
                if len(dim2range) != 2:
                    raise ValueError('dim2range must have two elements!')
            else: raise TypeError('dim2range must be a tuple or list!')

        if debug:
            print('%s range = %f, %f' % (dims[0], dim1range[0], dim1range[1]))
            print('%s range = %f, %f' % (dims[1], dim2range[0], dim2range[1]))

        # Now, Regrid.
        grid1 = np.arange(dim1range[0], dim1range[1]+cellsize, cellsize)
        grid2 = np.arange(dim2range[0], dim2range[1]+cellsize, cellsize)

        for key in list(self.data.keys()):
            self.data[key] = griddata(self.grid[dims[0]], self.grid[dims[1]],
                                      self.data[key], grid1, grid2)

        # Change grid, gridtype, gridsize, and npoints to match new layout.
        self.gridtype = 'Regular'
        self.grid[dims[0]] = grid1
        self.grid[dims[1]] = grid2
        self.gridsize = [len(grid1), len(grid2)]
        self.npoints = len(grid1) * len(grid2)


    ######################
    # TRACING TOOLS
    ######################
    def get_stream(self, x, y, xvar, yvar, method='rk4', style='mag'):
        '''
        Trace a 2D streamline through the domain, returning a Stream
        object to the caller.

        x and y set the starting point for the tracing.

        xvar and yvar are string keys to self.data that define the
        vector field through which this function traces.  

        The method kwarg sets the numerical method to use for the
        tracing.  Default is Runge-Kutta 4 (rk4).
        '''

        stream = Stream(self, x, y, xvar, yvar, style=style)

        return stream

    ######################
    # VISUALIZATION TOOLS
    ######################
    def planet_add_closed(self, ax, style='mag', DoImf=False, DoOpen=False):
        '''
        Create an array of field lines closed to the central body in the
        domain.  Add these lines to axes "ax".  Default line style is
        white solid.

        Note that this should currently only be used for GSM y=0 cuts
        of the magnetosphere.

        Method:
        First, the title angle is approximated by tracing a dipole-like
        field line and finding the point of maximum radial distance on
        that line.  This is used as the magnetic equator.  From this 
        equator, many lines are traced starting at the central body
        radius.  More lines are grouped together at higher magnetic
        latitude to give good coverage at larger L-shells.  Once an
        open field line is detected, no more lines are traced.  This
        is done on the day and night side of the central body.

        Because this method rarely captures the position of the night
        side x-line, more field lines are traced by marching radially
        outwards away from the furthest point from the last traced and
        closed field line.  This is repeated until open lines are found.
        '''
        
        from numpy import (arctan, cos, sin, where, pi, log, 
                           arange, sqrt, linspace)

        # Approximate the central body.
        stream = self.get_stream(3.0, 0, 'bx', 'bz')
        r = stream.x**2 + stream.y**2
        loc, = where(r==r.max())
        tilt = arctan(stream.y[loc[0]]/stream.x[loc[0]])
        
        # Initial values:
        daymax   = tilt + pi/2.0
        nightmax = tilt + 3.0*pi/2.0

        # Day side:
        n = arange(25)+1.0
        angle = tilt + 5.0*pi*log(n)/(12.0*log(n[-1]))
        for theta in angle:
            x = self.para['rbody'] * cos(theta)
            y = self.para['rbody'] * sin(theta)
            stream = self.get_stream(x,y,'bx','bz', style=style)
            if stream.y[0] > self.para['rbody']:
                daymax = theta
                break
            savestream = stream
            stream.plot(ax)

        # Add IMF field lines.
        if DoImf:
            stream = savestream
            r = sqrt(stream.x**2 + stream.y**2)
            loc, = where(r==r.max())
            x_mp = stream.x[loc[0]]+0.15
            y_mp = stream.y[loc[0]]
            delx = 2.0
            for i, x in enumerate(arange(x_mp, 15.0, delx)):
                # From dayside x-line out and up:
                y =y_mp-x_mp+x
                stream = self.get_stream(x, y, 'bx', 'bz', style=style)
                stream.plot(ax)
                # From top of magnetosphere down:
                y =x_mp+15.0-x+delx/3.0
                stream = self.get_stream(x-delx/3.0, y, 'bx', 
                                         'bz', style=style)
                stream.plot(ax)
                # From bottom of mag'sphere down:
                y =x_mp-10.0-x+2.0*delx/3.0
                stream = self.get_stream(x-2.0*delx/3.0, y, 'bx', 
                                         'bz', style=style)
                stream.plot(ax)

        # Night side:
        angle = pi + tilt + pi*log(n)/(2.5*log(n[-1]))
        for theta in angle:
            x = self.para['rbody'] * cos(theta)
            y = self.para['rbody'] * sin(theta)
            stream = self.get_stream(x,y,'bx','bz', style=style)
            if stream.open:
                nightmax = theta
                break
            savestream = stream
            stream.plot(ax)

        # March down tail.
        stream = savestream
        r = sqrt(stream.x**2 + stream.y**2)
        loc, = where(r==r.max())
        x1 = stream.x[loc[0]]
        y1 = stream.y[loc[0]]
        x = x1
        y = y1
        while (x-1.5)>min(self.grid['x']):
            #print "Closed extension at ", x-1.5, y
            #ax.plot(x-1.5, y, 'g^', ms=10)
            stream = self.get_stream(x-1.5, y, 'bx', 'bz', style=style)
            r = sqrt(stream.x**2 + stream.y**2)
            if stream.open:
                break
            stream.plot(ax)
            loc, = where(r==r.max())
            x = stream.x[loc[0]]
            y = stream.y[loc[0]]

        if x1 == x:
            stream = self.get_stream(x1+1.0, y1, 'bx', 'bz')
            r = sqrt(stream.x**2 + stream.y**2)
            loc, = where(r==r.max())
            x1 = stream.x[loc[0]]
            y1 = stream.y[loc[0]]

        # Add more along neutral sheet.
        m = (y-y1)/(x-x1)
        #print "Slope = ", m
        #print "From y, y1 = ", y, y1
        #print "and x, x1 = ", x, x1
        #ax.plot([x,x1],[y,y1], 'wo', ms=10)
        xmore = arange(x, self.grid['x'].min(), -1.5)
        ymore = m*(xmore-x)+y
        for x, y  in zip(xmore[1:], ymore[1:]):
            #print "Starting at x, y = ", x, y
            stream = self.get_stream(x, y, 'bx', 'bz', style=style)
            stream.plot(ax)
        #ax.plot(xmore, ymore, 'r+')

        # Add open field lines.
        if DoOpen:
            for theta in linspace(daymax,0.99*(2.0*(pi+tilt))-nightmax, 15):
                x = self.para['rbody'] * cos(theta)
                y = self.para['rbody'] * sin(theta)
                stream = self.get_stream(x,y,'bx','bz')
                if stream.open: stream.plot(ax)
                x = self.para['rbody'] * cos(theta+pi)
                y = self.para['rbody'] * sin(theta+pi)
                stream = self.get_stream(x,y,'bx','bz')
                if stream.open: stream.plot(ax)


    def add_body(self, ax=None, facecolor='k', **extra_kwargs):
        '''
        Creates a circle of radius=self.para['rbody'] and returns the
        MatPlotLib Ellipse patch object for plotting.  If an axis is specified
        using the "ax" keyword, the patch is added to the plot.
        Default color is black; extra keywords are handed to the Ellipse
        generator function.
        '''
        from matplotlib.patches import Ellipse

        if 'rbody' not in self.para:
            raise KeyError('rbody not found in self.para!')

        dbody = 2.0 * self.para['rbody']
        body = Ellipse((0,0),dbody,dbody,facecolor=facecolor, zorder=1000,
                       **extra_kwargs)
        if ax != None:
            ax.add_artist(body)
        
    def contourf(self, ax, dim1, dim2, data, *extra_args, **extra_kwargs):
        '''
        Create an axis object containing a contourf
        '''
        import matplotlib.pyplot as plt

        assert self.gridtype == 'Regular', \
            'Irregular grid cannot be contoured.  Use self.regrid() first.'
        
        if not isinstance(ax, plt.Axes):
            raise TypeError('ax must be an Axes instance!')
        contour = ax.contourf(self.grid[dim1],self.grid[dim2], 
                              self.data[data], *extra_args, **extra_kwargs)
        return contour
    
    def contouru(self, dim1, dim2, data, 
                 cellsize = 1.0, dim1range = [0,1], dim2range = [0,1]):
        '''
        For unstructured grids, matplotlib's contour function will
        not work.  It is critical to first interpolate onto a regular
        grid and then use the contour function.
        
        contouru acts as a wrapper for matplotlib's contour that first
        performs the interpolation before calling contour.  The value 
        returned is a matplotlib contour object.

        Contouru should be used for quick, first-look plots only.  The
        regridding of the data happens every time it is called, making it
        slow.  The result of the rebinning is not saved, so advanced 
        features such as streamline tracing will not be made available.
        
        Usage:

        kwargs:
        cellsize -- the new cell size resolution to which to interpolate.
        dim1range and dim2range -- two element lists that set the new range
        of the data.  Useful for reducing the domain in order to increase
        the speed of this function.

        '''
        import matplotlib.pyplot as plt
        from matplotlib.mlab import griddata

        # Set default range as entire domain.
        dim1range = [ min(self.grid[dim1]), max(self.grid[dim1]) ]
        dim2range = [ min(self.grid[dim2]), max(self.grid[dim2]) ]

        # Regrid data.
        grid1 = np.arange(dim1range[0], dim1range[1]+cellsize, cellsize)
        grid2 = np.arange(dim2range[0], dim2range[1]+cellsize, cellsize)
        newdat= griddata(self.grid[dim1], self.grid[dim2],
                         self.data[data], grid1, grid2)
        cont = plt.contourf(grid1, grid2, newdat)
