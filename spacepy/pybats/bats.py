'''
A PyBats module for handling input, output, and visualization of 
binary SWMF output files taylored to BATS-R-US-type data.
'''

import numpy as np
from spacepy.pybats import IdlBin, LogFile

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
        grid = bats['grid'].attrs['dims']

        # Trace forward
        x1, y1 = trc(bats[self.xvar], bats[self.yvar], 
                     self.xstart, self.ystart,
                     bats[grid[0]], bats[grid[1]])
        # Trace backward
        x2, y2 = trc(bats[self.xvar], bats[self.yvar], 
                     self.xstart, self.ystart,
                     bats[grid[0]], bats[grid[1]], ds=-0.1)
        # Join lines together such that point 0 is beginning of line
        # and point -1 is the end (when moving parallel to line.)
        self.x = array(x2[::-1].tolist() + x1[1:].tolist())
        self.y = array(y2[::-1].tolist() + y1[1:].tolist())

        # Check if line is closed to body.
        if bats.attrs.has_key('rbody'):
            r1 = sqrt(self.x[0]**2.0  + self.y[0]**2.0)
            r2 = sqrt(self.x[-1]**2.0 + self.y[-1]**2.0)
            if (r1 < bats.attrs['rbody']) and (r2 < bats.attrs['rbody']):
                self.open = False

    def plot(self, ax, *args, **kwargs):
        '''
        Add streamline to axes object "ax". 
        '''
        ax.plot(self.x, self.y, self.style, *args, **kwargs)

class Bats2d(IdlBin):
    '''
    An object class of parent pybats.idlbin taylored to BATS-R-US output.
    This function requires the Matplotlib griddata function.
    
    '''
    # Init by calling IdlBin init and then building qotree, etc.
    def __init__(self, filename):
        import spacepy.pybats.qotree as qo
        reload(qo)
        from numpy import array
        # Read file.
        IdlBin.__init__(self, filename)

        # Parse grid into quad tree.
        if self['grid'].attrs['gtype'] != 'Regular':
            xdim, ydim = self['grid'].attrs['dims'][0:2]
            self.tree=qo.QTree(array([self[xdim],self[ydim]]))
            

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
        if not self.has_key('rho'):
            raise AttributeError("Number density not found; required.")
        if not self.has_key('p'):
            raise AttributeError("Pressure not found.")
        
        # Calculate value then adjust units.
        self['t'] = self['p'] / self['rho']
        if units == 'ev':
            self['t'] = self['t'] * 6250.0
        elif units == 'kev':
            self['t'] = self['t'] * 6.250
        elif units == 'k':
            self['t'] = self['t'] * 72432275.82
        else:
            raise ValueError("Units of type %s not recognized." % units)
        self['t'].attrs={'units':units}

    def calc_b(self):
        '''
        Calculates total B-field strength using all three B components.
        Retains units of components.  Additionally, the unit vector
        b-hat is calculated as well.
        '''
        from numpy import sqrt

        self['b'] = sqrt(self['bx']**2.0 + self['by']**2.0 + self['bz']**2.0)
        self['b'].attrs={'units':self['bx'].attrs['units']}

        self['bx_hat'] = self['bx'] / self['b']
        self['by_hat'] = self['by'] / self['b']
        self['bz_hat'] = self['bz'] / self['b']

        self['bx_hat'].attrs={'units':'unitless'}
        self['by_hat'].attrs={'units':'unitless'}
        self['bz_hat'].attrs={'units':'unitless'}

    def calc_Econv(self):
        '''
        Calculates the convective electric field, -UxB.  Works for default
        MHD units of nT and km/s; if these units are not correct, an 
        exception will be raised.  Returns E_convective in mV/m.
        '''
        from copy import copy

        # Some quick declarations for more readable code.
        ux = self['ux']; uy=self['uy']; uz=self['uz']
        bx = self['bx']; by=self['by']; bz=self['bz']

        # Check units.  Should be nT(=Volt*s/m^2) and km/s.
        if (bx.attrs['units']!='nT') or (ux.attrs['units']!='km/s'):
            raise Exception('Incorrect units!  Should be km/s and nT.')

        # Calculate; return in millivolts per meter
        self['ex'] = -1.0*(uy*bz - uz*by) / 1000.0
        self['ey'] = -1.0*(uz*bx - ux*bz) / 1000.0
        self['ez'] = -1.0*(ux*by - uy*bx) / 1000.0
        self['ex'].attrs={'units':'mV/m'}
        self['ey'].attrs={'units':'mV/m'}
        self['ez'].attrs={'units':'mV/m'}
        

    def calc_beta(self):
        '''
        Calculates plasma beta (ratio of plasma to magnetic pressure, 
        indicative of who - plasma or B-field - is "in charge" with regards
        to flow.
        Assumes:
        -pressure in units of nPa
        -B in units of nT.
        Resulting value is unitless.
        Values where b_total = zero are set to -1.0 in the final array.
        '''
        from numpy import pi

        if not self.has_key('b'):
            self.calc_b()
        mu_naught = 4.0E2 * pi # Mu_0 x unit conversion (nPa->Pa, nT->T)
        temp_b = self['b']**2.0
        temp_b[temp_b<1E-8] =  -1.0*mu_naught*self['p'][temp_b==0.0]
        temp_beta=mu_naught*self['p']/temp_b
        #temp_beta[self['b']<1E-9] = -1.0
        self['beta']=temp_beta
        self['beta'].attrs={'units':'unitless'}

    def calc_jxb(self):
        '''
        Calculates the JxB force assuming:
        -Units of J are uA/m2, units of B are nT.
        Under these assumptions, the value returned is force density (N/m^3).
        '''
        from numpy import sqrt
        # Unit conversion (nT, uA, cm^-3 -> T, A, m^-3) to N/m^3.
        conv = 1E-15
        #conv = 1E-21 / self['rho'] <- old conversion to pure force (N).
        # Calculate curl, convert units.
        self['jbx']=(self['jy']*self['bz']-self['jz']*self['by'])*conv
        self['jby']=(self['jz']*self['bx']-self['jx']*self['bz'])*conv
        self['jbz']=(self['jx']*self['by']-self['jy']*self['bx'])*conv
        self['jb']=sqrt(self['jbx']**2 +
                        self['jby']**2 +
                        self['jbz']**2)

    def calc_alfven(self):
        '''
        Calculate the Alfven speed, B/(mu*rho)^(1/2) in km/s.
        The variable is saved under key "alfven" in self.data.
        '''
        from numpy import sqrt, pi
        if not self.has_key('b'):
            self.calc_b()
        #M_naught * conversion from #/cm^3 to kg/m^3
        mu_naught = 4.0E-7 * pi * 1.6726E-27 * 1.0E6
        # Alfven speed in km/s:
        self['alfven']= (self['b']*1E-9) / sqrt(mu_naught*self['rho']) /1000.
        self['alfven'].attrs={'units':'km/s'}

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
                except AttributeError, Error:
                    print('WARNING: Did not perform %s: %s' % (command, Error))

    #####################
    # Other calculations
    #####################
    def gradP_regular(self, cellsize=None, dim1range=-1, dim2range=-1):
        '''
        Calculate pressure gradient on a regular grid.
        Note that if the Bats2d object is not on a regular grid, one of 
        two things will happen.  If kwarg cellsize is set, the value of 
        cellsize will be used to call self.regrid and the object will
        be regridded using a cellsize of cellsize.  Kwargs dim1range and
        dim2range can be used in the same way they are used in self.regrid
        to restrict the regridding to a smaller domain.  If cellsize is
        not set and the object is on an irregular grid, an exception is
        raised.

        The gradient is calculated using numpy.gradient.  The output units
        are force density (N/m^3).  Three variables are added to self.data:
        gradp(dim1), gradp(dim2), gradp.  For example, if the object is an
        equatorial cut, the variables gradpx, gradpy, and gradp would be
        added representing the gradient in each direction and then the 
        magnitude of the vector.
        '''
        from numpy import gradient, sqrt

        if self.gridtype != 'Regular':
            if not cellsize:
                raise ValueError('Grid must be regular or ' + \
                    'cellsize must be given.')
            self.regrid(cellsize, dim1range=dim1range, dim2range=dim2range)

        # Order our dimensions alphabetically.
        newvars = []
        for key in sorted(self.grid.keys()):
            newvars.append('gradp'+key)

        dx=self.resolution*6378000.0 # RE to meters
        p =self['p']*10E-9        # nPa to Pa
        self[newvars[0]], self[newvars[1]]=gradient(p, dx, dx)
        self['gradp']=sqrt(self[newvars[0]]**2.0 + self[newvars[1]]**2.0)

    def regrid(self, cellsize=1.0, dim1range=-1, dim2range=-1, debug=False):
        '''
        Re-bin data to regular grid of spacing cellsize.  Action is 
        performed on all data entries in the bats2d object.

        '''
        from matplotlib.mlab import griddata

        if self['grid'].attrs['gtype'] == 'Regular': return
        
        # Order our dimensions alphabetically.
        dims = self['grid'].attrs['dims']
        if debug: print("Ordered dimensions: ", dims)

        # Check to see if dimranges are 2-element lists.
        # If not, either set defaults or raise exceptions.
        if dim1range == -1:
            dim1range = [self[dims[0]].min(),self[dims[0]].max()]
        else:
            if isinstance(dim1range, ( type(()), type([]) ) ):
                if len(dim1range) != 2:
                    raise ValueError('dim1range must have two elements!')
            else: raise TypeError('dim1range must be a tuple or list!')
        if dim2range == -1:
            dim2range = [self[dims[1]].min(),self[dims[1]].max()]
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

        for key in self.keys():
            # Skip grid-type entries.
            if key in (self['grid'].attrs['dims']+['grid']): continue
            self[key] = griddata(self[dims[0]], self[dims[1]],
                                 self[key], grid1, grid2)

        # Change grid, gridtype, gridsize, and npoints to match new layout.
        self['grid'].attrs['gtype'] = 'Regular'
        self['grid'].attrs['npoints'] = len(grid1) * len(grid2)
        self['grid'].attrs['resolution'] = cellsize
        self[dims[0]] = grid1; self[dims[1]] = grid2

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
    def add_grid_plot(self, target=None, loc=111, DoLabel=True, 
                      title='BATS-R-US Grid Layout'):
        '''
        Create a plot of the grid resolution by coloring regions of constant
        resolution.  Kwarg "target" specifies where to place the plot and can
        be a figure, an axis, or None.  If target is a figure, a new subplot
        is created.  The subplot location is set by the kwarg "loc", which
        defaults to 111.  If target is an axis, the plot is placed into that
        axis object.  If target is None, a new figure and axis are created
        and used to display the plot. 

        Resolution labels can be disabled by setting kwarg DoLabel to False.

        Plot title is set using the 'title' kwarg, defaults to 'BATS-R-US
        Grid Layout'.

        Note that if target is not an axis object, the axis will automatically
        flip the direction of positive X GSM and turn on equal aspect ratio.
        In other words, if you want a customized axis, it's best to create
        one yourself.

        Figure and axis, even if none given, are returned.
        '''
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.ticker import MultipleLocator
        from numpy import linspace

        # Get dimensions over which we shall plot.
        xdim, ydim = self['grid'].attrs['dims'][0:2]

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
            ax.set_aspect('equal')
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(loc)
            ax.set_aspect('equal')
        # Set plot range based on quadtree.
        ax.set_xlim([self.tree[1].lim[0],self.tree[1].lim[1]])
        ax.set_ylim([self.tree[1].lim[2],self.tree[1].lim[3]])
        # Plot.

        for key in self.tree.keys():
            self.tree[key].plotbox(ax)
        self.tree.plot_res(ax)

        # Customize plot.
        ax.set_xlabel('GSM %s' % xdim.upper())
        ax.set_ylabel('GSM %s' % ydim.upper())
        ax.set_title(title)
        if xdim=='x':
            ax.invert_xaxis()
        if ydim=='y':
            ax.invert_yaxis()
        self.add_body(ax)

        return fig, ax


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

#        # Add more along neutral sheet.
#        m = (y-y1)/(x-x1)
#        #print "Slope = ", m
#        #print "From y, y1 = ", y, y1
#        #print "and x, x1 = ", x, x1
#        #ax.plot([x,x1],[y,y1], 'wo', ms=10)
#        xmore = arange(x, self.grid['x'].min(), -1.5)
#        ymore = m*(xmore-x)+y
#        for x, y  in zip(xmore[1:], ymore[1:]):
#            #print "Starting at x, y = ", x, y
#            stream = self.get_stream(x, y, 'bx', 'bz', style=style)
#            stream.plot(ax)
#        #ax.plot(xmore, ymore, 'r+')

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


    def add_planet(self, ax=None, rad=1.0, ang=0.0, **extra_kwargs):
        '''
        Creates a circle of radius=self.para['rbody'] and returns the
        MatPlotLib Ellipse patch object for plotting.  If an axis is specified
        using the "ax" keyword, the patch is added to the plot.

        Unlike the add_body method, the circle is colored half white (dayside)
        and half black (nightside) to coincide with the direction of the 
        sun. Additionally, because the size of the planet is not intrinsically
        known to the MHD file, the kwarg "rad", defaulting to 1.0, sets the
        size of the planet.

        Extra keywords are handed to the Ellipse generator function.
        '''

        from matplotlib.patches import Circle, Wedge

        if not self.attrs.has_key('rbody'):
            raise KeyError('rbody not found in self.para!')

        body = Circle((0,0), rad, fc='w', zorder=1000, **extra_kwargs)
        arch = Wedge((0,0), rad, 90.+ang, -90.+ang, fc='k', 
                     zorder=1001, **extra_kwargs)
        
        if ax != None:
            ax.add_artist(body)
            ax.add_artist(arch)

        return body, arch

    def add_body(self, ax=None, facecolor='lightgrey', DoPlanet=True, ang=0.0, 
                 **extra_kwargs):
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

        if not self.attrs.has_key('rbody'):
            raise KeyError('rbody not found in self.para!')

        dbody = 2.0 * self.attrs['rbody']
        body = Ellipse((0,0),dbody,dbody,facecolor=facecolor, zorder=999,
                       **extra_kwargs)

        if DoPlanet:
            self.add_planet(ax, ang=ang)
        if ax != None:
            ax.add_artist(body)

    def add_pcolor(self, dim1, dim2, value, zlim=None, target=None, loc=111, 
                   title=None, xlabel=None, ylabel=None, dolabel=True,
                   ylim=None, xlim=None, add_cbar=False, clabel=None,
                   add_body=True, *args, **kwargs):
        '''doc forthcoming'''

        import matplotlib.pyplot as plt

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax  = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,10))
            ax  = fig.add_subplot(loc)

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=self[value].min(); zlim[1]=self[value].max()

        if self['grid'].attrs['gtype']=='Regular':
            pass
        else:
            # Indices corresponding to QTree dimensions:
            ix=self['grid'].attrs['dims'].index(dim1)
            iy=self['grid'].attrs['dims'].index(dim2)
            for k in self.tree.keys():
                # Plot only leafs of the tree.
                if not self.tree[k].isLeaf: continue
                leaf=self.tree[k]
                x=leaf.cells[ix]
                y=leaf.cells[iy]
                z=self[value][leaf.locs]
                pcol=ax.pcolormesh(x,y,z,vmin=zlim[0],vmax=zlim[1],**kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(pcol, pad=0.01)
            if clabel==None: 
                clabel="%s (%s)" % (value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        if ylabel==None: ylabel='%s ($R_{E}$)'%dim2.upper()
        if xlabel==None: xlabel='%s ($R_{E}$)'%dim1.upper()
        ax.set_ylabel(ylabel); ax.set_xlabel(xlabel)
        if type(xlim)==type([]) and len(xlim)==2:
            ax.set_xlim(xlim)
        if type(ylim)==type([]) and len(ylim)==2:
            ax.set_ylim(ylim)

        # Add body/planet.
        self.add_body(ax)

        return fig, ax, pcol, cbar                              

    def add_contour(self, dim1, dim2, value, levs=30, target=None, loc=111, 
                    title=None, xlabel=None, ylabel=None, dolabel=True,
                    ylim=None, xlim=None, add_cbar=False, clabel=None,
                    filled=True, add_body=True, dolog=False, logrange=[.01,50.],
                    *args, **kwargs):
        '''doc forthcoming.'''
        
        import matplotlib.pyplot as plt
        from matplotlib.colors import (LogNorm, Normalize)
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax  = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,10))
            ax  = fig.add_subplot(loc)

        # Set contour command based on grid type.
        if self['grid'].attrs['gtype'] != 'Regular':  # Non-uniform grids.
            if filled:
                contour=ax.tricontourf   
            else:
                contour=ax.tricontour   
        else:   # Uniform grids.
            if filled:
                contour=ax.contourf
            else:
                contour=ax.contour

        # If a log-scale plot, take log of values.
        if dolog:
            # Set cbar tick locators, contour norms.
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
            # Trim small values.
            z=np.where(self[value]>logrange[0], self[value], 1.01*logrange[0])
            # Set up levels if not explicitly set.
            if not isinstance(levs, (list, np.ndarray)):
                levs=np.power(10, np.linspace(np.log10(logrange[0]),
                                              np.log10(logrange[1]), levs))
        else:
            # Use defaults.  We will need to edit this to
            # not overwrite kwargs.
            norm=None
            ticks=None
            fmt=None
            z=self[value]

        # Plot contour.
        cont=contour(self[dim1],self[dim2],z,levs,*args, norm=norm, **kwargs)
        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel==None: 
                clabel="%s (%s)" % (value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        if ylabel==None: ylabel='%s ($R_{E}$)'%dim2.upper()
        if xlabel==None: xlabel='%s ($R_{E}$)'%dim1.upper()
        ax.set_ylabel(ylabel); ax.set_xlabel(xlabel)
        if type(xlim)==type([]) and len(xlim)==2:
            ax.set_xlim(xlim)
        if type(ylim)==type([]) and len(ylim)==2:
            ax.set_ylim(ylim)

        # Add body/planet.  Determine where the sun is first.
        if dim1.lower()=='x':
            ang=0.0
        elif dim2.lower()=='x':
            ang=90.0
        else: ang=0.0
        self.add_body(ax, ang=ang)

        return fig, ax, cont, cbar

class BatsMag(object):
    '''
    Holds data for a single magnetometer station.
    '''
    def __getitem__(self, key):
        return self.data[key]
    def __setitem__(self, key, value):
        self.data[key]=value
    def keys(self):
        return self.data.keys()

    def __init__(self, nlines, time, gmvars=(), ievars=()):
        from numpy import zeros

        self.time=time
        self.nlines=nlines
        
        self.data={}
        self.loc ={}
        self.loc['x']=np.zeros(nlines)
        self.loc['y']=np.zeros(nlines)
        self.loc['z']=np.zeros(nlines)

        # Create IE and GM specific containers.
        for key in gmvars:
            self['gm_'+key]=np.zeros(nlines)
        for key in ievars:
            self['ie_'+key]=np.zeros(nlines)
            
    def parse_gmline(self, i, line, namevar):
        '''
        Parse a single line from a GM_mag*.out file and put into
        the proper place in the magnetometer arrays.

        Usage: self.parse_gmline(i, line, namevar)
        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts=line.split()
        self.loc['x'][i]=float(parts[9])
        self.loc['y'][i]=float(parts[10])
        self.loc['z'][i]=float(parts[11])
        for j, key in enumerate(namevar):
            self['gm_'+key][i]=float(parts[j+12])

    def parse_ieline(self, i, line, namevar):
        '''
        Parse a single line from a IE_mag*.out file and put into
        the proper place in the magnetometer arrays.

        Usage: self.parse_gmline(i, line, namevar)
        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts=line.split()
        for j, key in enumerate(namevar):
            self['ie_'+key][i]=float(parts[j+11])

    def recalc(self):
        '''
        Calculate total dB from GM and IE; compute H-component.
        '''
        from numpy import sqrt, zeros

        # New containers:
        self['totaln']=np.zeros(self.nlines)
        self['totale']=np.zeros(self.nlines)
        self['totald']=np.zeros(self.nlines)

        for key in self.keys():
            if key[-2:]=='Bn':
                self['totaln']=self['totaln']+self[key]
            if key[-2:]=='Be':
                self['totale']=self['totale']+self[key]
            if key[-2:]=='Bd':
                self['totald']=self['totald']+self[key]
        #self['totalh']=sqrt(self['totaln']**2.0 + self['totale']**2.0)

    def plot(self, value, target=None, loc=111, label=None, **kwargs):
        '''
        Plot value 'value' against self.time onto target "target" given as a 
        kwarg.  Target may be a figure or axis but defaults to None, in 
        which a new figure and axis are created.  If target is a matplotlib
        figure, a new axis is created at subplot location 111 (can be changed
        using kwarg "loc").  If target is a matplotlib Axes object, the line
        is added to the plot as if Axes.plot() was used.  The line label 
        defaults to the value key but can be customized with the "label"
        kwarg.  All extra kwargs are handed to Axes.plot.

        Three values are returned: the Figure object, Axis object, and 
        newly created line object.  These can be used to further customize
        the figure, axis, and line as necessary.

        Example: Plot dBn onto an existing axis with line color blue, 
        line style dashed, and line label "Wow!":
        self.plot(ax, 'dBn', label='Wow!', lc='b')
        '''

        import matplotlib.pyplot as plt
        from spacepy.pybats import apply_smart_timeticks

        if not label:
            label=value

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
            line=ax.plot(self.time, self[value], label=label, **kwargs)
            apply_smart_timeticks(ax, self.time, dolabel=True)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
            line=ax.plot(self.time, self[value], label=label, **kwargs)
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(loc)
            line=ax.plot(self.time, self[value], label=label, **kwargs)
            apply_smart_timeticks(ax, self.time, dolabel=True)

        return fig, ax, line

class MagFile(object):
    '''
    BATS-R-US magnetometer files are powerful tools for both research and
    operations.  BatsMag objects open, parse, and visualize such output.

    The dB calculated by the SWMF requires two components: GM (BATSRUS)
    and IE (Ridley_serial).  The data is spread across two files: GM_mag*.dat
    and IE_mag*.dat.  The former contains dB caused by FACs and the changing
    global field.  The latter contains the dB caused by Pederson and Hall 
    currents in the ionosphere. 
    '''
    
    def __init__(self, filename, ie_name=None, find_ie=False):
        self.gmfile=filename
        self.iefile=ie_name
        self.readfiles()

    def __getitem__(self, key):
        'Grab item from self.mag[key]'
        return self.mag[key]

    def __setitem__(self, key, value):
        'Set item in self.mag[key].'
        self.mag[key]=value

    def readfiles(self):
        '''
        Read and parse GM file and IE file (if name given.)
        '''
        import datetime as dt
        from numpy import zeros

        # Slurp lines.
        infile = open(self.gmfile, 'r')
        lines=infile.readlines()
        infile.close()

        # Parse header.
        trash=lines.pop(0) # Get station names.
        nmags=int((trash.split(':')[0]).split()[0])
        names=trash.split(':')[1] # Remove number of magnetometers.
        self.namemag = names.split()
        # Check nmags vs number of mags in header.
        if nmags != len(self.namemag):
            raise BaseException('ERROR: GM file claims %i magnetomers, lists %i'
                % (nmags, len(self.namemag)))
        trash=lines.pop(0)
        self.gm_namevar = trash.split()[12:] #Vars skipping time, iter, and loc.
        self.nmag=len(self.namemag)
        nlines = len(lines)/self.nmag

        # If there is an IE file, Parse that header, too.
        if self.iefile:
            infile=open(self.iefile, 'r')
            ielns =infile.readlines()
            infile.close()
            trash=ielns.pop(0)
            nmags=int((trash.split(':')[0]).split()[0])
            iestats=(trash.split(':')[1]).split()
            # Check nmags vs number of mags in header.
            if nmags != len(iestats):
                raise BaseException('ERROR: IE file claims %i magnetomers, lists %i' % 
                                    (nmags, len(self.namemag)))
            if iestats != self.namemag:
                raise RuntimeError("Files do not have matching stations!")
            self.ie_namevar=ielns.pop(0).split()[11:]
            if (len(ielns)/self.nmag) != (nlines-1):
                #raise RuntimeError, "Files do nat have same number of lines!"
                # Adjust number of lines read to match (GM often has more.)
                print('Number of lines do not match: GM=%d, IE=%d!' % \
                    (nlines-1, len(ielns)/self.nmag))
                nlines=min(ielns, nlines-1)
        else:
            self.ie_namevar=()

        # Build containers.
        self.time=np.zeros(nlines, dtype=object)
        self.iter=np.zeros(nlines)
        self.mag={}
        for name in self.namemag:
            self[name]=BatsMag(nlines, self.time, 
                               self.gm_namevar, self.ie_namevar)

        # Read file data.
        for i in range(nlines):
            line = lines.pop(0)
            parts=line.split()
            self.iter[i]=int(parts[0])
            self.time[i]=dt.datetime(
                int(parts[1]), #year
                int(parts[2]), #month
                int(parts[3]), #day
                int(parts[4]), #hour
                int(parts[5]), #minute
                int(parts[6]), #second
                int(parts[7])*1000 #microsec
                )
            self[self.namemag[0]].parse_gmline(i, line, self.gm_namevar)
            for j in range(1, self.nmag):
                self[self.namemag[j]].parse_gmline(i, lines.pop(0), 
                                                   self.gm_namevar)
            if self.iefile and i>0:
                line=ielns.pop(0)
                self[self.namemag[0]].parse_ieline(i, line, self.ie_namevar)
                for j in range(1, self.nmag):
                    self[self.namemag[j]].parse_ieline(i, ielns.pop(0), 
                                                       self.ie_namevar)
        self.recalc()
        
        # Get time res.
        self.dt=(self.time[1]-self.time[0]).seconds/60.0

        
        
    def recalc(self):
        '''
        Recalculate the values for each magnetometer.  This involves:
        -summing GM and IE contributions.
        -Calculating the H component.
        '''
        for mag in self.namemag:
            self[mag].recalc()

class GeoIndFile(LogFile):
    '''
    Geomagnetic Index files are a specialized BATS-R-US output that contain
    geomagnetic indices calculated from simulated ground-based magnetometers.
    Currently, the only index instituted is Kp through the faKe_p setup.  
    Future work will expand the system to include Dst, AE, etc.

    GeoIndFiles are a specialized subclass of pybats.LogFile.  It includes
    additional methods to quickly visualize the output, perform data-model
    comparisons, and more.
    '''

    def kp_quicklook(self, target=None, pos=111, label='fa$K$e$_{P}$', 
                     **kwargs):
        '''
        Similar to "dst_quicklook"-type functions, this method fetches observed
        Kp from the web and plots it alongside the Kp read from the GeoInd file.
        Usage:
        >>> obj.kp_quicklook(target=SomeMplTarget)
        The target kwarg works like in other PyBats plot functions: it can be
        a figure, an axes, or None, and it determines where the plot is placed.

        Other kwargs customize the line.  Label defaults to fa$K$e$_{P}$, extra
        kwargs are passed to pyplot.plot.
        '''
        import matplotlib.pyplot as plt
        
        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(pos)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(111)
        
        ax.plot(self.time, self['Kp'], label=label,**kwargs)
        ax.set_ylabel('$K_{P}$')
        ax.set_xlabel('Time from '+ self.time[0].isoformat()+' UTC')
        apply_smart_timeticks(ax, self.time)

        try:
            import spacepy.pybats.kyoto as kt
        except ImportError:
            print("kyotodst package unavailable.")
            return fig, ax
        
        try:
            stime = self.time[0]; etime = self.time[-1]
            if not hasattr(self, 'obs_kp'):
                self.obs_kp = kt.kpwebfetch(stime.year, stime.month, 
                                            etime.year, etime.month)
        except BaseException, args:
            print('WARNING! Failed to fetch Kyoto Kp: ', args)
        else:
            self.obs_kp.histplot(target=ax, color='k', ls='--')
            ax.legend()
            apply_smart_timeticks(ax, self.time)

        return fig, ax
