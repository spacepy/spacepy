#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the CIMI Model
, one of the choices for the IM module in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
import spacepy.plot.apionly
from spacepy.plot import set_target
from spacepy.pybats import PbData, IdlFile
from spacepy.datamodel import dmarray




class CimiEq(PbData):
    '''
    Class for loading a CIMIeq.out file.
    At instantiation time, user may wish to set the start date and time of
    the simulation using the starttime kwarg.  If not given, start time
    will default to Jan. 1st, 2000, 00:00UT.
    '''
    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(CimiEq, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(starttime)


    def __repr__(self):
        return 'CimiEq single output file %s' % (self.attrs['file'])

    def _read(self, starttime):
        '''
        Read CimiEq file; should only be called upon instantiation.
        '''

        # Reads file assuming it is ascii
        # Slurp whole file.
        lines=open(self.attrs['file'], 'r').readlines()

        # Determine size of file.
        nTimes=lines.count(lines[0])
        nLat =int(lines[2].split()[0])
        nLon =int(lines[2].split()[1])


        self.attrs['nLat'] =nLat
        self.attrs['nLon'] =nLon
        self.attrs['nTime']=nTimes
        time = float((lines[1].split())[1])
        self.attrs['strtime'] = "{:04.0f}h{:02.0f}m{:06.3f}s".format(
            np.floor(time//3600), np.floor(time % 3600//60.0), time % 60.0)
        self.attrs['time'] = starttime+dt.timedelta(seconds=time)


        nDataLength = nLat*nLon

        # Start building time array.
        self['time']=dmarray(np.zeros(nTimes, dtype=object))
        self['time'].attrs['units'] = 's'
        # Get variable names; slice off g and rbody
        var=(lines[4].split())[0:-2]
        #Initialize grids to hold variables, strip unit infomation if it exits and store the units in the attrs dictionary of the variabls in var.
        #For CIMIEq_*.out type data file var names are likely written as 'VAR[UNIT]'.
        for i, v in enumerate(var):
            if '[' in v:
                v, unit = v[:-1].split('[')#v[:-1] cuts off the ], by spliting on [ we get list(VAR,UNIT).
                self[v]=dmarray(np.zeros((nTimes, nDataLength)))
                self[v].attrs['units']=unit
                var[i] = v
            elif v == 'x' or v == 'y':
                self[v]=dmarray(np.zeros((nTimes, nDataLength)))
                self[v].attrs['units'] = 'Re'
                var[i] = v
            else:
                self[v]=dmarray(np.zeros((nTimes, nDataLength)))
                self[v].attrs['units'] = None
                var[i] = v
        self._rawvar=var
        # initialize arrays to hold data

        # Loop through rest of data to fill arrays.
        for i in range(nTimes):
            t=float((lines[i*(nDataLength+5)+1].split())[1])
            self['time'][i]=starttime+dt.timedelta(seconds=t)
            for j, l in enumerate(lines[i*(nDataLength+5)+5:(i+1)*(nDataLength+5)]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][i,j]=float(parts[k])

    def m3_to_cc(self):
        """
        This function transforms density data from m3 to cm3.
        It does this by searching for self[VAR].attrs['units'] == '/m3'.
        If you're desnity values are not saved with this unit exactly it won't work!
        """
        for var in self.keys():
            if var.attrs['units'] == '/m3':
                self[var] = self[var] / 1e6
                self[var].attrs['units'] = '/cc'

###############################################################################
############################START PLOTTING FUNCTIONS#########################
###############################################################################
    def add_planet(self, ax=None, rad=1.0, ang=0.0, **extra_kwargs):
        '''
        Creates a circle of radius=self.attrs['rbody'] and returns the
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

        if 'rbody' not in self.attrs:
            raise KeyError('rbody not found in self.attrs!')

        body = Circle((0, 0), rad, fc='w', zorder=1000, **extra_kwargs)
        arch = Wedge((0, 0), rad, 90.+ang, -90.+ang, fc='k',
                     zorder=1001, **extra_kwargs)

        if ax is not None:
            ax.add_artist(body)
            ax.add_artist(arch)

        return body, arch

    def add_body(self, ax=None, facecolor='lightgrey', DoPlanet=True, ang=0.0,
                 **extra_kwargs):
        '''
        Creates a circle of radius=self.attrs['rbody'] and returns the
        MatPlotLib Ellipse patch object for plotting.  If an axis is specified
        using the "ax" keyword, the patch is added to the plot.
        Default color is light grey; extra keywords are handed to the Ellipse
        generator function.

        Because the body is rarely the size of the planet at the center of
        the modeling domain, add_planet is automatically called.  This can
        be negated by using the DoPlanet kwarg.
        '''
        from matplotlib.patches import Ellipse

        if 'rbody' not in self.attrs:
            raise KeyError('rbody not found in self.attrs!')

        dbody = 2.0 * self.attrs['rbody']
        body = Ellipse((0, 0), dbody, dbody, facecolor=facecolor, zorder=999,
                       **extra_kwargs)

        if DoPlanet:
            self.add_planet(ax, ang=ang)
        if ax is not None:
            ax.add_artist(body)

    @staticmethod
    def _adjust_dialplot(ax, rad, title='12', labelsize=15, add_body=True):
        '''
        Ram output is often visualized with equatorial dial plots.  This
        function quickly adjusts those plots to give a uniform, clean
        look and feel.
        '''
        from numpy import max, pi
        from matplotlib.ticker import MultipleLocator

        from spacepy.pybats.ram import add_body_polar

        # Constrain range of plot:
        ax.set_ylim([0, max(rad)])

        # Set MLT labels:
        lt_labels = ['06', title, '18', '00']
        xticks = [0, pi/2, pi, 3*pi/2]
        ax.set_xticks(xticks)
        ax.set_xticklabels(lt_labels)
        ax.tick_params('x', labelsize=14)

        # Set L labels and grid.  Turn off label at L=0.
        ax.yaxis.set_major_locator(MultipleLocator(2))
        ax.figure.canvas.draw()
        labs = [item.get_text() for item in ax.get_yticklabels()]
        if labelsize > 0:
            ax.set_yticklabels(labs, color='w', size=labelsize,
                            backgroundcolor='k')
            labels = ax.get_yticklabels()
            labels[0].set_visible(False)
        else:
            ax.set_yticklabels('')
        labels = ax.get_yticklabels()
        labels[0].set_visible(False)

        # Turn on grid.
        ax.grid(True, c='w', lw=1.5, ls=':')

        # Change background color so labels stand out.
        ax.set_facecolor('gray')
        if add_body:
            add_body_polar(ax)

    def add_dialplot(self, var, rmax=12, lon=None, nlev=51, zlim=None,
                      target=None, loc=111, title=None,
                      add_cbar=True, clabel=None,
                      dolog=False, add_body=True, filled=True,
                      *args, **kwargs):
        '''
        Create an equatorial "dial plot" of variable *var*: a polar contour
        plotted against L-shell/radius and magnetic local time (angle).

        Simple example:

        >>> self.add_dialplot('Plas', rmax=8, dolog=True, cmap='Greens_r')

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new *polar* axis is created to fill that figure at
        subplot location **loc**.  If **target** is a matplotlib Axes
        object, the plot is placed into that axis (it must already be
        a polar axis).

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *None* if not used.)

        Unlike :meth:`add_slice`, there is no *time* kwarg here -- a
        CimiEq object represents a single snapshot, so add_dialplot just
        plots whatever is currently stored on it (if the object somehow
        holds more than one time, the most recent is used).

        =========== ==========================================================
        Kwarg       Description
        =========== ==========================================================
        var         The variable to be plotted. (string)
        rmax        Maximum L-shell/radius shown on the dial. Defaults to 10.
        lon         Longitude/MLT grid, in radians, to plot against.
                    Defaults to the longitudes calculated from **x** and
                    **y**, which handles CIMI's grid shifting from one
                    snapshot to the next.
        nlev        Number of contour levels. Defaults to 51.
        zlim        Sets contour range. Defaults to variable max/min.
        target      Set plot destination. Defaults to new figure.
        loc         Set subplot location. Defaults to 111.
        title       Sets title of axes. Default is no title.
        add_cbar    Adds colorbar to plot. Defaults to *True*.
        clabel      Sets colorbar label. Defaults to **var** and units.
        dolog       Sets use of logarithmic scale/colorbar. Defaults to
                    *False*.
        add_body    Places planetary body at the center of the dial.
                    Defaults to *True*.
        filled      Use filled contours instead of contour lines. Defaults
                    to *True*.
        =========== ==========================================================

        Extra args and kwargs (e.g. *cmap*) are handed to the matplotlib
        contour command.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import LogLocator, LogFormatterMathtext

        fig, ax = set_target(target, figsize=(8, 8), loc=loc, polar=True)

        nLat, nLon = self.attrs['nLat'], self.attrs['nLon']

        def _grid(key):
            '''Pull a variable, drop any leading time axis down to the
            most recent snapshot, and reshape it onto CIMI's (nLat, nLon)
            grid.'''
            data = np.asarray(self[key])
            if data.ndim == 2:
                data = data[-1, :]
            return np.reshape(data, (nLat, nLon), order='F')

        value = _grid(var)

        # Longitude/MLT and radius (L-shell) grids.  CIMI's internal grid
        # isn't static through a run, so by default we calculate these
        # fresh from x & y (rather than assume a fixed set of longitudes)
        # to stay correct from snapshot to snapshot.
        x, y = _grid('x'), _grid('y')
        r = np.sqrt(x**2 + y**2)
        if lon is None:
            # +pi/2 lines up x-axis-sunward (noon) with the "12" position
            # used by _adjust_dialplot's MLT tick labels.
            lon = np.arctan2(y, x) + np.pi/2

        # Get max/min if none given.
        if zlim is None:
            zlim = [value.min(), value.max()]
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z = np.where(value > zlim[0], value, 1.01*zlim[0])
            z[z > zlim[1]] = zlim[1]
            norm = LogNorm()
            ticks = LogLocator()
            fmt = LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z = value
            norm = None
            ticks = None
            fmt = None

        # Repeat the first longitude column at +360deg so the contour
        # wraps cleanly all the way around the dial instead of leaving
        # a wedge-shaped gap.
        lon = np.concatenate([lon, lon[:, [0]] + 2*np.pi], axis=1)
        r = np.concatenate([r, r[:, [0]]], axis=1)
        z = np.concatenate([z, z[:, [0]]], axis=1)

        # Plot contour.
        contour = ax.contourf if filled else ax.contour
        cont = contour(lon, r, np.asarray(z), levs, *args, norm=norm,
                       **kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(cont, ax=ax, ticks=ticks, format=fmt,
                                pad=0.1, shrink=.9)
            if clabel is None:
                if 'units' in self[var].attrs:
                    clabel = f"{var} ({self[var].attrs['units']})"
                else:
                    clabel = f"{var}"
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        if title:
            ax.set_title(title)

        # Rotate so noon points to the right edge, then apply the
        # standard dial-plot look and feel (grid, MLT/L labels, body).
        ax.set_theta_zero_location('S')
        self._adjust_dialplot(ax, rmax, labelsize=14, add_body=add_body)

        return fig, ax, cont, cbar

    def add_slice(self, var, time, nlev=51, zlim=None, xlim=None, ylim=None,
                  target=None,
                  loc=111, title=None,
                  add_cbar=True, clabel=None,
                  show_pts=False, dolog=False, add_body=True,
                  figsize=(8.34,7),
                  *args, **kwargs):
        '''
        Create a plot of variable *var* at a given time.

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import (LogNorm, Normalize)
        from matplotlib.patches import Circle
        from matplotlib.ticker import (LogLocator, LogFormatter,
                                       LogFormatterMathtext, MultipleLocator)
        from matplotlib.patches import Wedge

        # Grab the slice of data that we want:
        if type(time) == type(self['time'][0]):
            if time not in self['time']:
                raise ValueError('Time not in object')
            time = np.arange(self.attrs['nTime'])[self['time']==time][0]

        fig, ax = set_target(target, loc=loc, figsize=figsize)
        ax.set_aspect('equal')

        # Grab values from correct time/location.
        x = self['x'][time, :]
        y = self['y'][time, :]
        value = self[var][time, :]

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=value.min(); zlim[1]=value.max()
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        #add points at "infinity" with values set zmin so that color
        #background fills in smoothly
        x=x.append(x,-15.0)
        x=x.append(x,-15.0)
        x=x.append(x,15.0)
        x=x.append(x,15.0)

        y=y.append(y,-15.0)
        y=y.append(y,15.0)
        y=y.append(y,-15.0)
        y=y.append(y,15.0)

        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z=np.where(value>zlim[0], value, 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z=value
            norm=None
            ticks=None
            fmt=None

        # Create contour plot.
        cont=ax.tricontourf(np.asarray(x), np.asarray(y), np.asarray(z), \
                            np.asarray(levs), *args, norm=norm, **kwargs)
        if show_pts:
            ax.plot(x, y, '+w')

        if xlim != None:
            ax.set_xlim(xlim[0],xlim[1])

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])

        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel==None:
                clabel="%s" % (var)
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.

        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        #ax.set_yticks([]), ax.set_xticks([])
        ax.set_xlabel('X [Re]')
        ax.set_ylabel('Y [Re]')

        if add_body:
            angle=-90.0
            radius=1.0
            colors=('w','k')
            theta1, theta2 = angle, angle + 180
            center = (0,0)
            w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], ec='k', **kwargs)
            w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], ec='k', **kwargs)
            for wedge in [w1,w2]:
                ax.add_artist(wedge)


        return fig, ax, cont, cbar

    def interpolate(self, var, xpt, ypt, time):
        '''
        extract solution at given point and time
        '''
        from scipy.interpolate import griddata

        #find closest time to do the interpolation
        dt=self['time']-time
        iTime=np.where(abs(dt) == np.min(abs(dt))) [0][0]

        #get number of points
        nPoints=len(self['x'][iTime,:])

        #create point pairs
        Points=np.zeros([nPoints,2])
        for i in range(nPoints):
            Points[i][0]=self['x'][iTime,i]
            Points[i][1]=self['y'][iTime,i]

        Values=self[var][iTime,:]

        InterpValue = griddata(Points, Values, (xpt, ypt),
                               method='linear',fill_value=0.0)

        return InterpValue

class CimiLog(PbData):
    '''
    A class for reading and plotting a complete set of CIMI log files.
    Open using a glob string that encompasses all of the log files that are
    intended to be read, e.g. 'CIMI*.log'.  Note that unix wildcards are
    accepted.
    '''

    def __init__(self, filepattern, starttime=None, *args, **kwargs):
        super(CimiLog, self).__init__(*args, **kwargs) # Init as PbData.
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(filepattern, starttime)

    def __repr__(self):
        return 'CIMI log files:', self['files']

    def _read(self, filepattern, starttime):
        '''
        Read all ascii line files; should only be called upon instantiation.
        '''

        from glob import glob

        self['files']=sorted(glob(filepattern))



        #get total length of combined data
        nDataLength = 0
        for file in self['files']:
            # Reads file assuming it is ascii
            # Slurp whole file.
            f=open(file, 'r')
            lines=f.readlines()
            f.close()
            nDataLength = nDataLength+len(lines)-2


        iCounter = 0
        offset = 0
        for file in self['files']:
            # Reads file assuming it is ascii
            # Slurp whole file.
            f=open(file, 'r')
            lines=f.readlines()
            f.close()

            if iCounter == 0:
                # Get variable names; pop g and rbody
                var=(lines[1].split())
                self._rawvar=var
                # initialize arrays to hold data
                for v in var:
                    self[v]=dmarray(np.zeros((nDataLength)))

            #read in data
            for j,l in enumerate(lines[2:]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][offset+j]=float(parts[k])

            iCounter=iCounter+1
            offset = offset +len(lines)-2

        # Start building time array.
        self['time']=np.zeros(nDataLength, dtype=object)

        for i in range(nDataLength):
            t= self['t'][i]
            self['time'][i]=starttime+dt.timedelta(seconds=t)

    def plotlog(self, vars, xlim=None, ylim=None,
                  target=None,
                  loc=111, title=None,
                  dolog=False,
                  figsize=(8.34,7),
                  *args, **kwargs):
        '''
        Create a plot of variable *var* at a given time.

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.
        '''

        import matplotlib.pyplot as plt


        fig, ax = set_target(target, loc=loc, figsize=figsize)
        #ax.set_aspect('equal')

        # Grab values from correct time/location.
        for var in vars:
            x = self['time'][:]
            y = self[var][:]
            ax.plot(x,y,label=var)


        # Create levels and set norm based on dolog.
        if dolog:
            plt.yscale('log')

        #add legend
        ax.legend()



        if xlim != None:
            ax.set_xlim(xlim[0],xlim[1])

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])



        return fig, ax
