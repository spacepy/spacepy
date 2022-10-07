#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Polar Wind Outflow Model
(PWOM), one of the choices for the PW module in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
from spacepy.plot import set_target
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

class Efile(PbData):
    '''
    Class for loading and parsing ionospheric electrodynamics ascii files.
    
    '''
    
    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(Efile, self).__init__(*args, **kwargs)
        self.attrs['file']=filename
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0)
        self._read(filename, starttime)

    def _read(self, filename, starttime):
        '''
        Read, parse, and load data into data structure.
        Should only be called by __init__.
        '''
        
        from re import match, search
        from spacepy.pybats import parse_tecvars

        # Save time (starttime + run time)
        # Run time must be found in file name.
        m=search('Time(\d+)\.dat', filename)
        if m:
            self.attrs['runtime']=int(m.group(1))
        else:
            self.attrs['runtime']=0
        self.attrs['time']=starttime + dt.timedelta(
            seconds=self.attrs['runtime'])

        # Open file
        f=open(filename, 'r')
        
        # Parse header.
        varlist = parse_tecvars(f.readline())
        m = search('.*I\=\s*(\d+)\,\s*J\=\s*(\d+)',f.readline()).groups()
        nLat, nLon = int(m[0]), int(m[1])
        self.attrs['nLat']=nLat
        self.attrs['nLon']=nLon

        # Create arrays.
        for k, u in varlist:
            self[k]=dmarray(np.zeros( (nLat*nLon) ), attrs={'units':u})

        for i in range(nLat*nLon):
            #for j in range(nLon):
            parts=f.readline().split()
            for k, (key, u) in enumerate(varlist):
                self[key][i]=float(parts[k])

        # Lat, CoLat, and Lon:
        #self['lon'] = dmarray(np.pi/2 - np.arctan(self['y']/self['x']),
        #                      attrs={'units','deg'})*180.0/np.pi
        #self['lon'][self['x']==0.0] = 0.0
        xy=np.sqrt(self['x']**2+self['y']**2)+0.0000001
        self['lon']=np.arcsin(self['y']/xy)
        self['lat']  =dmarray(np.arcsin(self['z']), {'units':'deg'})*180.0/np.pi
        self['colat']=90.0 - self['lat']

        f.close()

class Line(PbData):
    '''
    Class for loading a single field line output file.
    At instantiation time, user may wish to set the start date and time of 
    the simulation using the starttime kwarg.  If not given, start time
    will default to Jan. 1st, 2000, 00:00UT.
    '''
    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(Line, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(starttime)

    def __repr__(self):
        return 'PWOM single field line output file %s' % (self.attrs['file'])

    def _read(self, starttime):
        '''
        Read ascii line file; should only be called upon instantiation.
        '''
        
        # Slurp whole file.
        f=open(self.attrs['file'], 'r')
        lines=f.readlines()
        f.close()

        # Determine size of file.
        nTimes=lines.count(lines[0])
        nAlts =int(lines[2].strip())
        self.attrs['nAlt']=nAlts; self.attrs['nTime']=nTimes

        # Start building time array.
        self['time']=np.zeros(nTimes, dtype=object)

        # Get variable names; pop radius (altitude).
        var=(lines[4].split())[1:-1]
        self._rawvar=var

        # Get altitude, which is constant at all times.
        self['r']=dmarray(np.zeros(nAlts),{'units':'km'})
        for i, l in enumerate(lines[5:nAlts+5]):
            self['r'][i]=float(l.split()[0])

        # Create 2D arrays for data that is time and alt. dependent.
        # Set as many units as possible.
        for v in var:
            self[v]=dmarray(np.zeros((nTimes, nAlts)))
            if v=='Lat' or v=='Lon':
                self[v].attrs['units']='deg'
            elif v[0]=='u':
                self[v].attrs['units']='km/s'
            elif v[0:3]=='lgn':
                self[v].attrs['units']='log(cm-3)'
            elif v[0]=='T':
                self[v].attrs['units']='K'
            else:
                self[v].attrs['units']=None

        # Loop through rest of data to fill arrays.
        for i in range(nTimes):
            t=float((lines[i*(nAlts+5)+1].split())[1])
            self['time'][i]=starttime+dt.timedelta(seconds=t)
            for j, l in enumerate(lines[i*(nAlts+5)+5:(i+1)*(nAlts+5)]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][i,j]=float(parts[k+1])

class Lines(PbData):
    '''
    A class for reading and plotting a complete set of PWOM line output files.
    Open using a glob string that encompasses all of the lines that are
    intended to be read, e.g. 'north*.out'.  Note that unix wildcards are
    accepted.
    '''

    def __init__(self, lines, starttime=None, *args, **kwargs):
        super(Lines, self).__init__(*args, **kwargs) # Init as PbData.
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(lines, starttime)
        self.calc_flux()

    def __repr__(self):
        return 'PWOM field line group of %i separate line files.' %\
            (self.attrs['nFile'])

    def _read(self, lines, starttime):
        '''
        Read all ascii line files; should only be called upon instantiation.
        '''

        from glob import glob
        
        self['files']=glob(lines)
        
        # Read first file; use it to set up all arrays.
        l1=Line(self['files'][0], starttime=starttime)
        nA=l1.attrs['nAlt']
        nT=l1.attrs['nTime']
        nF=len(self['files'])
        self['r'] = l1['r']
        
        self.attrs['nAlt']=nA; self.attrs['nTime']=nT; self.attrs['nFile']=nF
        self['time'] = l1['time']

        # Create all arrays.
        for v in l1._rawvar:
            self[v]=dmarray(np.zeros((nF, nT, nA)))
            self[v][0,:,:]=l1[v]
            if v=='Lat' or v=='Lon':
                self[v].attrs['units']='deg'
            elif v[0]=='u':
                self[v].attrs['units']='km/s'
            elif v[0:3]=='lgn':
                self[v].attrs['units']='log(cm-3)'
            elif v[0]=='T':
                self[v].attrs['units']='K'
            else:
                self[v].attrs['units']=None

        # Place data into arrays.
        for i, f in enumerate(self['files'][1:]):
            l=Line(f)
            for v in l1._rawvar:
                self[v][i+1,:,:]=l[v]

        # Get non-log densities.
        logVars = []
        for v in self:
            if v[:3] == 'lgn': logVars.append(v)
        for v in logVars:
            self[v[2:]] = dmarray(10.**self[v], {'units':r'$cm^{-3}$'})


    def calc_flux(self):
        '''
        Calculate flux in units of #/cm2/s.  Variables saved as self[species+'Flux'].
        '''

        for s in ['H', 'O', 'He', 'e']:
            self[s+'Flux']=dmarray(1000*self['n'+s]*self['u'+s],  {'units':'$cm^{-2}s^{-1}$'})

    def _get_cartXY(self):
        '''
        For plotting, generate a set of x-y values such that plots can be
        produced on a cartesian grid using tricontourf.
        '''
        
        # Pass if values already exist.
        if hasattr(self, '_xGSM'):
            return

        # Unit conversions.
        r = self['r']/6371.0 + 1.0
        lat   = np.pi*self['Lat']/180.0
        colat = 90. - self['Lat']
        lon   = np.pi*self['Lon']/180.0

        # Values in a Cartesian plane.
        self._xGSM = r*np.cos(lat)*np.sin(lon)
        self._yGSM = r*np.cos(lat)*np.cos(lon)

        self._xLat = colat*np.sin(lon)
        self._yLat = colat*np.cos(lon)
        
    def add_slice(self, var, alt, time, nlev=31, zlim=None, target=None, 
                  loc=111, title=None, latoffset=1.05,
                  rlim=50., add_cbar=True, clabel=None,
                  show_pts=False, show_alt=True, dolog=False, 
                  lats=[75., 60.], colats=None, figsize=(8.34,7),
                  *args, **kwargs):
        '''
        Create a plot of variable *var* at altitude *alt*.

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

        # Grab the slice of data that we want:
        if type(time) == type(self['time'][0]):
            if time not in self['time']:
                raise ValueError('Time not in object')
            time = np.arange(self.attrs['nTime'])[self['time']==time][0]

        fig, ax = set_target(target, loc=loc, figsize=figsize)
        ax.set_aspect('equal')

        # Create values in Cartesian plane.
        self._get_cartXY()

        # Grab values from correct time/location.
        x = self._xLat[:,time, alt]
        y = self._yLat[:,time, alt]
        value = self[var][:, time, alt]

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=value.min(); zlim[1]=value.max()
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

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
        cont=ax.tricontourf(x, y, z, levs, *args, norm=norm, **kwargs)
        
        # Label altitude.
        if show_alt:
            ax.text(0.05, 0.05, r'Alt.={:.1f}$km/s$ ({:.2f}$R_E$)'.format(
                self['r'][alt], self['r'][alt]/6371.0+1), size=12, 
                    transform=ax.transAxes, bbox={'fc':'white', 'ec':'k'})

        if show_pts:
            ax.plot(x, y, '+w')

        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel==None: 
                clabel="%s (%s)" % (var, self[var].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        ax.set_yticks([]), ax.set_xticks([])
        ax.set_ylabel(r'Sun $\rightarrow$')
        colat_lim = 90.-rlim
        ax.set_xlim([-1*colat_lim, colat_lim])
        ax.set_ylim([-1*colat_lim, colat_lim])

        if colats:
            for colat in colats:
                r = latoffset*colat/np.sqrt(2)
                ax.add_patch(Circle([0,0],colat,fc='none',ec='k',ls='dashed'))
                ax.text(r, r, '{:.0f}'.format(colat)+r'$^{\circ}$',
                        rotation=315,ha='center', va='center')
        else:
            for lat in lats:
                colat = 90 - lat
                r = latoffset*colat/np.sqrt(2)
                ax.add_patch(Circle([0,0],colat,fc='none',ec='k',ls='dashed'))
                ax.text(r, r, '{:.0f}'.format(lat)+r'$^{\circ}$', 
                        rotation=315,ha='center', va='center')
                

        return fig, ax, cont, cbar
