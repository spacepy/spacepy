#!/usr/bin/env python

'''
The PyBats submodule for handling input and output for the Dynamic Global
Core Plasma Model (DGCPM), a plasmasphere module of the SWMF.
'''

import numpy as np
from spacepy.pybats import PbData

def _adjust_dialplot(ax, rad, title='12',labelsize=15):
    '''
    Ram output is often visualized with equatorial dial plots.  This
    function quickly adjusts those plots to give a uniform, clean
    look and feel.
    '''
    from numpy import max, pi
    from matplotlib.ticker import MultipleLocator

    from spacepy.pybats.ram import add_body_polar

    # Constrain range of plot:
    ax.set_ylim([0,max(rad)])

    # Set MLT labels:
    lt_labels = ['06', title, '18', '00']
    xticks    = [   0,   pi/2,   pi, 3*pi/2]
    ax.set_xticks(xticks)
    ax.set_xticklabels(lt_labels)
    ax.tick_params('x', labelsize=labelsize)

    # Set L labels and grid.  Turn off label at L=0.
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.figure.canvas.draw()
    labs = [item.get_text() for item in ax.get_yticklabels()]
    ax.set_yticklabels(labs, color='w', size=labelsize, backgroundcolor='k')
    labels=ax.get_yticklabels()
    labels[0].set_visible(False)

    # Turn on grid.
    ax.grid(True, c='w', lw=1.5, ls=':')

    # Change background color so labels stand out.
    ax.set_axis_bgcolor('gray')
    add_body_polar(ax)

class PlasmaFile(PbData):
    '''
    DGCPM's plasma files contain a plethora of values from the model's 2D plane.
    This class produces objects that load the data and sort it into a PyBats
    data object.
    '''

    def __init__(self, filename, *args, **kwargs):
        '''
        Reads the data; sorts into arrays.
        '''

        import datetime as dt
        from spacepy.datamodel import dmarray

        super(PlasmaFile, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file'] = filename

        # Some visualization defaults.
        self._default_cmaps={'n':'Greens_r', 'pot':'bwr'}
        self._default_zlims={'n':[1, 1000], 'pot':[-150,150]}

        # With no real header, these files require sufficient a priori
        # knowledge of the contents.  Let's hammer that out here.
        # Order is non-negotiable.
        # Some units are not what are actually in the file; see conversions
        # below (i.e. density defaults to m^-3, but that's unweildy.)
        varlist = [('n','cm^{-3}'), ('x','R_E'), ('y','R_E'), 
                   ('b_open','boolean'),
                   ('pot','kV',), ('corot','kV'), ('vr','m/s'), 
                   ('vphi','m/s'), ('source','m^-3'),('fluxphi','m^{-3}s^{-1}'),
                   ('fluxr','m^{-3}s^{-1}'), ('totaln','#'), ('vol','km^3')]

        # Read the header of the file.
        infile = open(filename, 'r')
        parts = infile.readline().split()
        nLat, nLon = int(parts[0]), int(parts[1])

        # Create arrays and fill 'em up.
        # First couple are 1D...
        self['lat'] = dmarray(np.zeros(nLat), {'units':'degrees'})
        self['lat'][:] = infile.readline().split()
        self['lon'] = dmarray(np.zeros(nLon), {'units':'degrees'})
        self['lon'][:] = infile.readline().split()

        # Rest are 2D.
        for v in varlist:
            self[v[0]] = dmarray(
                np.array(infile.readline().split(), dtype=float).reshape(
                    (nLat,nLon), order='F'), {'units':v[1]})

        # Final variables are time-related.
        self.attrs['time'] = dt.datetime(1965,1,1,0,0,0) + \
            dt.timedelta(seconds=float(infile.readline()))

        # That's it for our file.
        infile.close()

        # Calculate L-Shell for convenience.
        self['L'] = dmarray(np.cos(self['lat']*np.pi/180.0)**-2,{'units':'R_E'})

        # Sometimes, SI units aren't the best.
        self['n']    /= 100**3
        self['pot']  /= 1000.0
        self['corot']/= 1000.0

        # Some other useful attributes.
        self.attrs['dLat'] = self['lat'][1] - self['lat'][0]
        self.attrs['dLon'] = self['lon'][1] - self['lon'][0]
        self.attrs['dL']   = self['L'][1]   - self['L'][0]


    def add_pcolor(self, var, zlim=None, target=None, loc=111, title=None,
                   Lmax=None, add_cbar=False, clabel=None, dolog=False, 
                   **kwargs):
        
        from numpy import linspace, pi
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc, polar=True)
            doAdjust = True
        elif type(target).__base__ == plt.Axes:
            ax  = target
            fig = ax.figure
            doAdjust = False
        elif type(target).__base__ == plt.PolarAxes:
            ax  = target
            fig = ax.figure
            doAdjust = False
        else:
            fig = plt.figure(figsize=(10.5,8))
            ax  = fig.add_subplot(loc, polar=True)
            doAdjust = True

        # Get max/min if none given.
        if zlim==None:
            if var in self._default_zlims:
                zlim = self._default_zlims[var]
            else:
                zlim = [self[var].min(), self[var].max()]
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Logarithmic scale?
        if dolog:
            z=np.where(self[var]>zlim[0], self[var], 1.01*zlim[0])
            norm=LogNorm()
        else:
            z=self[var]
            norm=None

        # Set up proper meshes. 
        nL, nPhi = len(self['L']), len(self['lon'])
        dL, dPhi = self.attrs['dL'], self.attrs['dLon']*pi/180.0
        phi = linspace(-1.0*dPhi/2.0, 2.*pi-dPhi/2.0, nPhi+1)-pi/2.
        r   = linspace(self['L'][0]-dL/2.0, self['L'][-1]+dL/2.0, nL+1)

        # Set default color tables based on variable plotted.
        if ('cmap' not in kwargs) and var in self._default_cmaps:
            kwargs['cmap'] = self._default_cmaps[var]

        # -------Plot away------
        # Mesh:
        pcol = ax.pcolormesh(phi, r, z, vmin=zlim[0], vmax=zlim[1],
                             norm=norm, **kwargs)
        # Add cbar if requested:
        if add_cbar:
            cbar=plt.colorbar(pcol, pad=0.08, shrink=.8)
            if clabel==None: 
                clabel="%s ($%s$)" % (var, self[var].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None
        
        # Adjust plot appropriately.
        if doAdjust:
            if not Lmax:
                # Default to inside ghost cells.
                Lmax = self['L'][-3]
                if title:
                    ax.set_title(title+'\n'+self.attrs['time'].isoformat(), 
                                 position=(0,1), ha='left', size=14)
            _adjust_dialplot(ax, Lmax, labelsize=14)

        return fig, ax, pcol, cbar


    def add_contour(self, var, zlim=None, target=None, loc=111, title=None,
                    Lmax=None, add_cbar=False, clabel=None, dolog=False, 
                    filled=True, nLev=31, **kwargs):
        
        from numpy import linspace, pi
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import MultipleLocator, LogLocator

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc, polar=True)
            doAdjust = True
        elif type(target).__base__ == plt.Axes:
            ax  = target
            fig = ax.figure
            doAdjust = False
        elif type(target).__base__ == plt.PolarAxes:
            ax  = target
            fig = ax.figure
            doAdjust = False
        else:
            fig = plt.figure(figsize=(10.5,8))
            ax  = fig.add_subplot(loc, polar=True)
            doAdjust = True

        # Set function based on boolean "filled":
        if filled:
            contour = ax.contourf
        else:
            contour = ax.contour

        # Get max/min if none given.
        if zlim==None:
            if var in self._default_zlims:
                zlim = self._default_zlims[var]
            else:
                zlim = [self[var].min(), self[var].max()]
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]), 
                                            np.log10(zlim[1]), nLev))
            z=np.where(self[var]>zlim[0], self[var], 1.01*zlim[0])
            norm=LogNorm()
        else:
            levs = np.linspace(zlim[0], zlim[1], nLev)
            z=self[var]
            norm=None

        # Allow results to cross phi=360.
        phi = np.concatenate( (self['lon'], [360.0]) ) * np.pi/180. - np.pi/2.
        z = np.concatenate( (z, np.array([z[:,0]]).transpose()  ), 1)

        # Set default color tables based on variable plotted.
        if ('cmap' not in kwargs) and ('colors' not in kwargs) and \
                (var in self._default_cmaps):
            kwargs['cmap'] = self._default_cmaps[var]

        # -------Plot away------
        # Contour:
        cont = contour(phi, self['L'], z, levs, norm=norm, **kwargs)
        # Add cbar if requested:
        if add_cbar:
            cbar=plt.colorbar(cont, pad=0.08, shrink=.8)
            if clabel==None: 
                clabel="%s ($%s$)" % (var, self[var].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None
        
        # Adjust plot appropriately.
        if doAdjust:
            if not Lmax:
                # Default to inside ghost cells.
                Lmax = self['L'][-3]
                if title:
                    ax.set_title(title+'\n'+self.attrs['time'].isoformat(), 
                                 position=(0,1), ha='left', size=14)
            _adjust_dialplot(ax, Lmax, labelsize=14)

        return fig, ax, cont, cbar
