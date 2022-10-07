#!/usr/bin/env python

'''
The PyBats submodule for handling input and output for the Dynamic Global
Core Plasma Model (DGCPM), a plasmasphere module of the SWMF.
'''

import numpy as np
from spacepy.plot import applySmartTimeTicks, set_target
from spacepy.pybats import PbData


def _adjust_dialplot(ax, rad, title='12', labelsize=15):
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
    add_body_polar(ax)


def saturation(L):
    '''
    Return saturation density of a flux tube as a function of L-shell in units
    of :math:`#/cm^3`. Formula is from Carpenter and Anderson, JGR, 1992.
    '''

    return 10**(-0.3145*L + 3.9043)


def refill_flux(n, L, tau=1.5, chi='Auto'):
    '''
    Calculate the refill flux of a flux tube of density *n*, L-shell *L*,
    using a refill time constant of *tau*.  If kwarg *chi* is set to 'Auto',
    then zenith angle dependence is calculated using an assumed dipole field.
    To remove this dependence, set *chi* to a solar zenith angle in degrees.
    '''

    if chi == 'Auto':
        chi = np.pi/2. - np.arccos(1./np.sqrt(L))
    else:
        chi *= np.pi/180.0

    # Convenient values:
    vmin = 1879379699284659.0  # Taken right from pbo.f!
    Lmin = 1.29762044

    # Start with fmax.
    fmax = 2.*vmin*Lmin*saturation(Lmin)/(tau*24.*3600.0)

    # Calculate refill rate.
    nsat = saturation(L)
    flux = (nsat-n)/nsat * fmax * np.sin(chi) * L**-3

    return flux


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
        import gzip
        from spacepy.datamodel import dmarray

        super(PlasmaFile, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file'] = filename

        # Some visualization defaults.
        self._default_cmaps = {'n':'Greens_r', 'pot':'bwr'}
        self._default_zlims = {'n':[1, 1000], 'pot':[-150, 150]}

        # Set units and variable names.
        # Some units are not what are actually in the file; see conversions
        # below (i.e. density defaults to m^-3, but that's unweildy.)
        units = [('n', 'cm^{-3}'), ('x', 'R_E'), ('y', 'R_E'),
                 ('b_open', 'boolean'), ('theta', 'degrees'),
                 ('phi', 'degrees'), ('pot', 'kV',), ('corot', 'kV'),
                 ('vr', 'm/s'), ('vphi', 'm/s'), ('source', 'm^-3'),
                 ('fluxphi', 'm^{-3}s^{-1}'), ('fluxr', 'm^{-3}s^{-1}'),
                 ('totaln', '#'), ('vol', 'km^3')]
        units = dict(units)

        # Open file; unzip if necessary.
        if self.attrs['file'][-3:] == '.gz':
            infile = gzip.open(self.attrs['file'], 'rt')
        else:
            infile = open(self.attrs['file'], 'r')

        # Read the header of the file.
        parts = infile.readline().split()
        nLat, nLon = int(parts[-2]), int(parts[-1])
        self.attrs['time'] = dt.datetime.strptime(parts[0],
                                                  'T=%Y%m%d_%H%M%S_000')
        varlist = infile.readline().lower().split()[2:]

        # Create containers for data:
        for v in varlist:
            unit = '' if v not in units else units[v]
            self[v] = dmarray(np.zeros((nLat, nLon)), {'units': unit})

        # Read rest of file and sort data into arrays:
        for line in infile.readlines():
            parts = line.split()
            i, j = int(parts.pop(0))-1, int(parts.pop(0))-1
            for v, x in zip(varlist, parts):
                self[v][i, j] = x

        # That's it for our file.
        infile.close()

        # Create some "helper" variables:
        self['lat'] = dmarray(self['theta'][:, 0], {'units':'degrees'})
        self['lon'] = dmarray(self['phi'][0, :], {'units':'degrees'})
        self['mlt'] = dmarray(self['lon']/15.0, {'units':'hours'})

        # Calculate L-Shell for convenience.
        self['L'] = dmarray(np.cos(self['lat']*np.pi/180.)**-2,
                            {'units':'R_E'})

        # Sometimes, SI units aren't the best.
        if 'n' in self:
            self['n'] /= 100.0**3
        if 'pot' in self:
            self['pot'] /= 1000.0
        if 'corot' in self:
            self['corot'] /= 1000.0

        # Some other useful attributes.
        self.attrs['dLat'] = self['lat'][1] - self['lat'][0]
        self.attrs['dLon'] = self['lon'][1] - self['lon'][0]
        self.attrs['dL'] = self['L'][1] - self['L'][0]

    def calc_E(self):
        '''
        Differentiate the equatorial electric potential to arrive at electric field
        (stored using keys 'Er', 'Ephi', and 'E').
        '''

        from spacepy.datamodel import dmarray
        from spacepy.pybats.batsmath import d_dx, d_dy

        conv = 1.0E6/6371000.0  # Unit conversion: kV -> mV, 1/Re -> 1/m.
        x, L = np.meshgrid(np.zeros(self['n'].shape[-1]), np.array(self['L']))

        Er = d_dy(self['pot']*conv,   self.attrs['dL'])
        Ephi = d_dx(self['pot']*conv/L, self.attrs['dLon']*np.pi/180.0)
        E = np.sqrt(Er**2 + Ephi**2)

        self['E'] = dmarray(E, {'units':'mV/m'})
        self['Er'] = dmarray(Er, {'units':'mV/m'})
        self['Ephi'] = dmarray(Ephi, {'units':'mV/m'})

    def add_pcolor(self, var, zlim=None, target=None, loc=111, title=None,
                   Lmax=None, add_cbar=False, clabel=None, dolog=False,
                   labelsize=14, **kwargs):
        '''
        Create a polar pcolor plot of variable *var* and add it to *target*.
        Extra keyword arguments are handed to matplotlib's pcolor command.

        Parameters
        ==========
        var : string
             The variable within the object to plot.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object
        cont : matplotlib contour object
        cbar : matplotlib colorbar object

        Other Parameters
        ================
        zlim : two-element list or array
             Set the color bar range.  Some variables have default ranges.
             Others will use max and min of *var*.
        Lmax : real
             Set the radial extent of the plot.
        add_cbar : bool
             Set whether to add a color bar or not.  Defaults to **False**.
        dolog : bool
             If **True**, use a log scale to plot *var*.
             Defaults to **False**.
        labelsize : int
             Sets the font size of the labels.  Defaults to 14.
        title : string
             Sets the plot title.
             Defaults to 'auto', using the variable label.
        clabel : string
             Set label for the color bar.
             Defaults to *var* and associated units.
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

        '''

        from numpy import linspace, pi
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10.5, 8), loc=loc, polar=True)

        # Get max/min if none given.
        if zlim is None:
            if var in self._default_zlims:
                zlim = self._default_zlims[var]
            else:
                zlim = [self[var].min(), self[var].max()]
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Logarithmic scale?
        if dolog:
            z = np.where(self[var] > zlim[0], self[var], 1.01*zlim[0])
            norm = LogNorm()
        else:
            z = self[var]
            norm = None

        # Set up proper meshes.
        nL, nPhi = len(self['L']), len(self['lon'])
        dL, dPhi = self.attrs['dL'], self.attrs['dLon']*pi/180.0
        phi = linspace(-1.0*dPhi/2.0, 2.*pi-dPhi/2.0, nPhi+1)-pi/2.
        r = linspace(self['L'][0]-dL/2.0, self['L'][-1]+dL/2.0, nL+1)

        # Set default color tables based on variable plotted.
        if ('cmap' not in kwargs) and var in self._default_cmaps:
            kwargs['cmap'] = self._default_cmaps[var]

        # -------Plot away------
        # Mesh:
        pcol = ax.pcolormesh(phi, r, z, vmin=zlim[0], vmax=zlim[1],
                             norm=norm, **kwargs)
        # Add cbar if requested:
        if add_cbar:
            cbar = plt.colorbar(pcol, pad=0.08, shrink=.8, ax=ax)
            if clabel is None:
                clabel = "%s ($%s$)" % (var, self[var].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar = None

        # Adjust plot appropriately.
        if not Lmax:
            # Default to inside ghost cells.
            Lmax = self['L'][-3]
            if title:
                ax.set_title(title+'\n'+self.attrs['time'].isoformat(),
                             position=(0, 1), ha='left', size=14)
        _adjust_dialplot(ax, Lmax, labelsize=labelsize)

        return fig, ax, pcol, cbar

    def add_contour(self, var, zlim=None, target=None, loc=111, title=None,
                    Lmax=None, add_cbar=False, clabel=None, dolog=False,
                    filled=True, nLev=31, labelsize=14, **kwargs):
        '''
        Create a polar contour plot of variable *var* and add it to *target*.
        Extra keyword arguments are handed to matplotlib's contourf command.

        Parameters
        ==========
        var : string
             The variable within the object to plot.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object
        cont : matplotlib contour object
        cbar : matplotlib colorbar object

        Other Parameters
        ================
        nLev : int
             Sets the number of contour levels.  Defaults to 31.
        zlim : two-element list or array
             Set the color bar range.  Some variables have default ranges.  Others
             will use max and min of *var*.
        Lmax : real
             Set the radial extent of the plot.
        add_cbar : bool
             Set whether to add a color bar or not.  Defaults to **False**.
        dolog : bool
             If **True**, use a log scale to plot *var*.  Defaults to **False**.
        labelsize : int
             Sets the font size of the labels.  Defaults to 14.
        title : string
             Sets the plot title.  Defaults to 'auto', using the variable label.
        clabel : string
             Set label for the color bar.  Defaults to *var* and associated units.
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

        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import LogLocator

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10.5, 8), loc=loc, polar=True)

        # Set function based on boolean "filled":
        if filled:
            contour = ax.contourf
        else:
            contour = ax.contour

        # Get max/min if none given.
        if zlim is None:
            if var in self._default_zlims:
                zlim = self._default_zlims[var]
            else:
                zlim = [self[var].min(), self[var].max()]
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nLev))
            z = np.where(self[var] > zlim[0], self[var], 1.01*zlim[0])
            z[z > zlim[-1]] = zlim[-1]
            norm = LogNorm()
            lct = LogLocator()
        else:
            levs = np.linspace(zlim[0], zlim[1], nLev)
            z = self[var]
            norm = None
            lct = None

        # Allow results to cross phi=360.
        phi = np.concatenate((self['lon'], [360.0])) * np.pi/180. - np.pi/2.
        z = np.concatenate((z, np.array([z[:, 0]]).transpose()), 1)

        # Set default color tables based on variable plotted.
        if ('cmap' not in kwargs) and ('colors' not in kwargs) and \
                (var in self._default_cmaps):
            kwargs['cmap'] = self._default_cmaps[var]

        # -------Plot away------
        # Contour:
        cont = contour(phi, self['L'], z, levs, norm=norm, **kwargs)
        # Add cbar if requested:
        if add_cbar:
            cbar = plt.colorbar(cont, pad=0.08, shrink=.8, ticks=lct, ax=ax)
            if clabel is None:
                clabel = f"{var} (${self[var].attrs['units']}$)"
            cbar.set_label(clabel)
        else:
            cbar = None

        # Adjust plot appropriately.
        if not Lmax:
            # Default to inside ghost cells.
            Lmax = self['L'][-3]
            if title:
                ax.set_title(title+'\n'+self.attrs['time'].isoformat(),
                             position=(0, 1), ha='left', size=14)
        _adjust_dialplot(ax, Lmax, labelsize=labelsize)

        return fig, ax, cont, cbar

    def add_separatrix(self, target=None, loc=111, Lmax=None,
                       title=None, **kwargs):
        '''
        Attempts to locate the separatrix (separator between closed and
        open drift paths) by finding the minimum velocity in the domain and
        tracing the path of constant potential that passes through that point.
        The figure, axes, and contour object containing the separatrix line
        are returned to the user.

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        Kwargs that set line characteristics behave in the typical matplotlib
        manner (i.e. "colors" can be set to either a color name or a
        hexidecimal specifier.)

        =========== ==========================================================
        Kwarg       Description
        ----------- ----------------------------------------------------------
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        linewidths  Set width of plotted line.  Defaults to 3.0.
        colors      Set color of line.  Defaults to 'orange'.
        linestyles  Set line style.  Defaults to 'dashed'.
        ----------- ----------------------------------------------------------

        '''

        from numpy import pi

        # Set default line behavior.
        if 'linewidths' not in kwargs:
            kwargs['linewidths'] = 3.0
        if 'colors' not in kwargs:
            kwargs['colors'] = 'orange'
        if 'linestyles' not in kwargs:
            kwargs['linestyles'] = 'solid'

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10.5, 8), loc=loc, polar=True)
        doAdjust = not target == ax

        # Find the stagnation point by looking for E = 0.
        # Get the potential value at that spot.
        # Exclude near-boundary points and limit azimuthal extent.
        if 'E' not in self:
            self.calc_E()
        i = np.arange(self['L'].size)[self['L'] < 9.75]
        j = np.arange(self['lon'].size)[(self['lon'] > 250.) &
                                        (self['lon'] < 360.)]
        pmin = self['pot'][self['E'] == self['E'][i[0]:i[-1],
                                                  j[0]:j[-1]].min()][0]

        # Create a contour that only follows that potential curve.
        phi = np.concatenate((self['lon'], [360.0])) * pi/180. - pi/2.
        z = np.concatenate((self['pot'],
                            np.array([self['pot'][:, 0]]).transpose()), 1)
        cnt = ax.contour(phi, self['L'], z, levels=[pmin], **kwargs)

        # Adjust plot appropriately.
        if doAdjust:
            if not Lmax:
                # Default to inside ghost cells.
                Lmax = self['L'][-3]
                if title:
                    ax.set_title(title+'\n'+self.attrs['time'].isoformat(),
                                 position=(0, 1), ha='left', size=14)
            _adjust_dialplot(ax, Lmax, labelsize=14)

        # Return bits to user.
        return fig, ax, cnt


class MltSlice(PbData):
    '''
    Open and handle an MLT Slice output file.
    '''

    def __init__(self, filename, *args, **kwargs):
        '''
        Reads the data; sorts into arrays.
        '''

        import datetime as dt
        from spacepy.datamodel import dmarray
        from matplotlib.dates import date2num

        super(MltSlice, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file'] = filename

        f = open(filename, 'r')

        # Parse header.
        self.attrs['mlt'] = float(f.readline().split()[-1])

        line = f.readline().split('=')[-1]
        self['L'] = dmarray(np.array(line.split(), dtype=float),
                            {'units':'$R_E$'})

        # Parse remainder of file.
        lines = f.readlines()
        self['n'] = dmarray(np.zeros([len(lines), len(self['L'])]),
                            {'units':'cm^{-3}'})
        self['time'] = dmarray(np.zeros(len(lines), dtype=object))

        for i, l in enumerate(lines):
            p = l.split()
            self['time'][i] = dt.datetime(int(p[0]), int(p[1]), int(p[2]),
                                          int(p[3]), int(p[4]), int(p[5]),
                                          int(p[6])*1000)
            self['n'][i, :] = p[7:]

        # Some "hidden" variables for plotting.
        self._dtime = date2num(self['time'])
        self._dy = self['L'][1] - self['L'][0]

    def add_lut(self, target=None, loc=111, cmap='Greens_r', zlim=[1, 1000],
                add_cbar=True, clabel='Density $cm^{-3}$', xlabel='full',
                title=None, grid=True, ntick=5):
        '''
        Plot log(density) as a contour against L-Shell (y-axis) and
        universal time (x-axis) using the PyBats *target* method of other
        standard plotting methods.  Four items are returned: the Matplotlib
        Figure, Axes, Mesh, and ColorBar objects used (if cbar is set to
        **False**, the returned ColorBar object is simply set to **False**.)

        ========== =======================================================
        Kwarg      Description
        ---------- -------------------------------------------------------
        target     Select plot destination.  Defaults to new figure/axis.
        loc        The location of any generated subplots.  Default is 111.
        add_cbar   Toggles the automatic colorbar.  Default is**True**.
        cmap       Selects Matplotlib color table.  Defaults to *Greens_r*.
        zlim       Limits for z-axis.  Defaults to [0.1, 1000]
        clabel     Sets colorbar label.  Defaults to units.
        xlabel     Sets x-axis labels, use 'full', 'ticks', or **None**.
        title      Sets axis title; defaults to **None**.
        grid       Show white dotted grid?  Defaults to **True**
        ntick      Number of attempted cbar ticks.  Defaults to 5.
        ========== =======================================================

        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Set ax and fig based on given target.
        fig, ax = set_target(target, loc=loc)

        # Enforce values to be within limits.
        z = np.where(self['n'] > zlim[0], self['n'], 1.01*zlim[0])
        z[z > zlim[1]] = zlim[1]

        # Create plot:
        mesh = ax.pcolormesh(self._dtime, self['L'], z.transpose(),
                             cmap=plt.get_cmap(cmap), norm=LogNorm(),
                             vmin=zlim[0], vmax=zlim[-1])

        # Use LT ticks and markers on y-axis:
        ax.set_ylabel('L-Shell')

        # White ticks, slightly thicker:
        ax.tick_params(axis='both', which='both', color='w', width=1.2)

        # Grid marks:
        if grid:
            ax.grid(c='w')

        if title:
            ax.set_title(title)
        if xlabel == 'full':
            # Both ticks and label.
            applySmartTimeTicks(ax, self['time'], dolabel=True)
        elif xlabel == 'ticks':
            # Ticks, but no date label.
            applySmartTimeTicks(ax, self['time'], dolabel=False)
        else:
            # A blank x-axis is often useful.
            applySmartTimeTicks(ax, self['time'], dolabel=False)
            ax.set_xticklabels('')
        # Add cbar as necessary:
        if add_cbar:
            cbar = plt.colorbar(mesh, ax=ax, pad=0.01, shrink=0.85)
            cbar.set_label(clabel)
        else:
            cbar = None

        return fig, ax, mesh, cbar


class Lslice(PbData):
    '''
    Open an L-Slice output file.
    '''

    def __init__(self, filename, *args, **kwargs):
        '''
        Reads the data; sorts into arrays.
        '''

        import datetime as dt
        from spacepy.datamodel import dmarray
        from matplotlib.dates import date2num

        super(Lslice, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['file'] = filename

        f = open(filename, 'r')

        # Parse header.
        self.attrs['L'] = float(f.readline().split()[-1])
        self['mlt'] = dmarray(np.array(f.readline().split()[1:], dtype=float),
                              {'units':'Hours'})

        # Parse remainder of file.
        lines = f.readlines()
        self['n'] = dmarray(np.zeros([len(lines), len(self['mlt'])]),
                            {'units':'cm^{-3}'})
        self['time'] = dmarray(np.zeros(len(lines), dtype=object))

        for i, l in enumerate(lines):
            p = l.split()
            self['time'][i] = dt.datetime(int(p[0]), int(p[1]), int(p[2]),
                                          int(p[3]), int(p[4]), int(p[5]),
                                          int(p[6])*1000)
            self['n'][i, :] = p[7:]

        # Some "hidden" variables for plotting.
        self._dtime = date2num(self['time'])
        self._dy = self['mlt'][1] - self['mlt'][0]
        self._y = np.arange(0, self['mlt'][-1]+2*self._dy, self._dy)

    def add_ltut(self, target=None, loc=111, cmap='Greens_r', zlim=[1, 1000],
                 add_cbar=True, clabel='Density $cm^{-3}$', xlabel='full',
                 ylim=[4, 20], title=None, grid=True, ntick=5):
        '''
        Plot log(density) as a contour against local time (y-axis) and
        universal time (x-axis) using the PyBats *target* method of other
        standard plotting methods.  Four items are returned: the Matplotlib
        Figure, Axes, Mesh, and ColorBar objects used (if cbar is set to
        **False**, the returned ColorBar object is simply set to **False**.)

        ========== =======================================================
        Kwarg      Description
        ---------- -------------------------------------------------------
        target     Select plot destination.  Defaults to new figure/axis.
        loc        The location of any generated subplots.  Default is 111.
        add_cbar   Toggles the automatic colorbar.  Default is**True**.
        cmap       Selects Matplotlib color table.  Defaults to *Greens_r*.
        zlim       Limits for z-axis.  Defaults to [0.1, 1000]
        ylim       Sets the MLT range on the y-axis.  Defaults to [4,20].
        clabel     Sets colorbar label.  Defaults to units.
        xlabel     Sets x-axis labels, use 'full', 'ticks', or **None**.
        title      Sets axis title; defaults to **None**.
        grid       Show white dotted grid?  Defaults to **True**
        ntick      Number of attempted cbar ticks.  Defaults to 5.
        ========== =======================================================

        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Set ax and fig based on given target.
        fig, ax = set_target(target, loc=loc)
        # Enforce values to be within limits.
        z = np.where(self['n'] > zlim[0], self['n'], 1.01*zlim[0])
        z[z > zlim[1]] = zlim[1]

        # Create plot:
        mesh = ax.pcolormesh(self._dtime, self._y, z.transpose(),
                             cmap=plt.get_cmap(cmap), norm=LogNorm(),
                             vmin=zlim[0], vmax=zlim[-1])

        # Use LT ticks and markers on y-axis:
        ax.set_yticks([6, 12, 18])
        ax.set_yticklabels(['Dawn', 'Noon', 'Dusk'])
        ax.set_ylim(ylim)

        # White ticks, slightly thicker:
        ax.tick_params(axis='both', which='both', color='w', width=1.2)

        # Grid marks:
        if grid:
            ax.grid(c='w')

        if title:
            ax.set_title(title)
        if xlabel == 'full':
            # Both ticks and label.
            applySmartTimeTicks(ax, self['time'], dolabel=True)
        elif xlabel == 'ticks':
            # Ticks, but no date label.
            applySmartTimeTicks(ax, self['time'], dolabel=False)
        else:
            # A blank x-axis is often useful.
            applySmartTimeTicks(ax, self['time'], dolabel=False)
            ax.set_xticklabels('')
        # Add cbar as necessary:
        if add_cbar:
            cbar = plt.colorbar(mesh, ax=ax, pad=0.01, shrink=0.85)
            cbar.set_label(clabel)
        else:
            cbar = None

        return fig, ax, mesh, cbar
