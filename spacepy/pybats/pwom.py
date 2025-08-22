#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Polar Wind Outflow Model
(PWOM), one of the choices for the PW module in the SWMF.
'''

# Global imports:
import re

import numpy as np
import datetime as dt

from spacepy.plot import set_target
from spacepy.pybats import PbData, IdlFile
from spacepy.datamodel import dmarray


class Efile(PbData):
    '''
    Class for loading and parsing ionospheric electrodynamics ascii files.

    '''

    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(Efile, self).__init__(*args, **kwargs)
        self.attrs['file'] = filename
        if not starttime:
            starttime = dt.datetime(2000, 1, 1, 0, 0)
        self._read(filename, starttime)

    def _read(self, filename, starttime):
        '''
        Read, parse, and load data into data structure.
        Should only be called by __init__.
        '''

        from re import search
        from spacepy.pybats import parse_tecvars

        # Save time (starttime + run time)
        # Run time must be found in file name.
        m = search(r'Time(\d+)\.dat', filename)
        if m:
            self.attrs['runtime'] = int(m.group(1))
        else:
            self.attrs['runtime'] = 0
        self.attrs['time'] = starttime + dt.timedelta(
            seconds=self.attrs['runtime'])

        # Open file
        f = open(filename, 'r')

        # Parse header.
        varlist = parse_tecvars(f.readline())
        m = search(r'.*I\=\s*(\d+)\,\s*J\=\s*(\d+)', f.readline()).groups()
        nLat, nLon = int(m[0]), int(m[1])
        self.attrs['nLat'] = nLat
        self.attrs['nLon'] = nLon

        # Create arrays.
        for k, u in varlist:
            self[k] = dmarray(np.zeros((nLat*nLon)), attrs={'units': u})

        for i in range(nLat*nLon):
            parts = f.readline().split()
            for k, (key, u) in enumerate(varlist):
                self[key][i] = float(parts[k])

        xy = np.sqrt(self['x']**2+self['y']**2)+0.0000001
        self['lon'] = np.arcsin(self['y']/xy)
        self['lat'] = dmarray(np.arcsin(self['z']),
                              {'units': 'deg'})*180.0/np.pi
        self['colat'] = 90.0 - self['lat']

        f.close()


class Line(IdlFile):
    '''
    Class for loading a single field line output file.
    At instantiation time, user may wish to set the start date and time of
    the simulation using the starttime kwarg.  If not given, start time
    will default to Jan. 1st, 2000, 00:00UT.

    Parameters
    ----------
    filename : str
        Name of line file to open
    starttime : datetime.datetime
        Simulation time corresponding to line file.
    '''
    def __init__(self, filename, starttime=None, *args, **kwargs):
        # Initialize as an IdlFile. Filename handled by parent reader and
        # stashed as self.attrs['filename']
        super(Line, self).__init__(filename, header=None, *args, **kwargs)

        if not starttime:
            starttime = dt.datetime(2000, 1, 1, 0, 0, 0)

        # Parse PWOM's units from variable names.
        for v in self.keys():
            if v == 'grid':
                continue
            matched = re.match('(\w+)\[(.*)\]', v)
            if matched:
                self[matched.groups()[0]] = self.pop(v)
                self[matched.groups()[0]].attrs['units'] = matched.groups()[1]

    def __repr__(self):
        return f"PWOM single field line output file {self.attrs['file']}"


class Lines(PbData):
    '''
    A class for reading and plotting a complete set of PWOM line output files.
    Open using a glob string that encompasses all of the lines that are
    intended to be read, e.g. 'north*.out'.  Note that unix wildcards are
    accepted.
    '''

    def __init__(self, lines, starttime=None, *args, **kwargs):
        super(Lines, self).__init__(*args, **kwargs)  # Init as PbData.
        if not starttime:
            starttime = dt.datetime(2000, 1, 1, 0, 0, 0)
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

        self['files'] = glob(lines)

        # Read first file; use it to set up all arrays.
        l1 = Line(self['files'][0], starttime=starttime)
        nA = l1.attrs['nAlt']
        nT = l1.attrs['nTime']
        nF = len(self['files'])
        self['r'] = l1['r']

        self.attrs['nAlt'] = nA
        self.attrs['nTime'] = nT
        self.attrs['nFile'] = nF
        self['time'] = l1['time']

        # Create all arrays.
        for v in l1._rawvar:
            self[v] = dmarray(np.zeros((nF, nT, nA)))
            self[v][0, :, :] = l1[v]
            if v == 'Lat' or v == 'Lon':
                self[v].attrs['units'] = 'deg'
            elif v[0] == 'u':
                self[v].attrs['units'] = 'km/s'
            elif v[0:3] == 'lgn':
                self[v].attrs['units'] = 'log(cm-3)'
            elif v[0] == 'T':
                self[v].attrs['units'] = 'K'
            else:
                self[v].attrs['units'] = None

        # Place data into arrays.
        for i, f in enumerate(self['files'][1:]):
            ln = Line(f)
            for v in l1._rawvar:
                self[v][i+1, :, :] = ln[v]

        # Get non-log densities.
        logVars = [v for v in self if v[:3] == 'lgn']
        for v in logVars:
            self[v[2:]] = dmarray(10.**self[v], {'units': r'$cm^{-3}$'})

    def calc_flux(self):
        '''
        Calculate flux in units of #/cm2/s.

        Variables saved as self[species+'Flux'].
        '''

        for s in ['H', 'O', 'He', 'e']:
            self[s+'Flux'] = dmarray(1000*self['n'+s]*self['u'+s],
                                     {'units': '$cm^{-2}s^{-1}$'})

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
        lat = np.pi*self['Lat']/180.0
        colat = 90. - self['Lat']
        lon = np.pi*self['Lon']/180.0

        # Values in a Cartesian plane.
        self._xGSM = r*np.cos(lat)*np.sin(lon)
        self._yGSM = r*np.cos(lat)*np.cos(lon)

        self._xLat = colat*np.sin(lon)
        self._yLat = colat*np.cos(lon)

    def add_slice(self, var, alt, time, nlev=31, zlim=None, target=None,
                  loc=111, title=None, latoffset=1.05,
                  rlim=50., add_cbar=True, clabel=None,
                  show_pts=False, show_alt=True, dolog=False,
                  lats=[75., 60.], colats=None, figsize=(8.34, 7),
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
        from matplotlib.patches import Circle
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import (LogLocator, LogFormatterMathtext)

        # Grab the slice of data that we want:
        if time is self['time'][0]:
            if time not in self['time']:
                raise ValueError('Time not in object')
            time = np.arange(self.attrs['nTime'])[self['time'] == time][0]

        fig, ax = set_target(target, loc=loc, figsize=figsize)
        ax.set_aspect('equal')

        # Create values in Cartesian plane.
        self._get_cartXY()

        # Grab values from correct time/location.
        x = self._xLat[:, time, alt]
        y = self._yLat[:, time, alt]
        value = self[var][:, time, alt]

        # Get max/min if none given.
        if zlim is None:
            zlim = [0, 0]
            zlim[0] = value.min()
            zlim[1] = value.max()
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z = np.where(value > zlim[0], value, 1.01*zlim[0])
            norm = LogNorm()
            ticks = LogLocator()
            fmt = LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z = value
            norm = None
            ticks = None
            fmt = None

        # Create contour plot.
        cont = ax.tricontourf(x, y, z, levs, *args, norm=norm, **kwargs)

        # Label altitude.
        if show_alt:
            ax.text(0.05, 0.05, r'Alt.={:.1f}$km/s$ ({:.2f}$R_E$)'.format(
                self['r'][alt], self['r'][alt]/6371.0+1), size=12,
                    transform=ax.transAxes, bbox={'fc': 'white', 'ec': 'k'})

        if show_pts:
            ax.plot(x, y, '+w')

        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel is None:
                clabel = "%s (%s)" % (var, self[var].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Set title, labels, axis ranges (use defaults where applicable.)
        if title:
            ax.set_title(title)
        ax.set_yticks([]), ax.set_xticks([])
        ax.set_ylabel(r'Sun $\rightarrow$')
        colat_lim = 90.-rlim
        ax.set_xlim([-1*colat_lim, colat_lim])
        ax.set_ylim([-1*colat_lim, colat_lim])

        if colats:
            for colat in colats:
                r = latoffset*colat/np.sqrt(2)
                ax.add_patch(Circle([0, 0], colat, fc='none',
                                    ec='k', ls='dashed'))
                ax.text(r, r, '{:.0f}'.format(colat)+r'$^{\circ}$',
                        rotation=315, ha='center', va='center')
        else:
            for lat in lats:
                colat = 90 - lat
                r = latoffset*colat/np.sqrt(2)
                ax.add_patch(Circle([0, 0], colat, fc='none',
                                    ec='k', ls='dashed'))
                ax.text(r, r, '{:.0f}'.format(lat)+r'$^{\circ}$',
                        rotation=315, ha='center', va='center')

        return fig, ax, cont, cbar
