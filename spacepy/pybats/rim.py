#!/usr/bin/env python
'''
Classes, functions, and methods for reading, writing, and plotting output
from the Ridley Ionosphere Model (RIM) and the similar legacy code,
Ridley Serial.

Copyright 2010 Los Alamos National Security, LLC.

.. rubric:: Classes

.. autosummary::
    :toctree:
    :template: clean_class.rst

    Iono
    OvalDebugFile

.. rubric:: Functions

.. autosummary::
    :toctree:

    get_iono_cb
    tex_label
'''

import os
import numpy as np
from spacepy import deprecated
from spacepy.plot import set_target
from spacepy.pybats import PbData, dmarray

def get_iono_cb(ct_name='bwr'):
    '''
    Several custom colorbars used by RIM and AMIE have become standard when
    visualizing data from these models.  These are 'blue_white_red' and
    'white_red', used for data that have positive and negative values and
    for data that have only positive values, respectively.  This function
    builds and returns these colorbars when called with the initials of the
    color table name as the only argument.

    Other Parameters
    ================
    ct_name : str
       Select the color table.  Can be 'bwr' for blue-white-red or 'wr' for
       white-red.  Defaults to 'bwr'.

    Examples
    ========
    >>> bwr_map = get_iono_cb('bwr')
    >>> wr_map  = get_iono_cb('wr')

    '''

    from matplotlib.colors import LinearSegmentedColormap as lsc

    if ct_name=='bwr':
        table = {
            'red':  [(0.,0.,.0),(.34,0.,0.),(.5,1.,1.),(1.,1.,1.)],
            'green':[(0.,0.,0.),(.35,1.,1.),(.66,1.,1.),(1.,0.,0.)],
            'blue' :[(0.,1.,1.),(.5,1.,1.),(.66,0.,0.),(.85,0.,0.),(1.,.1,.1)]
            }
        cmap = lsc('blue_white_red',table)
    elif ct_name=='wr':
        table = {
            'red':  [(0.,1.,1.),(1.,1.,1.)],
            'green':[(0.,1.,1.),(1.,0.,0.)],
            'blue' :[(0.,1.,1.),(1.,.0,.0)]
            }
        cmap = lsc('white_red',table)

    return cmap

def tex_label(varname):
    '''
    Many variable names used in the Ridley Ionosphere Model look much better
    in LaTeX format with their proper Greek letters.  This function takes
    a variable name, and if it is recognized, returns a properly formatted
    string that uses MatPlotLib's MathText functionality to display the
    proper characters.  If it is not recognized, the varname is returned.

    Parameters
    ==========
    varname : string
       The variable to convert to a LaTeX label.

    Examples
    ========
    >>>tex_label('n_phi')
    '\\Phi_{Ionosphere}'
    >>>tex_label('Not Recognized')
    'Not Recognized'

    '''

    if varname[:2]=='n_' or varname[:2]=='s_':
        varname=varname[2:]

    known = {
        'phi': r'$\Phi_{Ionosphere}$',
        'sigmah':r'$\sigma_{Hall}$',
        'sigmap':r'$\sigma_{Peder}$',
        'jr':r'$J_{radial}$',
        'mA/m^2':r'$\mu A/m^{2}$',
        'W/m2':r'$W/m^2$',
        'eV':'$eV$'
        }


    if varname in known:
        label = known[varname]
    else:
        label = varname

    return label

class Iono(PbData):
    '''
    A class for handling 2D output from the Ridley Ionosphere Model.
    Instantiate an object as follows:

    >>> iono = rim.Iono('filename.idl')

    ...where filename.idl is the name of a RIM 2D output file.
    '''

    def __init__(self, infile, *args, **kwargs):
        super(Iono, self).__init__(*args, **kwargs)
        self.attrs['file']=infile
        self.readascii()

    def readascii(self):
        '''
        Read an ascii ".idl" output file and load the contents into the object.
        '''

        from sys import version_info
        import re
        import datetime as dt
        from numpy import zeros, reshape
        import gzip

        # slurp entire file.
        if self.attrs['file'][-3:]=='.gz':
            try:
                infile = gzip.open(self.attrs['file'],'rt')
            except ValueError: #Python2.7 (Windows) compatibility
                infile = gzip.open(self.attrs['file'])
        else:
            infile = open(self.attrs['file'], 'r')
        raw = infile.readlines()
        infile.close()

        # Parse header
        title = raw[raw.index('TITLE\n')+1]
        self.attrs['title'] = title[title.index('"')+1:title.rindex('"')]

        i = raw.index('NUMERICAL VALUES\n')
        self.attrs['nvars'] = int(raw[i+1].split()[0])
        self.attrs['ntheta']= int(raw[i+2].split()[0])
        self.attrs['nphi']  = int(raw[i+3].split()[0])

        # Convenience:
        nphi, ntheta = self.attrs['nphi'], self.attrs['ntheta']

        i = raw.index('TIME\n')
        self.attrs['time'] = dt.datetime(
            int(raw[i+1].split()[0]),      #year
            int(raw[i+2].split()[0]),      #month
            int(raw[i+3].split()[0]),      #day
            int(raw[i+4].split()[0]),      #hour
            int(raw[i+5].split()[0]),      #min
            int(raw[i+6].split()[0]),      #sec
            int(raw[i+7].split()[0])*1000  #microsec
            )

        i = raw.index('SIMULATION\n')
        self.attrs['iter']    =   int(raw[i+1].split()[0])
        self.attrs['simtime'] = float(raw[i+2].split()[0])

        i = raw.index('DIPOLE TILT\n')
        self.tilt = zeros(2)
        self.tilt[0] = float(raw[i+1].split()[0])
        self.tilt[1] = float(raw[i+2].split()[0])

        i = raw.index('VARIABLE LIST\n')
        namevar = []
        units   = {}
        for j in range(i+1,i+self.attrs['nvars']+1):
            match = re.match(r'\s*\d+\s+([\w\s\W]+)\[([\w\s\W]+)\]',raw[j])
            if match:
                name = (match.group(1).strip()).lower()
                namevar.append(name)
                units[name] = match.group(2).strip()
            else:
                raise ValueError('Could not parse %s' % raw[j])


        ### Read all data ###

        # Create data arrays
        nPts = self.attrs['ntheta']*self.attrs['nphi']
        for key in namevar:
            self['n_'+key] = dmarray(zeros(nPts), {'units':units[key]})
            self['s_'+key] = dmarray(zeros(nPts), {'units':units[key]})
        i = raw.index('BEGIN NORTHERN HEMISPHERE\n')+1

        # Some compilers insert line breaks automatically when fortran format
        # string is not adequately specified.  Let's see if that's the
        # case here: how many lines does it take to cover all variables?
        nvars, nvarline, nwrap = len(namevar), 0, 0
        while nvarline<nvars:
            nvarline += len(raw[i+nwrap].split())
            nwrap += 1

        # Fill data arrays:
        for j in range(nPts):
            # Build list of variables; accounting for line-wrapping:
            parts = []
            iLine = i + j*nwrap
            for iwrap in range(nwrap):
                parts += raw[iLine+iwrap].split()
            for k in range(self.attrs['nvars']):
                self['n_'+namevar[k]][j] = parts[k]
        i = raw.index('BEGIN SOUTHERN HEMISPHERE\n')+1
        for j in range(nPts):
            parts = []
            iLine = i + j*nwrap
            for iwrap in range(nwrap):
                parts += raw[iLine+iwrap].split()
            for k in range(self.attrs['nvars']):
                self['s_'+namevar[k]][j] = parts[k]

        # Create 2-D arrays.
        for key in namevar:
            nkey, skey = 'n_'+key, 's_'+key
            self[nkey] = reshape(self[nkey], (ntheta, nphi), 'F')
            self[skey] = reshape(self[skey], (ntheta, nphi), 'F')

        # Some extra grid info:
        self.dlon = self['n_psi'  ][0,3]-self['n_psi'  ][0,2]
        self.dlat = self['n_theta'][3,0]-self['n_theta'][2,0]

    def calc_j(self):
        '''
        Calculate total horizontal current as related values.  Each will be
        stored into *self* using the typical key-value approach.  Calculations
        are done for both the northern and southern hemisphere with the
        appropriate prefixes ('n_' and 's_') applied to each key.

        | key  | Description |
        |------|-------------|
        | j    | Total horizontal current, $sqrt(jx^2+jy^2+jz^2)$ |
        | jphi | Azimuthal current, positive values are eastward. |


        Parameters
        ==========

        Returns
        =======
        True

        Examples
        ========
        >>> a = rim.Iono('spacepy/tests/data/pybats_test/it000321_104510_000.idl.gz')
        >>> a.calc_j()
        >>> print(a['n_jphi'])

        '''

        # Loop over hemispheres:
        for h in ('n_', 's_'):
            # Calculate total horizontal current
            self[h+'j'] = np.sqrt(  self[h+'jx']**2
                                  + self[h+'jy']**2
                                  + self[h+'jz']**2 )

            # Calculate total azimuthal current (i.e., electrojets):
            self[h+'jphi'] = self[h+'jy'] * np.cos(np.pi/180. * self[h+'psi']) \
                - self[h+'jx'] * np.sin(np.pi/180. * self[h+'psi'])

        return True

    def calc_I(self):
        '''
        Integrate radial current to get the following values:
        I:     The total radial current over one hemisphere.
        Iup:   The total upward current over one hemisphere.
        Idown: The total downward current over one hemisphere.

        Values are stored in the object with the prefix 'n_' or 's_'
        to indicate the hemisphere and may be accessed via self['n_I'], etc.

        Values are stored in units of mega Amps.

        Parameters
        ==========


        Returns
        =======
        True


        Examples
        ========
        >>> a = rim.Iono('spacepy/tests/data/pybats_test/it000321_104510_000.idl.gz')
        >>> a.calc_I()
        >>> print(a['n_Iup'])

        '''

        # Calculate some physically meaningful values/units
        units = 1E-6*1E-6  # micro amps to amps, amps to MegaAmps
        R = (6371.0+110.0)*1000.0 # Radius of Earth + iono altitude
        dTheta = np.pi*self.dlat/180.
        dPhi   = np.pi*self.dlon/180.

        # -----NORTHERN HEMISPHERE-----
        # Get relevant values:
        colat     = self['n_theta']*np.pi/180.
        integrand = self['n_jr']*np.sin(colat)*dTheta*dPhi
        # Get locations of "up" and "down"
        loc_up = self['n_jr']>0
        loc_do = self['n_jr']<0
        self['n_I']     = units*R**2 * np.sum(integrand)
        self['n_Iup']   = units*R**2 * np.sum(integrand[loc_up])
        self['n_Idown'] = units*R**2 * np.sum(integrand[loc_do])

        # -----SOUTHERN HEMISPHERE-----
        # Get relevant values:
        colat     = self['s_theta']*np.pi/180.
        integrand = self['s_jr']*np.sin(colat)*dTheta*dPhi
        # Get locations of "up" and "down"
        loc_up = self['s_jr']>0
        loc_do = self['s_jr']<0
        self['s_I']     = units*R**2 * np.sum(integrand)
        self['s_Iup']   = units*R**2 * np.sum(integrand[loc_up])
        self['s_Idown'] = units*R**2 * np.sum(integrand[loc_do])

    def add_cont(self, var, target=None, n=50, maxz=False, lines=False,
                 cmap=False, add_cbar=False, label=None, loc=111,
                 xticksize=12, yticksize=12, max_colat=40, **kwargs):
        '''
        Create a polar contour of variable *var*.  Plot will be either drawn
        on a new matplotlib figure and axes, or you can specify a plot target
        using the *target* kwarg.

        Parameters
        ==========
        var : str
           The name of the variable to plot.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object
        cnt : matplotlib contour object
        cb  : matplotlib colorbar object

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
        n : int
            Set number of levels.  Default is 50.  If unfilled lines are
            added (see *lines* kwarg), multiples of three provide coherence
            between filled and unfilled contours.
        lines : bool
            Add unfilled black solid/dashed contours to plot for additional
            contrast.  Default is **True**.
        maxz : real
            Set the max/min value for the color bar.  Default is set by data.
        max_colat : real
            Set the co-latitude range of the plot, in degrees.  Defaults to 40.
        cmap : str
            Set the colormap.  Default is to use Matplotlib's "Reds" for data
            that is strictly positive and "Seismic" for diverging data.
            Alternatively, legacy Ridley Ionosphere Model color maps can be
            loaded using "l_wr" (white red) or "l_bwr" (blue-white-red),
            where the "l_" prefix indicates legacy and not Matplotlib color maps.
        add_cbar : bool
            Add colorbar to plot.  Default is **False** which will
            not add one to the plot.
        label : str or None
            Label to place at top of plot.  Defaults to variable name.
        xticksize : int
            Size of longitude markers in text points.  Defaults to 12.
        yticksize : int
            Size of latitude markers in text points.  Defaults to 12.

        Extra keywords are passed to the contourf routine.

        '''
        # Get only what we need to decrease runtime.
        from math import pi
        from numpy import linspace
        from matplotlib.colors import Normalize
        from matplotlib.ticker import MaxNLocator, MultipleLocator
        from matplotlib.pyplot import clabel, colorbar

        fig, ax = set_target(target, polar=True, loc=loc)

        hemi = var[:2]

        # user defined variables may not have hemisphere marking.
        # Assume nothern hemi in those cases.
        if hemi!='n_' and hemi!='s_': hemi = 'n_'

        # Set levels and ticks:
        if label==None:
            label=tex_label(var)
        lt_labels = ['06',    label, '18',   '00']
        xticks    = [   0,   pi/2,   pi, 3*pi/2]
        lct = MultipleLocator(10)
        minz = self[var].min()
        if minz < 0.0:
            if not maxz:
                maxz = max([abs(minz),self[var].max()])
            crange = Normalize(vmin=-1.*maxz, vmax=maxz)
            levs = linspace(-1.*maxz, maxz, n)
        else:
            if not maxz:
                maxz = self[var].max()
            crange = Normalize(vmin=0., vmax=maxz)
            levs = linspace(0., maxz, n)

        # Get color map if not given:
        if not cmap:
            if self[var].min() >= 0.0:
                cmap='Reds'
            else:
                cmap='seismic'
        # Search for legacy maps:
        elif 'l_' in cmap:
            cmap=get_iono_cb(cmap[2:])

        # Set the latitude based on hemisphere:
        theta = self[hemi+'theta']
        if 's_' in hemi: theta = 180-self[hemi+'theta']

        # Create contour:
        cnt1 = ax.contourf(self[hemi+'psi']*pi/180.0+pi/2., theta,
                           np.array(self[var]), levs, norm=crange, cmap=cmap,
                           **kwargs)
        # Set xtick label size, increase font of top label.
        labels = ax.get_xticklabels()
        for l in labels: l.set_size(xticksize)
        labels[1].set_size(xticksize*1.25)

        if lines:
            nk = int(round(n/3.0))
            cnt2 = ax.contour(self[hemi+'psi']*pi/180.0+pi/2., theta,
                              np.array(self[var]), nk, colors='k')
            #clabel(cnt2,fmt='%3i',fontsize=10)

        if add_cbar:
            cbarticks = MaxNLocator(7)
            cbar = colorbar(cnt1, ticks=cbarticks, shrink=0.75, pad=0.08, ax=ax)
            cbar.set_label(tex_label(self[var].attrs['units']))
        else:
            cbar=False
        ax.set_xticks(xticks)
        ax.set_xticklabels(lt_labels)
        ax.yaxis.set_major_locator(lct)
        ax.set_ylim([0,max_colat])

        # Use text function to manually add pretty ticks.
        ax.set_yticklabels('') # old ticks off.
        opts = {'size':yticksize, 'rotation':-45, 'ha':'center', 'va':'center'}
        for theta in [80.,70.,60.]:
            txt = '{:02.0f}'.format(theta)+r'$^{\circ}$'
            ax.text(pi/4., 90.-theta, txt, color='w', weight='heavy', **opts)
            ax.text(pi/4., 90.-theta, txt, color='k', weight='light', **opts)

        return fig, ax, cnt1, cbar

class OvalDebugFile(PbData):
    '''
    The auroral oval calculations in RIM may spit out special debug files that
    are extremely useful.  This class handles reading and plotting the data
    contained within those files.
    '''

    def __init__(self, infile, *args, **kwargs):
        super(OvalDebugFile, self).__init__(*args, **kwargs)
        self.attrs['file']=infile
        self._readascii()

    def _readascii(self):
        '''
        Loads data & populates object upon instantiation.
        '''

        import datetime as dt

        # Slurp in lines:
        f = open(self.attrs['file'])
        lines = f.readlines()
        f.close()

        # Parse header to get longitude in radians:
        l = lines.pop(0)
        self['lon'] = np.array(l.split('=')[-1].split(), dtype=float)* np.pi/180.
        l = lines.pop(0)

        # Some helper vars:
        nLons = self['lon'].size
        nLine = len(lines)

        # Create container arrays:
        self['time'] = np.zeros(nLine, dtype=object)
        self['oval'] = np.zeros( (nLine, nLons) )

        # Parse rest of file:
        for j,l in enumerate(lines):
            self['time'][j]   = dt.datetime.strptime(l[:19], '%Y %m %d %H %M %S')
            self['oval'][j,:] = l.split()[7:]


    def get_oval(self, time, interp=True):
        '''
        Given a datetime, *time*, interpolate the oval position to that time.
        If *time* is outside the range of self['time'], an exception will be
        raised.  Linear interpolation is used unless kwarg *interp* is False.
        In that case, the value of the oval nearest to *time* will be used
        without interpolation.

         Parameters
        ==========
        time : datetime.datetime
           The time at which the oval is requested.

        Returns
        =======
        oval : numpy.ndarray
           Oval location at time *time* in radians.  Zero corresponds to local
           noon.

        Other Parameters
        ================
        interp : bool
            Set interpolation behavior.  If **True**, linear interpolation is
            used.  If **False**, the oval location nearest to *time* is used
            without interpolation.

        Examples
        ========
        >>>
        '''

        from matplotlib.dates import date2num

        # If *time* is outside the bounds of self['time'], raise exception
        # to avoid extrapolation.
        if time<self['time'][0] or time>self['time'][-1]:
            raise ValueError('Given time outside object range ' +
                             'and requires extrapolation')

        # Turn datetimes into numbers.
        time     = date2num(time)
        ovaltime = date2num(self['time'])

        # Start by obtaining the indices of the time array that bound
        index = np.arange(self['time'].size)
        i1 = index[time>=ovaltime][-1] # Value before time
        i2 = index[time<=ovaltime][ 0] # Value after time

        # Check for exact matches:
        dT1 = np.abs(ovaltime[i1]-time)
        dT2 = np.abs(ovaltime[i2]-time)
        if dT1==0: return self['oval'][i1]
        if dT2==0: return self['oval'][i2]

        # If no interpolation, just send back nearest neighbor:
        if not interp:
            if dT1<dT2:
                return self['oval'][i1]
            else:
                return self['oval'][i2]

        # If interpolation, get slope and y-intercept vectors:
        dT = ovaltime[i2]-ovaltime[i1]
        m  = (self['oval'][i2]-self['oval'][i1])/dT
        b  = self['oval'][i2] - m*ovaltime[i2]

        # Interpolate and return:
        return m*time + b


    def add_oval_line(self, time, *args, **kwargs):
        '''
        Adds the location of the auroral oval at time *time* as a line plot on to
        *target*, which must be a Matplotlib figure/axes.

        If *target* not given, a new axes object is created.  If *target* is
        not an existing axes object, a new axes object is created and customized
        to be an ionosphere polar plot.

        Extra arguments/kwargs are sent to :meth:`matplotlib.pyplot.plot`.

        Parameters
        ==========
        time : integer or datetime
           Sets the time at which to plot the line.  If an integer is used,
           the integer is used to index the internal array naively.  If a
           datetime object is given and is within the bounds of self['time'],
           the oval location will be interpolated to *time*.  Control of the
           interpolation happens via the *interp* keyword.

        Returns
        =======
        fig  : matplotlib figure object
        ax   : matplotlib axes object
        line : matplotlib line object

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
        interp : bool
            Control the behavior of time interpolation when the *time* argument
            is given as a datetime object.  Defaults to **True**, meaning that
            oval location is interpolated in time.

        Examples
        ========
        >>>
        '''

        import datetime as dt

        # Handle kwargs not being handed to plot:
        target=None; loc=111; interp=True
        if 'target' in kwargs: target=kwargs.pop('target')
        if 'interp' in kwargs: interp=kwargs.pop('interp')
        if 'loc'    in kwargs: loc=kwargs.pop('loc')

        # Set plot targets:
        fig, ax = set_target(target, polar=True, loc=loc)

        # Get oval interpolated in time:
        if type(time) == dt.datetime:
            oval = self.get_oval(time, interp=interp)
        elif type(time) == int:
            oval = self['oval'][time]
        else:
            raise TypeError('Unrecognized type for *time*:'+type(time))

        # Plot oval:
        line = ax.plot(self['lon']+np.pi/2., oval, *args, **kwargs)

        return fig, ax, line

