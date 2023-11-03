'''
A PyBats module for handling input, output, and visualization of
binary SWMF output files taylored to BATS-R-US-type data.

.. currentmodule:: spacepy.pybats.bats

.. rubric:: Classes

.. autosummary::
    :template: clean_class.rst
    :toctree:

    BatsLog
    Stream
    Bats2d
    Mag
    MagFile
    GeoIndexFile
    VirtSat
'''
import sys

import numpy as np
from spacepy.pybats import PbData, IdlFile, LogFile, calc_wrapper
from spacepy.plot import set_target, applySmartTimeTicks
from spacepy.datamodel import dmarray

# Module-level variables:
# recognized species:
mass = {'hp': 1.0, 'op': 16.0, 'he': 4.0,
        'sw': 1.0, 'o': 16.0, 'h': 1.0, 'iono': 1.0, '': 1.0}

RE = 6371000  # Earth radius in meters.

# Module-level functions:


def _calc_ndens(obj):
    '''
    Given an object of subclass :class:`~spacepy.pybats.PbData` that uses
    standard BATS-R-US variable names, calculate the number density from
    all mass density variables (i.e., all variables that end in *rho*).
    New variables are added to the object.  Variable names are constructed
    by taking the mass-density variable, *Speciesrho*, and replacing *rho*
    with *N*.  Total number density is also saved as *N*.

    Composition information is also saved by taking each species and
    calculating the percent of the total number density via
    *fracspecies* = 100% x *speciesN*/*N*.

    Many species that are implemented in multi-species and multi-fluid
    BATS-R-US equation files have known mass.  For example, *oprho* is known
    a-priori to be singly ionized oxygen mass density, so *opN* is created
    by dividing *oprho* by 16.  If the fluid/species is not recognized, it is
    assumed to be hyrogen.  The single atom/molecule mass is saved in the
    attributes of the new variable.

    This function should be called by object methods of the same name.
    It is placed at the module level because it is used by many different
    classes.

    Parameters
    ==========
    obj : :class:`~spacepy.pybats.PbData` object
       The object on which to act.

    Other Parameters
    ================

    Returns
    =======
    True

    Examples
    ========
    >>> import spacepy.pybats.bats as pbs
    >>> mhd = pbs.Bats2d('spacepy-code/spacepy/pybats/slice2d_species.out')
    >>> pbs._calc_ndens(mhd)

    '''

    species = []
    names = []

    # Get name of Rho (case sensitive check):
    rho = 'Rho'*('Rho' in obj) + 'rho'*('rho' in obj)

    # Find all species: the variable names end or begin with "rho".
    # Take care not to double-count.
    for k in obj:
        # Ends with rho?
        if (k[-3:] == rho) and (k != rho) and (k[:-3] + 'N' not in obj):
            species.append(k)
            names.append(k[:-3])
        # Begins with rho?
        if (k[:3] == rho) and (k != rho) and (k[3:] + 'N' not in obj):
            species.append(k)
            names.append(k[3:])

    # Individual number density
    for s, n in zip(species, names):
        if n.lower() in mass:
            m = mass[n.lower()]
        else:
            m = 1.0
        obj[n+'N'] = dmarray(obj[s]/m,
                             attrs={'units': '$cm^{-3}$', 'amu mass': m})

    # Total N is sum of individual  number densities.
    obj['N'] = dmarray(np.zeros(obj[rho].shape),
                       attrs={'units': '$cm^{-3}$'})
    if species:
        # Total number density:
        for n in names:
            obj['N'] += obj[n+'N']
        # Composition as fraction of total per species:
        for n in names:
            obj[n+'Frac'] = dmarray(100.*obj[n+'N']/obj['N'],
                                    {'units': 'Percent'})
    else:
        # No individual species => no composition, simple ndens.
        obj['N'] += dmarray(obj[rho], attrs={'units': '$cm^{-3}$'})


# Classes:

class BatsLog(LogFile):
    '''
    A specialized version of :class:`~spacepy.pybats.LogFile` that includes
    special methods for plotting common BATS-R-US log file values, such as
    D$_{ST}$.

    .. autosummary::

        ~BatsLog.add_dst_quicklook

    .. automethod:: add_dst_quicklook
    '''

    def fetch_obs_dst(self):
        '''
        Fetch the observed Dst index for the time period covered in the
        logfile.  Return *True* on success.

        Observed Dst is automatically fetched from the Kyoto World Data Center
        via the :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoDst` object, which holds the observed
        Dst, is stored as *self.obs_dst* for future use.
        '''

        import spacepy.pybats.kyoto as kt

        # Return if already obtained:
        if hasattr(self, 'obs_dst'):
            return True

        # Start and end time to collect observations:
        stime = self['time'][0]
        etime = self['time'][-1]

        # Attempt to fetch from Kyoto website:
        try:
            self.obs_dst = kt.fetch('dst', stime, etime)
        # Warn on failure:
        except BaseException as args:
            raise Warning('Failed to fetch Kyoto Dst: ', args)
            return False

        return True

    def fetch_obs_sym(self):
        '''
        Fetch the observed SYM-H index for the time period covered in the
        logfile.  Return *True* on success.

        Observed SYM-H is automatically fetched from the Kyoto World Data
        Center via the :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoSym` object, which holds the observed
        Dst, is stored as *self.obs_sym* for future use.
        '''

        import spacepy.pybats.kyoto as kt

        # Return if already obtained:
        if hasattr(self, 'obs_sym'):
            return True

        # Start and end time to collect observations:
        stime = self['time'][0]
        etime = self['time'][-1]

        # Attempt to fetch from Kyoto website:
        try:
            self.obs_sym = kt.fetch('sym', stime, etime)
        # Warn on failure:
        except BaseException as args:
            raise Warning('Failed to fetch Kyoto Dst: ', args)
            return False

        return True

    def add_dst_quicklook(self, target=None, loc=111, plot_obs=False,
                          epoch=None, add_legend=True, plot_sym=False,
                          dstvar=None, lw=2.0, obs_kwargs={'ls': '-.'},
                          sym_kwargs={'ls': '-'}, **kwargs):
        '''
        Create a quick-look plot of Dst (if variable present in file)
        and compare against observations.

        Like all *add_\\* * methods in Pybats, the *target* kwarg determines
        where to place the plot.
        If kwarg *target* is **None** (default), a new figure is
        generated from scratch.  If *target* is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot location
        *loc* (defaults to 111).  If target is a matplotlib Axes object,
        the plot is placed into that axis at subplot location *loc*.

        With newer versions of BATS-R-US, new dst-like variables are included,
        named 'dst', 'dst-sm', 'dstflx', etc.  This subroutine will attempt
        to first use 'dst-sm' as it is calculated consistently with
        observations.  If not found, 'dst' is used.  Users may choose which
        value to use via the *dstvar* kwarg.

        Observed Dst and SYM-H is automatically fetched from the Kyoto World
        Data Center via the :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoDst` or
        :class:`spacepy.pybats.kyoto.KyotoSym` object, which holds the observed
        Dst/SYM-H, is stored as *self.obs_dst* for future use.
        The observed line can be customized via the *obs_kwargs* kwarg, which
        is a dictionary of plotting keyword arguments.

        If kwarg *epoch* is set to a datetime object, a vertical dashed line
        will be placed at that time.

        The figure and axes objects are returned to the user.

        .. versionchanged:: 0.5.0
            Separated kwargs for Sym-H and Dst lines.
        '''

        from datetime import datetime

        # Set the correct dst value to be used.
        if not dstvar:
            if 'dst_sm' in self:
                dstvar = 'dst_sm'
            else:
                dstvar = 'dst'

        if dstvar not in self:
            return None, None

        fig, ax = set_target(target, figsize=(10, 4), loc=loc)

        if 'label' not in kwargs:
            kwargs['label'] = 'BATS-R-US $D_{ST}$ (Biot-Savart)'

        if 'label' not in obs_kwargs:
            obs_kwargs['label'] = 'Obs. Dst'
        if 'label' not in sym_kwargs:
            sym_kwargs['label'] = 'Obs. SYM-H'

        ax.plot(self['time'], self[dstvar], lw=lw, **kwargs)
        ax.hlines(0.0, self['time'][0], self['time'][-1],
                  'k', ':', label='_nolegend_')
        applySmartTimeTicks(ax, self['time'])
        ax.set_ylabel('D$_{ST}$ ($nT$)')
        ax.set_xlabel('Time from ' + self['time'][0].isoformat()+' UTC')

        # Add observations (Dst and/or SYM-H):
        if plot_obs:
            # Attempt to fetch observations, plot if success.
            if self.fetch_obs_dst():
                ax.plot(self.obs_dst['time'], self.obs_dst['dst'],
                        **obs_kwargs)
                applySmartTimeTicks(ax, self['time'])
        if plot_sym:
            # Attempt to fetch SYM-h observations, plot if success.
            if self.fetch_obs_sym():
                ax.plot(self.obs_sym['time'], self.obs_sym['sym-h'],
                        **sym_kwargs)
                applySmartTimeTicks(ax, self['time'])

        # Place vertical line at epoch:
        if type(epoch) == datetime:
            yrange = ax.get_ylim()
            ax.vlines(epoch, yrange[0], yrange[1], linestyles='dashed',
                      colors='k', linewidths=1.5)
            ax.set_ylim(yrange)

        # Apply legend
        if add_legend:
            ax.legend(loc='best')
        if target is None:
            fig.tight_layout()

        return fig, ax


class Extraction(PbData):
    '''
    A class for creating and visualizing extractions from other
    :class:`~spacepy.pybats.PbData` 2D data sets.  Bilinear interpolation is
    used to obtain data between points.  At present, this class only
    works with :class:`~spacepy.pybats.bats.Bats2d` objects, but will be
    generalized in the future.  Though it can be instantiated in typical
    fashion, It is best to create these objects via other object methods (e.g.,
    :func:`~spacepy.pybats.bats.Bats2d.extract`)

    Parameters
    ----------
    x : float or sequence
        X value(s) for points at which to extract data values.
    y : float or sequence
        Y value(s) for points at which to extract data values.
    dataset : Bats
        :class:`~spacepy.pybats.bats.Bats2d` object from which to extract
        values.

    Other Parameters
    ----------------
    var_list : string or sequence of strings
         List of values to extract from *dataset*.  Defaults to 'all', for all
         values within *dataset*.

    '''

    def __init__(self, x, y, dataset, var_list='all', *args, **kwargs):
        # Initialize as a PbData object.
        super(Extraction, self).__init__(*args, **kwargs)

        # Store x, y internally:
        # If our x, y locations are not numpy arrays, fix that.
        # Additionally, convert scalars to 1-element vectors.
        x, y = np.array(x), np.array(y)
        if not x.shape:
            x = x.reshape(1)
        if not y.shape:
            y = y.reshape(1)

        # Default: all variables are extracted except coordinates.
        if var_list == 'all':
            var_list = list(dataset.keys())
            for v in ('x', 'y', 'z', 'grid'):
                if v in var_list:
                    var_list.remove(v)
        elif type(var_list) is str:
            var_list = var_list.split()
        elif type(var_list) != list:
            var_list = []
        else:
            raise TypeError('Kwarg var_list must be string, list, or None.')

        # Stash var list, dataset internally:
        self._var_list = var_list
        self._dataset = dataset

        # Extract along trace:
        self.extract(x, y)

    def extract(self, xpts, ypts):
        # Perform actual extraction:
        from spacepy.pybats.batsmath import interp_2d_reg

        # Some convenience variables:
        x, y = xpts, ypts
        data = self._dataset

        for v, value in zip(data['grid'].attrs['dims'], (x, y)):
            self[v] = dmarray(value, attrs=data[v].attrs)
        # Create data object for holding extracted values.
        # One vector for each value w/ same units as parent object.
        for v in self._var_list:
            self[v] = dmarray(np.zeros(x.shape), attrs=data[v].attrs)

        # Some helpers:
        xAll = data[data['grid'].attrs['dims'][0]]
        yAll = data[data['grid'].attrs['dims'][1]]

        # Navigate the MHD quad tree, interpolating as we go.
        for k in data.qtree:
            # Only leafs, of course!
            if not data.qtree[k].isLeaf:
                continue
            # Find all points that lie in this block.
            # If there are none, just keep going.
            pts = \
                (x >= data.qtree[k].lim[0]) & (x <= data.qtree[k].lim[1]) & \
                (y >= data.qtree[k].lim[2]) & (y <= data.qtree[k].lim[3])
            if not pts.any():
                continue
            locs = data.qtree[k].locs
            for v in self._var_list:
                self[v][pts] = interp_2d_reg(x[pts], y[pts], xAll[locs],
                                             yAll[locs], data[v][locs])


class Stream(Extraction):
    '''
    A class for streamlines.  Contains all of the information about
    the streamline, including extracted variables.

    Upon instantiation, the object will trace through the vector
    field determined by the "[x/y]field" values and the Bats object
    "bats".

    Parameters
    ----------
    bats : Bats
        :class:`~spacepy.pybats.bats.Bats2d` object through which to trace.
    xstart : float
        X value of location to start the trace.
    ystart : float
        Y value of location to start the trace.
    xfield : str
        Name of variable in ``bats`` which contains X values of the field
    yfield : str
        Name of variable in ``bats`` which contains Y values of the field

    Other Parameters
    ----------------
    style : str
        Sets line style, including colors. See :meth:`set_style` for details.
        (Default 'mag')
    type : str
        (Default 'streamline')
    method : str
        Integration method. The default is Runge-Kutta 4 ('rk4') which gives
        a good blend of speed and accuracy. See the test functions in
        :mod:`~spacepy.pybats.trace2d` for more info.  The other option is
        a simple Euler's method approach ('eul'). (Default 'rk4')
    extract : bool
        (Default: False) Extract variables along stream trace and save within
        object.
    maxPoints : int
        (Default : 20000) Maximum number of integration steps to take.
    var_list : string or sequence of strings
        (Default : 'all') List of values to extract from *dataset*.
        Defaults to 'all', for all values within *bats*.

    Notes
    -----
    .. Not really "notes" but need to keep this section from being parsed
       as parameters

    .. rubric:: Methods

    .. autosummary::

        ~Stream.set_style
        ~Stream.treetrace
        ~Stream.trace
        ~Stream.plot

    .. automethod:: set_style
    .. automethod:: treetrace
    .. automethod:: trace
    .. automethod:: plot
    '''

    def __init__(self, bats, xstart, ystart, xfield, yfield, style='mag',
                 type='streamline', method='rk4', var_list='all',
                 extract=False, maxPoints=20000, *args, **kwargs):

        # Key values:
        self.xstart = xstart  # X and Y starting
        self.ystart = ystart  # points in the field.
        self.xvar = xfield  # Strings that list the variables
        self.yvar = yfield  # that will be used for tracing.

        # Descriptors:
        self.open = True
        self.status = 'open'
        self.method = method

        # Do tracing if a tracing method has been set.
        self.treetrace(bats, maxPoints=maxPoints)

        # Now, initialize as Extraction using self.x, self.y.
        # This handles the extraction as well as converting self.x, y into
        # dmarrays with the correct names.
        super(Stream, self).__init__(self.x, self.y, bats,
                                     var_list=var_list*extract)

        # Place parameters into attributes:
        self.attrs['start'] = [xstart, ystart]
        self.attrs['trace_var'] = [xfield, yfield]
        self.attrs['method'] = method

        # Set style
        self.set_style(style)

    def set_style(self, style):
        '''
        Set the line style either using a simple matplotlib-type style
        string or using a preset style type.  Current types include:

        'mag' : treat line as a magnetic field line.  Closed lines are
                white, other lines are black.
        '''

        import re

        # Here, we can set line color and style based on
        # line characteristics.  Right now, only one preset is
        # available.
        if style == 'mag':
            if self.open:
                self.style = 'k-'
                self.color = 'k'
            else:
                self.style = 'w-'
                self.color = 'w'
        else:
            self.style = style
            col = re.match('([bgrcmykw])', style)
            if col:
                self.color = col.groups()[0]
            else:
                self.color = 'k'

    def treetrace(self, bats, maxPoints=20000):
        '''
        Trace through the vector field using the quad tree.
        '''
        if self.method == 'euler' or self.method == 'eul':
            from spacepy.pybats.trace2d import trace2d_eul as trc
        elif self.method == 'rk4':
            from spacepy.pybats.trace2d import trace2d_rk4 as trc
        else:
            raise ValueError('Tracing method {} not recognized!'.format(
                self.method))

        # Get name of dimensions in order.
        grid = bats['grid'].attrs['dims']

        # Find starting block and set starting locations.
        block = bats.find_block(self.xstart, self.ystart)
        xfwd = [self.xstart]
        yfwd = [self.ystart]
        xnow, ynow = self.xstart, self.ystart

        # Trace forwards.
        while block:
            # Grab indices of all values inside current block.
            # Ghost cell check is for experimental testing:
            if hasattr(bats.qtree[block], 'ghost'):
                loc = bats.qtree[block].ghost
            # Typical approach is no ghost cell info:
            else:
                loc = bats.qtree[block].locs
            # Trace through this block.
            x, y = trc(bats[self.xvar][loc], bats[self.yvar][loc],
                       xnow, ynow, bats[grid[0]][loc][0, :],
                       bats[grid[1]][loc][:, 0], ds=0.01)
            # Update location and block:
            xnow, ynow = x[-1], y[-1]
            newblock = bats.find_block(xnow, ynow)
            # If we didn't leave the block, stop tracing.
            # Additionally, if inside rBody, stop.
            radnow = xnow**2+ynow**2
            if (block == newblock) or radnow < bats.attrs['rbody']*0.8:
                block = False
            elif newblock:
                block = newblock
            else:
                block = False
            # Append to full trace vectors.
            xfwd = np.append(xfwd, x[1:])
            yfwd = np.append(yfwd, y[1:])
            del x, y
            # It's possible to get stuck swirling around across
            # a few blocks.  If we spend a lot of time tracing,
            # call it quits.
            if xfwd.size > maxPoints:
                block = False

        # Trace backwards.  Same Procedure as above.
        block = bats.find_block(self.xstart, self.ystart)
        xbwd = [self.xstart]
        ybwd = [self.ystart]
        xnow, ynow = self.xstart, self.ystart
        while block:
            if hasattr(bats.qtree[block], 'ghost'):
                loc = bats.qtree[block].ghost
            else:
                loc = bats.qtree[block].locs
            x, y = trc(bats[self.xvar][loc], bats[self.yvar][loc],
                       xnow, ynow, bats[grid[0]][loc][0, :],
                       bats[grid[1]][loc][:, 0], ds=-0.01)
            xnow, ynow = x[-1], y[-1]
            newblock = bats.find_block(xnow, ynow)
            radnow = xnow**2+ynow**2
            if (block == newblock) or radnow < bats.attrs['rbody'] * .8:
                block = False
            elif newblock:
                block = newblock
            else:
                block = False
            # Append to full trace vectors.
            xbwd = np.append(x[::-1], xbwd)
            ybwd = np.append(y[::-1], ybwd)
            if xbwd.size > maxPoints:
                block = False

        # Trim duplicate points created when backwards tracing.
        # There's always at least 1 duplicate.
        for i in range(1, xbwd.size):
            if xbwd[-(i+1)]-xfwd[0] != 0:
                break

        # Combine foward and backward traces.
        self.x = np.append(xbwd[:-i], xfwd)
        self.y = np.append(ybwd[:-i], yfwd)

        # If planetary run w/ body:
        # 1) Check if line is closed to body.
        # 2) Trim points within body.
        if 'rbody' in bats.attrs:
            # Radial distance:
            r = np.sqrt(self.x**2.0 + self.y**2.0)
            # Closed field line?  Lobe line?  Set status:
            if (r[0] < bats.attrs['rbody']) and (r[-1] < bats.attrs['rbody']):
                self.open = False
                self.status = 'closed'
            elif (r[0] > bats.attrs['rbody']) and \
                 (r[-1] < bats.attrs['rbody']):
                self.open = True
                self.status = 'north lobe'
            elif (r[0] < bats.attrs['rbody']) and \
                 (r[-1] > bats.attrs['rbody']):
                self.open = True
                self.status = 'south lobe'
            # Trim the fat!
            limit = bats.attrs['rbody']*.8
            self.x, self.y = self.x[r > limit], self.y[r > limit]

    def trace(self, bats):
        '''
        Trace through the vector field.
        '''
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
        self.x = np.array(x2[::-1].tolist() + x1[1:].tolist())
        self.y = np.array(y2[::-1].tolist() + y1[1:].tolist())

        # Check if line is closed to body.
        if 'rbody' in bats.attrs:
            r1 = np.sqrt(self.x[0]**2.0 + self.y[0]**2.0)
            r2 = np.sqrt(self.x[-1]**2.0 + self.y[-1]**2.0)
            if (r1 < bats.attrs['rbody']) and (r2 < bats.attrs['rbody']):
                self.open = False

    def plot(self, ax, *args, **kwargs):
        '''
        Add streamline to axes object "ax".
        '''
        ax.plot(self.x, self.y, self.style, *args, **kwargs)


class Bats2d(IdlFile):
    '''
    A child class of :class:`~pybats.IdlFile` tailored to 2D BATS-R-US output.

    Calculations
    ------------
    New values can be added via the addition of new keys.  For example,
    a user could add radial distance to an equatorial Bats2d object as follows:

    >>> import numpy as np
    >>> from spacepy.pybats import bats
    >>> mhd = bats.Bats2d('tests/data/pybats_test/z=0_sine.out')
    >>> mhd['rad'] = np.sqrt( mhd['x']**2 + mhd['y']**2 )

    Note, however, that if the user switches the data frame in a *.outs file
    to access data from a different epoch, these values will need to be
    updated.

    An exception to this is built-in `calc_*` methods, which perform common
    MHD/fluid dynamic calculations (i.e., Alfven wave speed, vorticity, and
    more.)  These values are updated when the data frame is switched (see the
    `switch_frame` method).

    For example, to calculate number density and the fraction of the total
    density for each species in a multi-fluid simulation,

    >>> mhd = bats.Bats2d('tests/data/pybats_test/cut_multispecies.out')
    >>> mhd.calc_ndens()

    In an interactive session, users can tab-complete to explore the possible
    built-in calculation methods.

    Plotting
    --------
    While users can employ Matplotlib to plot values, a set of built-in
    methods are available to expedite plotting.  These are the
    `add_<plot type>` methods.  These methods always have the following
    keyword arguments that allow users to optionally build more complicated
    plots: *target* and *loc*.  The *target* kwarg tells the plotting method
    where to place the plot and can either be a Matplotlib figure or axes
    object.  If it's an axes object, *loc* sets the subplot location using
    the typical matplotlib syntax (e.g., `loc=121`).  The default behavior is
    to create a new figure and axes object.

    This approach allows a user to over-plot contours, field lines, and
    other plot artists as well as combine different subplots onto a single
    figure.  Continuing with our example above, let's plot the grid layout
    for our file as well as equatorial pressure and flow streamlines:

    >>> import matplotlib.pyplot as plt
    >>> fig = plt.Figure(figsize=(8,6))
    >>> mhd.add_grid_plot(target=fig, loc=121)
    >>> mhd.add_contour('x','y','p', target=fig, loc=122)
    >>> mhd.add_stream_scatter('ux', 'uy', target=fig, loc=122)

    Useful plotting methods include the following:

    | Plot Method        | Description                                    |
    | ------------------ | ---------------------------------------------- |
    | add_grid_plot      | Create a quick-look diagram of the grid layout |
    | add_contour        | Create a contour plot of a given variable      |
    | add_pcolor         | Add a p-color (no-interpolation contour) plot  |
    | add_stream_scatter | Scatter stream traces (any vector field)       |
    | add_b_magsphere    | Add magnetic field lines for X-Z plane cuts    |
    | add_planet         | Add a simple black/white planet at the origin  |
    | add_body           | Add an inner boundary at the origin            |

    Extracting and Stream Tracing
    -----------------------------
    Extracting values via interpolation to arbitrary points and creating
    stream traces through any vector field (e.g., velocity or magnetic field)
    are aided via the use of the following object methods:

    | Method     | Description                                        |
    | ---------- | -------------------------------------------------- |
    | extract    | Interpolate to arbitrary points and extract values |
    | get_stream | Integrate stream lines through vector fields       |


    Be sure to read the docstring information of :class:`~pybats.IdlFile` to
    see how to handle multi-frame files (*.outs) and for a list of critical
    attributes.

    .. versionchanged:: 0.5.0

       Unstructured data are now presented as in the files. See
       `~pybats.IdlFile` for details.
    '''
    # Init by calling IdlFile init and then building qotree, etc.
    def __init__(self, filename, *args, blocksize=8, **kwargs):

        # Stash blocksize as an attribute:
        self.blocksize = blocksize

        # Create quad tree object attribute:
        self._qtree = None

        # Read file.
        IdlFile.__init__(self, filename, keep_case=False, *args,
                         **kwargs)

        # Behavior of output files changed Jan. 2017:
        # Check for 'r' instead of 'rbody' in attrs.
        if 'r' in self.attrs and 'rbody' not in self.attrs:
            self.attrs['rbody'] = self.attrs['r']

    def switch_frame(self, *args, **kwargs):
        '''
        For files that have more than one data frame (i.e., `*.outs` files),
        load data from the *iframe*-th frame into the object replacing what is
        currently loaded.
        '''

        # Reset our quad tree before reading our file:
        self._qtree = None

        # Switch frames using parent method:
        super(Bats2d, self).switch_frame(*args, **kwargs)

    @property
    def qtree(self):
        if self._qtree is None:

            import spacepy.pybats.qotree as qo

            # Parse grid into quad tree.
            if self['grid'].attrs['gtype'] != 'Regular':
                xdim, ydim = self['grid'].attrs['dims'][0:2]
                try:
                    points = np.array([self[xdim], self[ydim]])
                    self._qtree = qo.QTree(points, blocksize=self.blocksize)
                except Exception:
                    from traceback import print_exc
                    print_exc()
                    self._qtree = False
                    self.find_block = lambda: False
            else:
                self._qtree = False
                self.find_block = lambda: False

        return self._qtree

    @property
    def find_block(self):
        if self.qtree:
            return self._qtree.find_leaf
        else:
            return lambda: False

    ####################
    # CALCULATIONS
    ####################

    @calc_wrapper
    def calc_temp(self, units='eV'):
        '''
        Calculate plasma temperature for each fluid.  Number density is
        calculated using *calc_ndens* if it hasn't been done so already.

        Temperature is obtained via density and pressure through the simple
        relationship P=nkT.

        Use the units kwarg to set output units.  Current choices are
        KeV, eV, and K.  Default is eV.
        '''

        from spacepy.datamodel import dmarray

        units = units.lower()

        # Create dictionary of unit conversions.
        conv = {'ev': 6241.50935,   # nPa/cm^3 --> eV.
                'kev': 6.24150935,  # nPa/cm^3 --> KeV.
                'k': 72429626.47}   # nPa/cm^3 --> K.

        # Calculate number density if not done already.
        if 'N' not in self:
            self.calc_ndens()

        # Find all number density variables.
        for key in list(self.keys()):
            # Next variable if not number density:
            if key[-1] != 'N':
                continue
            # Next variable if no matching pressure:
            if not key[:-1]+'p' in self:
                continue
            self[key[:-1]+'t'] = dmarray(
                conv[units] * self[key[:-1]+'p']/self[key],
                attrs={'units': units})

    @calc_wrapper
    def calc_b(self):
        '''
        Calculates total B-field strength using all three B components.
        Retains units of components.  Additionally, the unit vector
        b-hat is calculated as well.
        '''

        if 'b' in self:
            return

        self['b'] = np.sqrt(self['bx']**2 + self['by']**2 + self['bz']**2)
        self['b'].attrs = {'units': self['bx'].attrs['units']}

        self['bx_hat'] = self['bx'] / self['b']
        self['by_hat'] = self['by'] / self['b']
        self['bz_hat'] = self['bz'] / self['b']

        self['bx_hat'].attrs = {'units': 'unitless'}
        self['by_hat'].attrs = {'units': 'unitless'}
        self['bz_hat'].attrs = {'units': 'unitless'}

    @calc_wrapper
    def calc_j(self):
        '''
        Calculates total current density strength using all three J components.
        Retains units of components, stores in self['j']
        '''
        self['j'] = np.sqrt(self['jx']**2 + self['jy']**2.0 + self['jz']**2)
        self['j'].attrs = {'units': self['jx'].attrs['units']}

    @calc_wrapper
    def calc_uperp(self):
        '''
        Calculate the magnitude of the velocity perpendicular to the
        magnetic field: $\\vec{U} \\times \\hat{b}$.  Result maintains units
        of velocity.

        Values are calculated for each fluid and stored as self['u_perp']
        (using species prefixes as needed).

        Notes
        =====
        .. versionadded:: 0.5.0
        '''

        # Ensure b_hat is calculated:
        self.calc_b()

        # Get all fluid species variables.  Save in new list.
        species = []
        for k in self:
            if (k[-2:]) == 'ux':
                species.append(k[:-2])

        # Calculate perpendicular velocity.  Separate step to avoid
        # changing dictionary while looping over keys.
        for s in species:
            # Build components:
            ux = self[s+'uy']*self['bz_hat']-self[s+'uz']*self['by_hat']
            uy = self[s+'uz']*self['bx_hat']-self[s+'ux']*self['bz_hat']
            uz = self[s+'ux']*self['by_hat']-self[s+'uy']*self['bx_hat']

            # Get magnitude:
            self[s+'u_perp'] = np.sqrt(ux**2+uy**2+uz**2)

    @calc_wrapper
    def calc_upar(self):
        '''
        Calculate $\\vec{U} \\cdot \\vec{B}$ and store as 'upar' in *self*.
        Result maintains units of velocity.  Values are calculated for each
        fluid.

        Values arestored as self['u_par']; species prefixes are used for
        multi-species/fluid files.

        Notes
        =====
        .. versionadded:: 0.5.0
        '''

        # Ensure b_hat is calculated:
        self.calc_b()

        # Get all fluid species variables.  Save in new list.
        species = []
        for k in self:
            if (k[-2:]) == 'ux':
                species.append(k[:-2])

        # Calculate U dot B.  Separate step to avoid
        # changing dictionary while looping over keys.
        for s in species:
            self[s+'u_par'] = \
                self[s+'ux']*self['bx_hat']   \
                + self[s+'uy']*self['by_hat'] \
                + self[s+'uz']*self['bz_hat']

    @calc_wrapper
    def calc_E(self):
        '''
        Calculates the MHD electric field, -UxB.  Works for default
        MHD units of nT and km/s; if these units are not correct, an
        exception will be raised.  Stores E in mV/m.

        Values are saved as self['Ex'], self['Ey'], self['Ez'], and
        self['E'].
        '''

        # Some quick declarations for more readable code.
        ux, uy, uz = self['ux'], self['uy'], self['uz']
        bx, by, bz = self['bx'], self['by'], self['bz']

        # Check units.  Should be nT(=Volt*s/m^2) and km/s.
        if (bx.attrs['units'] != 'nT') or (ux.attrs['units'] != 'km/s'):
            raise Exception('Incorrect units!  Should be km/s and nT.')

        # Calculate; return in millivolts per meter
        self['Ex'] = -1.0*(uy*bz - uz*by) / 1000.0
        self['Ey'] = -1.0*(uz*bx - ux*bz) / 1000.0
        self['Ez'] = -1.0*(ux*by - uy*bx) / 1000.0
        self['Ex'].attrs = {'units': 'mV/m'}
        self['Ey'].attrs = {'units': 'mV/m'}
        self['Ez'].attrs = {'units': 'mV/m'}

        # Total magnitude.
        self['E'] = np.sqrt(self['Ex']**2+self['Ey']**2+self['Ez']**2)

    @calc_wrapper
    def calc_ndens(self):
        '''
        Calculate number densities for each fluid.  Species mass is ascertained
        via recognition of fluid name (e.g. OpRho is clearly oxygen).  A full
        list of recognized fluids/species can be found by exploring the
        dictionary *mass* found in :mod:`~spacepy.pybats.bats`.  Composition is
        also calculated as percent of total number density.

        New values are saved using the keys *speciesN* (e.g. *opN*) and
        *speciesFrac* (e.g. *opFrac*).
        '''

        # Use shared function.
        _calc_ndens(self)

    @calc_wrapper
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

        Values are stored as self['beta'].
        '''

        if 'b' not in self:
            self.calc_b()
        mu_naught = 4.0E2 * np.pi  # Mu_0 x unit conversion (nPa->Pa, nT->T)
        temp_b = self['b']**2.0
        temp_b[temp_b < 1E-8] = -1.0*mu_naught*self['p'][temp_b == 0.0]
        temp_beta = mu_naught*self['p']/temp_b
        self['beta'] = temp_beta
        self['beta'].attrs = {'units': 'unitless'}

    @calc_wrapper
    def calc_jxb(self):
        '''
        Calculates the JxB force assuming:
        -Units of J are uA/m2, units of B are nT.
        Under these assumptions, the value calculated is force density with
        units of nN/m^3.

        Values are stored as 'jb' (magnitude) and 'jbx', 'jby', 'jbz'.
        '''
        from spacepy.datamodel import dmarray

        # Unit conversion (nT, uA, cm^-3 -> nT, A, m^-3) to nN/m^3.
        conv = 1E-6
        # Calculate cross product, convert units.
        self['jbx'] = dmarray((self['jy']*self['bz'] -
                               self['jz']*self['by'])*conv,
                              {'units': 'nN/m^3'})
        self['jby'] = dmarray((self['jz']*self['bx'] -
                               self['jx']*self['bz'])*conv,
                              {'units': 'nN/m^3'})
        self['jbz'] = dmarray((self['jx']*self['by'] -
                               self['jy']*self['bx'])*conv,
                              {'units': 'nN/m^3'})
        self['jb'] = dmarray(np.sqrt(self['jbx']**2 +
                                     self['jby']**2 +
                                     self['jbz']**2), {'units': 'nN/m^3'})

    @calc_wrapper
    def calc_alfven(self):
        '''
        Calculate the Alfven speed, B/(mu*rho)^(1/2) in km/s.  This is
        performed for each species and the total fluid.
        The variable is saved under key "alfven" in self.data.
        '''
        from spacepy.datamodel import dmarray

        if 'b' not in self:
            self.calc_b()
        # M_naught * conversion from #/cm^3 to kg/m^3
        mu_naught = 4.0E-7 * np.pi * 1.6726E-27 * 1.0E6

        # Get all rho-like variables.  Save in new list.
        rho_names = []
        for k in self:
            if (k[-3:]) == 'rho':
                rho_names.append(k)

        # Calculate Alfven speed in km/s.  Separate step to avoid
        # changing dictionary while looping over keys.
        for k in rho_names:
            self[k[:-3]+'alfven'] = dmarray(self['b']*1E-12 /
                                            np.sqrt(mu_naught*self[k]),
                                            attrs={'units': 'km/s'})

    @calc_wrapper
    def _calc_divmomen(self):
        '''
        Calculate the divergence of momentum, i.e.
        $\\rho(u \\dot \\nabla)u$.
        This is currently exploratory.
        '''

        from spacepy.datamodel import dmarray
        from spacepy.pybats.batsmath import d_dx, d_dy

        if self.qtree is False:
            raise ValueError('calc_divmomen requires a valid qtree')

        # Create empty arrays to hold new values.
        size = self['ux'].shape
        self['divmomx'] = dmarray(np.zeros(size), {'units': 'nN/m3'})
        self['divmomz'] = dmarray(np.zeros(size), {'units': 'nN/m3'})

        # Units!
        c1 = 1000. / 6371.0  # km2/Re/s2 -> m/s2
        c2 = 1.6726E-21      # AMU/cm3 -> kg/m3
        c3 = 1E9             # N/m3 -> nN/m3.

        for k in self.qtree:
            # Calculate only on leafs of quadtree.
            if not self.qtree[k].isLeaf:
                continue

            # Extract values from current leaf.
            leaf = self.qtree[k]
            ux = self['ux'][leaf.locs]
            uz = self['uz'][leaf.locs]

            self['divmomx'][leaf.locs] = ux*d_dx(ux, leaf.dx) + \
                uz*d_dy(ux, leaf.dx)
            self['divmomz'][leaf.locs] = ux*d_dx(uz, leaf.dx) + \
                uz*d_dy(uz, leaf.dx)

        # Unit conversion.
        self['divmomx'] *= self['rho']*c1*c2*c3
        self['divmomz'] *= self['rho']*c1*c2*c3

    @calc_wrapper
    def calc_vort(self, conv=1000./RE):
        '''
        Calculate the vorticity (curl of bulk velocity) for the direction
        orthogonal to the cut plane.  For example, if output file is
        a cut in the equatorial plane (GSM X-Y plane), only the z-component
        of the curl is calculated.

        Output is saved as self['wD'], where D is the resulting dimension
        (e.g., 'wz' for the z-component of vorticity in the X-Y plane).

        Vorticity is calculated for total fluid and by-species for
        multi-fluid output files.

        Parameters
        ==========
        None

        Other Parameters
        ================
        conv : float
           Required unit conversion such that output units are 1/s.
           Defaults to 1/RE (in km), which assumes grid is in RE and
           velocity is in km/s.
        '''

        from spacepy.pybats.batsmath import d_dx, d_dy

        if self.qtree is False:
            raise ValueError('calc_vort requires a valid qtree')

        # Determine which direction to calculate based on what direction
        # is not present.  Save appropriate derivative operators, order
        # useful directions, create new variable name.
        dims = self['grid'].attrs['dims']
        if 'x' not in dims:
            w = 'wx'
            dim1, dim2 = 'z', 'y'
            dx1, dx2 = d_dy, d_dx
        elif 'y' not in dims:
            w = 'wy'
            dim1, dim2 = 'x', 'z'
            dx1, dx2 = d_dy, d_dx
        else:
            w = 'wz'
            dim1, dim2 = 'x', 'y'
            dx1, dx2 = d_dy, d_dx

        # Get all fluid species variables.  Save in new list.
        species = []
        for k in self:
            if (k[-2:]) == 'ux':
                species.append(k[:-2])

        # Calculate vorticity by species.
        for s in species:
            # Create new arrays to hold curl.
            size = self[s+'ux'].shape
            self[s+w] = dmarray(np.zeros(size), {'units': '1/s'})

            # Navigate quad tree, calculate curl at every leaf.
            for k in self.qtree:
                # Plot only leafs of the tree.
                if not self.qtree[k].isLeaf:
                    continue

                # Get location of points and extract velocity:
                leaf = self.qtree[k]
                u1 = self[s+'u'+dim1][leaf.locs]
                u2 = self[s+'u'+dim2][leaf.locs]

                # Calculate curl
                self[s+w][leaf.locs] = conv * (dx1(u1, leaf.dx) -
                                               dx2(u2, leaf.dx))

    @calc_wrapper
    def calc_gradP(self):
        '''
        Calculate the pressure gradient force.
        '''

        from spacepy.datamodel import dmarray
        from spacepy.pybats.batsmath import d_dx, d_dy

        if self.qtree is False:
            raise ValueError('calc_gradP requires a valid qtree')

        if 'p' not in self:
            raise KeyError('Pressure not found in object')

        # Create new arrays to hold pressure.
        dims = self['grid'].attrs['dims']
        size = self['p'].shape
        self['gradP'] = dmarray(np.zeros(size), {'units': 'nN/cm^3'})
        for d in dims:
            self['gradP_'+d] = dmarray(np.zeros(size), {'units': 'nN/m^3'})

        for k in self.qtree:
            # Plot only leafs of the tree.
            if not self.qtree[k].isLeaf:
                continue

            # Extract leaf; place pressure into 2D array.
            leaf = self.qtree[k]
            z = self['p'][leaf.locs]

            # Calculate derivatives; place into new dmarrays.
            # Unit conversion: P in nPa => gradP in nN/m3, dx=Re
            # Convert to nN/m3 by multiplying by 6378000m**-1.
            # Clever indexing maps values back to "unordered" array so that
            # values can be plotted like all others.
            conv = 1.0/-6378000.0
            self['gradP_'+dims[0]][leaf.locs] = d_dx(z, leaf.dx)*conv
            self['gradP_'+dims[1]][leaf.locs] = d_dy(z, leaf.dx)*conv

        # Scalar magnitude:
        for d in dims:
            self['gradP'] += self['gradP_'+d]**2
        self['gradP'] = np.sqrt(self['gradP'])

    @calc_wrapper
    def calc_utotal(self):
        '''
        Calculate bulk velocity magnitude: $u^2 = u_X^2 + u_Y^2 + u_Z^2$.
        This is done on a per-fluid basis.
        '''

        from spacepy.datamodel import dmarray

        species = []

        # Find all species, the variable names end in "ux".
        for k in self:
            if (k[-2:] == 'ux') and (k[:-2]+'u' not in self):
                species.append(k[:-2])

        units = self['ux'].attrs['units']

        for s in species:
            self[s+'u'] = dmarray(np.sqrt(self[s+'ux']**2 +
                                          self[s+'uy']**2 +
                                          self[s+'uz']**2),
                                  attrs={'units': units})

    @calc_wrapper
    def _calc_Ekin(self, units='eV'):
        '''
        Calculate average kinetic energy per particle using
        $E=\frac{1}{2}mv^2$.  Note that this is not the same as energy
        density.  Units are $eV$.
        '''

        from spacepy.datamodel import dmarray

        raise Warning("This calculation is unverified.")

        conv = 0.5 * 0.0103783625  # km^2-->m^2, amu-->kg, J-->eV.
        if units.lower == 'kev':
            conv = conv/1000.0

        species = []

        # Find all species, the variable names end in "rho".
        for k in self:
            if (k[-3:] == 'rho') and (k[:-3]+'Ekin' not in self):
                species.append(k[:-3])

        for s in species:
            # THIS IS WRONG HERE: 1/2mV**2?  Notsomuch.
            self[s+'Ekin'] = dmarray(np.sqrt(self[s+'ux']**2 +
                                             self[s+'uy']**2 +
                                             self[s+'uz']**2)
                                     * conv * mass[s.lower()],
                                     attrs={'units': units})

    def calc_all(self, exclude=[]):
        '''
        Perform all variable calculations (e.g. calculations that
        begin with 'calc_').  Any exceptions raised by functions that
        could not be peformed (typicaly from missing variables) are
        discarded.
        '''
        for command in dir(self):
            if (command[0:5] == 'calc_') and (command != 'calc_all') \
                 and (command not in exclude):
                try:
                    eval('self.'+command+'()')
                except AttributeError:
                    raise Warning('Did not perform {0}: {1}'.format(
                        command, sys.exc_info()[0]))

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

        if self.gridtype != 'Regular':
            if not cellsize:
                raise ValueError('Grid must be regular or ' +
                                 'cellsize must be given.')
            self.regrid(cellsize, dim1range=dim1range, dim2range=dim2range)

        # Order our dimensions alphabetically.
        newvars = []
        for key in sorted(self.grid.keys()):
            newvars.append('gradp'+key)

        dx = self.resolution * 6378000.0  # RE to meters
        p = self['p']*10E-9               # nPa to Pa
        self[newvars[0]], self[newvars[1]] = np.gradient(p, dx, dx)
        self['gradp'] = np.sqrt(self[newvars[0]]**2.0 + self[newvars[1]]**2.0)

    def cfl(self, dt, xcoord='x', ycoord='z'):
        """
        Calculate the CFL number in each cell given time step dt.

        Result is stored in self['cfl'].
        """

        if 'u' not in self:
            self.calc_utotal()

        cfl = np.zeros(self['u'].shape)

        for key in list(self.qtree.keys()):
            child = self.qtree[key]
            if not child.isLeaf:
                continue
            pts = np.where(np.logical_and.reduce((self[xcoord] >= child.lim[0],
                                                  self[xcoord] <= child.lim[1],
                                                  self[ycoord] >= child.lim[2],
                                                  self[ycoord] <= child.lim[3])
                                                 ))

            cfl[pts] = self['u'][pts]*dt/child.dx

        self['cfl'] = dmarray(cfl, attrs={'units': ''})

    def vth(self, m_avg=3.1):
        """
        Calculate the thermal velocity. m_avg denotes the average ion
        mass in AMU.

        Result is stored in self['vth'].
        """

        m_avg_kg = m_avg*1.6276e-27
        ndensity = self['rho']/m_avg*1e6
        self['vth'] = dmarray(np.sqrt(self['p']*1e-9/ndensity/(m_avg_kg))
                              / 1000, attrs={'units': 'km/s'})

    def gyroradius(self, velocities=('u', 'vth'), m_avg=3.1):
        """
        Calculate the ion gyroradius in each cell.

        velocities to use in calculating the gyroradius are listed by name in
        the sequence argument velocities. If more than one variable is given,
        they are summed in quadrature.

        m_avg denotes the average ion mass in AMU.

        Result is stored in self['gyroradius']
        """

        if 'vth' in velocities:
            try:
                self['vth']
            except KeyError:
                self.vth()

        if 'u' in velocities and 'u' not in self:
            self.calc_utotal()

        velocities_squared_sum = self[velocities[0]]**2

        for vname in velocities[1:]:
            velocities_squared_sum += self[vname]**2

        v = np.sqrt(velocities_squared_sum)*1000

        if 'b' not in self:
            self.calc_b()

        B = self['b']*1e-9

        m_avg_kg = m_avg*1.6276e-27

        q = 1.6022e-19

        self['gyroradius'] = dmarray(m_avg_kg*v/(q*B)/6378000,
                                     attrs={'units': 'Re'})

    def plasma_freq(self, m_avg=3.1):
        """
        Calculate the ion plasma frequency.

        m_avg denotes the average ion mass in AMU.

        Result is stored in self['plasm_freq'].
        """

        m_avg_kg = m_avg*1.6276e-27
        ndensity = self['rho']/m_avg*1e6
        q = 1.6022e-19
        self['plasma_freq'] = dmarray(np.sqrt(4*np.pi*ndensity*q**2/m_avg_kg),
                                      attrs={'units': 'rad/s'})

    def inertial_length(self, m_avg=3.1):
        """
        Calculate the ion inertial length.

        m_avg denotes the average ion mass in AMU.

        Result is stored in self['inertial_length'].
        """

        try:
            self['plasma_freq']
        except KeyError:
            self.plasma_freq(m_avg=m_avg)

        if 'alfven' not in self:
            self.calc_alfven()

        self['inertial_length'] = dmarray(self['alfven']/self['plasma_freq']
                                          / 6378000, attrs={'units': 'Re'})

    def regrid(self, cellsize=1.0, dim1range=-1, dim2range=-1, debug=False):
        '''
        Re-bin data to regular grid of spacing cellsize.  Action is
        performed on all data entries in the bats2d object.

        '''
        from matplotlib.mlab import griddata

        if self['grid'].attrs['gtype'] == 'Regular':
            return

        # Order our dimensions alphabetically.
        dims = self['grid'].attrs['dims']
        if debug:
            print("Ordered dimensions: ", dims)

        # Check to see if dimranges are 2-element lists.
        # If not, either set defaults or raise exceptions.
        if dim1range == -1:
            dim1range = [self[dims[0]].min(), self[dims[0]].max()]
        else:
            if isinstance(dim1range, (type(()), type([]))):
                if len(dim1range) != 2:
                    raise ValueError('dim1range must have two elements!')
            else:
                raise TypeError('dim1range must be a tuple or list!')
        if dim2range == -1:
            dim2range = [self[dims[1]].min(), self[dims[1]].max()]
        else:
            if isinstance(dim2range, (type(()), type([]))):
                if len(dim2range) != 2:
                    raise ValueError('dim2range must have two elements!')
            else:
                raise TypeError('dim2range must be a tuple or list!')

        if debug:
            print('%s range = %f, %f' % (dims[0], dim1range[0], dim1range[1]))
            print('%s range = %f, %f' % (dims[1], dim2range[0], dim2range[1]))

        # Now, Regrid.
        grid1 = np.arange(dim1range[0], dim1range[1]+cellsize, cellsize)
        grid2 = np.arange(dim2range[0], dim2range[1]+cellsize, cellsize)

        for key in self:
            # Skip grid-type entries.
            if key in (self['grid'].attrs['dims']+['grid']):
                continue
            self[key] = griddata(self[dims[0]], self[dims[1]],
                                 self[key], grid1, grid2)

        # Change grid, gridtype, gridsize, and npoints to match new layout.
        self['grid'].attrs['gtype'] = 'Regular'
        self['grid'].attrs['npoints'] = len(grid1) * len(grid2)
        self['grid'].attrs['resolution'] = cellsize
        self[dims[0]] = grid1
        self[dims[1]] = grid2

    #############################
    # EXTRACTION/INTERPOLATION
    #############################
    def extract(self, x, y, **kwargs):
        '''
        For x, y of a 1D curve, extract values along that curve
        and return slice as a new :class:`~spacepy.pybats.bats.Extraction`
        object.  Valid keyword arguments are the same as for
        :class:`~spacepy.pybats.bats.Extraction`.
        '''

        return Extraction(x, y, self, **kwargs)

    ######################
    # TRACING TOOLS
    ######################
    def get_stream(self, x, y, xvar, yvar, method='rk4', style='mag',
                   maxPoints=40000, extract=False):
        '''
        Trace a 2D streamline through the domain, returning a Stream
        object to the caller.

        x and y set the starting point for the tracing.

        xvar and yvar are string keys to self.data that define the
        vector field through which this function traces.

        The method kwarg sets the numerical method to use for the
        tracing.  Default is Runge-Kutta 4 (rk4).
        '''
        startvals = [x, y]
        dims = self['grid'].attrs['dims']
        for v, d in zip(startvals, dims):
            if v < self[d].min() or v > self[d].max():
                raise ValueError('Start value {} out of range for variable {}.'
                                 .format(v, d))
        stream = Stream(self, x, y, xvar, yvar, style=style,
                        maxPoints=maxPoints, method=method, extract=extract)

        return stream

    ######################
    # VISUALIZATION TOOLS
    ######################
    def add_grid_plot(self, target=None, loc=111, do_label=True,
                      do_fill=True, show_nums=False, show_borders=True,
                      cmap='jet_r', title='BATS-R-US Grid Layout'):
        '''
        Create a plot of the grid resolution by coloring regions of constant
        resolution.  Kwarg "target" specifies where to place the plot and can
        be a figure, an axis, or None.  If target is a figure, a new subplot
        is created.  The subplot location is set by the kwarg "loc", which
        defaults to 111.  If target is an axis, the plot is placed into that
        axis object.  If target is None, a new figure and axis are created
        and used to display the plot.

        Resolution labels can be disabled by setting kwarg do_label to False.

        Plot title is set using the 'title' kwarg, defaults to 'BATS-R-US
        Grid Layout'.

        Note that if target is not an axis object, the axis will automatically
        flip the direction of positive X GSM and turn on equal aspect ratio.
        In other words, if you want a customized axis, it's best to create
        one yourself.

        Figure and axis, even if none given, are returned.

         Parameters
        ----------
        target : Matplotlib Figure or Axes object
           Set plot destination.  Defaults to new figure.
        loc : 3-digit integer
           Set subplot location.  Defaults to 111.
        do_label : bool
           Adds resolution legend to righthand margin.  Defaults to True.
        do_fill : bool
           If true, fill blocks with color indicating grid resolution.
        show_nums : bool
           Adds quadtree values to each region.  Useful for debugging.
           Defaults to False.
        show_borders : True
           Show black borders around each quadtree leaf.  Defaults to True.
        cmap : string
           Sets the color map used to color each region.  Must be a Matplotlib
           named colormap.  Defaults to 'jet_r'.
        title : string
           Sets the title at the top of the plot.  Defaults to
           'BATS-R-US Grid Layout'.
        '''

        if self['grid'].attrs['gtype'] == 'Regular':
            raise ValueError('Function not compatable with regular grids')

        # Get dimensions over which we shall plot.
        xdim, ydim = self['grid'].attrs['dims'][0:2]

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)
        ax.set_aspect('equal')

        # Set plot range based on quadtree.
        ax.set_xlim([self.qtree[1].lim[0], self.qtree[1].lim[1]])
        ax.set_ylim([self.qtree[1].lim[2], self.qtree[1].lim[3]])
        # Plot.

        if show_borders:
            for key in list(self.qtree.keys()):
                self.qtree[key].plotbox(ax)
        self.qtree.plot_res(ax, tag_leafs=show_nums, do_label=do_label,
                            cmap=cmap, do_fill=do_fill)

        # Customize plot.
        ax.set_xlabel('GSM %s' % xdim.upper())
        ax.set_ylabel('GSM %s' % ydim.upper())
        ax.set_title(title)
        if 'rbody' in self.attrs or 'body' in self.attrs:
            self.add_body(ax)

        return fig, ax

    def add_stream_scatter(self, xcomp, ycomp, nlines=100, target=None,
                           loc=111, method='rk4', maxPoints=1E6, xlim=None,
                           ylim=None, narrow=0, arrsize=12, arrstyle='->',
                           start_points=None, **kwargs):
        '''
        Add a set of stream traces to a figure or axes that are distributed
        evenly but randomly throughout the plot domain.

        Lines will be seeded randomly over a given spatial range given by
        `xlim` and `ylim` **OR** the range of the axes (if `target` is set to
        a non-empty axes object) **OR** over the entire object domain (in that
        order).

        Extra keyword args are handed to matplotlib's LineCollection object:
        :class:`matplotlib.collections.LineCollection`.

        Parameters
        ----------
        xcomp : string
             The first component of the vector field to trace (e.g., 'bx').
        ycomp : string
             The second component of the vector field to trace (e.g., 'bz').

        Other Parameters
        ----------------
        target : Matplotlib Figure or Axes object
            Set plot destination.  Defaults to new figure.
        loc : 3-digit integer
            Set subplot location.  Defaults to 111.
        xlim : Two-element list/tuple
            Set the range in 1st dimension over which lines will be seeded.
        ylim : Two-element list/tuple
            Set the range in 2nd dimension over which lines will be seeded.
        nlines : int
            Number of stream lines to create; default is 100.
        start_points : nlinesx2 array
            Set start_points to define starting location of traces instead of
            using random points.  This is useful for creating timeseries of
            plots.
        maxPoints : int
            Set the maximum number of points in a single trace.
            Defaults to 1E6.
        narrow : int
            Add "n" arrows to each line to indicate direction.  Default is
            zero, or no lines.  If narrow=1, arrows will be placed at
            *start_points*.
        arrstyle : string
            Set the arrow style in the same manner as Matplotlib's
            annotate function.  Default is '->'.
        arrsize : int
            Set the size, in points, of each directional arrow.  Default is 12.

        Returns:
        --------
        fig : matplotlib Figure object
        ax  : matplotlib Axes object
        collect : matplotlib Collection object of trace results
        start_points : nlines x 2 array of line starting points
        '''

        from numpy.random import sample
        from matplotlib.collections import LineCollection
        from spacepy.plot import add_arrows

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)

        # Try to determine the order of the dimensions used:
        dims = self['grid'].attrs['dims']  # Default to standard order.
        flip = False
        for letters in zip(xcomp, ycomp):  # Check for reverse order.
            if letters == dims[::-1]:
                flip = True  # flip it!

        # Set if using the axes limits is a viable option for setting limits
        # for region over which to seed lines:
        use_ax_lims = False
        if type(target) == type(ax):
            # Is this a new, empty axes object? Check for any relevant artists:
            use_ax_lims = bool(ax.texts) or bool(ax.containers) \
                or bool(ax.lines) or bool(ax.patches) or bool(ax.images) \
                or bool(ax.artists)

        # Set range over which to place lines.  Use keyword values if provided
        # OR subset of axes ranges that fit in domain (if axes are reasonable)
        if xlim is None:
            xlim = [self[dims[0]].min(), self[dims[0]].max()]
            if use_ax_lims:
                axlim = ax.get_xlim()
                xlim = [max(xlim[0], axlim[0]), min(xlim[1], axlim[1])]
        if ylim is None:
            ylim = [self[dims[1]].min(), self[dims[1]].max()]
            if use_ax_lims:
                axlim = ax.get_ylim()
                ylim = [max(ylim[0], axlim[0]), min(ylim[1], axlim[1])]

        # If initial source points not given, create a random set:
        if start_points is None:
            # Get random points.
            start_points = sample([nlines, 2])
            # Scale to limits:
            start_points[:, 0] = (xlim[1]-xlim[0])*start_points[:, 0]+xlim[0]
            start_points[:, 1] = (ylim[1]-ylim[0])*start_points[:, 1]+ylim[0]
        else:
            nlines = start_points.shape[-1]

        # Extract stream traces, organize x and y coords:
        lines = []
        for xstart, ystart in start_points:
            # Some index errors crop up from time to time.
            # While a better solution is dug up, we use a try
            # block for the time being.
            try:
                stream = self.get_stream(xstart, ystart, xcomp, ycomp,
                                         method=method, maxPoints=maxPoints)
            except IndexError:
                continue

            iflip = 1-2*flip
            lines.append(np.array([stream.x, stream.y][::iflip]).transpose())

        # Create line collection & plot.
        collect = LineCollection(lines, **kwargs)
        ax.add_collection(collect)

        # Set the plot limits to match
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        # Add arrows if requested.  If one point per line, use starting points.
        # Otherwise, distribute arrows along lines.
        if narrow == 1:
            add_arrows(collect, n=narrow, size=arrsize,
                       positions=start_points, style=arrstyle)
        elif narrow > 1:
            add_arrows(collect, n=narrow, size=arrsize, style=arrstyle)

        return fig, ax, collect, start_points

    def find_earth_lastclosed(self, tol=np.pi/360., method='rk4',
                              max_iter=100, debug=False):
        '''
        For Y=0 cuts, attempt to locate the last-closed magnetic field line
        for both day- and night-sides.  This is done using a bisection
        approach to precisely locate the transition between open and closed
        geometries.  The method stops once this transition is found within
        a latitudinal tolerance of *tol*, which defaults to $\\pi/360.$, or
        one-half degree.  The tracing *method* can be set via keyword and
        defaults to 'rk4' (4th order Runge Kutta, see
        :class:`~spacepy.pybats.bats.Stream` for more information).
        The maximum number of iterations the algorithm will take is set
        by *max_iter*, which defaults to 100.  Latitudinal footprints of the
        last closed field lines at the inner boundary (not the ionosphere!)
        are also returned.

        This method returns 5 objects:

        * The dipole tilt in radians
        * A tuple of the northern/southern hemisphere polar angle of
          footpoints for the dayside last-closed field line.
        * A tuple of the northern/southern hemisphere polar angle of
          footpoints for the nightside last-closed field line.
        * The dayside last-closed field line as a
          :class:`~spacepy.pybats.bats.Stream` object.
        * The nightside last-closed field line as a
          :class:`~spacepy.pybats.bats.Stream` object.

        In each case, the angle is defined as elevation from the positive
        x-axis, in radians.

        '''

        # Get the dipole tilt by tracing a field line near the inner
        # boundary.  Find the max radial distance; tilt angle == angle off
        # equator of point of min R (~=max |B|).
        x_small = self.attrs['rbody']*-1.2  # check nightside.
        stream = self.get_stream(x_small, 0, 'bx', 'bz', method=method)
        r = stream.x**2 + stream.y**2
        loc = r == r.max()
        tilt = np.arctan(stream.y[loc]/stream.x[loc])[0]

        if debug:
            print('Dipole is tilted {} degress above the z=0 plane.'.format(
                tilt*180./np.pi))

        # Dayside- start by tracing from plane of min |B| and perp. to that:
        R = self.attrs['rbody']*1.15
        s1 = self.get_stream(R*np.cos(tilt), R*np.sin(tilt), 'bx', 'bz',
                             method=method)

        # Get initial angle and step.
        theta = tilt
        dTheta = np.pi/4.  # Initially, search 90 degrees.
        nIter = 0
        while (dTheta > tol) or (s1.open):
            nIter += 1

            # Are we closed or open?  Day or nightside?
            closed = not s1.open     # open or closed?
            isNig = s1.x.mean() < 0  # line on day or night side?
            isDay = not isNig

            # Adjust the angle towards the open-closed boundary.
            theta += (closed and isDay) * dTheta  # adjust nightwards.
            theta -= (s1.open or isNig) * dTheta  # adjust daywards.
            # Trace at the new theta to further restrict angular range:
            s1 = self.get_stream(R*np.cos(theta), R*np.sin(theta), 'bx', 'bz',
                                 method=method)
            # Reduce angular step:
            dTheta /= 2.
            if nIter > max_iter:
                if debug:
                    print('Did not converge before reaching max_iter')
                break

        # Possible to land on open or nightside line.
        # If this happens, inch back to dayside.
        isNig = s1.x.mean() < 0
        while (s1.open or isNig):
            theta -= tol / 2  # inch daywards.
            s1 = self.get_stream(R*np.cos(theta), R*np.sin(theta), 'bx', 'bz',
                                 method=method)
            isNig = s1.x.mean() < 0

        # Use last line to get southern hemisphere theta:
        npts = int(s1.x.size/2)
        r = np.sqrt(s1.x**2+s1.y**2)  # Distance from origin.
        # This loc finds the point(s) nearest to Rbody.
        loc = np.abs(r-self.attrs['rbody']) == np.min(np.abs(
            r[:npts]-self.attrs['rbody']))  # point closest to IB.
        xSouth, ySouth = s1.x[loc], s1.y[loc]
        # "+ 0" syntax is to quick-copy object.
        theta_day = [theta+0, 2*np.pi+np.arctan(ySouth/xSouth)[0]+0]
        day = s1

        # Nightside: Use more points in tracing (lines are long!)
        theta += tol/2.0  # Nudge nightwards.

        # Set dTheta to half way between equator and dayside last-closed:
        dTheta = (np.pi+tilt-theta)/2.

        nIter = 0
        while (dTheta > tol) or (s1.open):
            nIter += 1

            s1 = self.get_stream(R*np.cos(theta), R*np.sin(theta), 'bx', 'bz',
                                 method=method, maxPoints=1E6)
            # Closed?  Nightside?
            closed = not s1.open
            isNig = s1.x.mean() < 0
            isDay = not isNig

            theta -= (closed and isNig) * dTheta  # closed? move poleward.
            theta += (s1.open or isDay) * dTheta  # open?   move equatorward.

            # Don't cross over into dayside territory.
            if theta < theta_day[0]:
                theta = theta_day[0]+tol
                s1 = self.get_stream(R*np.cos(theta), R*np.sin(theta),
                                     'bx', 'bz', method=method, maxPoints=1E6)
                if debug:
                    print('No open flux over polar cap.')
                break

            dTheta /= 2.
            if nIter > max_iter:
                if debug:
                    print('Did not converge before reaching max_iter')
                break

        # Use last line to get southern hemisphere theta:
        npts = int(s1.x.size/2)  # Similar to above for dayside.
        r = np.sqrt(s1.x**2+s1.y**2)
        loc = np.abs(r-self.attrs['rbody']) == np.min(np.abs(
            r[:npts]-self.attrs['rbody']))
        xSouth, ySouth = s1.x[loc], s1.y[loc]
        theta_night = [theta+0, 2*np.pi+np.arctan2(ySouth, xSouth)[0]+0]
        night = s1
        # plt.plot(s1.x, s1.y, 'r-')

        return tilt, theta_day, theta_night, day, night

    def add_b_magsphere(self, target=None, loc=111,  style='mag',
                        DoLast=True, DoOpen=True,
                        compX='bx', compY='bz', narrow=0, arrsize=12,
                        method='rk4', tol=np.pi/720., DoClosed=True,
                        colors=None, linestyles=None, maxPoints=1E6,
                        nOpen=5, nClosed=15, arrstyle='->',
                        add_body=False, **kwargs):
        '''
        Create an array of field lines closed to the central body in the
        domain.  Add these lines to Matplotlib target object *target*.
        If no *target* is specified, a new figure and axis are created.

        Note that this should currently only be used for GSM y=0 cuts
        of the magnetosphere.

        A tuple containing the figure, axes, and LineCollection object
        is returned.

        Basic styling (color and linestyle) can be handled with the
        *style*, *colors*, and *linestyles* kwargs.  *style* can accept
        style names as defined in :class:`~spacepy.pybats.bats.Stream`, which
        colors and styles lines based on characteristics (e.g., open, closed).
        The default is 'mag', which colors open lines black and closed lines
        white.  Alternatively, this kwarg works in a similar manner as
        it does in :function:`~matplotlib.pyplot.plot`,
        i.e., a string code such as "b-" (a solid blue line) or 'r:' (a
        dotted red line), etc.  Both *colors* and *linestyles* work much
        as they do for :class:`~matplotlib.collections.LineCollection`, but
        only a single value (not a list or tuple) should be provided.
        *colors* can be a CSS4 color name, an RGB tuple, or a string hex code.
        *linestyles* can be the name of the style (e.g., "dashed") or a
        shortcut compatable with the *style* kwarg (e.g., "--").  See the
        documentation for the associated Matplotlib classes & functions to
        see all options.  Note that *linestyles* and *colors* override
        *style*.

        If the styling kwargs are used, they will set the colors for all
        lines except last-closed boundaries.  Users may control groups
        individually using multiple calls and plotting one group at a time.
        Note that *colors* and *linestyles* kwargs will override *style*;
        *colors* allows for more flexibility concerning color choice.

        Algorithm:  This method, unlike its predecessor, starts by finding
        the last closed field lines via
        :func:`~spacepy.pybats.bats.Bats2d.find_earth_lastclosed`.  It then
        fills the regions between the open and closed regions.  Currently, it
        does not treat purely IMF field lines.

        ========== ===========================================================
        Kwarg      Description
        ========== ===========================================================
        target     The figure or axes to place the resulting lines.
        style      The color coding system for field lines.  Defaults to 'mag'.
                   See :class:`spacepy.pybats.bats.Stream`.  Because lines are
                   added as a :class:`~matplotlib.collections.LineCollection`,
                   only certain styles are allowed (i.e., line styles only,
                   no marker styles).
        loc        The location of the subplot on which to place the lines.
        DoLast     Plot last-closed lines as red lines.  Defaults to **True**.
        DoOpen     Plot open field lines.  Defaults to **True**.
        DoClosed   Plot closed field lines.  Defaults to **True**.
        nOpen      Number of closed field lines to trace per hemisphere.
                   Defaults to 5.
        nClosed    Number of open field lines to trace per hemisphere.
                   Defaults to 15.
        narrow     Add "n" arrows to each line to indicate direction.
                   Default is zero, or no arrows.
        arrstyle   Set the arrow style in the same manner as Matplotlib's
                   annotate function.  Default is '->'.
        arrsize    Set the size, in points, of each directional arrow.
                   Default is 12.
        method     The tracing method; defaults to 'rk4'.   See
                   :class:`spacepy.pybats.bats.Stream`.
        tol        Tolerance for finding open-closed boundary; see
                   :func:`~spacepy.pybats.bats.Bats2d.find_earth_lastclosed`.
        compX      Name of x-variable through which to trace, defaults to 'bx'.
        compY      Name of y-variable through which to trace, defaults to 'bz'.
        colors     Matplotlib-compatable color name (single) to apply to lines.
        maxPoints  Set the maximum number of points in an field line
                   integration.  Defaults to one million.
        linestyles A single line style indicator, defaults to '-';
                   see :class:`~matplotlib.collections.LineCollection` for
                   possible options.
        add_body   Add an planet/body to the figure. Defaults to False.
        ========== ===========================================================

        Extra kwargs are passed to Matplotlib's LineCollection class as
        described above.

        Returns
        =======
        fig : matplotlib Figure object
        ax  : matplotlib Axes object
        collect : matplotlib Collection object of trace results

        Examples
        ========
        >>> import matplotlib.pyplot as plt
        >>> from spacepy.pybats import bats
        >>> # Open a 2D slice, add a pressure contour.
        >>> # Example file in spacepy/tests/data/pybats_test/:
        >>> mhd = bats.Bats2d('./y0_binary.out')
        >>> mhd.add_contour('x','z','p')
        >>> # Add field lines using default styling:
        >>> mhd.add_b_magsphere(target=plt.gca())
        >>> # Add a subset of lines using custom styling:
        >>> mhd.add_b_magsphere(target=plt.gca(), DoLast=False, DoOpen=False, style='g--')

        '''
        import re
        from matplotlib.collections import LineCollection
        from spacepy.plot import add_arrows

        # Set ax and fig based on given target.
        adj_lims = not target  # If no target set, adjust axes limits.
        fig, ax = set_target(target, figsize=(10, 10), loc=111)
        if add_body:
            self.add_body(ax)

        # Lines, colors, and styles:
        lines = []
        cols = []

        # Try to get line style from "style" string.
        # Default to regular line if not successful.
        if not linestyles:
            lstyle = re.sub(r'\w', '', style)
            if not lstyle:
                linestyles = '-'
            else:
                linestyles = lstyle

        # Start by finding open/closed boundary.
        tilt, thetaD, thetaN, last1, last2 = self.find_earth_lastclosed(
            method=method, tol=tol)

        # Useful parameters for the following traces:
        R = self.attrs['rbody']
        dTheta = 1.5*np.pi/180.
        dThetaN = .05*np.abs(thetaN[0]-thetaD[0])
        dThetaS = .05*np.abs(thetaN[1]-thetaD[1])

        # Do closed field lines #
        if DoClosed:
            for tDay, tNit in zip(
                    np.linspace(0,     thetaD[0]-dTheta, nClosed),
                    np.linspace(np.pi, thetaN[1]-dTheta, nClosed)):
                x, y = R*np.cos(tDay), R*np.sin(tDay)
                sD = self.get_stream(x, y, compX, compY, method=method,
                                     maxPoints=maxPoints, style=style)
                x, y = R*np.cos(tNit), R*np.sin(tNit)
                sN = self.get_stream(x, y, compX, compY, method=method,
                                     maxPoints=maxPoints, style=style)
                # Append to lines, colors.
                lines.append(np.array([sD.x, sD.y]).transpose())
                lines.append(np.array([sN.x, sN.y]).transpose())
                cols.append(sD.color)
                cols.append(sN.color)

        # Do open field lines
        if DoOpen:
            for tNorth, tSouth in zip(
                    np.linspace(thetaD[0]+dThetaN, thetaN[0]-dThetaN, nOpen),
                    np.linspace(thetaN[1]+dThetaS, thetaD[1]-dThetaS, nOpen)):
                x, y = R*np.cos(tNorth), R*np.sin(tNorth)
                sD = self.get_stream(x, y, compX, compY, method=method,
                                     maxPoints=maxPoints, style=style)
                x, y = R*np.cos(tSouth), R*np.sin(tSouth)
                sN = self.get_stream(x, y, compX, compY, method=method,
                                     maxPoints=maxPoints, style=style)
                # Append to lines, colors.
                lines.append(np.array([sD.x, sD.y]).transpose())
                lines.append(np.array([sN.x, sN.y]).transpose())
                cols.append(sD.color)
                cols.append(sN.color)

        # Finalize Collection
        # If colors is given, replace what is given from
        # individual lines.  Keep the list-approach, however.
        if colors:
            cols = [colors]*len(cols)

        # Add last-closed field lines at end so they are plotted "on top".
        if DoLast:
            lines += [np.array([last1.x, last1.y]).transpose(),
                      np.array([last2.x, last2.y]).transpose()]
            cols += 2 * ['r']

        # Create line collection & plot.
        collect = LineCollection(lines, colors=cols, linestyles=linestyles,
                                 **kwargs)
        ax.add_collection(collect)

        # Add lines if required:
        if narrow > 0:
            add_arrows(collect, n=narrow, size=arrsize, style=arrstyle)

        # On fresh axes, adjust limits from default ([0,1]):
        if adj_lims:
            # Set defaults:
            xlim, ylim = [0, 1], [0, 1]

            # Get x,y locations along each line:
            points = [path.vertices for path in collect.get_paths()]

            # Get max/min from each line, update lims:
            for p in points:
                xlim = min(xlim[0], p.min(0)[0]), max(xlim[1], p.max(0)[0])
                ylim = min(ylim[0], p.min(0)[1]), max(ylim[1], p.max(0)[1])

            # Convert to arrays for element arithmatic:
            xlim, ylim = np.array(xlim), np.array(ylim)

            # Add a buffer:
            dX, dY = min(5, xlim[1]-xlim[0]), min(5, ylim[1]-ylim[0])
            xlim += (-dX, dX)
            ylim += (-dY, dY)

            # Set new axes limits:
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        return fig, ax, collect

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

    def add_pcolor(self, dim1, dim2, value, zlim=None, target=None, loc=111,
                   title=None, xlabel=None, ylabel=None,
                   ylim=None, xlim=None, add_cbar=False, clabel=None,
                   add_body=True, dolog=False, *args, **kwargs):
        '''
        Create a pcolor plot of variable **value** against **dim1** on the
        x-axis and **dim2** on the y-axis.  Pcolor plots shade each
        computational cell with the value at the cell center.  Because no
        interpolation or smoothing is used in the visualization, pcolor plots
        are excellent for examining the raw output.

        Simple example:

        >>> self.add_pcolor('x', 'y', 'rho')

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        =========== ==========================================================
        Kwarg       Description
        =========== ==========================================================
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        title       Sets title of axes.  Default is no title.
        xlabel      Sets x label of axes.  Defaults to **dim1** and units.
        ylabel      Sets y label of axes.  Defaults to **dim2** and units.
        xlim        Sets limits of x-axes.  Defaults to whole domain.
        ylim        Sets limits of y-axes.  Defaults to whole domain.
        zlim        Sets color bar range.  Defaults to variable max/min.
        add_cbar    Adds colorbar to plot.  Defaults to *False*.
        clabel      Sets colorbar label.  Defaults to **var** and units.
        add_body    Places planetary body in plot.  Defaults to **True**.
        dolog       Sets use of logarithmic scale.  Defaults to **False**.
        =========== ==========================================================
        '''

        import numbers
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, Normalize

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)

        # Get max/min if none given.
        if zlim is None:
            zlim = [0, 0]
            zlim[0] = self[value].min()
            zlim[1] = self[value].max()
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Logarithmic scale?
        if dolog:
            z = np.where(self[value] > zlim[0], self[value], 1.01*zlim[0])
            norm = LogNorm(vmin=zlim[0], vmax=zlim[1])
        else:
            z = self[value]
            norm = Normalize(vmin=zlim[0], vmax=zlim[1])

        if self['grid'].attrs['gtype'] == 'Regular':
            pass
        else:
            # Indices corresponding to QTree dimensions:
            ix = self['grid'].attrs['dims'].index(dim1)
            iy = self['grid'].attrs['dims'].index(dim2)
            for k in self.qtree:
                # Plot only leafs of the tree.
                if not self.qtree[k].isLeaf:
                    continue
                leaf = self.qtree[k]
                x = leaf.cells[ix]
                y = leaf.cells[iy]
                z_local = z[leaf.locs]
                pcol = ax.pcolormesh(x, y, z_local, norm=norm, **kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(pcol, ax=ax, pad=0.01)
            if clabel is None:
                if 'units' in self[value].attrs:
                    clabel = f"{value}, ({self[value].attrs['units']})"
                else:
                    clabel = f"{value}"
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Set title, labels, axis ranges (use defaults where applicable.)
        if title:
            ax.set_title(title)
        if ylabel is None:
            ylabel = '%s ($R_{E}$)' % dim2.upper()
        if xlabel is None:
            xlabel = '%s ($R_{E}$)' % dim1.upper()
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        try:
            assert len(xlim) == 2
            assert isinstance(xlim[0], numbers.Number)
            assert isinstance(xlim[1], numbers.Number)
        except (TypeError, AssertionError):
            if xlim is not None:
                raise ValueError('add_pcolor: xlim must be list- ' +
                                 'or array-like and have 2 elements')
        else:
            ax.set_xlim(xlim)
        try:
            assert len(ylim) == 2
            assert isinstance(ylim[0], numbers.Number)
            assert isinstance(ylim[1], numbers.Number)
        except (TypeError, AssertionError):
            if ylim is not None:
                raise ValueError('add_pcolor: ylim must be list- ' +
                                 'or array-like and have 2 elements')
        else:
            ax.set_ylim(ylim)

        # Add body/planet.  Determine where the sun is first.
        if dim1 == 'x':
            ang = 0.0
        elif dim2 == 'x':
            ang = 90.0
        else:
            ang = 0.0
        if add_body and 'rbody' in self.attrs:
            self.add_body(ax, ang=ang)

        return fig, ax, pcol, cbar

    def add_contour(self, dim1, dim2, value, nlev=30, target=None, loc=111,
                    title=None, xlabel=None, ylabel=None,
                    ylim=None, xlim=None, add_cbar=False, clabel=None,
                    filled=True, add_body=True, dolog=False, zlim=None,
                    *args, **kwargs):
        '''
        Create a contour plot of variable **value** against **dim1** on the
        x-axis and **dim2** on the y-axis.

        Simple example:

        >>> self.add_contour('x', 'y', 'rho')

        If kwarg **target** is None (default), a new figure is
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object,
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        =========== ==========================================================
        Kwarg       Description
        =========== ==========================================================
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        nlev        Number of contour levels.  Defaults to 30.
        title       Sets title of axes.  Default is no title.
        xlabel      Sets x label of axes.  Defaults to **dim1** and units.
        ylabel      Sets y label of axes.  Defaults to **dim2** and units.
        xlim        Sets limits of x-axes.  Defaults to whole domain.
        ylim        Sets limits of y-axes.  Defaults to whole domain.
        zlim        Sets contour range.  Defaults to variable max/min.
        add_cbar    Adds colorbar to plot.  Defaults to *False*.
        clabel      Sets colorbar label.  Defaults to **var** and units.
        add_body    Places planetary body in plot.  Defaults to **True**.
        dolog       Sets use of logarithmic scale.  Defaults to **False**.
        =========== ==========================================================
        '''

        import numbers
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import LogLocator, LogFormatterMathtext

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10, 10), loc=loc)

        # Get max/min if none given.
        if zlim is None:
            zlim = [0, 0]
            zlim[0] = self[value].min()
            zlim[1] = self[value].max()
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Set contour command based on grid type.
        if self['grid'].attrs['gtype'] != 'Regular':  # Non-uniform grids.
            if filled:
                contour = ax.tricontourf
            else:
                contour = ax.tricontour
        else:   # Uniform grids.
            if filled:
                contour = ax.contourf
            else:
                contour = ax.contour

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z = np.where(self[value] > zlim[0], self[value], 1.01*zlim[0])
            norm = LogNorm()
            ticks = LogLocator()
            fmt = LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z = self[value]
            norm = None
            ticks = None
            fmt = None

        # Plot contour.
        cont = contour(self[dim1], self[dim2], np.array(z),
                       levs, *args, norm=norm, **kwargs)
        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(cont, ax=ax, ticks=ticks, format=fmt, pad=0.01)
            if clabel is None:
                if 'units' in self[value].attrs:
                    clabel = f"{value}, ({self[value].attrs['units']})"
                else:
                    clabel = f"{value}"
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Set title, labels, axis ranges (use defaults where applicable.)
        if title:
            ax.set_title(title)
        if ylabel is None:
            ylabel = '%s ($R_{E}$)' % dim2.upper()
        if xlabel is None:
            xlabel = '%s ($R_{E}$)' % dim1.upper()
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        try:
            assert len(xlim) == 2
            assert isinstance(xlim[0], numbers.Number)
            assert isinstance(xlim[1], numbers.Number)
        except (TypeError, AssertionError):
            if xlim is not None:
                raise ValueError('add_contour: xlim must be list- ' +
                                 'or array-like and have 2 elements')
        else:
            ax.set_xlim(xlim)
        try:
            assert len(ylim) == 2
            assert isinstance(ylim[0], numbers.Number)
            assert isinstance(ylim[1], numbers.Number)
        except (TypeError, AssertionError):
            if ylim is not None:
                raise ValueError('add_contour: ylim must be list- or ' +
                                 'array-like and have 2 elements')
        else:
            ax.set_ylim(ylim)

        # Add body/planet.  Determine where the sun is first.
        if dim1 == 'x':
            ang = 0.0
        elif dim2 == 'x':
            ang = 90.0
        else:
            ang = 0.0
        if add_body and 'rbody' in self.attrs:
            self.add_body(ax, ang=ang)

        return fig, ax, cont, cbar


class ShellSlice(IdlFile):
    '''
    Shell slices are special MHD outputs where the domain is interpolated
    onto a spherical slice in 1, 2, or 3 dimensions.  Some examples
    include radial or azimuthal lines, spherical shells, or 3D wedges.

    The *Shell* class reads and handles these output types.
    '''

    def __init__(self, filename, format='binary', *args, **kwargs):

        from spacepy.pybats import parse_filename_time

        IdlFile.__init__(self, filename, format=format, keep_case=False)

        # Extract time from file name:
        i_iter, runtime, time = parse_filename_time(self.attrs['file'])
        if 'time' not in self.attrs:
            self.attrs['time'] = time
        if 'iter' not in self.attrs:
            self.attrs['iter'] = i_iter

        # Create some helper variables for plotting and calculations
        d2r = np.pi/180.  # Convert degrees to radians.

        # Get grid spacing.  If npoints ==1, set to 1 to avoid math errors.
        self.drad = (self['r'][-1] - self['r'][0])/max(self['grid'][0]-1, 1)
        self.dlon = (self['lon'][-1]-self['lon'][0])/max(self['grid'][1]-1, 1)
        self.dlat = (self['lat'][-1]-self['lat'][0])/max(self['grid'][2]-1, 1)

        self.dphi = d2r*self.dlon
        self.dtheta = d2r*self.dlat

        # Get spherical, uniform grid in units of r_body/radians:
        self.lon, self.r, self.lat = np.meshgrid(
            np.array(self['lon']), np.array(self['r']), np.array(self['lat']))
        self.phi = d2r*self.lon
        self.theta = d2r*(90-self.lat)

    @calc_wrapper
    def calc_urad(self):
        '''
        Calculate radial velocity.
        '''
        ur = self['ux']*np.sin(self.theta)*np.cos(self.phi) + \
            self['uy']*np.sin(self.theta)*np.sin(self.phi) + \
            self['uz']*np.cos(self.theta)

        self['ur'] = dmarray(ur, {'units': self['ux'].attrs['units']})

    @calc_wrapper
    def calc_radflux(self, var, conv=1000. * (100.0)**3):
        '''
        For variable *var*, calculate the radial flux of *var* through each
        grid point in the slice as self[var]*self['ur'].
        Resulting value stored as "var_rflx".
        '''

        if var+'_rflx' in self:
            return

        # Make sure we have radial velocity.
        if 'ur' not in self:
            self.calc_urad()

        # Calc flux:
        self[var+'_rflx'] = self[var] * self['ur'] * conv

    @calc_wrapper
    def calc_radflu(self, var):
        '''
        For variable *var*, calculate the radial fluence, or the
        spatially integrated radial flux through 2D surfaces of
        constant radius.

        Resulting variable stored as "var_rflu".  Result will be an array
        with one value for each radial distance within the object.
        '''

        # Need at least 2D in angle space:
        if self.dphi == 0 or self.dtheta == 0:
            raise ValueError('Fluence can be calculated for >=2D surfaces.')

        # Trim flux off of val name:
        if '_rflx' in var:
            var = var[:-5]

        # Convenience:
        flux = var + '_rflx'
        flu = var + '_rflu'

        # Make sure flux exists:
        if flux not in self:
            self.calc_radflux(var)
        if flu in self:
            return

        # Create output container, one point per radial distance:
        self[flu] = np.zeros(self['grid'][0])

        # Integrate over all radii.
        # Units: convert R to km and cm-3 to km.
        for i, R in enumerate(self['r']):
            self[flu] = np.sum(
                (R * 6371.0)**2 * self[flux][i, :, :] *
                np.sin(self.theta[i, :, :]) *
                self.dtheta * self.dphi * 1000.**2)

    def add_cont_shell(self, value, irad=0, target=None, loc=111,
                       zlim=None, dolabel=True, add_cbar=False,
                       dofill=True, nlev=51, colat_max=90,
                       rotate=np.pi/2, dolog=False, clabel=None,
                       latticks=15., yticksize=14, extend='both', **kwargs):
        '''
        For slices that cover a full hemisphere or more, create a polar
        plot of variable *value*.

        Extra keywords are sent to the matplotlib contour command.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import (LogLocator, LogFormatterMathtext,
                                       MultipleLocator)

        fig, ax = set_target(target, figsize=(10, 10), loc=loc, polar=True)

        # Get max/min if none given.
        if zlim is None:
            zlim = [0, 0]
            zlim[0] = self[value][irad, :, :].min()
            zlim[1] = self[value][irad, :, :].max()

            # For log space, no negative zlimits.
            if dolog and zlim[0] <= 0:
                zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Create levels and set norm based on dolog.
        if dolog:  # Log space!
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            norm = LogNorm()
            ticks = LogLocator()
            fmt = LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            norm = None
            ticks = None
            fmt = None

        # Select proper contour function based on fill/don't fill.
        if dofill:
            func = ax.contourf
        else:
            func = ax.contour

        # Plot result.  Rotate "rotate" radians to get sun in right spot.
        # Plot against colatitude to arrange results correctly.
        cnt = func(self.phi[irad, :, :]+rotate, 90-self.lat[irad, :, :],
                   np.array(self[value][irad, :, :]), levs, norm=norm,
                   extend=extend, **kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(cnt, ax=ax, ticks=ticks, format=fmt,
                                shrink=.85)
            if clabel is None:
                if 'units' in self[value].attrs:
                    clabel = f"{value}, ({self[value].attrs['units']})"
                else:
                    clabel = f"{value}"
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Adjust latitude
        ax.set_ylim([0, colat_max])

        # Adjust atitude and add better labels:
        ax.yaxis.set_major_locator(MultipleLocator(latticks))
        ax.set_ylim([0, colat_max])
        ax.set_yticklabels('')
        opts = {'size': yticksize, 'rotation': -45,
                'ha': 'center', 'va': 'center'}
        for theta in np.arange(90-latticks, 90-colat_max, -latticks):
            txt = '{:02.0f}'.format(theta)+r'$^{\circ}$'
            ax.text(np.pi/4., 90.-theta, txt, color='w', weight='extra bold',
                    **opts)
            ax.text(np.pi/4., 90.-theta, txt, color='k', weight='light',
                    **opts)

        # Use MLT-type labels.
        lt_labels = ['Noon', '18', '00',   '06']
        xticks = [0, np.pi/2,   np.pi, 3*np.pi/2]
        xticks = np.array(xticks) + rotate

        # Apply x-labels:
        ax.set_xticks(xticks)
        ax.set_xticklabels(lt_labels)

        return fig, ax, cnt, cbar


class Mag(PbData):
    '''
    A container for data from a single BATS-R-US virtual magnetometer.  These
    work just like a typical :class:`spacepy.pybats.PbData` object.  Beyond
    raw magnetometer data, additional values are calculated and stored,
    including total pertubations (the sum of all global and ionospheric
    pertubations as measured by the magnetometer).  Users will be interested
    in methods :meth:`~spacepy.pybats.bats.Mag.add_comp_plot` and
    :meth:`~spacepy.pybats.bats.Mag.calc_dbdt`.

    Instantiation is best done through :class: `spacepy.pybats.MagFile`
    objects, which load and parse organize many virtual magnetometers from a
    single output file into a single object.  However, they can be created
    manually, though painfully.  Users must instantiate by handing the
    new object the number of lines that will be parsed (rather, the number
    of data points that will be needed), a time vector, and (optionally)
    the list of variables coming from the GM and IE module.  While the
    latter two are keyword arguments, at least one should be provided.
    Next, the arrays whose keys were given by the *gmvars* and *ievars*
    keyword arguments in the instantiation step can either be filled manually
    or by using the :meth:`~spacepy.pybats.bats.Mag.parse_gmline` and
    :meth:`~spacepy.pybats.bats.Mag.parse_ieline` methods to parse lines of
    ascii data from a magnetometer output file.  Finally, the
    :meth:`~spacepy.pybats.bats.Mag.recalc` method should be called to
    calculate total perturbation.
    '''
    def __init__(self, nlines, time, gmvars=(), ievars=(), *args, **kwargs):
        super(Mag, self).__init__(*args, **kwargs)  # Init as PbData.

        self['time'] = time
        self.attrs['nlines'] = nlines

        self['x'] = np.zeros(nlines)
        self['y'] = np.zeros(nlines)
        self['z'] = np.zeros(nlines)

        # Create IE and GM specific containers.
        for key in gmvars:
            self[key] = np.zeros(nlines)
        for key in ievars:
            self['ie_'+key] = np.zeros(nlines)

    def parse_gmline(self, i, line, namevar):
        '''
        Parse a single line from a GM_mag*.out file and put into
        the proper place in the magnetometer arrays.  The line should
        have the same number of variables as was initially given to the
        :class:`~spacepy.pybats.bats.Mag` object.  This method is best
        used through the :class:`~spacepy.pybats.bats.MagFile` class interface.

        Usage:

        >>> self.parse_gmline(i, line, namevar)

        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts = line.split()
        self['x'][i] = float(parts[9])
        self['y'][i] = float(parts[10])
        self['z'][i] = float(parts[11])
        for j, key in enumerate(namevar):
            self[key][i] = float(parts[j+12])

    def parse_ieline(self, i, line, namevar):
        '''
        Parse a single line from a IE_mag*.out file and put into
        the proper place in the magnetometer arrays.  The line should
        have the same number of variables as was initially given to the
        :class:`~spacepy.pybats.bats.Mag` object.  This method is best
        used through the :class:`~spacepy.pybats.bats.MagFile` class interface.


        Usage:

        >>> self.parse_gmline(i, line, namevar)

        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts = line.split()
        for j, key in enumerate(namevar):
            self['ie_'+key][i] = float(parts[j+11])

    def _recalc(self):
        '''
        Calculate total :math:`\\Delta B` from GM and IE; store under object
        keys *totaln*, *totale*, and *totald* (one for each component of the
        HEZ coordinate system).

        This function should only be called to correct legacy versions of
        magnetometer files.
        '''

        # If values already exist, do not overwrite.
        if 'dBn' in self:
            return

        # New containers:
        self['totaln'] = np.zeros(self.attrs['nlines'])
        self['totale'] = np.zeros(self.attrs['nlines'])
        self['totald'] = np.zeros(self.attrs['nlines'])

        for key in list(self.keys()):
            if key[-2:] == 'Bn':
                self['totaln'] = self['totaln']+self[key]
            if key[-2:] == 'Be':
                self['totale'] = self['totale']+self[key]
            if key[-2:] == 'Bd':
                self['totald'] = self['totald']+self[key]

        # Old names -> new names:
        varmap = {'totaln': 'dBn', 'totale': 'dBe', 'totald': 'dBd',
                  'gm_dBn': 'dBnMhd', 'gm_dBe': 'dBeMhd', 'gm_dBd': 'dBdMhd',
                  'gm_facdBn': 'dBnFac', 'gm_facdBe': 'dBeFac',
                  'gm_facdBd': 'dBdFac', 'ie_JhdBn': 'dBnHal',
                  'ie_JhdBe': 'dBeHal', 'ie_JhdBd': 'dBdHal',
                  'ie_JpBn': 'dBnPed',  'ie_JpBe': 'dBePed',
                  'ie_JpBd': 'dBdPed'}

        # Replace variable names.
        for key in list(self.keys()):
            if key in varmap:
                self[varmap[key]] = self.pop(key)

    def calc_h(self):
        '''
        Calculate the total horizontal perturbation, 'H', using the pythagorean
        sum of the two horizontal components (north-south and east-west
        components):

        $\\Delta B_H = \\sqrt{\\Delta B_N^2 + \\Delta B_E^2}$
        '''

        allvars = list(self.keys())

        for v in allvars:
            # Find all dB-north variables:
            if v[:3] == 'dBn':
                v_east = v.replace('dBn', 'dBe')
                self[v.replace('dBn', 'dBh')] = dmarray(
                    np.sqrt(self[v]**2+self[v_east]**2), {'units': 'nT'})

    def calc_dbdt(self):
        '''
        Calculate the time derivative of all dB-like variables and save as
        'dBdt[direction][component].  For example, the time derivative of
        dBeMhd will be saved as dBdteMhd.

        |dB/dt|_h is also calculated following the convention of
        Pulkkinen et al, 2013:
        $|dB/dt|_H = \\sqrt{(\\dB_N/dt)^2 + (dB_E/dt)^2}$

        A 2nd-order accurate centeral difference method is used to
        calculate the time derivative.  For the first and last points,
        2nd-order accurate forward and backward differences are taken,
        respectively.
        '''

        # Do not calculate twice.
        if 'dBdtn' in self:
            return

        # Get dt values:
        dt = np.array([x.total_seconds() for x in np.diff(self['time'])])

        # Loop through variables:
        oldvars = list(self.keys())
        for k in oldvars:
            if 'dB' not in k:
                continue

            # Create new variable name and container:
            new = k.replace('dB', 'dBdt')
            self[new] = dmarray(np.zeros(self.attrs['nlines']),
                                {'units': 'nT/s'})

            # Central diff:
            self[new][1:-1] = (self[k][2:]-self[k][:-2])/(dt[1:]+dt[:-1])

            # Forward diff:
            self[new][0] = (-self[k][2] + 4*self[k][1] - 3*self[k][0]) \
                / (dt[1]+dt[0])

            # Backward diff:
            self[new][-1] = (3*self[k][-1] - 4*self[k][-2] + self[k][-3]) \
                / (dt[-1]+dt[-2])

        self['dBdth'] = np.sqrt(self['dBdtn']**2+self['dBdte']**2)

    def add_plot(self, value, style='-', target=None, loc=111, label=None,
                 **kwargs):
        '''
        Plot **value**, which should be a key corresponding to a data vector
        stored in the :class:`~spacepy.pybats.bats.Mag` object,
        against the object's *time*.  The **target** kwarg specifies the
        destination of the plot.  If not set, **target** defaults to None and
        a new figure and axis will be created.  If target is a matplotlib
        figure, a new axis is created at subplot location 111 (which can be
        changed
        using kwarg **loc**).  If target is a matplotlib Axes object, the line
        is added to the plot as if ``Axes.plot()`` was used.  The line label,
        which is used for setting labels on figure legends,
        defaults to the value key but can be customized with the "label"
        kwarg.  The **style** keyword accepts basic Matplotlib line style
        strings such as '--r' or '+g'.  This string will be passed on to the
        plot command to customize the line.

        All extra kwargs are handed to ``Axes.plot``, allowing the user to set
        any additional options (e.g., line color and style, etc.).

        Three values are returned: the Figure object, Axis object, and
        newly created line object.  These can be used to further customize
        the figure, axis, and line as necessary.

        Example: Plot total :math:`\\Delta B_N` onto an existing axis with line
        color blue, line style dashed, and line label "Wow!":

        >>> self.plot('dBn', target='ax', label='Wow!', lc='b', ls='--')

        Example: Plot total :math:`\\Delta B_N` on a new figure, save returned
        values and overplot additional values on the returned axis.  Default
        labels and line styles are used in this example.

        >>> fig, ax, line = self.plot('n')
        >>> self.plot('dBe', target = ax)

        '''

        if not label:
            label = value

        # Set figure and axes based on target:
        fig, ax = set_target(target, figsize=(10, 4), loc=loc)

        ax.plot(self['time'], self[value], style, label=label, **kwargs)
        applySmartTimeTicks(ax, self['time'], dolabel=True)

        return fig, ax

    def add_comp_plot(self, direc, target=None, add_legend=True,
                      loc=111, lw=2.0):
        '''
        Create a plot with,  or add to an existing plot, an illustration of
        how the
        separate components sum together to make the total disturbance in a
        given orthongal direction (arg **direc**).  The three possible
        components are 'n' (northwards, towards the magnetic pole), 'e'
        (eastwards), or 'd' (downwards towards the center of the Earth.)  The
        components of the total disturbance in any on direction are
        magnetospheric currents ('gm_dB'), gap-region field-aligned currents
        ('gm_facdB'), and ionospheric Hall and Pederson currents ('ie_Jp' and
        'ie_Jh').

        Example usage:

        >>> self.add_comp_plot('n')

        This will create a new plot with the total disturbance in the 'n'
        direction along with line plots of each component that builds this
        total.  This method uses the familiar PyBats **target** kwarg system
        to allow users to add these plots to existing figures or axes.

        Parameters
        ==========
        direc : string
           Indicate the direction to plot: either 'n', 'e', 'd', or
           'h' if calculated.

        Other Parameters
        ================
        target : Matplotlib Figure or Axes object
           Set plot destination.  Defaults to new figure.
        loc : 3-digit integer
           Set subplot location.  Defaults to 111.
        add_legend : bool
           Add legend to plot.  Defaults to True.
        lw : float
           Set the width of the lines.  Defaults to 2.0  Total field
           is always 1.5 times thicker.
        '''

        # Set plot targets.
        fig, ax = set_target(target, figsize=(10, 4), loc=loc)

        prefix = 'dB'+direc

        # Use a dictionary to assign line styles, widths.
        styles = {prefix+'Mhd': '--', prefix+'Fac': '--',
                  prefix+'Hal': '-.', prefix+'Ped': '-.',
                  prefix: '-'}
        widths = {prefix+'Mhd': lw, prefix+'Fac': lw,
                  prefix+'Hal': lw, prefix+'Ped': lw,
                  prefix: 1.5*lw}
        colors = {prefix+'Mhd': '#FF6600', prefix+'Fac': 'r',
                  prefix+'Hal': 'b',       prefix+'Ped': 'c',
                  prefix: 'k'}

        # Labels:
        labels = {prefix+'Mhd': r'$J_{Mag}$',  prefix+'Fac': r'$J_{Gap}$',
                  prefix+'Hal': r'$J_{Hall}$', prefix+'Ped': r'$J_{Peder}$',
                  prefix: r'Total $\Delta B'+'_{}$'.format(direc)}

        # Plot.
        for k in sorted(self):
            if ('dB'+direc not in k) or (k == 'time') or (k[:4] == 'dBdt'):
                continue
            ax.plot(self['time'], self[k], label=labels[k],
                    lw=widths[k], c=colors[k])  # ,ls=styles[k]

        # Ticks, zero-line, and legend:
        applySmartTimeTicks(ax, self['time'], True, True)
        ax.hlines(0.0, self['time'][0], self['time'][-1],
                  linestyles=':', lw=2.0, colors='k')
        if add_legend:
            ax.legend(ncol=3, loc='best')

        # Axis labels:
        ax.set_ylabel(r'$\Delta B_{%s}$ ($nT$)' % (direc.upper()))

        if target is None:
            fig.tight_layout()

        return fig, ax


class MagFile(PbData):
    '''
    BATS-R-US magnetometer files are powerful tools for both research and
    operations.  :class:`~spacepy.pybats.bats.MagFile` objects open, parse,
    and visualize such output.

    The $\\delta B$ calculated by the SWMF requires two components: GM (BATSRUS)
    and IE (Ridley_serial).  The data is spread across two files: GM_mag*.dat
    and IE_mag*.dat.  The former contains $\\delta B$ caused by gap-region
    (i.e., inside the inner boundary) FACs and the changing global field.
    The latter contains the $\\delta B$ caused by Pederson and Hall
    currents in the ionosphere.  :class:`~spacepy.pybats.bats.MagFile` objects
    can open one or both of these files at a time; when both are opened, the
    total $\\delta B$ is calculated and made available to the user.

    Usage:

    >>> # Open up the GM magnetometer file only.
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag')
    >>>
    >>> # Open up both the GM and IE file [LEGACY SWMF ONLY]
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag', 'IE_file.mag')
    >>>
    >>> # Open up the GM magnetometer file; search for the IE file.
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag', find_ie=True)

    Note that the **find_ie** kwarg uses a simple search assuming the data
    remain in a typical SWMF-output organizational tree (i.e., if the results
    of a simulation are in folder *results*, the GM magnetometer file can be
    found in *results/GM/* or *results/GM/IO2/* while the IE file can be found
    in *results/IE/* or *results/IE/ionosphere/*).  It will also search the
    present working directory.  This method is not robust; the user must take
    care to ensure that the two files correspond to each other.
    '''
    def __init__(self, filename, ie_name=None, find_ie=False, *args, **kwargs):

        from glob import glob

        super(MagFile, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['gmfile'] = filename

        # Try to find the IE file based on our current location.
        if find_ie and not ie_name:
            basedir = filename[0:filename.rfind('/')+1]
            if glob(basedir + '../IE/IE_mag_*.mag'):
                self.attrs['iefile'] = glob(basedir + '../IE/IE_mag_*.mag')[-1]
            elif glob(basedir + '../../IE/ionosphere/IE_mag_*.mag'):
                self.attrs['iefile'] = \
                    glob(basedir + '../../IE/ionosphere/IE_mag_*.mag')[-1]
            elif glob(basedir + '/IE_mag_*.mag'):
                self.attrs['iefile'] = glob(basedir + '/IE_mag_*.mag')[-1]
            else:
                self.attrs['iefile'] = None
        else:
            self.attrs['iefile'] = ie_name

        # Set legacy mode to handle old variable names:
        self.legacy = find_ie or bool(ie_name)

        self.readfiles()

    def readfiles(self):
        '''
        Read and parse GM file and IE file (if name given.)
        '''
        import datetime as dt

        # Slurp lines.
        infile = open(self.attrs['gmfile'], 'r')
        lines = infile.readlines()
        infile.close()

        # Parse header:
        # Get number of stations.
        nmags = int((lines[0].split(':')[0]).split()[0])

        # Get station names.
        names = lines[0].split(':')[1]
        namemag = names.split()
        self.attrs['namemag'] = namemag

        # Check nmags vs number of mags in header.
        if nmags != len(namemag):
            raise BaseException(
                'ERROR: GM file claims %i magnetomers, lists %i'
                % (nmags, len(namemag)))

        # Grab variable names.  Use legacy mode if necessary:
        prefix = 'gm_' * self.legacy
        # skip time, iter, and loc; add prefix to var names:
        gm_namevar = [prefix + x for x in lines[1].split()[12:]]

        # Set number of mags and records.
        self.attrs['nmag'] = len(namemag)
        nrecords = (len(lines)-2)//nmags

        # If there is an IE file, Parse that header, too.
        if self.attrs['iefile']:
            infile = open(self.attrs['iefile'], 'r')
            ielns = infile.readlines()
            infile.close()
            nmags = int((ielns[0].split(':')[0]).split()[0])
            iestats = (ielns[0].split(':')[1]).split()
            # Check nmags vs number of mags in header.
            if nmags != len(iestats):
                raise BaseException(
                    'ERROR: IE file claims %i magnetomers, lists %i'
                    % (nmags, len(namemag)))
            if iestats != self.attrs['namemag']:
                raise RuntimeError("Files do not have matching stations.")
            ie_namevar = ielns[1].split()[11:]
            self.attrs['ie_namevar'] = ie_namevar
            if (len(ielns)/self.attrs['nmag']) != (nrecords-1):
                print('Number of lines do not match: GM=%d, IE=%d!' %
                      (nrecords-1, len(ielns)/self.attrs['nmag']))
                nrecords = min(ielns, nrecords-1)
        else:
            ie_namevar = ()
            self.attrs['ie_namevar'] = ()

        # Build containers.
        self['time'] = np.zeros(nrecords, dtype=object)
        self['iter'] = np.zeros(nrecords, dtype=float)
        for name in namemag:
            self[name] = Mag(nrecords, self['time'], gm_namevar, ie_namevar)

        data_buffer = np.zeros((nrecords, nmags, (len(gm_namevar)+3)))

        # Read file data.
        for i in range(nrecords):
            line = lines[i*nmags+2]
            parts = line.split()
            self['iter'][i] = parts[0]
            self['time'][i] = dt.datetime(
                int(parts[1]),  # year
                int(parts[2]),  # month
                int(parts[3]),  # day
                int(parts[4]),  # hour
                int(parts[5]),  # minute
                int(parts[6]),  # second
                int(parts[7])*1000  # microsec
                )
            for j in range(nmags):
                line = lines[i*nmags+j+2]
                if j > 0:
                    parts = line.split()
                values = [float(part) for part in parts[9:]]
                data_buffer[i, j] = values

            if self.attrs['iefile'] and i > 0:
                line = ielns[i*nmags+2]
                self[namemag[0]].parse_ieline(i, line, ie_namevar)
                for j in range(1, nmags):
                    self[namemag[j]].parse_ieline(i, ielns[i*nmags+j+2],
                                                  ie_namevar)

        for j in range(nmags):
            mag = self[namemag[j]]
            mag['x'] = data_buffer[:, j, 0]
            mag['y'] = data_buffer[:, j, 1]
            mag['z'] = data_buffer[:, j, 2]
            for k, key in enumerate(gm_namevar):
                mag[key] = data_buffer[:, j, k+3]

        # Sum up IE/GM components if necessary (legacy only):
        if self.legacy:
            self._recalc()

        # Get time res.
        self.attrs['dt'] = (self['time'][1] - self['time'][0]).seconds / 60.0

    def _recalc(self):
        '''
        Old magnetometer files had different variable names and did not
        contain the total perturbation.  This function updates variable names
        and sums all contributions from all models/regions to get total
        :math:`\\Delta B`.

        This function is only required for legacy results.  New versions of
        the SWMF include both GM, IE, and total perturbations in a single
        file.
        '''
        for mag in self.attrs['namemag']:
            self[mag]._recalc()

    def calc_h(self):
        '''
        For each magnetometer object, calculate the horizontal component of
        the perturbations using the pythagorean sum of the two horizontal
        components (north-south and east-west components):

        $\\Delta B_H = \\sqrt{\\Delta B_N^2 + \\Delta B_E^2}$
        '''
        for k in self:
            if k == 'time' or k == 'iter':
                continue
            self[k].calc_h()

    def calc_dbdt(self):
        '''
        For each magnetometer object, calculate the horizontal component of
        the perturbations.

        |dB/dt|_h is also calculated following the convention of
        Pulkkinen et al, 2013:
        $|dB/dt|_H = \\sqrt{(\\dB_N/dt)^2 + (dB_E/dt)^2}$
        '''
        for k in self:
            if k == 'time' or k == 'iter':
                continue
            self[k].calc_dbdt()


class MagGridFile(IdlFile):
    '''
    Magnetometer grids are a recent addition to BATS-R-US: instead of
    specifying a small set of individual stations, the user can specify a
    grid of many stations spanning a latitude/longitude range.  The files
    are output in the usual :class:`spacepy.pybats.IdlFile` format.  This
    class handles the reading, manipulating, and visualization of these files.
    '''

    def __init__(self, *args, **kwargs):
        import re
        from spacepy.pybats import parse_filename_time

        # Initialize as an IdlFile.
        super(MagGridFile, self).__init__(header=None, *args, **kwargs)

        # Additional header parsing:
        head = self.attrs['header']
        match = re.search(r'\((\w+)\).*\[(\w+)\].*\[(\w+)\]', head)
        coord, unit1, unit2 = match.groups()

        self['grid'].attrs['coord'] = coord

        # Extract time from file name:
        i_iter, runtime, time = parse_filename_time(self.attrs['file'])
        if 'time' not in self.attrs:
            self.attrs['time'] = time
        if 'iter' not in self.attrs:
            self.attrs['iter'] = i_iter

        # Set units based on header parsing:
        for v in self:
            if v == 'grid':
                continue
            elif v in self['grid'].attrs['dims']:
                self[v].attrs['units'] = unit1
            else:
                self[v].attrs['units'] = unit2

    @calc_wrapper
    def calc_h(self):
        '''
        Calculate the total horizontal perturbation, 'h', using the pythagorean
        sum of the two horizontal components (east and west components).
        '''

        allvars = list(self.keys())

        for v in allvars:
            # Find all dB-north variables:
            if v[:3] == 'dBn':
                v_east = v.replace('dBn', 'dBe')
                self[v.replace('dBn', 'dBh')] = dmarray(
                    np.sqrt(self[v]**2 + self[v_east]**2), {'units': 'nT'})

    def add_contour(self, value, nlev=30, target=None, loc=111,
                    title=None, xlabel=None, ylabel=None,
                    ylim=None, xlim=None, add_cbar=False, clabel=None,
                    filled=True, dolog=False, zlim=None,
                    coords=None, add_conts=False,
                    *args, **kwargs):
        '''
        Add a lat-lon contour plot of variable *value*.
        '''

        from spacepy.pybats import mhdname_to_tex

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import (LogLocator, FuncFormatter,
                                       LogFormatterMathtext, MultipleLocator)

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10, 7), loc=loc)

        # Get max/min if none given.
        if zlim is None:
            zlim = [self[value].min(), self[value].max()]
            if zlim[1]-zlim[0] > np.max(zlim):
                maxval = max(abs(zlim[0]), abs(zlim[1]))
                zlim = [-maxval, maxval]

        # No zero-level for log scales:
        if dolog and zlim[0] <= 0:
            zlim[0] = np.min([0.0001, zlim[1]/1000.0])

        # Better default color maps:
        if 'cmap' not in kwargs:
            # If zlim spans positive and negative:
            if zlim[1]-zlim[0] > np.max(zlim):
                kwargs['cmap'] = 'bwr'
            else:
                kwargs['cmap'] = 'Reds'

        # Set contour command based on filled/unfilled contours:
        if filled:
            contour = ax.contourf
        else:
            contour = ax.contour

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]),
                                            np.log10(zlim[1]), nlev))
            z = np.where(self[value] > zlim[0], self[value], 1.01*zlim[0])
            norm = LogNorm()
            ticks = LogLocator()
            fmt = LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z = self[value]
            norm = None
            ticks = MultipleLocator((zlim[1]-zlim[0])/10)  # fix this
            fmt = None

        # Add Contour to plot:
        cont = contour(self['Lon'], self['Lat'], np.array(np.transpose(z)),
                       levs, *args, norm=norm, **kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar = plt.colorbar(cont, ax=ax, ticks=ticks, format=fmt, pad=0.01)
            if clabel is None:
                varname = mhdname_to_tex(value)
                if 'units' in self[value].attrs:
                    units = mhdname_to_tex(self[value].attrs['units'])
                else:
                    units = 'unitless'
                clabel = f"{varname} ({units})"
            cbar.set_label(clabel)
        else:
            cbar = None  # Need to return something, even if none.

        # Set title, labels, axis ranges (use defaults where applicable.)
        if title:
            ax.set_title(title)
        coord_sys = self['grid'].attrs['coord']
        if ylabel is None:
            ylabel = 'Latitude ({})'.format(coord_sys)
        if xlabel is None:
            xlabel = 'Longitude ({})'.format(coord_sys)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        if type(xlim) is not None:
            ax.set_xlim(xlim)
        if type(ylim) is not None:
            ax.set_ylim(ylim)

        # If a brand-new axes was created, use custom ticks:
        if not issubclass(type(target), plt.Axes):
            fmttr = FuncFormatter(lambda x, pos:
                                  '{:4.1f}'.format(x)+r'$^{\circ}$')
            ax.xaxis.set_major_formatter(fmttr)
            ax.yaxis.set_major_formatter(fmttr)

        # If a brand-new figure was created, use tight-layout.
        if target is None:
            fig.tight_layout()

        return fig, ax, cont, cbar

    def interp(self, var, lons, lats):
        '''
        Interpolate variable *var* to point(s) *lats*, *lons*.  *lats* and
        *lons* can either be scalar or a numpy-like array.  Values interpolated
        to these positions is returned to the caller either as a scalar or
        array.

        Simple bilinear interpolation is used.
        '''

        from spacepy.pybats.batsmath import interp_2d_reg as intp

        return intp(lats, lons, self['Lat'], self['Lon'], self[var])


class GeoIndexFile(LogFile):
    '''
    Geomagnetic Index files are a specialized BATS-R-US output that contain
    geomagnetic indices calculated from simulated ground-based magnetometers.
    Currently, the only index instituted is Kp through the faKe_p setup.
    Future work will expand the system to include Dst, AE, etc.

    GeoIndFiles are a specialized subclass of pybats.LogFile.  It includes
    additional methods to quickly visualize the output, perform data-model
    comparisons, and more.
    '''

    def __repr__(self):
        return 'GeoIndexFile object at %s' % (self.attrs['file'])

    def __init__(self, filename, keep_case=True, *args, **kwargs):
        '''
        Load ascii file located at self.attrs['file'].
        '''
        # Call super __init__:
        super(GeoIndexFile, self).__init__(filename, *args, **kwargs)

        # Re-parse the file header; look for key info.
        head = (self.attrs['descrip']).replace('=', ' ')
        parts = head.split()
        if 'DtOutput=' in head:
            self.attrs['dt'] = float(parts[parts.index('DtOutput=') + 1])
        if 'SizeKpWindow' in head:
            self.attrs['window'] = float(parts[parts.index(
                        'SizeKpWindow(Mins)')+1])
        if 'Lat' in head:
            self.attrs['lat'] = float(parts[parts.index('Lat') + 1])
        if 'K9' in head:
            self.attrs['k9'] = float(parts[parts.index('K9') + 1])

    def fetch_obs_kp(self):
        '''
        Fetch the observed Kp index for the time period covered in the
        logfile.  Return *True* on success.

        Observed Kp is automatically fetched from the Kyoto World Data Center
        via the :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoKp` object, which holds the observed
        kp, is stored as *self.obs_kp* for future use.
        '''

        import spacepy.pybats.kyoto as kt

        # Return if already obtained:
        if hasattr(self, 'obs_kp'):
            return True

        # Start and end time to collect observations:
        stime = self['time'][0]
        etime = self['time'][-1]

        # Attempt to fetch from Kyoto website:
        try:
            self.obs_kp = kt.fetch('kp', stime, etime)
        # Warn on failure:
        except BaseException as args:
            raise Warning('Failed to fetch Kyoto Kp: ', args)
            return False

        return True

    def fetch_obs_ae(self):
        '''
        Fetch the observed AE index for the time period covered in the
        logfile.  Return *True* on success.

        Observed AE is automatically fetched from the Kyoto World Data Center
        via the :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoAe` object, which holds the observed
        AE, is stored as *self.obs_ae* for future use.
        '''

        import spacepy.pybats.kyoto as kt

        # Return if already obtained:
        if hasattr(self, 'obs_ae'):
            return True

        # Start and end time to collect observations:
        stime = self['time'][0]
        etime = self['time'][-1]

        # Attempt to fetch from Kyoto website:
        try:
            self.obs_ae = kt.fetch('ae', stime, etime)
        # Warn on failure:
        except BaseException as args:
            raise Warning('Failed to fetch Kyoto AE: ', args)
            return False

        return True

    def add_kp_quicklook(self, target=None, loc=111, label=None,
                         plot_obs=False, add_legend=True,
                         obs_kwargs={'c': 'k', 'ls': '--', 'lw': 2}, **kwargs):
        '''
        Similar to "dst_quicklook"-type functions, this method fetches observed
        Kp from the web and plots it alongside the Kp read from the GeoIndex
        file.

        Usage:
        >>> obj.kp_quicklook(target=SomeMplTarget)
        The target kwarg works like in other PyBats plot functions: it can be
        a figure, an axes, or None, and it determines where the plot is placed.

        Other kwargs customize the line.  Label defaults to fa$K$e$_{P}$, extra
        kwargs are passed to pyplot.plot.

        Observed Kp can be added via the *plot_obs* kwarg.  Kp is automatically
        fetched from the Kyoto World Data Center via the
        :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoKp` object, which holds the observed
        Kp, is stored as *self.obs_kp* for future use.
        The observed line can be customized via the *obs_kwargs* kwarg, which
        is a dictionary of plotting keyword arguments.
        '''

        # Set up plot target.
        fig, ax = set_target(target, figsize=(10, 4), loc=loc)

        # Create label:
        if not label:
            label = 'fa$K$e$_{P}$'
            if ('lat' in self.attrs) and ('k9' in self.attrs):
                label += ' (Lat=%04.1f$^{\\circ}$, K9=%03i)' % \
                    (self.attrs['lat'], self.attrs['k9'])

        # Sometimes, the "Kp" varname is caps, sometimes not.
        kp = 'Kp'
        if kp not in self.keys():
            kp = 'kp'

        ax.plot(self['time'], self['Kp'], label=label, **kwargs)
        ax.set_ylabel('$K_{P}$')
        ax.set_xlabel('Time from ' + self['time'][0].isoformat()+' UTC')
        applySmartTimeTicks(ax, self['time'])

        if target is None:
            fig.tight_layout()

        if plot_obs:
            # Attempt to fetch Kp:
            if self.fetch_obs_kp():
                # If successful, add to plot:
                self.obs_kp.add_histplot(target=ax, **obs_kwargs)
                applySmartTimeTicks(ax, self['time'])

        if add_legend:
            ax.legend(loc='best')
        return fig, ax

    def add_ae_quicklook(self, target=None, loc=111, label=None,
                         plot_obs=False, val='AE', add_legend=True,
                         obs_kwargs={'c': 'k', 'ls': '--', 'lw': 1.5},
                         **kwargs):
        '''
        Similar to "dst_quicklook"-type functions, this method fetches observed
        AE indices from the web and plots it alongside the corresponding
        AE read from the GeoInd file.  Because there are four AE-like indices
        (AL, AU, AE, and AO), the kwarg *val* specifies which to plot
        (default is AE).

        Usage:
        >>> obj.ae_quicklook(target=SomeMplTarget)
        The target kwarg works like in other PyBats plot functions: it can be
        a figure, an axes, or None, and it determines where the plot is placed.

        Other kwargs customize the line.  Extra kwargs are passed to
        pyplot.plot

        Observed AE can be added via the *plot_obs* kwarg.  AE is automatically
        fetched from the Kyoto World Data Center via the
        :mod:`spacepy.pybats.kyoto` module.  The associated
        :class:`spacepy.pybats.kyoto.KyotoAe` object, which holds the observed
        AE, is stored as *self.obs_kp* for future use.
        The observed line can be customized via the *obs_kwargs* kwarg, which
        is a dictionary of plotting keyword arguments.
        '''

        # Set up plot target.
        fig, ax = set_target(target, figsize=(10, 4), loc=loc)

        if not label:
            label = 'Virtual {}'.format(val)

        ax.plot(self['time'], self[val], label=label, **kwargs)
        ax.set_ylabel('{} ($nT$)'.format(val))
        ax.set_xlabel('Time from ' + self['time'][0].isoformat() + ' UTC')
        applySmartTimeTicks(ax, self['time'])

        if target is None:
            fig.tight_layout()

        if plot_obs:
            # Check for label in obs. kwargs:
            if 'label' not in obs_kwargs:
                obs_kwargs['label'] = 'Obs. ' + val

            if self.fetch_obs_ae():
                ax.plot(self.obs_ae['time'], self.obs_ae[val.lower()],
                        **obs_kwargs)
                applySmartTimeTicks(ax, self['time'])

        if add_legend:
            ax.legend(loc='best')

        return fig, ax


class VirtSat(LogFile):
    '''
    A :class:`spacepy.pybats.LogFile` object tailored to virtual satellite
    output; includes special satellite-specific plotting methods.
    '''

    def __init__(self, *args, **kwargs):

        from re import findall
        from scipy.interpolate import interp1d
        from matplotlib.dates import date2num

        super(VirtSat, self).__init__(keep_case=False, *args, **kwargs)
        # Attempt to extract the satellite's name and save it.
        try:
            s = self.attrs['file']
            searchstr = r'sat_([\-\w]+)_\d+_n\d+\.sat|sat_(\w+)_n\d+\.sat'
            a = findall(searchstr, s)[0]
            name = list(filter(None, a))[0]
        except IndexError:
            name = None
        self.attrs['name'] = name

        # Create interpolation functions for position.
        self._interp = {}
        tnum = date2num(self['time'])
        for x in 'xyz':
            self._interp[x] = interp1d(tnum, self[x])

    def __repr__(self):
        return 'Satellite object {} at {}'.format(self.attrs['name'],
                                                  self.attrs['file'])

    def calc_ndens(self):
        '''
        Calculate number densities for each fluid.  Species mass is ascertained
        via recognition of fluid name (e.g. OpRho is clearly oxygen).  A full
        list of recognized fluids/species can be found by exploring the
        dictionary *mass* found in :mod:`~spacepy.pybats.bats`.  Composition is
        also calculated as percent of total number density.

        New values are saved using the keys *speciesN* (e.g. *opN*) and
        *speciesFrac* (e.g. *opFrac*).
        '''

        _calc_ndens(self)

    def calc_temp(self, units='eV'):
        '''
        Calculate plasma temperature for each fluid.  Number density is
        calculated using *calc_ndens* if it hasn't been done so already.

        Temperature is obtained via density and pressure through the simple
        relationship P=nkT.

        Use the units kwarg to set output units.  Current choices are
        KeV, eV, and K.  Default is eV.
        '''

        from spacepy.datamodel import dmarray

        units = units.lower()

        # Create dictionary of unit conversions.
        conv = {'ev': 6241.50935,   # nPa/cm^3 --> eV.
                'kev': 6.24150935,  # nPa/cm^3 --> KeV.
                'k': 72429626.47}   # nPa/cm^3 --> K.

        # Calculate number density if not done already.
        if 'N' not in self:
            self.calc_ndens()

        # Find all number density variables.
        for key in list(self.keys()):
            # Next variable if not number density:
            if key[-1] != 'N':
                continue
            # Next variable if no matching pressure:
            if not key[:-1] + 'p' in self:
                continue
            self[key[:-1] + 't'] = dmarray(
                conv[units] * self[key[:-1] + 'p']/self[key],
                attrs={'units': units})

    def calc_bmag(self):
        '''
        Calculates total magnetic field magnitude via:
        $|B|=\\sqrt{B_X^2+B_Y^2+B_Z^2}$.  Results stored as variable name *b*.
        '''

        if 'b' not in self:
            self['b'] = np.sqrt(self['bx']**2 +
                                self['by']**2 +
                                self['bz']**2)

        return True

    def calc_magincl(self, units='deg'):
        '''
        Magnetic inclination angle (a.k.a. inclination angle) is defined:
        $\\Theta = sin^{-1}(B_Z/B)$.  It is a crucial value when examining
        magnetic dynamics about geosychronous orbit, the tail, and other
        locations.

        This function calculates the magnetic inclination for the Virtual
        Satellite object and saves it as *b_incl*.  Units default to degrees;
        the keyword **units** can be changed to 'rad' to change this.
        '''

        if 'b' not in self:
            self.calc_bmag()

        incl = np.arcsin(self['bz'] / self['b'])

        if units == 'deg':
            self['b_incl'] = dmarray(incl * 180. / np.pi, {'units': 'degrees'})
        elif units == 'rad':
            self['b_incl'] = dmarray(incl, {'units': 'radians'})
        else:
            raise ValueError('Unrecognized units.  Use "deg" or "rad"')

        return True

    def get_position(self, time):
        '''
        For an arbitrary time, *time*, return a tuple of coordinates for
        the satellite's position at that time in GSM coordinates.

        *time* should be type datetime.datetime, matplotlib.dates.date2num
        output (i.e., number of days (fraction part represents hours,
        minutes, seconds) since 0001-01-01 00:00:00 UTC, *plus* *one*), or
        a sequence of either.
        The satellite's position is interpolated (linearly) to *time*.
        '''

        from matplotlib.dates import date2num

        # Test if sequence.
        if isinstance(time, (list, np.ndarray)):
            testval = time[0]
        else:
            testval = time

        # Test if datetime or not.
        if isinstance(testval, type(self['time'][0])):
            time = date2num(time)

        # Interpolate, using "try" as to not pass time limits and extrapolate.
        try:
            loc = (self._interp['x'](time),
                   self._interp['y'](time),
                   self._interp['z'](time))
        except ValueError:
            loc = (None, None, None)

        return loc

    def add_sat_loc(self, time, target, plane='XY', dobox=False,
                    dolabel=False, size=12, c='k', **kwargs):
        '''
        For a given axes, *target*, add the location of the satellite at
        time *time* as a circle.  If kwarg *dolabel* is True, the satellite's
        name will be used to label the dot.  Optional kwargs are any accepted
        by matplotlib.axes.Axes.plot.  The kwarg *plane* specifies the plane
        of the plot, e.g., 'XY', 'YZ', etc.
        '''

        # Get Axes' limits:
        xlim = target.get_xlim()
        ylim = target.get_ylim()

        plane = plane.lower()
        loc = self.get_position(time)
        if None in loc:
            return
        x = loc['xyz'.index(plane[0])]
        y = loc['xyz'.index(plane[1])]

        # Do not label satellite if outside axes bounds.
        if (x < min(xlim)) or (x > max(xlim)) or \
           (y < min(ylim)) or (y > max(ylim)):
            dolabel = False

        target.plot(x, y, 'o', **kwargs)
        if dolabel:
            xoff = 0.03 * (xlim[1] - xlim[0])
            if dobox:
                target.text(x + xoff, y, self.attrs['name'],
                            bbox={'fc': 'w', 'ec': 'k'},
                            size=size, va='center')
            else:
                target.text(x + xoff, y, self.attrs['name'],
                            size=size, va='center', color=c)

        # Restore Axes' limits.
        target.set_ylim(ylim)
        target.set_xlim(xlim)

    def add_orbit_plot(self, plane='XY', target=None, loc=111, rbody=1.0,
                       title=None, trange=None, add_grid=True, style='g.',
                       adjust_axes=True, add_arrow=True,
                       arrow_kwargs={'color': 'g', 'width': .05, 'ec': 'k'},
                       **kwargs):
        '''
        Create a 2D orbit plot in the given plane (e.g. 'XY' or 'ZY').
        Extra kwargs are handed to the plot function.


        Parameters
        ==========
        None

        Other Parameters
        ================
        plane : string
           Set the plane in which to plot (XY, XZ, etc.)  Defaults to XY.
        target : Matplotlib Figure or Axes object
           Set plot destination.  Defaults to new figure.
        loc : 3-digit integer
           Set subplot location.  Defaults to 111.
        adjust_axes : bool
           If True, axes will be customized to best display orbit (equal
           aspect ratio, grid on, planet drawn, etc.).  Defaults to True.
        style: string
           A matplotlib line style specifier.  Defaults to 'g.'.
        title : string
           Set title of axes.
        rbody : real
           Set radius of model inner boundary.  Defaults to 1.0
        trange : list of datetimes
           Set the time range to plot.  Defaults to None, or entire dataset.
        add_grid : bool
           Turn on or off grid style of axes.  Default is True.
        add_arrow : bool
           Add arrow at end of orbit path to indicate direction.
           Default is True.
        arrow_kwargs : dict
           Dictionary of arrow kwargs.  Defaults to match line style.
        '''

        from spacepy.pybats import add_body
        from spacepy.pybats.ram import grid_zeros

        fig, ax = set_target(target, figsize=(5, 5), loc=loc)

        plane = plane.upper()

        # Set time range of plot.
        if not trange:
            trange = [self['time'].min(), self['time'].max()]
        tloc = (self['time'] >= trange[0]) & (self['time'] <= trange[-1])

        # Extract orbit X, Y, or Z.
        plane = plane.lower()
        if plane[0] in ['x', 'y', 'z']:
            x = self[plane[0]][tloc]
        else:
            raise ValueError('Bad dimension specifier: ' + plane[0])
        if (plane[1] in ['x', 'y', 'z']) and (plane[0] != plane[1]):
            y = self[plane[1]][tloc]
        else:
            raise ValueError('Bad dimension specifier: ' + plane[1])

        # Actually plot:
        ax.plot(x, y, style, **kwargs)

        # Add arrow to plot:
        if add_arrow:
            x_arr, y_arr = x[-1], y[-1]
            dx, dy = x[-1]-x[-2], y[-1] - y[-2]
            ax.arrow(x_arr, y_arr, dx, dy, **arrow_kwargs)

        # Finish customizing axis.
        if adjust_axes:
            ax.axis('equal')
            ax.set_xlabel('GSM %s' % (plane[0].upper()))
            ax.set_ylabel('GSM %s' % (plane[1].upper()))
            if title:
                ax.set_title(title)
            if add_grid:
                ax.grid()
            grid_zeros(ax)
            add_body(ax, rad=rbody, add_night=('X' in plane.upper()))

        return fig, ax
