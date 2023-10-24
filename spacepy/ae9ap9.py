#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module for reading and dealing with AE9AP9 data files.

See https://www.vdl.afrl.af.mil/programs/ae9ap9/ to download the model.
This is not a AE9AP9 runner.

Authors: Brian Larsen, Steve Morley
Institution: Los Alamos National Laboratory
Contact: balarsen@lanl.gov

Copyright 2015 Los Alamos National Security, LLC.

This module provides a convenient class for handling data from
AE9/AP9 (and legacy models provided by the software).

.. rubric:: Class
.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    Ae9Data

Though the class is derived from SpacePy's SpaceData, the class also provides several methods
targeted at the AE9/AP9 output. Additional functions for working with the data are provided.

.. rubric:: Functions
.. autosummary::
    :template: clean_function.rst
    :toctree: autosummary

    readFile
    parseHeader

"""

import datetime as dt
import gzip
import os
import re
import warnings

from dateutil import relativedelta
import numpy as np

import spacepy.coordinates as spc
import spacepy.datamodel as dm
import spacepy.plot as splot
import spacepy.time as spt
import spacepy.toolbox as tb

__contact__ = 'Brian Larsen, balarsen@lanl.gov'


class Ae9Data(dm.SpaceData):
    """Dictionary-like container for AE9/AP9/SPM data, derived from SpacePy's datamodel

    To inspect the variables within this class, use the tree method.
    To export the data to a CDF, HDF5 or JSON-headed ASCII file use the relevant "to" method
    (toCDF, toHDF5, toJSONheadedASCII).

    .. autosummary::

        ~Ae9Data.getLm
        ~Ae9Data.plotOrbit
        ~Ae9Data.plotSummary
        ~Ae9Data.plotSpectrogram
        ~Ae9Data.setUnits

    .. automethod:: getLm
    .. automethod:: plotOrbit
    .. automethod:: plotSummary
    .. automethod:: plotSpectrogram
    .. automethod:: setUnits
    """

    def setUnits(self, per=None):
        """Set units of energy and flux/fluence

        If keyword 'per' is set to None, this method reports the units currently set.
        To set energy in MeV and flux/fluence in 'per MeV', set 'per=MeV'. Valid options are
        'eV', 'keV', 'Mev' and 'GeV'.

        Other Parameters
        ----------------
        per : string (optional)
            Energy units for both energy and flux/fluence
        """

        curr = self['Energy'].attrs['UNITS']
        particle_var = self.attrs['varname']
        if not per:
            print('Energy units: {0}'.format(curr))
            print('{0} units: {1}'.format(particle_var, self[particle_var].attrs['UNITS']))
            return
        unitlist = ['eV', 'keV', 'MeV', 'GeV']
        faclist = [1, 1e3, 1e6, 1e9]
        if per not in unitlist:
            raise ValueError("Units of {0} are not supported: Valid options are {1}".format(per, unitlist))
        else:
            if per == curr: return
            unitidx_from = unitlist.index(curr)
            unitidx_to = unitlist.index(per)
            self[particle_var] /= faclist[unitidx_from]  # convert to eV
            self[particle_var] *= faclist[unitidx_to]  # convert to target units
            self[particle_var].attrs['UNITS'] = self[particle_var].attrs['UNITS'].replace(curr, per)
            self['Energy'] *= faclist[unitidx_from]  # convert to eV
            self['Energy'] /= faclist[unitidx_to]  # convert to target units
            self['Energy'].attrs['UNITS'] = self['Energy'].attrs['UNITS'].replace(curr, per)

    def getLm(self, alpha=[90], model='T89'):
        """Calculate McIlwain L for the imported AE9/AP9 run and add to object
        """
        import spacepy.irbempy as ib
        ticks = spt.Ticktock(self['Epoch'])
        loci = spc.Coords(self['Coords'], self['Coords'].attrs['COORD_SYS'], 'car')
        loci.ticks = ticks
        retvals = ib.get_Lm(ticks, loci, alpha, extMag=model)
        self['Lm'] = dm.dmarray(retvals['Lm'].squeeze(), attrs={ 'MODEL': model })

    def _makeOrbitAxis(self, ser1, ser2, ax_target):
        """Helper function. Not intended for direct use.
        
        cx is the set of coordinates for plotting, ax_target is an axes object"""
        l1 = ax_target.plot(ser1, ser2)
        if np.abs((ax_target.get_xlim()[1] - ax_target.get_xlim()[0])) < 1:
            refpt = np.abs(np.max(ax_target.get_xlim()))
            ax_target.set_xlim([-1.25 * refpt, 1.25 * refpt])
            ax_target.set_ylim(ax_target.get_xlim())
            l1[0].set_marker('o')
            c1 = splot.plt.Circle([0, 0], radius=1.0, fc='none', ec='k')
            ax_target.add_artist(c1)
        else:
            splot.dual_half_circle(ax=ax_target)
        return ax_target

    def plotSummary(self, timerange=None, coord_sys=None, fig_target=None, spec=False, orbit_params=(False, True),
                    **kwargs):
        """Generate summary plot of AE9/AP9/SPM data loaded
        
        spec : if True, plot spectrogram instead of flux/fluence lineplot, requires 'ecol' keyword
        """
        if timerange:
            if isinstance(timerange, spt.Ticktock):
                t1 = timerange.UTC[0]
                t2 = timerange.UTC[-1]
            elif isinstance(timerange[0], dt.datetime):
                t1 = timerange[0]
                t2 = timerange[-1]
            else:
                raise TypeError('Incorrect data type provided for timerange')
            # now select subset
            i_use = tb.tOverlapHalf([t1, t2], self['Epoch'])
            t_use = self['Epoch'][i_use]
            c_use = self['Coords'][i_use]
            f_use = self[self.attrs['varname']][i_use, ...]
        else:
            t_use = self['Epoch']
            c_use = self['Coords']
            f_use = self[self.attrs['varname']]

        if coord_sys and (coord_sys.upper() != c_use.attrs['COORD_SYS']):
            # TODO: We assume cartesian, make flexible so can take spherical
            cx = spc.Coords(c_use, c_use.attrs['COORD_SYS'], 'car')
            cx.ticks = spt.Ticktock(t_use)
            cx = cx.convert(coord_sys.upper(), 'car').data
        else:
            coord_sys = c_use.attrs['COORD_SYS']
            cx = c_use
        sys_subs = r'$_{' + coord_sys + r'}$'

        splot.style('spacepy')
        if orbit_params[0]:
            landscape = orbit_params[1]
            locs = [121, 122] if landscape else [211, 212]
        else:
            locs = [223, 224]
            landscape = True
        if fig_target:
            fig, ax1 = splot.set_target(fig_target, loc=locs[0])
        else:
            fig, ax1 = splot.set_target(fig_target, figsize=(8, 8), loc=locs[0])
        fig, ax2 = splot.set_target(fig, loc=locs[1])

        ax1 = self._makeOrbitAxis(cx[:, 0], cx[:, 1], ax1)
        ax2 = self._makeOrbitAxis(cx[:, 0], cx[:, 2], ax2)
        if landscape: ax1.set_xlabel('X{0} [{1}]'.format(sys_subs, self['Coords'].attrs['UNITS']))
        ax2.set_xlabel('X{0} [{1}]'.format(sys_subs, self['Coords'].attrs['UNITS']))
        ax1.set_ylabel('Y{0} [{1}]'.format(sys_subs, self['Coords'].attrs['UNITS']))
        ax2.set_ylabel('Z{0} [{1}]'.format(sys_subs, self['Coords'].attrs['UNITS']))
        ax1xl, ax1yl, ax2xl, ax2yl = ax1.get_xlim(), ax1.get_ylim(), ax2.get_xlim(), ax2.get_ylim()
        maxabslim = np.max(np.abs([ax1xl, ax1yl, ax2xl, ax2yl]))
        if np.abs((ax1.get_xlim()[1] - ax1.get_xlim()[0])) < 1:
            refpt = maxabslim
            ax1.set_xlim([-1.25 * refpt, 1.25 * refpt])
            ax1.set_ylim(ax1.get_xlim())
            refpt = ax2.get_xlim()[0]
            ax2.set_xlim([-1.25 * refpt, 1.25 * refpt])
            ax2.set_ylim(ax1.get_xlim())
        else:
            ax1.set_xlim([-maxabslim, maxabslim])
            ax1.set_ylim([-maxabslim, maxabslim])
            ax2.set_xlim([-maxabslim, maxabslim])
        ax2.set_ylim(ax2.get_xlim())
        ax1.invert_yaxis()
        ax1.invert_xaxis()
        ax2.invert_xaxis()
        ax1.set_aspect('equal')
        ax2.set_aspect('equal')

        if not orbit_params[0]:
            ax3 = splot.plt.subplot2grid((2, 2), (0, 0), colspan=2)
            if not spec:
                l3 = ax3.semilogy(t_use, f_use)
                ylab = '{0} ['.format(self.attrs['varname']) + re.sub(r'(\^[\d|-]*)+', _grp2mathmode,
                                                                      self[self.attrs['varname']].attrs['UNITS']) + ']'
                ax3.set_ylabel(ylab)
                for ll, nn in zip(l3, self['Energy']):
                    ll.set_label('{0} {1}'.format(nn, self['Energy'].attrs['UNITS']))
                ncol = len(self['Energy']) // 2 if len(self['Energy']) <= 6 else 3
                leg = ax3.legend(loc='center', bbox_to_anchor=(0.5, 1), ncol=ncol,
                                 frameon=True, framealpha=0.5)
                lims3 = ax3.get_ylim()
                newupper = 10 ** (np.log10(lims3[0]) + (np.log10(lims3[1] / lims3[0]) * 1.125))
                ax3.set_ylim([lims3[0], newupper])
                splot.applySmartTimeTicks(ax3, t_use)
                fig.tight_layout()
                pos3 = ax3.get_position()
                ax3.set_position([pos3.x0, pos3.y0, pos3.width, pos3.height * 0.8])
                # fig.suptitle('{model_type}\n'.format(**self.attrs) +
                #             '{0} - {1}'.format(t_use[0].isoformat()[:19], t_use[-1].isoformat()[:19]))
            else:
                if timerange: raise NotImplementedError('Time range selection not yet implemented for spectrograms')
                pos3 = ax3.get_position()
                ax3.set_position([pos3.x0, pos3.y0, pos3.width, pos3.height * 0.9])
                ecol = kwargs['ecol'] if 'ecol' in kwargs else 0
                ax3 = self.plotSpectrogram(target=ax3, ecol=ecol)
                splot.plt.subplots_adjust(wspace=0.3)
                pos1 = ax1.get_position()
                ax1.set_position([pos3.x0, pos1.y0, pos1.width, pos1.height])

        splot.revert_style()
        return fig

    def plotOrbit(self, timerange=None, coord_sys=None, landscape=True, fig_target=None):
        """
        Plot X-Y and X-Z projections of satellite orbit in requested coordinate system
        """
        fig = self.plotSummary(timerange=timerange, coord_sys=coord_sys,
                               fig_target=fig_target, orbit_params=(True, landscape))
        return fig

    def plotSpectrogram(self, ecol=0, pvars=None, **kwargs):
        '''
        Plot a spectrogram of the flux along the requested orbit, as a function of Lm and time

        Other Parameters
        ----------------
        zlim : list
            2-element list with upper and lower bounds for color scale
        colorbar_label : string
            text to appear next to colorbar (default is 'Flux' plus the units)
        ylabel : string
            text to label y-axis (default is 'Lm' plus the field model name)
        title : string
            text to appear above spectrogram (default is climatology model name, data type and energy)
        pvars : list
            list of plotting variable names in order [Epoch-like (X axis), Flux-like (Z axis), Energy (Index var for Flux-like)]
        ylim : list
            2-element list with upper and lower bounds for y axis
        '''
        import spacepy.plot as splot
        if 'Lm' not in self:
            self.getLm()
        sd = dm.SpaceData()

        if pvars is None:
            varname = self.attrs['varname']
            enname = 'Energy'
            epvar = 'Epoch'
        else:
            epvar = pvars[0]
            varname = pvars[1]
            enname = pvars[2]
        if len(self['Lm']) != len(self[epvar]): #Lm needs interpolating to new timebase
            import matplotlib.dates as mpd
            tmptimetarg, tmptimesrc = mpd.date2num(self[epvar]), mpd.date2num(self['Epoch'])
            sd['Lm'] = dm.dmarray(np.interp(tmptimetarg, tmptimesrc, self['Lm'], left=np.nan, right=np.nan))
        else:
            sd['Lm'] = self['Lm']
        #filter any bad Lm
        goodidx = sd['Lm']>1
        sd['Lm'] = sd['Lm'][goodidx]
        if 'ylim' in kwargs:
            Lm_lim = kwargs['ylim']
            del kwargs['ylim']
        else:
            Lm_lim = [2.0,8.0]
        #TODO: allow user-definition of bins in time and Lm
        sd['Epoch'] = dm.dmcopy(self[epvar])[goodidx] #TODO: assumes 1 pitch angle, generalize
        try:
            sd['1D_dataset'] = self[varname][goodidx,ecol] #TODO: assumes 1 pitch angle, generalize
        except IndexError: #1-D
            sd['1D_dataset'] = self[varname][goodidx]
        bins = [[],[]]
        if 'tbins' not in kwargs:
            bins[0] = spt.tickrange(self[epvar][0], self[epvar][-1], 3/24.).UTC
        else:
            bins[0] = kwargs['tbins']
            del kwargs['tbins']
        if 'Lbins' not in kwargs:
            bins[1] = np.arange(Lm_lim[0], Lm_lim[1], 1./4.)
        else:
            bins[1] = kwargs['Lbins']
            del kwargs['Lbins']
        spec = splot.spectrogram(sd, variables=['Epoch', 'Lm', '1D_dataset'], ylim=Lm_lim, bins=bins)

        if 'zlim' not in kwargs:
            zmax = 10 ** (int(np.log10(max(sd['1D_dataset']))) + 1)
            idx = np.logical_and(sd['1D_dataset'] > 0, sd['Lm'] > Lm_lim[0])
            idx = np.logical_and(idx, sd['Lm'] <= Lm_lim[1])
            zmin = 10 ** int(np.log10(min(sd['1D_dataset'][idx])))
            kwargs['zlim'] = [zmin, zmax]
        if 'colorbar_label' not in kwargs:
            flux_units = self[varname].attrs['UNITS']
            kwargs['colorbar_label'] = '{0} ['.format(varname) + re.sub(r'(\^[\d|-]*)+', _grp2mathmode, flux_units) + ']'
        if 'ylabel' not in kwargs:
            kwargs['ylabel'] = 'L$_M$' + ' ' + '[{0}]'.format(self['Lm'].attrs['MODEL'])
        if 'title' not in kwargs:
            kwargs['title'] = '{model_type} {varname}: '.format(**self.attrs) + \
                              '{0:5.2f} {1}'.format(self[enname][ecol], self[enname].attrs['UNITS'])
        reset_shrink = splot.mpl.mathtext.SHRINK_FACTOR
        splot.mpl.mathtext.SHRINK_FACTOR = 0.85
        splot.mpl.mathtext.GROW_FACTOR = 1 / 0.85
        ax = spec.plot(cmap='plasma', **kwargs)
        splot.mpl.mathtext.SHRINK_FACTOR = reset_shrink
        splot.mpl.mathtext.GROW_FACTOR = 1 / reset_shrink
        return ax


def readFile(fname, comments='#'):
    """
    read a model generated file into a datamodel.SpaceData object

    Parameters
    ==========
    fname : str
        filename of the file

    Other Parameters
    ================
    comments : str (optional)
        String that is the comments in the data file

    Returns
    =======
    out : :mod:`~spacepy.datamodel.SpaceData`
        Data contained in the  file

    Examples
    ========
    >>> from spacepy import ae9ap9
    >>> ae9ap9.readFile('ephem_sat.dat').tree(verbose=1)
    +
    |____Epoch (spacepy.datamodel.dmarray (121,))
    |____Coords (spacepy.datamodel.dmarray (121, 3))
    |____MJD (spacepy.datamodel.dmarray (121,))
    |____posComp (spacepy.datamodel.dmarray (3,))
    """
    if not os.path.isfile(fname):
        raise (ValueError("File {0} not found".format(fname)))
    # get the header information first
    header = parseHeader(fname)
    # and read in all the data
    data = np.loadtxt(fname, delimiter=header['delimiter'], comments=comments)
    # now parse the data    
    ans = Ae9Data()
    # parse the datetime if it is there (it is always first)
    if 'datetime' in header['columns'][0]:
        if header['time_format'] == 'eDOY':  # have to massage the data first
            year = data[:, 0].astype(int)
            frac = data[:, 1]
            time = spt.Ticktock([dt.datetime(y, 1, 1) + relativedelta.relativedelta(days=v)
                                 for y, v in zip(year, frac)], 'UTC')
            ans[header['time_format']] = dm.dmarray(data[:, 0:2])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            ans[header['time_format']].attrs['LABL_PTR_1'] = 'timeComp'
            ans['timeComp'] = dm.dmarray(['Year', 'eDOY'])
            ans['timeComp'].attrs['VAR_TYPE'] = 'metadata'
            del header['columns'][0:2]
            data = data[:, 2:]
        elif header['time_format'] == 'UTC':  # have to massage the data first
            tm = data[:, 0:6].astype(int)
            time = spt.Ticktock([dt.datetime(*v) for v in tm], 'UTC')
            ans[header['time_format']] = dm.dmarray(data[:, 0:6])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            ans[header['time_format']].attrs['LABL_PTR_1'] = 'timeComp'
            ans['timeComp'] = dm.dmarray(['Year', 'Month', 'Day', 'Hour', 'Minute', 'Second'])
            ans['timeComp'].attrs['VAR_TYPE'] = 'metadata'
            del header['columns'][0:6]
            data = data[:, 6:]
        else:
            time = spt.Ticktock(data[:, 0], header['time_format'])
            ans[header['time_format']] = dm.dmarray(data[:, 0])
            ans[header['time_format']].attrs['VAR_TYPE'] = 'support_data'
            del header['columns'][0]
            data = data[:, 1:]

        ans['Epoch'] = dm.dmarray(time.UTC)
        ans['Epoch'].attrs['VAR_TYPE'] = 'support_data'
    # parse the position, it is always next
    if 'posx' in header['columns'][0]:
        varname = 'Coords'
        pos = dm.dmarray(data[:, 0:3])
        ans[varname] = pos
        ans[varname].attrs['UNITS'] = header['coord_system'][1]
        ans[varname].attrs['FIELDNAM'] = varname
        ans[varname].attrs['LABLAXIS'] = 'Position'
        ans[varname].attrs['DEPEND_0'] = 'Epoch'
        ans[varname].attrs['LABL_PTR_1'] = 'posComp'
        ans[varname].attrs['VAR_TYPE'] = 'data'
        ans[varname].attrs['COORD_SYS'] = header['coord_system'][0]
        ans['posComp'] = dm.dmarray(['X', 'Y', 'Z'])
        ans['posComp'].attrs['VAR_TYPE'] = 'metadata'
        del header['columns'][0:3]
        data = data[:, 3:]

    # if there are PA they are now
    if header['columns'] and ('pitch' in header['columns'][0]):
        raise (NotImplementedError("Sorry have not implemented pitch angle resolved files yet"))
    ## This is an issue and the files have repeadedlines at the same time for each pitch
    ##  angle. They would all need to be read in them combines into a multi-d array

    # now parse fluence, flux, or doserate
    # do all the rest of the column headers match?
    col_arr = np.asarray(header['columns'])
    if header['columns'] and (col_arr == col_arr[0]).all():
        varname = col_arr[0][0].title()
        ans[varname] = dm.dmarray(data)
        ans[varname].attrs['UNITS'] = str(col_arr[0][1])
        ans[varname].attrs['FIELDNAM'] = varname
        ans[varname].attrs['LABLAXIS'] = varname
        ans[varname].attrs['VAR_TYPE'] = 'data'
        ans[varname].attrs['SCALETYP'] = 'log'
        header['varname'] = varname
        ans[varname].attrs[
            'Description'] = '{flux_direction} {varname} based on {flux_type} from the {model_type} model'.format(
            **header)
        ans[varname].attrs['DEPEND_0'] = 'Epoch'
        if ans[varname].shape[1] == header['energy'][0].shape[0]:
            ans[varname].attrs['DEPEND_1'] = 'Energy'
        if 'percentile' in header:
            ans[varname].attrs['TITLE'] = '{0} percentile'.format(header['percentile'])

    # create the energy variable
    if 'energy' in header:
        ans['Energy'] = dm.dmarray(header['energy'][0])
        ans['Energy'].attrs['UNITS'] = header['energy'][1]
        ans['Energy'].attrs['FIELDNAM'] = 'Energy'
        ans['Energy'].attrs['LABLAXIS'] = 'Energy'
        ans['Energy'].attrs['VAR_TYPE'] = 'data'
        ans['Energy'].attrs['SCALETYP'] = 'log'
        ans['Energy'].attrs['VAR_TYPE'] = 'support_data'

    skip_header = ['columns', 'coord_system', 'delimiter', 'energy', 'flux_direction']
    for k in header:
        if k not in skip_header:
            ans.attrs[k] = header[k]
    return ans


def _parseInfo(header):
    """
    given a header parse and return the common information in all headers

    
    # Time format:              Modified Julian Date
    # Coordinate system:        GEI (Geocentric Equatorial Inertial) Cartesian in Earth radii
    # Data Delimiter:           comma
    """
    # get the time format, it is in all the files
    ans = { }
    for ind, val in enumerate(header):
        if "Time format:" in val:
            if "Modified Julian Date" in val:
                ans['time_format'] = 'MJD'
            elif "Year, day_of_year.frac" in val:
                ans['time_format'] = 'eDOY'
            elif "Year, Month, Day, Hour, Minute, Seconds" in val:
                ans['time_format'] = 'UTC'
            else:
                raise (NotImplementedError("Sorry can't read that time format yet: {0}".format(val)))
        elif "Coordinate system:" in val:
            coord_sys = val.split(':')[1].strip().split()[0]
            if "in Earth radii" in val:
                ans['coord_system'] = (coord_sys, 'Re')
            else:
                ans['coord_system'] = (coord_sys, '')
        elif "Data Delimiter:" in val:
            if "comma" in val:
                ans['delimiter'] = ','
            else:
                raise (NotImplementedError("Sorry can't read that delimiter yet"))
        elif "Ae9Ap9 Software Version:" in val:
            match = re.match(r'^Ae9Ap9 Software Version:.*(\d\.\d\d\.\d\d\d)$', val)
            version = match.group(1).strip()
            ans['software_version'] = version
        elif "Model type:" in val:
            match = re.match(r'^Model type:(.*)$', val)
            mtype = match.group(1).strip()
            ans['model_type'] = mtype
        elif "CRRESELE Software Version:" in val:
            match = re.match(r'^CRRESELE Software Version:.*(\d\.\d\d\.\d\d\d)$', val)
            version = match.group(1).strip()
            ans['software_version'] = version
            ans['model_type'] = 'CRRESELE/PRO'
        elif "Flux type:" in val:
            match = re.match(r'^Flux type:(.*)$', val)
            ftype = match.group(1).strip()
            ans['flux_type'] = ftype
        elif "Flux direction:" in val:
            if '=' not in val:
                match = re.match(r'^Flux direction:(.*)$', val)
                fdir = match.group(1).strip()
            else:
                pa = val.split('=')[1].strip()
                fdir = np.asarray(pa.split())
            ans['flux_direction'] = fdir
        elif "Energy levels" in val:
            match = re.match(r'^Energy levels.*\((.*)\):(.*)$', val)
            ans['energy'] = (np.asarray(match.group(2).strip().split()).astype(float), match.group(1).strip())
        # Get the orbital element propagator, two versions based on AE9 model changes
        # new format data file
        elif "Propagator" in val:  # New format
            match = re.search(r'Propagator:\ (.*)$', val)
            ans['propagator'] = match.group(1).strip()
        elif "generated from specified elements" in val:  # In both old and new
            match = re.search(r'^generated from specified elements.*:\ (.*)$', val)
            if match:  # But old has propagator on this line; process and warn
                warnings.warn(
                    "Support for orbit files from AE9AP9 model <1.5 is deprecated; please update to model 1.5 or later.",
                    DeprecationWarning)
                ans['propagator'] = match.group(1).strip()
    return ans


def _readHeader(fname):
    """
    read only the header from a Ae9Ap9 file
    """
    dat = []
    if fname.endswith('.gz'):
        try:
            fp = gzip.open(fname, 'rt')
        except ValueError: #Python 2.7 (Windows) compatibility
            fp = gzip.open(fname)
    else:
        fp = open(fname, 'rt')
    while True:
        tmp = fp.readline()
        if tmp[0] not in ('#', 35):
            break
        dat.append(tmp.strip())
    dat = [v[1:].strip() for v in dat]
    fp.close()
    return dat


def parseHeader(fname):
    """
    given an AE9AP9 output test file parse the header and return the information in a
    dictionary

        .. versionchanged:: 0.3.0

        The underlying AE9AP9 model changed the ephem file format and this
        reader was updated to match. Reading the old format will issue
        DeprecationWarning.

    Parameters
    ==========
    fname : str
        filename of the file

    Returns
    =======
    out : dict
        Dictionary of the header information in the file
    """
    ans = { }
    ans['filename'] = os.path.abspath(fname)

    header = _readHeader(fname)
    # get the information
    ans.update(_parseInfo(header))

    columns = header[-1].split(ans['delimiter'])
    columns_unq = _unique_elements_order(columns)
    len_columns_unq = len(columns)

    # pull the units off the column headings
    cols = []
    while columns_unq:
        tmp = columns_unq.pop(0)
        umatch = re.match(r'(.*)\((.*)\).*', tmp)
        if umatch:
            cols.append((umatch.group(1), umatch.group(2)))
        else:
            cols.append(tmp)
    # now go through and see if any non-tuple entry has the same name as a tuple entry
    to_remove = []
    for ii, c1 in enumerate(cols):
        if not isinstance(c1, tuple):
            for jj, c2 in enumerate(cols[ii:]):
                if c2[0] == c1:
                    to_remove.append(ii)
    if to_remove:
        for v in to_remove:
            del cols[v]
    # and now replicate the last entry enough times to make it the right length
    cols.extend([cols[-1]] * (len_columns_unq - len(cols)))

    ans['columns'] = cols

    if 'pctile' in os.path.basename(fname):  # capture the percentile
        match = re.match(r'.*_(\d\d).txt', os.path.basename(fname))
        ans['percentile'] = float(match.group(1))
    return ans


def _unique_elements_order(seq):
    """
    see http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
    """
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def combinePercentiles(files, timeframe='all', verbose=True):
    """
    combine files at different percentiles into one file with the spectra at different
    percentiles for easy plotting and analysis

    NOTE: Requires pandas for any timeframe other than 'all'

    Parameters
    ==========
    files : str
        Filenames of the percentile files

    Other Parameters
    ================
    timeframe : str
        Timeframe to average the input spectra over (either 'all' or a pandas understoop resample() time
    verbose : boolean
        Print out information while reading the files
    
    Returns
    =======
    out : :mod:`~spacepy.datamodel.SpaceData`
        Combined output spectra of the file
    """
    if not files:
        raise (ValueError("Must input files"))
    data = { }
    for fname in files:
        if verbose: print("Reading: {0}".format(fname))
        tmp = readFile(fname)
        if 'percentile' not in tmp.attrs:
            raise (ValueError("File {0} does not have a percentile key".format(fname)))
        data[tmp.attrs['percentile']] = tmp

    varnames = ['Flux', 'Fluence', 'Dose', None]
    for varname in varnames:
        if varname in data[list(data.keys())[0]]:
            break  # this leaves it where it was found
    if varname is None:
        raise (ValueError("Did not understand file type could not find one of {0}".format(varnames)))

    # now that they are all read in, we collapse them down to different time bases
    if timeframe == 'all':
        # this is just a full collapse
        # find the var to collapse
        d_comp = np.asarray([])
        ps = sorted(data.keys())[::-1]
        for p in ps:
            if data[p][varname].attrs['DEPEND_0'] == 'Epoch':
                d_comp = np.concatenate((d_comp, data[p][varname].mean(axis=0)))  # full time collapse
            else:
                d_comp = np.concatenate((d_comp, data[p][varname]))  # full time collapse
        d_comp = d_comp.reshape(len(ps), -1).T
    else:
        try:
            import pandas as pd
        except ImportError:
            raise (ImportError("panads is required for timeframe other than 'all'"))
        raise (NotImplementedError("Timeframes other than 'all' not yet implemented"))
        # make a dataframe
        # resample to timeframe
        # join everything
        # meet up at same place for output
    # now that we have all the data prep it to a spacedata
    ans = dm.SpaceData()
    ans[varname] = dm.dmarray(d_comp)
    if timeframe == 'all':
        ans[varname].attrs['DEPEND_0'] = 'Energy'
        ans[varname].attrs['DEPEND_1'] = 'Percentile'
        ans[varname].attrs['DISPLAY_TYPE'] = 'time_series'
    for k in data[p][varname].attrs:
        if 'DEPEND' in k:
            continue
        else:
            ans[varname].attrs[k] = data[p][varname].attrs[k]
    ans['Energy'] = data[p]['Energy']
    ans['Percentile'] = dm.dmarray(ps)
    ans['Percentile'].attrs['FIELDNAM'] = 'Percentile'
    ans['Percentile'].attrs['LABLAXIS'] = 'Percentile'
    ans['Percentile'].attrs['VALIDMIN'] = 0.0
    ans['Percentile'].attrs['VALIDMAX'] = 100.0
    ans['Percentile'].attrs['SCALETYP'] = 'support_data'
    return ans


def _getData(fnames):
    """
    helper routine for dm.toCDF and dm.toHDF5 and dm.toJSONheadedASCII since the data
    prep is all the same
    """
    if hasattr(fnames, 'strip'):  # it is a string
        return readFile(fnames)
    else:
        return combinePercentiles(fnames)


def _grp2mathmode(matchobj):
    """Helper function to latex-ify the units for plotting"""
    unitstr = matchobj.group()
    unitstr = unitstr[0] + '{' + unitstr[1:]
    return '$' + unitstr + '}$'
