#!/usr/bin/env python3

'''
Classes, functions, and methods for reading, writing, and plotting output
from the version of the Rice Convection Model embedded into the SWMF.

RCM follows *flux tube content*, $\eta_s$, at different *energy invariants*,
$\lambda_s$, across a grid of energy invariants and for three different species
(electrons, protons, and oxygen ions). The grid is ionospheric, however, the
GSM-coordinates of the equatorial crossing of field lines threading each
ionospheric grid point are given.

The density and energy can be obtained from the content and energy invariants
using the flux tube volume ($V$):

Energy:  $W_s = \lambda_s / V^{2/3}$
Density: $n_s = \eta_s / V$

These values are automatically generated upon object instantiation.
'''

import re

import numpy as np
import matplotlib.pyplot as plt

# from spacepy import deprecated
# from spacepy.plot import set_target
from spacepy.pybats import PbData, dmarray

filename = 'IM/plots/3d__max_t00000100.dat'


class Rcm3d(PbData):
    '''
    A class for handling 3D output from the RCM.
    Instantiate an object as follows:

    >>> iono = rcm.Rcm3d('filename.dat')

    ...where filename.dat is the name of a TecPlot-formatted ASCII output file.
    '''

    def __init__(self, filename, *args, **kwargs):
        super(Rcm3d, self).__init__(*args, **kwargs)
        self.attrs['file'] = filename
        self._read_tec()

    def _read_tec(self):
        '''
        Load a 3D tecplot IM file.
        '''

        nI, nJ, nK = 0, 0, 0

        # Open file and parse header:
        with open(self.attrs['file'], 'r') as f:
            for i in range(6):
                raw = f.readline()
                # Title line.
                if 'TITLE' in raw:
                    self.attrs['title'] = re.match('TITLE="(.*)"',
                                                   raw).groups()[0]
                # Get the size of the grid.
                if 'ZONE' in raw:
                    mstr = "ZONE.*I=\s+(\d+),\s+J=\s+(\d+),\s+K=\s+(\d+),"
                    match = re.match(mstr, raw)
                    nI, nJ, nK = [int(s) for s in match.groups()]
                # Parse variable names and units; skip IJK.
                if 'VARIABLES' in raw:
                    raw = raw.replace('"', '')  # Get rid of quotes.
                    parts = raw[10:].split(',')  # Break into parts.

                    # Break out into lists.
                    varnames = []
                    units = []
                    for p in parts[3:]:
                        # Look for matching units (e.g., [keV])
                        match = re.match('(.+)\[(.+)\]', p)
                        if match:
                            varnames.append(match.group(1).strip())
                            units.append(match.group(2).strip())
                        else:
                            # Default to blank units.
                            varnames.append(p.strip())
                            units.append('')

            # Slurp remaining contents.
            lines = f.readlines()

        # Build appropriate container:
        self.attrs['nI'], self.attrs['nJ'], self.attrs['nK'] = nI, nJ, nK
        for v, u in zip(varnames, units):
            self[v] = dmarray(np.zeros([nI, nJ, nK]), {'units': u})

        # Parse the rest of the file.
        for line in lines:
            # Scrub for bad entries.
            cleanline = re.sub('(\d)([+-]\d)', r'\1E\2', line)

            parts = cleanline.split()
            # Get indices of line:
            i, j, k = [int(s) - 1 for s in parts[:3]]
            # Parse values:
            for p, v in zip(parts[3:], varnames):
                self[v][i, j, k] = p

    def calc_E(self):
        '''
        Calculate energy using flux tube volume and energy invariant value.
        '''
        # Return if already calculated.
        if 'W' in self:
            return
        pass

    def calc_n(self):
        '''
        Calculate number density, in units of $cm^-3$, using flux tube content
        and flux tube volume.
        '''
        self['n'] = dmarray(self['EETA'] / self['V'], {'units': 'cm^-3'})
        pass