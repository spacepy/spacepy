"""International Geomagnetic Reference Field model

This module is intended primarily to support :mod:`~spacepy.coordinates`
rather than for direct end use, and the interface is subject to change.

Classes
-------
.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    IGRFCoefficients
    IGRF
"""


# std. lib.
import os
import datetime as dt
import warnings

# scientific stack
import numpy as np

# third-party or current
from spacepy import __path__ as basepath
from spacepy import DOT_FLN
import spacepy.time as spt


class IGRFCoefficients():
    """Read and store IGRF coefficients from data file

    Other Parameters
    ----------------
    fname : str, optional
        Filename to read from; defaults from the .spacepy data directory.
    """
    def __init__(self, fname=None):
        # Store coefficients and SV in nested arrays
        # This is triangular, e.g.,
        # g[2] has 3 elements (g[2][0:2])
        # h[5] has 5 elements (h[5][0:5])
        # Store top level as object arrays so the arrays inside
        # can have staggered lengths
        self.coeffs = {'g': np.empty(14, dtype=object),
                       'h': np.empty(14, dtype=object),
                       'g_SV': np.empty(14, dtype=object),
                       'h_SV': np.empty(14, dtype=object),
                       }

        # Open IGRF coefficients file...
        if fname is None:
            # Updates shoud be downloaded by toolbox.update and
            # put in the .spacepy/data directory
            # TODO: write getter for toolbox.update
            fname = os.path.join(DOT_FLN, 'data', 'igrfcoeffs.txt')
            if not os.path.exists(fname):
                # Fall back to IGRF13 coefficients
                fname = os.path.join('{0}'.format(basepath[0]), 'data', 'igrf13coeffs.txt')
        with open(fname) as fh:
            header = [fh.readline() for i in range(4)]
            data = fh.readlines()
        epochs = [float(yy) for yy in header[-1].strip().split()[3:-1]]
        self.epochs = np.array(epochs)
        self.datelow = dt.datetime(int(self.epochs[0]), 1, 1)
        self.lastepoch = dt.datetime(int(self.epochs[-1]), 1, 1)
        self.datehigh = dt.datetime(int(self.epochs[-1])+5, 1, 1)
        n_epochs = len(epochs)

        # Now we know how many epochs are in the current IGRF file,
        # we can construct the arrays to store the coefficients
        for n_idx in range(0, 14):
            self.coeffs['g'][n_idx] = np.empty((n_idx+1, n_epochs))
            self.coeffs['h'][n_idx] = np.empty((n_idx+1, n_epochs))
            self.coeffs['g_SV'][n_idx] = np.empty((n_idx+1))
            self.coeffs['h_SV'][n_idx] = np.empty((n_idx+1))

        # Parse the file and populat the coefficient arrays
        for line in data:
            line = line.strip().split()
            # Read g and h coefficients into tables
            vals = [float(v) for v in line[3:]]
            n_idx = int(line[1])
            m_idx = int(line[2])
            self.coeffs[line[0]][n_idx][m_idx, ...] = vals[:-1]
            svtag = 'g_SV' if line[0] == 'g' else 'h_SV'
            self.coeffs[svtag][n_idx][m_idx] = vals[-1]

        # set arrays to read only
        for key in ['g', 'h', 'g_SV', 'h_SV']:
            for arr in self.coeffs[key]:
                arr.setflags(write=0)


# load coefficients so it doesn't get done on instantiation of IGRF class
igrfcoeffs = IGRFCoefficients()


class IGRF():
    """International Geomagnetic Reference Field model

    Notes
    -----

    .. versionadded:: 0.3.0

    .. rubric:: Methods

    .. autosummary::

        ~IGRF.calcDipoleAxis
        ~IGRF.initialize

    .. rubric:: Data

    .. autosummary::

        dipole
        moment

    .. automethod:: calcDipoleAxis
    .. automethod:: initialize

    .. attribute:: dipole

        Characteristics of dipole (`dict`).

    .. attribute:: moment

        Dipole moments (`dict`).
    """
    dipole = {}
    """Characteristics of dipole (`dict`)."""

    moment = {}
    """Dipole moments (`dict`)."""

    def __init__(self):
        self.__status = {'coeffs': False,
                         'init': False,
                         'time_set': False,
                         }
        self.__coeffs = igrfcoeffs.coeffs
        self.__status['coeffs'] = True

    def initialize(self, time, limits='warn'):
        """Initialize model state to a particular time.

        Parameters
        ----------
        time : `~datetime.datetime`
            Time for which to initialize the model

        Other Parameters
        ----------------
        limits : `str`, optional
            Set to ``warn`` to warn about out-of-range times (default);
            any other value to error.
        """
        errmsg = 'IGRF: Requested time is outside valid range.\n'
        errmsg += 'Valid range is [{0}, {1}]'.format(igrfcoeffs.datelow, igrfcoeffs.datehigh)
        self.time = time
        if time < igrfcoeffs.datelow or time > igrfcoeffs.datehigh:
            if limits.lower() == 'warn':
                errmsg += '\nProceeding using effective date {0}'.format(igrfcoeffs.datehigh)
                warnings.warn(errmsg)
                self.time = igrfcoeffs.datehigh
            else:
                errmsg += '''\nUse "limits='warn'" to force out-of-range times '''
                errmsg += 'to use the nearest valid endpoint.'
                raise ValueError(errmsg)
        self.__valid_extrap = True if (self.time <= igrfcoeffs.datehigh
                                       and self.time > igrfcoeffs.lastepoch) else False
        self.__status['time_set'] = False
        self.__status['init'] = True

        self.calcDipoleAxis()
        self.__status['time_set'] = True

    def calcDipoleAxis(self):
        """Calculates dipole axis for initialized time.

        Populates :data:`moment` and :data:`dipole`.
        """
        if self.__status['time_set']:
            return
        if not self.__status['init']:
            warnings.warn('IGRF: Initialize for a time before invoking this method.')
            return

        # Find the indices for the epochs surrounding input time
        utc = self.time
        tstamp = spt.Ticktock(utc, 'UTC').MJD[0]

        # Compute the various IGRF dependent things like position of CD
        # TODO: So that full IGRF can be used eventually, for each in
        # __coeffs, interp to time and save as '[g|h][n][m]' e.g., 'h32' (entries in self.coeffs)
        if not self.__valid_extrap:
            # Interpolate coefficients linearly
            numepochs = [spt.Ticktock(dt.datetime(int(val), 1, 1)).MJD[0] for val in igrfcoeffs.epochs]
            g10 = np.interp(tstamp, numepochs, self.__coeffs['g'][1][0])
            g11 = np.interp(tstamp, numepochs, self.__coeffs['g'][1][1])
            h11 = np.interp(tstamp, numepochs, self.__coeffs['h'][1][1])
        else:
            # Extrapolate using secular variation
            ref_ep = spt.Ticktock(dt.datetime(int(igrfcoeffs.epochs[-1]), 1, 1)).MJD[0]
            diff_years = ((tstamp % 1 + tstamp//1) - (ref_ep % 1 + ref_ep//1))/365.25 # in years
            g10 = self.__coeffs['g_SV'][1][0]*diff_years + self.__coeffs['g'][1][0][-1]
            g11 = self.__coeffs['g_SV'][1][1]*diff_years + self.__coeffs['g'][1][1][-1]
            h11 = self.__coeffs['h_SV'][1][1]*diff_years + self.__coeffs['h'][1][1][-1]

        mom_sq = g10*g10 + g11*g11 + h11*h11
        mom = np.sqrt(mom_sq)

        # Compute dipole moments.
        self.moment = {'cd': mom,
                       'cd_McIlwain': 31165.3,
                       'cd_2010': 29950.126496,
                       }

        self.dipole = dict()
        self.dipole['cd_gcolat_rad'] = np.pi - np.arccos(g10/mom)
        self.dipole['cd_glon_rad'] = np.arctan(h11/g11)
        self.dipole['cd_gcolat'] = np.rad2deg(self.dipole['cd_gcolat_rad'])
        self.dipole['cd_glon'] = np.rad2deg(self.dipole['cd_glon_rad'])
