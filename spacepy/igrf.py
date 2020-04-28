# std. lib.
import os
import datetime as dt
import warnings

# scientific stack
import numpy as np

# third-party or current
from spacepy import __path__ as basepath
from spacepy import DOT_FLN

class IGRF():
    def __init__(self):
        self.__status = {'coeffs': False,
                         'init': False,
                         'time_set': False,
                         }
        self._readCoefficients()

    def _readCoefficients(self, fname=None):
        '''Read IGRF coefficients from data file
        '''
        # Store coefficients and SV in nested arrays
        # This is triangular, e.g.,
        # g[2] has 3 elements (g[2][0:2])
        # h[5] has 5 elements (h[5][0:5])
        # Store top level as object arrays so the arrays inside
        # can have staggered lengths
        self.__coeffs = {'g': np.empty(14, dtype=object),
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
        self._epochs = np.array(epochs)
        self._datelow = dt.datetime(int(self._epochs[0]), 1, 1)
        self._lastepoch = dt.datetime(int(self._epochs[-1]), 1, 1)
        self._datehigh = dt.datetime(int(self._epochs[-1])+5, 1, 1)
        n_epochs = len(epochs)

        # Now we know how many epochs are in the current IGRF file,
        # we can construct the arrays to store the coefficients
        for n_idx in range(0, 14):
            self.__coeffs['g'][n_idx] = np.empty((n_idx+1, n_epochs))
            self.__coeffs['h'][n_idx] = np.empty((n_idx+1, n_epochs))
            self.__coeffs['g_SV'][n_idx] = np.empty((n_idx+1))
            self.__coeffs['h_SV'][n_idx] = np.empty((n_idx+1))

        # Parse the file and populat the coefficient arrays
        for line in data:
            line = line.strip().split()
            # Read g and h coefficients into tables
            vals = [float(v) for v in line[3:]]
            n_idx = int(line[1])
            m_idx = int(line[2])
            self.__coeffs[line[0]][n_idx][m_idx, ...] = vals[:-1]
            svtag = 'g_SV' if line[0] == 'g' else 'h_SV'
            self.__coeffs[svtag][n_idx][m_idx] = vals[-1]

        self.__status['coeffs'] = True

    def initialize(self, time, limits='warn'):
        errmsg = 'IGRF: Requested time is outside valid range.\n'
        errmsg += 'Valid range is [{0}, {1}]'.format(self._datelow, self._datehigh)
        self.time = time
        if time < self._datelow or time > self._datehigh:
            if limits.lower() == 'warn':
                errmsg += '\nProceeding using effective date {0}'.format(self._datehigh)
                warnings.warn(errmsg)
                self.time = self._datehigh
            else:
                errmsg += '''\nUse "limits='warn'" to force out-of-range times '''
                errmsg += 'to use the nearest valid endpoint.'
                raise ValueError(errmsg)
        self.__valid_extrap = True if (self.time <= self._datehigh
                                       and self.time > self._lastepoch) else False
        self.__status['time_set'] = False
        self.__status['init'] = True

        self.calcDipoleAxis()
        self.__status['time_set'] = True

    def calcDipoleAxis(self):
        if self.__status['time_set']:
            return
        if not self.__status['init']:
            warnings.warn('IGRF: Initialize for a time before invoking this method.')
            return

        # Find the indices for the epochs surrounding input time
        utc = self.time
        # TODO: sort out numeric times... something like TAI in Ticktock
        tstamp = utc.timestamp()

        # Compute the various IGRF dependent things like position of CD
        # TODO: So that full IGRF can be used eventually, for each in
        # __coeffs, interp to time and save as '[g|h][n][m]' e.g., 'h32' (entries in self.coeffs)
        if not self.__valid_extrap:
            # Interpolate coefficients linearly
            numepochs = [dt.datetime(int(val), 1, 1).timestamp() for val in self._epochs]
            g10 = np.interp(tstamp, numepochs, self.__coeffs['g'][1][0])
            g11 = np.interp(tstamp, numepochs, self.__coeffs['g'][1][1])
            h11 = np.interp(tstamp, numepochs, self.__coeffs['h'][1][1])
        else:
            # Extrapolate using secular variation
            ref_ep = dt.datetime(int(self._epochs[-1]), 1, 1).timestamp()
            diff_years = (tstamp - ref_ep)/(86400*365.25)  # in years
            g10 = self.__coeffs['g_SV'][1][0]*diff_years + self.__coeffs['g'][1][0][-1]
            g11 = self.__coeffs['g_SV'][1][1]*diff_years + self.__coeffs['g'][1][1][-1]
            h11 = self.__coeffs['h_SV'][1][1]*diff_years + self.__coeffs['h'][1][1][-1]

        mom_sq = g10*g10 + g11*g11 + h11*h11
        mom = np.sqrt(mom_sq)

        # Compute dipole moments.
        self.moment = {'cd': mom,
                       'cd_McIlwain': 31165.3,
                       'cd_2010': 29950.126496041714,
                       }

        self.dipole = dict()
        self.dipole['cd_gcolat_rad'] = np.pi - np.arccos(g10/mom)
        self.dipole['cd_glon_rad'] = np.arctan(h11/g11)
        self.dipole['cd_gcolat'] = np.rad2deg(self.dipole['cd_gcolat_rad'])
        self.dipole['cd_glon'] = np.rad2deg(self.dipole['cd_glon_rad'])
