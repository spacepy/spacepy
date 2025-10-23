#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Global 
Ionosphere-Thermosphere Model (GITM), one of the choices for the UA module
in the SWMF.
'''

# Global imports:
import datetime
import glob
import numpy as np
import struct

from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

class GitmBin(PbData):
    '''
    Object to open, manipulate and visualize 1 and 3 dimensional GITM output
    stored in binary format.  Object inherits from spacepy.pybats.PbData; see
    that documentation for general information on how these objects work.

    GITM index ordering is [lon, lat, altitude].  For 2D or 1D files, any
    dimension that has only one value will have that dimension removed from the
    final arrays (e.g., a 2D cut through a single altitude will have the 
    altitude dimension removed; the final arrays will be lon X lat only.)
    '''

    def __init__(self, filenames, varlist=None, degrees=False, *args, **kwargs):
        """
        Initializes the holder for GITM outputs & reads the data.

        Parameters
        ==========
            filenames (str/list-like): Path(s) of GITM files to read. Works for one or
                multiple. Can also be a globbed path, like
                '/path/to/run/data/3DALL*.bin'. Cannot mix output types!
            varlist (list-like): Indices of the variables to read. Must be int's.
                Default=None, and all variables are read.

        Examples
        ========
            (to read all variables):
            
            > data = gitm.GitmBin('path/to/file/3DALL_t123_000.bin')

            (to read select variables):
            
            > data = gitm.GitmBin('path/to/file/3DALL_*.bin', varlist=[0,1,2,34])


        """
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.

        # Can accept single/multiple files. sanitize this:
        if isinstance(filenames, str):
            filenames = sorted(glob.glob(filenames))

        self.attrs['files'] = filenames

        # If varlist is just an int, we can handle it
        if isinstance(varlist, int):
            varlist = [varlist]

        # Read first header. Puts attrs into self
        self._read_header(varlist=varlist)

        # read in all the data
        self._readbin(varlist=varlist)
        if degrees:
            self.calc_deg()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _readbin(self, varlist: list[int]=None):
        '''
        Read binary file; should only be called upon instantiation.
        '''

        # pre-allocate the data arrays:
        # Shape is (nTimes, nLons, nLats, nAlts), these are all present in the header
        for var in self.attrs['var_idxs'].keys():
            self[var] = dmarray(np.empty([len(self.attrs['files']),
                                          self.attrs['nLon'],
                                          self.attrs['nLat'], 
                                          self.attrs['nAlt']]), #attrs
                                          )
            
        
        # skip the header: grid, nvars, version, variable names, time
        # all are 8-bytes followed by 4-byte startstop sentinel
        HeaderLength = 52 + 48 * self.attrs["nVarsTotal"] + 32
        nPtsTotal = self.attrs["nLon"]*self.attrs["nLat"]*self.attrs["nAlt"]
        iDataLength = nPtsTotal*8 + 4 + 4

        # Open, read, and parse the file into numpy arrays.
        # Note that Fortran writes integer buffers around records, so
        # we must parse those as well.
        isFirstTime = True
        for iFile, filename in enumerate(self.attrs['files']):
            with open(filename, 'rb') as file:
                # On the first try, we may fail because of wrong-endianess.
                # If that is the case, swap that endian and try again.
                endChar = '>'
                rawRecLen = file.read(4)
                rec_len = (struct.unpack(endChar + 'l', rawRecLen))[0]
                if rec_len > 10000 or rec_len < 0:
                    # Ridiculous record length implies wrong endian.
                    endChar='<'

                # Only read time if _read_header has not (it reads the 1st time)
                if not isFirstTime:
                    # Extract time, stored as 7 consecutive integers.
                    # time is placed after variable names, so we skip:
                    # 52(grid etc.) + ( 40[nvars] + 8 [head/foot] ) *nVarsTotal
                    file.seek(52 + 48 * self.attrs['nVarsTotal'])
                    (yy,mm,dd,hh,mn,ss,ms,startstop) = struct.unpack(endChar+'llllllli',file.read(32))
                    self.attrs['time'][iFile]=datetime.datetime(yy,mm,dd,hh,mn,ss,ms*1000)
                isFirstTime = False

                for varname, iVar in self.attrs['var_idxs'].items():
                    # Get to the right location in file
                    file.seek(HeaderLength+iVar*iDataLength)
                    # Pull out the data
                    s=struct.unpack(endChar+'l', file.read(4))[0]
                    self[varname][iFile, ...] = \
                        np.array(struct.unpack(endChar+'%id'%(nPtsTotal),file.read(s))
                                 ).reshape((self.attrs['nLon'],
                                            self.attrs['nLat'],
                                            self.attrs['nAlt']), order='F')

        # Finally we reduce dimensionality:
        for v in self.attrs['var_idxs']:
            self[v] = np.squeeze(self[v])

    
    def _read_header(self, varlist=None):
        """
        This opens up a single GITM output ans extracts some data from the header.
        Will, by default, only work on the first file in self.attrs['files']
        
        This will let us pre-allocate the arrays to hold the data

        Parameters
        ==========
            varlist (list): list of indices of variables to read. Default=None, 
                which reads all variables
        
        """
    
        with open(self.attrs['files'][0], 'rb') as file:
            # determine endianess 
            endChar='>'
            rawRecLen=file.read(4)
            recLen=(struct.unpack(endChar+'l',rawRecLen))[0]
            if (recLen>10000)or(recLen<0):
                # Ridiculous record length implies wrong endian.
                endChar='<'
                recLen=(struct.unpack(endChar+'l',rawRecLen))[0]

            # Read version; read fortran footer+header.
            self.attrs["version"] = struct.unpack(endChar+'d',file.read(recLen))[0]

            (_, recLen)=struct.unpack(endChar+'2l',file.read(8))

            # Read grid size information.
            (self.attrs["nLon"],self.attrs["nLat"],self.attrs["nAlt"]) = \
                struct.unpack(endChar+'lll',file.read(recLen))
            (_, recLen)=struct.unpack(endChar+'2l',file.read(8))

            # Read number of variables.
            self.attrs["nVarsTotal"]=struct.unpack(endChar+'l',file.read(recLen))[0]
            (_, recLen)=struct.unpack(endChar+'2l',file.read(8))

            # Collect variable names & indices:
            # This is going into a dict of {varname: index}
            self.attrs['var_idxs'] = {}

            for i in range(self.attrs["nVarsTotal"]):
                v = struct.unpack(endChar+'%is'%(recLen),file.read(recLen))[0]

                nVarsRead = 0
                # Only save it if it is in varlist, or varlist is None (read all vars)
                if varlist is None or i in [0, 1, 2] or i in varlist:
                    # All 3 coords (lon, lat, alt) are present in all outputs
                    self.attrs['var_idxs'][v.decode('utf-8').replace(" ","")] = i
                    nVarsRead += 1
                (oldLen, recLen)=struct.unpack(endChar+'2l',file.read(8))
                self.attrs['nVars'] = nVarsRead

            if varlist is not None and max(varlist) > self.attrs['nVarsTotal']:
                raise IndexError(
                    "Variable out of range!\n"
                    f" -> Only {self.attrs['nVarsTotal']} variables are available\n"
                    f" --> You provided: varlist={varlist}")

            # Extract time. 
            (yy,mm,dd,hh,mn,ss,ms) = struct.unpack(endChar+'lllllll',file.read(recLen))
            self.attrs["time"] = dmarray(np.zeros(len(self.attrs['files']), 
                                                  dtype='object'),
                                                  dtype='object') # is this necessary?
            self.attrs['time'][0] = datetime.datetime(yy,mm,dd,hh,mn,ss,ms*1000)

        return 

    def calc_deg(self):
        '''
        Gitm defaults to radians for lat and lon, and m for altitude. This is
        sometimes difficult to use. This method creates *dLat* and *dLon*,
        which is lat and lon in degrees, and converts *Altitude* to km.
        '''
        if 'Latitude' in self:
            self['dLat'] = dmarray(np.rad2deg(self['Latitude']), 
                                   attrs={'units':'degrees'})
        if 'Longitude' in self:
            self['dLon'] = dmarray(np.rad2deg(self['Longitude']), 
                                   attrs={'units':'degrees'})
        if 'Altitude' in self:
            self['Altitude'] = dmarray(self['Altitude']*1000, 
                                       attrs={'units':'kilometers'})
            
