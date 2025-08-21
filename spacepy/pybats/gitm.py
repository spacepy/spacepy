#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Global 
Ionosphere-Thermosphere Model (GITM), one of the choices for the UA module
in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
from struct import unpack
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray
import glob

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

        Args:
            filenames (str/list-like): Path(s) of GITM files to read. Works for one or
                multiple. Can also be a globbed path, like
                '/path/to/run/data/3DALL*.bin'. Cannot mix output types!
            varlist (list-like): Indices of the variables to read. Must be int's.
                Default=None, and all variables are read.

        Examples:
            (to read all variables):
            
            > data = gitm.GitmBin('path/to/file/3DALL_t123_000.bin')

            (to read select variables):
            
            > data = gitm.GitmBin('path/to/file/3DALL_*.bin', varlist=[0,1,2,34])


        """
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.
        
        #TODO: As written, this cannot read multiple output types in one call.
        # -> Either fix this (desc. how at '^[outputtypes]'), or check inputs
        # - Inputs can be checked by ensuring f.split('/')[-1][:5] is same for f in filelist

        # Can accept single/multiple files. sanitize this:
        if isinstance(filenames, str):
            if '*' in filenames: # We have to glob the path
                filenames = sorted(glob.glob(filenames))
            else:
                filenames = [filenames]
        self.attrs['files'] = filenames

        # TODO: add method/escape to just print the variable list (or all attrs).
        # Would need a call to _get_header() and then a print of self.attrs['vars']
        # -> ideally we do not read in any bin data to do this!

        # Read first header. Puts attrs into self
        self._read_header(varlist=varlist)

        # read in all the data
        self._readbin(varlist=varlist)
        if degrees:
            self.calc_deg()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _readbin(self, varlist: list=None):
        '''
        Read binary file; should only be called upon instantiation.
        '''

        # pre-allocate the data arrays:
        # Shape is (nTimes, nLons, nLats, nAlts), these are all present in the header
        for var in self.attrs['var_idxs'].keys():
            #TODO: Add a call to the variable parsing function described in _read_header()
            self[var] = dmarray(np.empty([len(self.attrs['files']),
                                          self.attrs['nLons'],
                                          self.attrs['nLats'], 
                                          self.attrs['nAlts']]), #attrs
                                          )
            
        
        # skip the version, recLen, shape, dimensions, etc. + start/stop byte
        HeaderLength = 84 + self.attrs["nVarsTotal"] * 48
        nPtsTotal = self.attrs["nLons"]*self.attrs["nLats"]*self.attrs["nAlts"]
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
                rec_len = (unpack(endChar + 'l', rawRecLen))[0]
                if rec_len > 10000 or rec_len < 0:
                    # Ridiculous record length implies wrong endian.
                    endChar='<'

                # Only read time if _read_header has not (it reads the 1st time)
                if not isFirstTime:
                    # Extract time, stored as 7 consecutive integers.
                    # time is placed after variable names, so we skip:
                    # 64(header) + ( 40[nvars] + 8 [head/foot] ) *nVarsTotal
                    file.seek(64 + 48 * self.attrs['nVarsTotal'])
                    (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',file.read(28))
                    self.attrs['time'][iFile]=dt.datetime(yy,mm,dd,hh,mn,ss,ms//1000)
                    isFirstTime = False

                for varname, iVar in self.attrs['var_idxs'].items():
                    # Get to the right location in file
                    file.seek(HeaderLength+iVar*iDataLength)
                    # Pull out the data
                    s=unpack(endChar+'l', file.read(4))[0]
                    self[varname][iFile, ...] = \
                        np.array(unpack(endChar+'%id'%(nPtsTotal),file.read(s))
                                 ).reshape((self.attrs['nLons'],
                                            self.attrs['nLats'],
                                            self.attrs['nAlts']), order='F')

        # Finally we reduce dimensionality:
        for v in self.attrs['var_idxs'].keys():
            self[v] = np.squeeze(self[v])

    
    def _read_header(self, varlist=None):
        """
        This opens up a single GITM output ans extracts some data from the header.
        Will, by default, only work on the first file in self.attrs['files']
        
        This will let us pre-allocate the arrays to hold the data

        Inputs
        ------
            varlist (list): list of indices of variables to read. Default=None, 
                which reads all variables
        
        """
    
        with open(self.attrs['files'][0], 'rb') as file:
            # determine endianess 
            endChar='>'
            rawRecLen=file.read(4)
            recLen=(unpack(endChar+'l',rawRecLen))[0]
            if (recLen>10000)or(recLen<0):
                # Ridiculous record length implies wrong endian.
                endChar='<'
                recLen=(unpack(endChar+'l',rawRecLen))[0]

            # Read version; read fortran footer+header.
            self.attrs["version"] = unpack(endChar+'d',file.read(recLen))[0]

            (_, recLen)=unpack(endChar+'2l',file.read(8))

            # Read grid size information.
            (self.attrs["nLons"],self.attrs["nLats"],self.attrs["nAlts"]) = \
                unpack(endChar+'lll',file.read(recLen))
            (_, recLen)=unpack(endChar+'2l',file.read(8))

            # Read number of variables.
            self.attrs["nVarsTotal"]=unpack(endChar+'l',file.read(recLen))[0]
            (_, recLen)=unpack(endChar+'2l',file.read(8))

            # Collect variable names & indices:
            # This is going into a dict of {varname: index}
            self.attrs['var_idxs'] = {}

            for i in range(self.attrs["nVarsTotal"]):
                v = unpack(endChar+'%is'%(recLen),file.read(recLen))[0]
                # TODO: Here we can add a call to a lookup table with variable info
                # - units!
                # - pretty name for plots
                # - user-friendly name
                # - etc.

                nVarsRead = 0
                # Only save it if it is in varlist!
                if varlist is None or i in [0, 1, 2]:
                    # NeedsReview, this reads coord info even if user didn't request it.
                    # All 3 coords (lon, lat, alt) are present in all outputs
                    self.attrs['var_idxs'][v.decode('utf-8').replace(" ","")] = i
                    nVarsRead += 1
                (oldLen, recLen)=unpack(endChar+'2l',file.read(8))
                self.attrs['nVars'] = nVarsRead

            # Extract time. 
            (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',file.read(recLen))
            self.attrs["time"] = dmarray(np.zeros(len(self.attrs['files']), 
                                                  dtype='object'),
                                                  dtype='object') # is this necessary?
            self.attrs['time'][0] = dt.datetime(yy,mm,dd,hh,mn,ss,ms*1000)

        return 

    def calc_deg(self):
        '''
        Gitm defaults to radians for lat and lon, and m for altitude. This is
        sometimes difficult to use.  This method creates *dLat* and *dLon*,
        which is lat and lon in degrees.
        '''
        from numpy import pi
        if 'Latitude' in self:
            self['dLat'] = dmarray(self['Latitude']*180.0/pi, 
                                   attrs={'units':'degrees'})
        if 'Longitude' in self:
            self['dLon'] = dmarray(self['Longitude']*180.0/pi, 
                                   attrs={'units':'degrees'})
        if 'Altitude' in self:
            self['dLon'] = dmarray(self['Altitude']*1000, 
                                   attrs={'units':'kilometers'})
            
