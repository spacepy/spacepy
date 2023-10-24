#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the Global 
Ionosphere-Thermosphere Model (GITM), one of the choices for the UA module
in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
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

    def __init__(self, filename, *args, **kwargs):
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        self._read()
        self.calc_deg()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _read(self):
        '''
        Read binary file; should only be called upon instantiation.
        '''

        from re import sub
        from spacepy.pybats import readarray
        
        # Open, read, and parse the file into numpy arrays.
        # Note that Fortran writes integer buffers around records, so
        # we must parse those as well.
        with open(self.attrs['file'], 'rb') as infile:
            # On the first try, we may fail because of wrong-endianess.
            # If that is the case, swap that endian and try again.
            endian='little'

            inttype=np.dtype(np.int32)
            EndChar='<'
            inttype.newbyteorder(EndChar)

            # detect double-precision file by checking record length
            # of version (stored as float). Simultaneously, test for
            # byte ordering:
            try:
                RecLen=np.fromfile(infile,dtype=inttype,count=1)
            except (ValueError,EOFError):
                endian='big'
                EndChar='>'
                inttype.newbyteorder(EndChar)
                infile.seek(0)
                RecLen=np.fromfile(infile,dtype=inttype,count=1)
            infile.seek(0) # return to beginning of file.
            self.attrs['endian']  = endian # Store endian order.
            
            # Set data types based on record length from above:
            if RecLen > 4:
                floattype=np.dtype(np.float64)
            else:
                floattype=np.dtype(np.float32)
            floattype.newbyteorder(EndChar)

            # Get code version:
            self.attrs['version']=readarray(infile,floattype,inttype)[0]

            # Read grid size information.  Because these three values
            # are written at the same time, they get "packed" together and
            # must be read together.
            header_fields_dtype=np.dtype([
                ('nLon',np.int32), ('nLat',np.int32), ('nAlt',np.int32)])
            header_fields_dtype.newbyteorder(EndChar)
            self.attrs['nLon'], self.attrs['nLat'], self.attrs['nAlt'] = \
                readarray(infile,header_fields_dtype,inttype)[0]
                
            # Read number of variables:
            self.attrs['nVars'] = readarray(infile,inttype,inttype)[0]

            # Collect variable names; decode into usable strings.
            var=[]
            for i in range(self.attrs['nVars']):
                var.append(readarray(infile,str,inttype).decode('utf-8'))

            # Extract time, stored as 7 consecutive integers.
            (yy,mm,dd,hh,mn,ss,ms)=readarray(infile,inttype,inttype)
            self['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms//1000)

            # Read the rest of the data.
            nTotal=self.attrs['nLon']*self.attrs['nLat']*self.attrs['nAlt']
            for val in var:
                # Trim variable names.
                v=sub(r'\[|\]', '', val).strip()
                self[v] = readarray(infile,floattype,inttype)
                # Reshape arrays, note that ordering in file is Fortran-like.
                self[v]=self[v].reshape( 
                    (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt']),
                    order='F')
                # Reduce dimensionality:
                self[v] = np.squeeze(self[v])

    def calc_deg(self):
        '''
        Gitm defaults to radians for lat and lon, which is sometimes difficult
        to use.  This method creates *dLat* and *dLon*, which is lat and lon
        in degrees.
        '''
        from numpy import pi
        if 'Latitude' in self:
            self['dLat'] = dmarray(self['Latitude']*180.0/pi, 
                                   attrs={'units':'degrees'})
        if 'Longitude' in self:
            self['dLon'] = dmarray(self['Longitude']*180.0/pi, 
                                   attrs={'units':'degrees'})
