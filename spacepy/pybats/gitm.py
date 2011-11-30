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

    GITM index ordering is [lon, lat, altitude]; data arrays read from file
    will always be of the same shape and size.
    '''

    def __init__(self, filename, *args, **kwargs):
        super(GitmBin, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        self._read()

    def __repr__(self):
        return 'GITM binary output file %s' % (self.attrs['file'])

    def _read(self):
        '''
        Read binary file; should only be called upon instantiation.
        '''

        from re import sub
        from struct import unpack
        
        # Read data and header info
        f=open(self.attrs['file'], 'rb')

        # Using the first FORTRAN header, determine endian.
        # Default is little.
        self.attrs['endian']='little'
        endChar='>'
        rawRecLen=f.read(4)
        recLen=(unpack(endChar+'l',rawRecLen))[0]
        if (recLen>10000)or(recLen<0):
            # Ridiculous record length implies wrong endian.
            self.attrs['endian']='big'
            endChar='<'
            recLen=(unpack(endChar+'l',rawRecLen))[0]

        # Read version; read fortran footer+header.
        self.attrs['version']=unpack(endChar+'d',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read grid size information.
        (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt'])=\
            unpack(endChar+'lll',f.read(recLen))
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Read number of variables.
        self.attrs['nVars']=unpack(endChar+'l',f.read(recLen))[0]
        (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Collect variable names.
        var=[]
        for i in range(self.attrs['nVars']):
            var.append(unpack(endChar+'%is'%(recLen),f.read(recLen))[0])
            (oldLen, recLen)=unpack(endChar+'2l',f.read(8))

        # Extract time.
        (yy,mm,dd,hh,mn,ss,ms)=unpack(endChar+'lllllll',f.read(recLen))
        self['time']=dt.datetime(yy,mm,dd,hh,mn,ss,ms/1000)
        (oldLen)=unpack(endChar+'l',f.read(4))


        # Read the rest of the data.
        nTotal=self.attrs['nLon']*self.attrs['nLat']*self.attrs['nAlt']
        for val in var:
            # Trim variable names.
            v=sub('\[|\]', '', val).strip()
            s=unpack(endChar+'l',f.read(4))[0]
            self[v]=dmarray(np.array(unpack(
                        endChar+'%id'%(nTotal),f.read(s))))
            # Reshape arrays, note that ordering in file is Fortran-like.
            self[v]=self[v].reshape( 
                (self.attrs['nLon'],self.attrs['nLat'],self.attrs['nAlt']),
                order='fortran')
            f.read(4)


#    def add_alt_slice(alt, var, target=None):
#        '''
#        Add a contour plot of variable var against lat/lon for a given
#        constant altitude slice, alt.
#        '''
