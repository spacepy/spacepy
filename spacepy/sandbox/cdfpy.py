#!/usr/bin/python
# -*- coding: utf-8 -*-

"""CDFPy -- Interface to NASA CDF in Python.


--++-- By Steve Morley and Brian Larsen--++--

smorley@lanl.gov/morley_steve@hotmail.com,
balarsen@lanl.gov
Los Alamos National Laboratory, ISR-1, MS D466
PO Box 1663, Los Alamos, NM 87545
"""

import cdfpy_interface as pcdfi #imports swig-wrapped shared library

__license__ = """Python Interface to NASA CDF, Copyright (C) 2010  Steven Morley

This interface to the NASA CDF library is released under GPL v3.0:
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__licence__ = __license__ #for those who speak English, rather than an odd dialect

class CDF(object):
    """CDF file class
    
    Not nearly complete...
    To Do: import (some) more functionality, including actually reading in values!

    """
    def __init__(self, verbose=True):
        self.data = {'GlobalAttr': {}, 'LocalAttr': {}, 'zVars': {}, 'rVars': {}}
        self.verbose = verbose
        self.returncodes = {
            'Encoding': ['Network', 'Sun', 'VAX', 'Dec', 'SGi', 'IBMPC',
                'IBMRS', 'Host', 'PPC', 'HP', 'NeXT', 'AlphaOSF1', 'AlphaVMSd'
                , 'AlphaVMSg', 'AlphaVMSi'],
            'Majority': ['Row', 'Column'],
            'Scope': ['Global', 'Variable']}
    
    def __str__(self):
        """Define String Representation of Sea object"""
        try:
            return """CDF Object - Filename:  %s 
            For file contents use the index method""" % (self.fname)
        except:
            return """Object not in use"""
    
    __repr__ = __str__
    
    def open(self, fname=None, verbose=False, load=True):
        """Opens and reads a CDF file into an object
        
        Optional arguments:
        fname - input filename. If not supplied will read from fname attribute (if set).
        load - (Default: True) Load all attributes and variables on open
        """
        if not fname:
            try:
                fname = self.fname
            except AttributeError:
                return 'Specify a filename'
        else:
            self.fname = fname
        try:
            try:
                assert self._id
                return 'File already open. Create a new CDF object...'
            except:
                pass
            n_id = pcdfi.new_CDFidp()
            a = pcdfi.PyCDFopen(fname, n_id)
            _id = pcdfi.CDFidp_value(n_id)
            self._id = _id
        except:
            return 'File open failed'
        
        if a==0: #to be replaced by CDFstatus handler
            print "File opened successfully. "
        else:
            return 'File open failed for some reason...'
        return None

    def close(self):
        """Closes a previously opened CDF file"""
        try:
            dum = self._id
        except AttributeError:
            return 'File not open'
        pcdfi.PyCDFclose(self._id)
        self._id = None #del pointer ref.
        del self.fname
        return None

    
    def index(self):
        """Prints variables in CDF data store"""
        try:
            return spacepy.toolbox.dictree(self.data, verbose=True)
        except AttributeError:
            return 'No data'
    
        return None


    def listzVars(self):
        """List the z-variables in the CDF file """
        


#End class. Begin namespace accessible functions

