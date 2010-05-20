#!/usr/bin/python
# -*- coding: utf-8 -*-

"""CDFPy -- Interface to NASA CDF in Python.


--++-- By Steve Morley --++--

smorley@lanl.gov/morley_steve@hotmail.com,
Los Alamos National Laboratory, ISR-1,
PO Box 1663, Los Alamos, NM 87545
"""

from spacepy import help
import cdfpy_interface as pcdfi #imports swig-wrapped shared library
from cdfpy_interface import new_longp, new_charp, longp_value
import spacepy

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
    
    def open(self, fname=None, load=True):
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
            self._id = pcdfi.CDFidp_value(n_id)
        except:
            return 'File open failed'
        
        if a==0: #to be replaced by CDFstatus handler
            print "File opened successfully. Querying file..."
        else:
            return 'File open failed for some reason...'
        #stuff to read attrs and vars goes here
        #put data from vars into attribute data (a dictionary)
        numDims, dimSizes = new_longp(), new_longp()
        encoding, majority = new_longp(), new_longp()
        maxrRec, numrVars = new_longp(), new_longp()
        maxzRec, numzVars  = new_longp(), new_longp()
        numAttrs = new_longp()
        a = pcdfi.PyCDFinquire(self._id, numDims, dimSizes, encoding, majority, \
            maxrRec, numrVars, maxzRec, numzVars, numAttrs)
        if a!=0: #to be replaced by CDFstatus handler
            return 'File inquiry failed... CDFstatus code: %d' % a
        if self.verbose:
            print 'Majority: %s' % self.returncodes['Majority'][longp_value(majority)-1]
            print 'Number of dimensions: %d' % longp_value(numDims)
            print 'Number of Attributes: %d' % longp_value(numAttrs)
            print 'Number of zVariables: %d' % longp_value(numzVars)
            print 'Number of rVariables: %d' % longp_value(numrVars)
            print 'Maximum r-record: %d' % longp_value(maxrRec)

        #get attribute info
        for i in range(longp_value(numAttrs)):
            attrScope, maxEntry = new_longp(), new_longp()
            a, attrName = pcdfi.PyCDFinquireAttr(self._id, long(i), attrScope, maxEntry)
            print 'Max Entry: %d' % longp_value(maxEntry)
            if pcdfi.longp_value(attrScope) == 1:
                self.data['GlobalAttr'][attrName]=[]
            else:
                self.data['LocalAttr'][attrName]=[]

            if load: #load data
                #if i==0:
                    #attrEntry = new_charp()
                attrEntry = pcdfi.PyCDFgetAttr(self._id, long(i), long(1))
                print ''
                if pcdfi.longp_value(attrScope) == 1:
                    self.data['GlobalAttr'][attrName].append(pcdfi.charp_value(attrEntry))
                else:
                    self.data['LocalAttr'][attrName].append(pcdfi.charp_value(attrEntry))
                #else:
                    #pass

        #get z-variable info
        for i in range(longp_value(numzVars)):
            dataType, numElements = new_longp(), new_longp()
            numDims, dimSizes = new_longp(), new_longp()
            recVary, dimVarys = new_longp(), new_longp()
            a, varName = pcdfi.PyCDFinquirezVar(self._id, long(i), dataType, numElements, \
                numDims, dimSizes, recVary, dimVarys)
            self.data['zVars'][varName]=[]
            if load:
                pass
        #get r-variable info
        for i in range(longp_value(numrVars)):
            dataType, numElements = new_longp(), new_longp()
            recVary, dimVarys = new_longp(), new_longp()
            a, varName = pcdfi.PyCDFinquirerVar(self._id, long(i), dataType, numElements, \
                recVary, dimVarys)
            print 'r-variable name: %s - num elements: %d' % (varName, pcdfi.longp_value(numElements))
            print 'r-variable dimVarys: %d - recVary: %d' % (pcdfi.longp_value(dimVarys), pcdfi.longp_value(recVary))
            self.data['rVars'][varName]=[]
            if load:
                pass

        print 'Data Structure'
        spacepy.toolbox.dictree(self.data)

        pcdfi.delete_CDFidp(n_id) #necessary?
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


#End class. Begin namespace accessible functions
def computeEpoch(year, month, day, hour=0, minute=0, second=0, msec=0):
    """Find CDF epoch to nearest millisecond"""
    #set up inputs (must be long-int)
    if year % int(year) != 0 or month % int(month) != 0 or day % int(day) != 0:
        return 'Year, month and day must be integer'
    if hour != 0 or minute != 0:
        if hour % int(hour) != 0 or minute % int(minute) != 0:
            return 'Hour and minute must be integer'
    if second != 0 and second % int(second) != 0:
        msec = round((second % int(second))*1000)
    year, month, day, hour, minute, second, msec = \
        int(year), int(month), int(day), int(hour), int(minute), int(second), int(msec)
    CDFepoch = pcdfi.computeEPOCH(year, month, day, hour, minute, second, msec)
    return CDFepoch

def epochBreakdown(CDFepoch):
    #use spacetime module rather than call CDF library
    import spacepy.spacetime as spt
    dum = TickTock(CDFepoch, 'CDF')
    dtobj = dum.UTC
    #dtobj = spacepy.getUTC(CDFepoch, 'CDF') #replace this line when ticktock is renamed/fixed
    
    return dtobj