import cdfpy_interface as pcdfi
import spacepy
import numpy as np

class cdfpy(object):
    """
    I need docs
    """
    # a few constants needed for the code
    

    def epoch2Ticktock(self, epoch):
        """
        I need docs
        I change the epoch to Ticktock
        """
        from spacepy.time import Ticktock
        from datetime import datetime
        year = pcdfi.new_longArray(1)
        month = pcdfi.new_longArray(1)
        day = pcdfi.new_longArray(1)
        hour = pcdfi.new_longArray(1)
        minute = pcdfi.new_longArray(1)
        second = pcdfi.new_longArray(1)
        msec = pcdfi.new_longArray(1)
        dt = np.array([])
        epoch = self._data['zVars']['Epoch']['data']
        for val in epoch:
            pcdfi.EPOCHbreakdown(val, year, month, day, hour, minute, second, msec)
            dt = np.append(dt, datetime(pcdfi.longArray_getitem(year, 0),
                                        pcdfi.longArray_getitem(month, 0),
                                        pcdfi.longArray_getitem(day, 0),
                                        pcdfi.longArray_getitem(hour, 0), 
                                        pcdfi.longArray_getitem(minute, 0),
                                        pcdfi.longArray_getitem(second, 0),
                                        pcdfi.longArray_getitem(msec, 0)*1000))
        # free these pointers!!!! 
        return Ticktock(dt, 'UTC')

    def __init__(self, filename = None, verbose = True):
        _data = {'GlobalAttr': {}, 'zVars': {}, 'rVars': {}}
        self._data = _data

        self.__CDF_MAX_DIMS = 10

        self.__datatypes = { '1':'CDF_INT1',
                             '2':'CDF_INT2',
                             '4':'CDF_INT4',
                             '11':'CDF_UINT1',
                             '12':'CDF_UINT2',
                             '14':'CDF_UINT4',
                             '21':'CDF_REAL4',
                             '22':'CDF_REAL8',
                             '31':'CDF_EPOCH',
                             '32':'CDF_EPOCH16',
                             '41':'CDF_BYTE',
                             '44':'CDF_FLOAT',
                             '45':'CDF_DOUBLE',
                             '51':'CDF_CHAR',
                             '52':'CDF_UCHAR' }

        self.__datatypesReverse = dict((v,k) for k, v in self.__datatypes.iteritems())
        
        self.__returncodes = {
            'Encoding': ['Network', 'Sun', 'VAX', 'Dec', 'SGi', 'IBMPC',
                         'IBMRS', 'Host', 'PPC', 'HP', 'NeXT', 'AlphaOSF1', 'AlphaVMSd'
                         , 'AlphaVMSg', 'AlphaVMSi'],
            'Majority': ['Row', 'Column'],
            'Scope': ['Global', 'Variable']}

        self.__formatCodes = {
            'CDF_INT2':'h',
            'CDF_INT4':'i',
            'CDF_UINT1':'B',
            'CDF_UINT2':'H',
            'CDF_UINT4':'L',
            'CDF_REAL4':'f',
            'CDF_REAL8':'d',
            'CDF_EPOCH':'d',
            'CDF_EPOCH16':None,  # need to figure this out
            'CDF_BYTE':'b',
            'CDF_FLOAT':'f',
            'CDF_DOUBLE':'d',
            'CDF_CHAR':'c',
            'CDF_UCHAR':'B' }


        self.__dtypeSize = {
            '1':1,
            '2':2,
            '4':4,
            '11':1,
            '12':2,
            '14':4,
            '21':4,
            '22':8,
            '31':8,
            '32':16,
            '41':1,
            '44':4,
            '45':8,
            '51':1,
            '52':1,
            'CDF_INT2':2,
            'CDF_INT4':4,
            'CDF_UINT1':1,
            'CDF_UINT2':2,
            'CDF_UINT4':4,
            'CDF_REAL4':4,
            'CDF_REAL8':8,
            'CDF_EPOCH':8,
            'CDF_EPOCH16':16,
            'CDF_BYTE':1,
            'CDF_FLOAT':4,
            'CDF_DOUBLE':8,
            'CDF_CHAR':1,
            'CDF_UCHAR':1 }


       
        

        self.filename = filename
        

    def __repr__(self):
        st = 'CDF: ' + self.filename
        return st

    def _errorHandle(self, code, verbose=False):
        """
        I need docs
        I need code
        """
        # Error codes < CDF WARN < Warning codes < CDF OK < Informational codes
        status_codes = {
            # informational codes
            1001: 'VIRTUAL_RECORD_DATA',
            1002: 'DID_NOT_COMPRESS',
            1003: 'VAR_ALREADY_CLOSED',
            1004: 'SINGLE_FILE_FORMAT',
            1005: 'NO_PADVALUE_SPECIFIED',   
            1006: 'NO_VARS_IN_CDF',      
            1007: 'MULTI_FILE_FORMAT',
            1008: 'SOME_ALREADY_ALLOCATED',
            1009: 'PRECEEDING_RECORDS_ALLOCATED',
            0: 'CDF_OK',                       
            # warning codes
            -1001: 'ATTR_NAME_TRUNC',         
            -1002: 'CDF_NAME_TRUNC',          
            -1003: 'VAR_NAME_TRUNC',          
            -1004: 'NEGATIVE_FP_ZERO',
            -1006: 'FORCED_PARAMETER',
            -1007: 'NA_FOR_VARIABLE',
            -2000: 'CDF_WARN',
            # Error codes
            -2001: 'ATTR_EXISTS',
            -2002: 'BAD_CDF_ID', 
            -2003: 'BAD_DATA_TYPE',
            -2004: 'BAD_DIM_SIZE', 
            -2005: 'BAD_DIM_INDEX',
            -2006: 'BAD_ENCODING', 
            -2007: 'BAD_MAJORITY', 
            -2008: 'BAD_NUM_DIMS', 
            -2009: 'BAD_REC_NUM', 
            -2010: 'BAD_SCOPE',  
            -2011: 'BAD_NUM_ELEMS',
            -2012: 'CDF_OPEN_ERROR',
            -2013: 'CDF_EXISTS',  
            -2014: 'BAD_FORMAT',  
            -2015: 'BAD_ALLOCATE_RECS',
            -2016: 'BAD_CDF_EXTENSION',
            -2017: 'NO_SUCH_ATTR', 
            -2018: 'NO_SUCH_ENTRY',
            -2019: 'NO_SUCH_VAR', 
            -2020: 'VAR_READ_ERROR',
            -2021: 'VAR_WRITE_ERROR',
            -2022: 'BAD_ARGUMENT',  
            -2023: 'IBM_PC_OVERFLOW',
            -2024: 'TOO_MANY_VARS', 
            -2025: 'VAR_EXISTS',  
            -2026: 'BAD_MALLOC', 
            -2027: 'NOT_A_CDF', 
            -2028: 'CORRUPTED_V2_CDF',
            -2029: 'VAR_OPEN_ERROR', 
            -2030: 'BAD_INITIAL_RECS',
            -2031: 'BAD_BLOCKING_FACTOR',
            -2032: 'END_OF_VAR',
            -2034: 'BAD_CDFSTATUS', 
            -2035: 'CDF_INTERNAL_ERROR',
            -2036: 'BAD_NUM_VARS',
            -2037: 'BAD_REC_COUNT',
            -2038: 'BAD_REC_INTERVAL',
            -2039: 'BAD_DIM_COUNT', 
            -2040: 'BAD_DIM_INTERVAL',
            -2041: 'BAD_VAR_NUM',
            -2042: 'BAD_ATTR_NUM',
            -2043: 'BAD_ENTRY_NUM',
            -2044: 'BAD_ATTR_NAME', 
            -2045: 'BAD_VAR_NAME',
            -2046: 'NO_ATTR_SELECTED',
            -2047: 'NO_ENTRY_SELECTED', 
            -2048: 'NO_VAR_SELECTED', 
            -2049: 'BAD_CDF_NAME',
            -2051: 'CANNOT_CHANGE',
            -2052: 'NO_STATUS_SELECTED',
            -2053: 'NO_CDF_SELECTED',
            -2054: 'READ_ONLY_DISTRIBUTION',
            -2055: 'CDF_CLOSE_ERROR',
            -2056: 'VAR_CLOSE_ERROR',
            -2058: 'BAD_FNC_OR_ITEM',
            -2060: 'ILLEGAL_ON_V1_CDF',
            -2063: 'BAD_CACHE_SIZE', 
            -2066: 'CDF_CREATE_ERROR', 
            -2067: 'NO_SUCH_CDF',
            -2068: 'VAR_CREATE_ERROR', 
            -2070: 'READ_ONLY_MODE',
            -2071: 'ILLEGAL_IN_zMODE', 
            -2072: 'BAD_zMODE',
            -2073: 'BAD_READONLY_MODE',
            -2074: 'CDF_READ_ERROR',
            -2075: 'CDF_WRITE_ERROR',
            -2076: 'ILLEGAL_FOR_SCOPE',
            -2077: 'NO_MORE_ACCESS', 
            -2079: 'BAD_DECODING',
            -2081: 'BAD_NEGtoPOSfp0_MODE',
            -2082: 'UNSUPPORTED_OPERATION',
            -2083: 'CDF_SAVE_ERROR',
            -2084: 'VAR_SAVE_ERROR',
            -2086: 'NO_WRITE_ACCESS',
            -2087: 'NO_DELETE_ACCESS',
            -2088: 'CDF_DELETE_ERROR',
            -2089: 'VAR_DELETE_ERROR',
            -2090: 'UNKNOWN_COMPRESSION',
            -2091: 'CANNOT_COMPRESS',
            -2092: 'DECOMPRESSION_ERROR',
            -2093: 'COMPRESSION_ERROR',
            -2096: 'EMPTY_COMPRESSED_CDF',
            -2097: 'BAD_COMPRESSION_PARM',
            -2098: 'UNKNOWN_SPARSENESS',
            -2099: 'CANNOT_SPARSERECORDS',
            -2100: 'CANNOT_SPARSEARRAYS',
            -2101: 'TOO_MANY_PARMS',
            -2102: 'NO_SUCH_RECORD',
            -2103: 'CANNOT_ALLOCATE_RECORDS',
            -2106: 'SCRATCH_DELETE_ERROR',
            -2107: 'SCRATCH_CREATE_ERROR',
            -2108: 'SCRATCH_READ_ERROR',
            -2109: 'SCRATCH_WRITE_ERROR',
            -2110: 'BAD_SPARSEARRAYS_PARM',
            -2111: 'BAD_SCRATCH_DIR',
            -2113: 'NOT_A_CDF_OR_NOT_SUPPORTED',
            -2123: 'CORRUPTED_V3_CDF',
            -2124: 'ILLEGAL_EPOCH_FIELD',
            -2125: 'BAD_CHECKSUM',
            -2126: 'CHECKSUM_ERROR',
            -2127: 'CHECKSUM_NOT_ALLOWED',
            -2128: 'ZLIB_DECOMPRESSION_ERROR' }
            
        if code == 0:
            if verbose:
                print(status_codes[code])
            return
        else:
            print("Error code %d: %s" % (code, status_codes[code]))
            return

    def _open(self, mode='r', verbose=False):
        """
        open a cdf file for reading

        @keyword mode: mode the open the file, r or w
        @type mode: string
        @keyword verbose: be verbose to the screen
        @type verbose: bool

        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)
        @version: V2: 26-Jul-2010 (BAL) added mode keyword and 'w' functionality

        """
        if self.filename == None:
            print("must specify a filename")
            return False
        try: self.n_id
        except AttributeError:
            n_id = pcdfi.new_CDFidp()
            if mode.lower() == 'r':
                a = pcdfi.PyCDFopen(self.filename, n_id)
            elif mode.lower() == 'w':
                a = pcdfi.PyCDFcreate(self.filename, n_id)
            else:
                print("did not understand CDF file mode (r/w)")
                return False
            self._errorHandle(a, verbose=verbose)
            _id = pcdfi.CDFidp_value(n_id)
            self._id = _id
            self.n_id = n_id
            return True
        print("CDF already opened, not reopening")
        return True

    def _readCDFInfo(self, verbose=False):
        """
        read info about the open CDF file
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        try: self.n_id
        except:
            print("CDF is not open")
            return False
        #stuff to read attrs and vars goes here
        #put data from vars into attribute data (a dictionary)
        numDims = pcdfi.new_longArray(1)
        dimSizes = pcdfi.new_longArray(1)
        encoding = pcdfi.new_longArray(1)
        majority = pcdfi.new_longArray(1)
        maxrRec = pcdfi.new_longArray(1)
        numrVars = pcdfi.new_longArray(1)
        maxzRec = pcdfi.new_longArray(1)
        numzVars = pcdfi.new_longArray(1)
        numAttrs = pcdfi.new_longArray(1)
        a = pcdfi.PyCDFinquireCDF(self._id, numDims, dimSizes, encoding, majority, \
                                  maxrRec, numrVars, maxzRec, numzVars, numAttrs)
        self._errorHandle(a)
        if verbose:
            print 'Majority: %s' % self.__returncodes['Majority'][pcdfi.longArray_getitem(majority, 0)-1]
            print 'Number of dimensions: %d' % pcdfi.longArray_getitem(numDims, 0)
            print 'Number of Attributes: %d' % pcdfi.longArray_getitem(numAttrs, 0)
            print 'Number of zVariables: %d' % pcdfi.longArray_getitem(numzVars, 0)
            print 'Number of rVariables: %d' % pcdfi.longArray_getitem(numrVars, 0)
            print 'Maximum r-record: %d' % pcdfi.longArray_getitem(maxrRec, 0)
        if pcdfi.longArray_getitem(numrVars, 0) != 0:
            raise(Exception("**** CDF has R-variables, they are not yet supported ****"))
        # lets add all these to the dict in 'GlobalAttr'
        self._data['GlobalAttr']['Majority'] = self.__returncodes['Majority'][pcdfi.longArray_getitem(majority, 0)-1]
        self._data['GlobalAttr']['numDims'] = pcdfi.longArray_getitem(numDims, 0)
        self._data['GlobalAttr']['numAttrs'] = pcdfi.longArray_getitem(numAttrs, 0)
        self._data['GlobalAttr']['numzVars'] = pcdfi.longArray_getitem(numzVars, 0)
        self._data['GlobalAttr']['numrVars'] = pcdfi.longArray_getitem(numrVars, 0)
        self._data['GlobalAttr']['maxrRec'] = pcdfi.longArray_getitem(maxrRec, 0)

        # free them
        pcdfi.delete_longArray(numDims)
        pcdfi.delete_longArray(dimSizes)
        pcdfi.delete_longArray(encoding)
        pcdfi.delete_longArray(majority)
        pcdfi.delete_longArray(maxrRec) 
        pcdfi.delete_longArray(numrVars)
        pcdfi.delete_longArray(maxzRec) 
        pcdfi.delete_longArray(numzVars)
        pcdfi.delete_longArray(numAttrs)
        
        return True 

            
    def _readGlobalAttributes(self, verbose=False):
        """
        read the global attributes from the open CDF file
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        try: self.n_id
        except:
            print("CDF is not open")
            return False

        try: self._data['GlobalAttr']['Majority']
        except (NameError, KeyError): self._readCDFInfo(verbose=verbose)

        #get global attribute info
        if verbose:
            print("Getting global attributes")
        attrScope = pcdfi.new_longArray(1)
        maxgEntry = pcdfi.new_longArray(1)
        maxrEntry = pcdfi.new_longArray(1)
        maxzEntry = pcdfi.new_longArray(1)
        dataType = pcdfi.new_longArray(1)
        numElements = pcdfi.new_longArray(1)
        attrName = ''

        for i in xrange(self._data['GlobalAttr']['numAttrs']): 
            a, attrName = pcdfi.PyCDFinquireAttr(self._id,
                                                 long(i),
                                                 attrScope,
                                                 maxgEntry,
                                                 maxrEntry,
                                                 maxzEntry )
            self._errorHandle(a)
            if verbose:
                print("\tFound attr: %d - %s it has scope %s" % (i, attrName,  self.__returncodes['Scope'][pcdfi.longArray_getitem(attrScope, 0)-1]))
                print("\t\tg:%d r:%d z:%d"  % (pcdfi.longArray_getitem(maxgEntry,0), pcdfi.longArray_getitem(maxrEntry,0), pcdfi.longArray_getitem(maxzEntry,0)))

            if verbose:
                print("\t\tFound %d entries" % (pcdfi.longArray_getitem(maxgEntry, 0)+1))

            if a ==0: ## this means that it gpt a global, if gets a local a = -2076
                self._data['GlobalAttr'][attrName] = {}
                self._data['GlobalAttr'][attrName]['value'] = []

            for ii in xrange(pcdfi.longArray_getitem(maxgEntry, 0)+1):
            # for ii in xrange(0+1):

                a = pcdfi.PyCDFinquireAttrgEntry (self._id,
                                                  long(i),
                                                  long(ii), 
                                                  dataType,
                                                  numElements )
                if verbose:
                    self._errorHandle(a)
            
                if a == 0:  ## this means that it gpt a global, if gets a local a = -2076
                    self._data['GlobalAttr'][attrName]['attrNum'] = i
                    self._data['GlobalAttr'][attrName]['dataType'] = self.__datatypes[str(pcdfi.longArray_getitem(dataType,0))]
                    self._data['GlobalAttr'][attrName]['numElements'] = pcdfi.longArray_getitem(numElements, 0)
                
                    value = ' '*  (pcdfi.longArray_getitem(numElements, 0) + 1) # need the +1 here 

    
                    a = pcdfi.PyCDFgetAttrgEntry(self._id,
                                                 long(i),
                                                 long(ii),
                                                 value)  # this only gets global vars
                    self._errorHandle(a)
                    self._data['GlobalAttr'][attrName]['value'].append(value)
                    self._data['GlobalAttr'][attrName]['attrNum'] = i
                    self._data['GlobalAttr'][attrName]['attrScope'] = pcdfi.longArray_getitem(attrScope, 0)

                    if verbose:
                        print("\t\t\t" + attrName + ' = ' +  self._data['GlobalAttr'][attrName]['value'][ii])
        pcdfi.delete_longArray(maxgEntry)
        pcdfi.delete_longArray(maxrEntry)
        pcdfi.delete_longArray(maxzEntry)
        pcdfi.delete_longArray(dataType) 
        pcdfi.delete_longArray(numElements)
        return True

            
## get and add compression info to the data['GlobalAttr']
        ## malloc error for some reason
## compressionType = pcdfi.longArray(1)m
## compressionParams = pcdfi.longArray(pcdfi.CDF_MAX_PARMS)
## compressionPercentage = pcdfi.longArray(1)
## a = pcdfi.PyCDFgetCompression(_id, compressionType, compressionParams, compressionPercentage)

    def _getzVarInfo(self, verbose=False):
        """
        read info about the zVars in the open CDF file
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        try: self.n_id
        except:
            print("CDF is not open")
            return False

        try: self._data['GlobalAttr']['numzVars']
        except  (NameError, KeyError): self._readGlobalAttributes(verbose=verbose)
        #get z-variable info
        if verbose:
            print("Getting zVariable info")
        zvarnum = {}
        # step through each zVar
        dimSizes = pcdfi.new_longArray(self.__CDF_MAX_DIMS)
        recVary  = pcdfi.new_longArray(1) 
        dimVarys = pcdfi.new_longArray(self.__CDF_MAX_DIMS)
        dataType = pcdfi.new_longArray(1)
        numElements = pcdfi.new_longArray(1)
        numDims = pcdfi.new_longArray(1)

        for i in xrange( self._data['GlobalAttr']['numzVars']):
            # get its name
            a, varName = pcdfi.PyCDFinquirezVar(self._id,
                                                long(i),
                                                dataType,
                                                numElements, 
                                                numDims,
                                                dimSizes,
                                                recVary,
                                                dimVarys)
            if verbose:
                print("Found zVar number %d: %s "  % (i, varName))
            # and create a dict with that name 
            try: self._data['zVars'][varName]
            except: self._data['zVars'][varName]={}
            self._data['zVars'][varName]['dataType'] = self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))]
            self._data['zVars'][varName]['numElements'] = pcdfi.longArray_getitem(numElements, 0)
            self._data['zVars'][varName]['numDims'] = pcdfi.longArray_getitem(numDims, 0)
            self._data['zVars'][varName]['dimSizes'] = np.array([], dtype=long)
            for ii in xrange(self._data['zVars'][varName]['numDims']):
                self._data['zVars'][varName]['dimSizes'] = np.append(self._data['zVars'][varName]['dimSizes'], pcdfi.longArray_getitem(dimSizes, ii))
            self._data['zVars'][varName]['dimVarys'] = np.array([], dtype=long)
            for ii in xrange(self._data['zVars'][varName]['numDims']):
                self._data['zVars'][varName]['dimVarys'] = np.append(self._data['zVars'][varName]['dimVarys'], pcdfi.longArray_getitem(dimVarys, ii ))
            self._data['zVars'][varName]['recVary'] = pcdfi.longArray_getitem(recVary, 0)
            self._data['zVars'][varName]['zVarnum'] = i
            self._data['zVars'][varName]['formatCode'] = self.__formatCodes[self._data['zVars'][varName]['dataType'] ]
    
            # how many records are there for each variable?
        numRecs =  pcdfi.new_longArray(1)
        for val in self._data['zVars']:
            a = pcdfi.PyCDFgetzVarAllocRecords(self._id,
                                               self._data['zVars'][val]['zVarnum'],
                                               numRecs)
            self._data['zVars'][val]['numRecs'] = pcdfi.longArray_getitem(numRecs, 0)
        pcdfi.delete_longArray(dimSizes)
        pcdfi.delete_longArray(recVary)                                                            
        pcdfi.delete_longArray(dimVarys)
        pcdfi.delete_longArray(dataType)
        pcdfi.delete_longArray(numElements)
        pcdfi.delete_longArray(numDims)
        pcdfi.delete_longArray(numRecs)                                            
        return True # add some error checking

    
    def _readzVar(self, zname_in, verbose=False):
        """
        read in the data and remaining metadata about a zVar

        @param zname_in: zVar name or number (or a list/tiple of names and/or numbers)
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        @todo: test more throurghly on row majority CDF files
        """
        import struct
        try: self.n_id
        except:
            print("CDF is not open")
            return False

        # quick and dirty hack here
        if not isinstance(zname_in, list) and not isinstance(zname_in, tuple):
            zname_in = [zname_in]
        for zname in zname_in:
            # this works for zvar names or numbers (based on string or not)
            if not isinstance(zname, str):
                # we were handed a zvar number get the associated name
                # this can probably be cleaned nicely, TODO
                found = False
                for zvarval in self._data['zVars']:
                    if self._data['zVars'][zvarval]['zVarnum'] == zname:
                        zname = zvarval
                        found = True
                        # they are unique so we can move on
                        break
            
            # TODO need some error checkoing here 

            # if we already have the data dont read it again
            try:
                self._data['zVars'][zvarval]['data']
                if verbose:
                    print("zVar %s already read, skipping" %(zvarval))
            except:
                # this works for Column majority, will need testing for row (but should work)
                if verbose:
                    print("\tReading data for: %s" % (zname))
                numRecs = self._data['zVars'][zname]['numRecs']
                dimSizes = self._data['zVars'][zname]['dimSizes']
                numDims = self._data['zVars'][zname]['numDims']
                dimSizes2 = 1
                indices = pcdfi.new_longArray(numDims)
                counts = pcdfi.new_longArray(numDims)
                intervals = pcdfi.new_longArray(numDims)
        
                for valii in dimSizes:
                    dimSizes2 = dimSizes2 * valii

                if  self._data['zVars'][zname]['dataType'] == 'CDF_CHAR' or \
                       self._data['zVars'][zname]['dataType'] == 'CDF_UCHAR':
                    buf = ' '*numRecs*dimSizes2
                    dataChar = True
                else:
                    # create the right kind of array
                    # buf = eval('pcdfi.new_%sArray(%d)' % ( self._data['zVars'][zname]['dataType'], numRecs*dimSizes2))
                    memSize = long(numRecs*dimSizes2*self.__dtypeSize[self._data['zVars'][zname]['dataType']])
                    buf = pcdfi.calloc_void(memSize)                    
                    dataChar = False

                                            

                    
                for ii in xrange(numDims):
                    pcdfi.longArray_setitem(indices, ii, 0) # start reading at the start of the variable
                    pcdfi.longArray_setitem(counts, ii, long(dimSizes[ii]))  # get the sizes
                    pcdfi.longArray_setitem(intervals, ii, 1)  ## read every entry

                a = pcdfi.PyCDFhyperGetzVarData(self._id,
                                                self._data['zVars'][zname]['zVarnum'],
                                                0,
                                                self._data['zVars'][zname]['numRecs'],
                                                1,
                                                indices,
                                                counts,
                                                intervals,
                                                buf)
                self._errorHandle(a)


                tmp = pcdfi.cdata(buf, memSize)
                        
                            
                if self._data['zVars'][zname]['dataType'] != 'CDF_CHAR' and \
                       self._data['zVars'][zname]['dataType'] != 'CDF_UCHAR' :
                    tmp = struct.unpack( self.__formatCodes[ \
                        self._data['zVars'][zname]['dataType']]*
                                         numRecs*dimSizes2,
                                         tmp)
                value = np.array(tmp)


            
                ## value = np.array([eval('pcdfi.%sArray_getitem(buf, %d)' % ( self._data['zVars'][zname]['dataType'], ic)) \
                ##                   for ic in xrange(numRecs*dimSizes2)])

                if self._data['GlobalAttr']['Majority'] == 'Column':
                    order = 'F'   # Fortran is Colum major
                    newShape = dimSizes
                    newShape = np.append(newShape, -1)
                else:
                    order = 'C'   # C is Row major
                    newShape = [-1]
                    newShape.extend(dimSizes)
                value = np.reshape(value, newShape, order=order)
                self._data['zVars'][zname]['data'] = value

                if verbose:
                    print("\t\tGetting local attributes")

                attrScope = pcdfi.new_longArray(1)
                maxgEntry = pcdfi.new_longArray(1)
                maxrEntry = pcdfi.new_longArray(1)
                maxzEntry = pcdfi.new_longArray(1)
                dataType = pcdfi.new_longArray(1)
                numElements = pcdfi.new_longArray(1)

                for i in xrange(self._data['GlobalAttr']['numAttrs']):
                    a, attrName = pcdfi.PyCDFinquireAttr(self._id,
                                                         long(i),
                                                         attrScope,
                                                         maxgEntry,
                                                         maxrEntry,
                                                         maxzEntry)
                    a = pcdfi.PyCDFinquireAttrzEntry (self._id,
                                                      i,
                                                      self._data['zVars'][zname]['zVarnum'],
                                                      dataType,
                                                      numElements)
                    if a != 0:  ## this means that it got a local, if gets a global a = -2076
                        continue
                    else:
                        # create a place to put the answer
                        memSize = pcdfi.longArray_getitem(numElements, 0)*self.__dtypeSize[str(pcdfi.longArray_getitem(dataType, 0))]
                        value = pcdfi.calloc_void(memSize)
                        a = pcdfi.PyCDFgetAttrzEntry(self._id,
                                                     long(i),
                                                     self._data['zVars'][zname]['zVarnum'],
                                                     value)
                        ## self._data['zVars'][zname][attrName] = pcdfi.cdata(value, memSize)
                        ## if self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] == 'CDF_REAL4' or \
                        ##       self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] == 'CDF_FLOAT' :
                        ##     self._data['zVars'][zname][attrName] = struct.unpack('f'* pcdfi.longArray_getitem(numElements, 0),
                        ##                                                         self._data['zVars'][zname][attrName])
                        ## if self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] == 'CDF_REAL8' or \
                        ##    self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] == 'CDF_DOUBLE' or \
                        ##    self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] == 'CDF_EPOCH':
                        ##     self._data['zVars'][zname][attrName] = struct.unpack('d'* pcdfi.longArray_getitem(numElements, 0),
                        ##                                                         self._data['zVars'][zname][attrName])

                        self._data['zVars'][zname][attrName] = pcdfi.cdata(value, memSize)
                        
                        if verbose:
                            print("\t\t\t%s %s %s %s %d" % (zname,
                                                            attrName,
                                                            pcdfi.cdata(value, memSize),
                                                            self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))],
                                                            pcdfi.longArray_getitem(numElements, 0)  ))
                            
                        if self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] != 'CDF_CHAR' and \
                        self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))] != 'CDF_UCHAR' :
                            self._data['zVars'][zname][attrName] = struct.unpack( self.__formatCodes[ \
                                self.__datatypes[str(pcdfi.longArray_getitem(dataType, 0))]]*
                                                                                 pcdfi.longArray_getitem(numElements, 0),
                                                                                 self._data['zVars'][zname][attrName])




                            
                    
                try: # dont delete them if they dont exist
                    attrScope
        
                    pcdfi.delete_longArray(attrScope)
                    pcdfi.delete_longArray(maxgEntry)
                    pcdfi.delete_longArray(maxrEntry)
                    pcdfi.delete_longArray(maxzEntry)
                    pcdfi.delete_longArray(dataType)
                    pcdfi.delete_longArray(numElements)
                    pcdfi.free_void(value)
                    if not dataChar: pcdfi.free_void(buf)

                except:
                    pass
        return True  # iguess it worked if we get here

    def _readAllzVars(self, verbose=False):
        """
        read data and metadata from all zVars in the open CDF file
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        try: self.n_id
        except:
            print("CDF is not open")
            return False

        for val in xrange(self._data['GlobalAttr']['numzVars']):
            print('Reading zVar %d of %d' % (val, self._data['GlobalAttr']['numzVars']-1))
            self._readzVar(val, verbose=verbose)
        return True # TODO needs some error checking


    def _close(self, verbose=False):
        """
        close the open CDF file
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        try: self._id
        except AttributeError:
            print("No cdf file open, cannot close")
            return False
        pcdfi.PyCDFcloseCDF(self._id)
        # clean up a few variables  so what we know the file is closed
        del self._id
        del self.n_id
        return True


    def _getzVarDepend(self, var_in, verbose=False):
        """
        read in the zVar data and metadata for the specified zVar's dependencies

        @param var_in: variable name or list/tuple of variable names
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        """
        if not isinstance(var_in, list) and not isinstance(var_in, tuple):
            var_in = [var_in]
        for val in var_in:
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_0'], verbose=verbose)
                if verbose: print("found a DEPEND_0 = %s" % (self._data['zVars'][val]['DEPEND_0']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_1'])
                if verbose: print("found a DEPEND_1 = %s" % (self._data['zVars'][val]['DEPEND_1']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_2'])
                if verbose: print("found a DEPEND_2 = %s" % (self._data['zVars'][val]['DEPEND_2']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_3'])
                if verbose: print("found a DEPEND_3 = %s" % (self._data['zVars'][val]['DEPEND_3']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_4'])
                if verbose: print("found a DEPEND_4 = %s" % (self._data['zVars'][val]['DEPEND_4']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_5'])
                if verbose: print("found a DEPEND_5 = %s" % (self._data['zVars'][val]['DEPEND_5']))
            except:
                pass
            try:
                self.getzVar(self._data['zVars'][val]['DEPEND_6'])
                if verbose: print("found a DEPEND_6 = %s" % (self._data['zVars'][val]['DEPEND_6']))
            except:
                pass
        return True
            

    def _createAttr(self, attrName, scope=pcdfi.GLOBAL_SCOPE, valNum=None, verbose=False):
        """
        I need docs
        """
        attrNum = pcdfi.new_longArray(1)
        a = pcdfi.PyCDFcreateAttr(self._id, attrName, scope, attrNum)
        self._errorHandle(a)  # maybe a strickerror in _errorHandle
        self._data['GlobalAttr'][attrName] = {}
        self._data['GlobalAttr'][attrName]['attrNum'] = pcdfi.longArray_getitem(attrNum, 0)
        self._data['GlobalAttr'][attrName]['value'] = None

        pcdfi.delete_longArray(attrNum)
        
        return True # add some error checking here

    def _putAttrVal(self, attrNameIN, attrValue, valNum=None, verbose=False):
        """
        I need docs
        """
        ## self._readGlobalAttributes(verbose=verbose)
        ## try:
        ##     self._data['GlobalAttr'][attrName]
        ## except:
        ##     print("attrName %s doesnt exist, must create it first with _createAttr()" % (attrName))
        ##     return False
        if isinstance(attrValue, str):
            dataType = self.__datatypesReverse['CDF_CHAR']
            if valNum == None:
                valNum = 0 # this si placeholder for search for the right answer
            # allocate memory for the void * input
            tmp = pcdfi.calloc_void(len(attrValue))  # do I need a +1 for a trailing null?
            # move the data to the location of the pointer
            pcdfi.memmove(tmp, attrValue)
            # make the call
            a = pcdfi.PyCDFputAttrgEntry(self._id,
                                         self._data['GlobalAttr'][attrNameIN]['attrNum'],
                                         0,
                                         long(dataType),
                                         len(attrValue),
                                         tmp)
        # add some elif here
        else:
            print("Only string supported now, sorry")
            return False

        self._errorHandle(a)
        pcdfi.free_void(tmp)

        # add the new info to the dict
        self._data['GlobalAttr'][attrNameIN]['numElements'] = len(attrValue)
        self._data['GlobalAttr'][attrNameIN]['value'] = attrValue
        
        return True
        
        
    def getGlobalAttr(self, verbose=False):
        """
        read Global attributes from the CDF file, they are saved in cdf.GlobalAttr
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        >>> cdf = cdfpy.cdfpy('l4_h0_mpa_20011111_v02.cdf')
        >>> cdf.getGlobalAttr()
        >>> cdf.GlobalAttr
        {'ADID_ref': {'attrNum': 8,
        'attrScope': 1,
        'dataType': 'CDF_CHAR',
        'numElements': 8,
        'value': 'NSSD0099'},
        'Data_type': {'attrNum': 5,
        'attrScope': 1,
        'dataType': 'CDF_CHAR',
        'numElements': 18,
        'value': 'H0>High Resolution'},
        'Data_version': {'attrNum': 6,
        'attrScope': 1,
        'dataType': 'CDF_CHAR',
        'numElements': 1,
        'value': ' '},
        
        """
        self._open(verbose=verbose)
        self._readCDFInfo(verbose=verbose)
        self._readGlobalAttributes(verbose=verbose)
        self._close(verbose=verbose)
        self.GlobalAttr = self._data['GlobalAttr']
        return True

    def listzVars(self, verbose=False):
        """
        list the zVars contained in the CDF file, they are saved in cdf.GlobalAttr
        
        @keyword verbose: be verbose to the screen
        @type verbose: bool
 
        @return: list of zVars
        @rtype: list
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        >>> cdf = cdfpy.cdfpy('l4_h0_mpa_20011111_v02.cdf')
        >>> cdf.listzVars() 
        ['etimesg',
        'Theta_l',
        'Pcounts',
        'Ecounts',
        'Theta_u',
        'thetaphi_hip',
        'Poten',
        'backgrd',
        'vel_lop',
        'dens_hip',
        'ece',
        'temp_e',
        'effecp',
        'Flags',
        'ecp',
        'temp_hip',
        'xyzgeo',
        'effece',
        'thetaphi_e',
        'Azanglp',
        'Epoch',
        'Uthours',
        'dens_lop',
        'dens_e',
        'Azangle',
        'Inval_mf']

        
        """
        self._open(verbose=verbose)
        self._readCDFInfo(verbose=verbose)
        self._readGlobalAttributes(verbose=verbose)
        self._getzVarInfo(verbose=verbose)
        self._close(verbose=verbose)
        zvars = [val for val in self._data['zVars']]
        self.zVars = zvars
        return self.zVars

    def getzVar(self, var_in, retval=False, depend=False, keepEpoch=False, verbose=False):
        """
        get the contents of a zVar (or list of zVars) from the CDF file, they are saved in cdf.data

        @param var_in: name of varialbe to read (or list/tuple of names)
        @type var_in: str
        @keyword retval: return the data dont just add it to the dict
        @type retval: bool
        @keyword depend: also get the dependent zVars
        @type depend: bool
        @keyword keepEpoch: keep the cdf epoch, dont change it to Ticktock
        @type keepEpoch: bool
        @keyword verbose: be verbose to the screen
        @type verbose: bool
        @return: True for success
        @rtype: bool
        
        @author: Brian Larsen
        @organization: Los Alamos National Lab
        @contact: balarsen@lanl.gov
        
        @version: V1: 20-Jul-2010 (BAL)

        >>> cdf = cdfpy.cdfpy('l4_h0_mpa_20011111_v02.cdf')
        >>> cdf.getzVar('temp_e', depend=True)

        
        """
        self._open(verbose=verbose)
        self._readCDFInfo(verbose=verbose)
        self._readGlobalAttributes(verbose=verbose)
        self._getzVarInfo(verbose=verbose)
        self._readzVar(var_in, verbose=verbose)
        self._close(verbose=verbose)
        # TODO fix this up nbut I cant figure it out to work woth both
        # var names and numbers
        if not isinstance(var_in, list) and not isinstance(var_in, tuple):
            var_in = [var_in]
        for val in var_in:
            try: self.data
            except: self.data = {}
            self.data[val] = self._data['zVars'][val]
            if depend:
                self._getzVarDepend(val, verbose=verbose)
            if val == 'Epoch':
                tt = self.epoch2Ticktock(self.data['Epoch']['data'])
                if keepEpoch:
                    self.data['Epoch']['Ticktock'] = tt
                else:
                    self.data['Epoch']['data'] = tt
                self.epoch = tt
        if retval:
            return self._data['zVars'][val]
        else:
            return True
    
            
        
        
        
                                   
        

if __name__ == "__main__":


    import cdfpy
    fname = 'l4_h0_mpa_20011111_v02.cdf'
    # initialize the class
    cdf = cdfpy.cdfpy(fname)

    # list the zVars so we know which to grab
    cdf.listzVars()

    
    # get dens_e and its dependencies
    cdf.getzVar('dens_e', depend=True)
    
    # epoch was also grabbed as this was a cdf that was put together right

    # also get the global attr so we have more info for the plot
    cdf.getGlobalAttr()



    # plot dens_e versus time
    import matplotlib.pyplot as plt
    fig=plt.figure(1)
    ax=fig.add_subplot(111)
    pp = ax.plot(cdf.epoch.UTC, cdf.data['dens_e']['data'], lw=2)


    ax.set_xlabel(cdf.data['Epoch']['LABLAXIS'] )
    ax.set_ylabel(cdf.data['dens_e']['LABLAXIS'] + ' (' + cdf.data['dens_e']['UNITS'] + ')')
    ax.set_title(cdf._data['GlobalAttr']['TITLE']['value'] + '\n' + cdf.epoch[0].ISO[0] )
    Mt, mt, fmt = spacepy.toolbox.smartTimeTicks(cdf.epoch.UTC)

    ax.xaxis.set_major_locator(Mt)
    ax.xaxis.set_minor_locator(mt)
    ax.xaxis.set_major_formatter(fmt)

    plt.draw()




    
