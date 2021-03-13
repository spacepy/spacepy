#!/usr/bin/env python

"""Unit testing for istp support"""

import datetime
import inspect
import os.path
import shutil
import sys
import tempfile
import unittest
import warnings

import numpy
import numpy.testing
import spacepy_testing
import spacepy
import spacepy.pycdf
import spacepy.pycdf.const
import spacepy.pycdf.istp


class ISTPTestsBase(unittest.TestCase):
    """Base class for ISTP check functions, provide a blank CDF to modify"""

    def setUp(self):
        """Setup: make an empty, open, writeable CDF"""
        self.tempdir = tempfile.mkdtemp()
        # We know what the backward-compatible default is, suppress it.
        warnings.filterwarnings(
            'ignore',
            message=r'^spacepy\.pycdf\.lib\.set_backward not called.*$',
            category=DeprecationWarning,
            module='^spacepy.pycdf$')
        try:
            self.cdf = spacepy.pycdf.CDF(os.path.join(
                self.tempdir, 'source_descriptor_datatype_19990101_v00.cdf'),
                                         create=True)
        finally:
            del warnings.filters[0]

    def tearDown(self):
        """Delete the empty cdf"""
        try:
            self.cdf.close()
        except:
            pass
        shutil.rmtree(self.tempdir)


class VariablesTests(ISTPTestsBase):
    """Tests of variable-checking functions"""
    longMessage = True

    def testAllVarFailure(self):
        """Call variable checks with a known bad one"""
        class BadTestClass(spacepy.pycdf.istp.VariableChecks):
            @classmethod
            def varraiseserror(cls, v):
                raise RuntimeError('Bad')
        data = spacepy.dmarray([1, 2, 3], dtype=numpy.int8, attrs={
            'FIELDNAM': 'var1',
            'FILLVAL': -128,
            })
        var = self.cdf.new('var1', data=data)
        errs = BadTestClass.all(var, catch=True)
        self.assertEqual(
            ['Test varraiseserror did not complete.'],
            errs)

    def testEmptyEntries(self):
        """Are there any CHAR entries of empty string"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf['var1'].attrs['attribute'] = ''
        errs = spacepy.pycdf.istp.VariableChecks.empty_entry(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual('Empty CHAR entry for attribute attribute.',
                         errs[0])

    def testDepends(self):
        """Dependencies, valid and not"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf['var1'].attrs['DEPEND_0'] = 'var2'
        errs = spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual('DEPEND_0 variable var2 missing.', errs[0])
        self.cdf['var2'] = [1, 2, 3]
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])))
        self.cdf['var1'].attrs['DELTA_PLUS_VAR'] = 'foobar'
        errs = spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual('DELTA_PLUS_VAR variable foobar missing.', errs[0])

    def testDeltas(self):
        """DELTA variables"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf['var1'].attrs['DELTA_PLUS_VAR'] = 'var2'
        self.cdf['var2'] = [1., 1, 1]
        errs = spacepy.pycdf.istp.VariableChecks.deltas(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual([
            'DELTA_PLUS_VAR type CDF_FLOAT does not match variable type '
            'CDF_BYTE.'], errs)
        del self.cdf['var2']
        self.cdf['var2'] = [1, 1, 1]
        self.cdf['var1'].attrs['UNITS'] = 'smoots'
        errs = spacepy.pycdf.istp.VariableChecks.deltas(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual([
            'DELTA_PLUS_VAR units do not match variable units.'], errs)
        self.cdf['var2'].attrs['UNITS'] = 'smoots'
        errs = spacepy.pycdf.istp.VariableChecks.deltas(self.cdf['var1'])
        self.assertEqual(0, len(errs))

    def testDeltasMisdimmed(self):
        """DELTA variable dimension size mismatch"""
        self.cdf['var1'] = [[1, 2, 3], [4, 5, 6]]
        self.cdf['var1'].attrs['DELTA_PLUS_VAR'] = 'var2'
        self.cdf.new('var2', data=[1, 1], recVary=False)
        errs = spacepy.pycdf.istp.VariableChecks.deltas(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual(
            'DELTA_PLUS_VAR shape (2,) does not match variable shape (3,).',
            errs[0])
        del self.cdf['var2']
        self.cdf.new('var2', data=[[1, 1, 1]])
        errs = spacepy.pycdf.istp.VariableChecks.deltas(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual(
            'DELTA_PLUS_VAR record count 1 does not match variable record'
            ' count 2.',
            errs[0])

    def testValidRangeDimensioned(self):
        """Validmin/validmax with multiple elements"""
        v = self.cdf.new('var1', data=[[1, 10], [2, 20], [3, 30]])
        v.attrs['VALIDMIN'] = [1, 20]
        v.attrs['VALIDMAX'] = [3, 30]
        v.attrs['FILLVAL'] = -100
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 10 at index [0 1] under VALIDMIN [ 1 20].',
                         errs[0])
        v.attrs['VALIDMIN'] = [1, 10]
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(0, len(errs))
        v[0, 0] = -100
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(0, len(errs))

    def testValidRangeDimensionMismatch(self):
        """Validmin/validmax with something wrong in dimensionality"""
        v = self.cdf.new('var1', data=[[1, 10], [2, 20], [3, 30]])
        v.attrs['VALIDMIN'] = [1, 10, 100]
        v.attrs['VALIDMAX'] = [3, 30, 127]
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(2, len(errs))
        self.assertEqual('VALIDMIN element count 3 does not match '
                         'first data dimension size 2.', errs[0])
        self.assertEqual('VALIDMAX element count 3 does not match '
                         'first data dimension size 2.', errs[1])

    def testValidRangeHighDimension(self):
        """Validmin/validmax with high-dimension variables"""
        v = self.cdf.new('var1',
                         data=numpy.reshape(numpy.arange(27.), (3, 3, 3,)))
        v.attrs['VALIDMIN'] = [1, 10, 100]
        v.attrs['VALIDMAX'] = [3, 30, 300]
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(2, len(errs))
        self.assertEqual('Multi-element VALIDMIN only valid with 1D variable.',
                         errs[0])
        self.assertEqual('Multi-element VALIDMAX only valid with 1D variable.',
                         errs[1])

    def testValidRangeWrongType(self):
        """Validmin/validmax not matching variable type"""
        v = self.cdf.new('var1', data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_INT4)
        v.attrs.new('VALIDMIN', data=1, type=spacepy.pycdf.const.CDF_INT2)
        v.attrs.new('VALIDMAX', data=3, type=spacepy.pycdf.const.CDF_INT2)
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        errs.sort()
        self.assertEqual(2, len(errs))
        self.assertEqual(
            ['VALIDMAX type CDF_INT2 does not match variable type CDF_INT4.',
             'VALIDMIN type CDF_INT2 does not match variable type CDF_INT4.'],
            errs)

    def testValidRangeIncompatibleType(self):
        """Validmin/validmax can't be compared to variable type"""
        v = self.cdf.new('var1', data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_INT4)
        v.attrs.new('VALIDMIN', data='2')
        v.attrs.new('VALIDMAX', data='5')
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        errs.sort()
        self.assertEqual(4, len(errs))
        self.assertEqual(
            ['VALIDMAX type CDF_CHAR does not match variable type CDF_INT4.',
             'VALIDMAX type CDF_CHAR not comparable to variable type CDF_INT4.',
             'VALIDMIN type CDF_CHAR does not match variable type CDF_INT4.',
             'VALIDMIN type CDF_CHAR not comparable to variable type CDF_INT4.'
            ],
            errs)

    def testValidRangeNRV(self):
        """Validmin/validmax"""
        v = self.cdf.new('var1', recVary=False, data=[1, 2, 3])
        v.attrs['VALIDMIN'] = 1
        v.attrs['VALIDMAX'] = 3
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))
        v.attrs['VALIDMIN'] = 2
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 1 at index 0 under VALIDMIN 2.', errs[0])
        v.attrs['VALIDMAX'] = 2
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(2, len(errs))
        self.assertEqual('Value 3 at index 2 over VALIDMAX 2.', errs[1])

    def testValidRangeNRVFillval(self):
        """Validmin/validmax with fillval set"""
        v = self.cdf.new('var1', recVary=False, data=[1, 2, 3])
        v.attrs['VALIDMIN'] = 1
        v.attrs['VALIDMAX'] = 3
        v.attrs['FILLVAL'] = 99
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))
        
        v.attrs['VALIDMIN'] = 2
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 1 at index 0 under VALIDMIN 2.', errs[0])
        
        v.attrs['VALIDMAX'] = 2
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(2, len(errs))
        self.assertEqual('Value 3 at index 2 over VALIDMAX 2.', errs[1])
        
        v.attrs['FILLVAL'] = 3
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 1 at index 0 under VALIDMIN 2.', errs[0])

        v.attrs['FILLVAL'] = 1
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 3 at index 2 over VALIDMAX 2.', errs[0])

    def testValidRangeFillvalFloat(self):
        """Validmin/validmax with fillval set, floating-point"""
        v = self.cdf.new('var1', recVary=False, data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_DOUBLE)
        v.attrs['VALIDMIN'] = 0
        v.attrs['VALIDMAX'] = 10
        #This is a bit contrived to force a difference between attribute
        #and value that's only the precision of the float
        v.attrs.new('FILLVAL', -1e31, type=spacepy.pycdf.const.CDF_FLOAT)
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))

        v[0] = -100
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value -100.0 at index 0 under VALIDMIN 0.0.', errs[0])

        v[0] = -1e31
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))

    def testValidRangeFillvalFloatWrongType(self):
        """Validmin/validmax with fillval, floating-point, but fillval string"""
        v = self.cdf.new('var1', recVary=False, data=[-1e31, 2, 3],
                         type=spacepy.pycdf.const.CDF_DOUBLE)
        v.attrs['VALIDMIN'] = 0
        v.attrs['VALIDMAX'] = 10
        v.attrs.new('FILLVAL', b'badstuff', type=spacepy.pycdf.const.CDF_CHAR)
        expected = ['Value -1e+31 at index 0 under VALIDMIN 0.0.']
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(len(expected), len(errs))
        for a, e in zip(sorted(errs), sorted(expected)):
            self.assertEqual(e, a)

    def testValidRangeFillvalDatetime(self):
        """Validmin/validmax with fillval set, Epoch var"""
        v = self.cdf.new(
            'var1', data=[datetime.datetime(2010, 1, i) for i in range(1, 6)],
            type=spacepy.pycdf.const.CDF_EPOCH)
        v.attrs['VALIDMIN'] = datetime.datetime(2010, 1, 1)
        v.attrs['VALIDMAX'] = datetime.datetime(2010, 1, 31)
        v.attrs['FILLVAL'] = datetime.datetime(9999, 12, 31, 23, 59, 59, 999000)
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))

        v[-1] = datetime.datetime(2010, 2, 1)
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Value 2010-02-01 00:00:00 at index 4 over VALIDMAX '
                         '2010-01-31 00:00:00.', errs[0])

        v[-1] = datetime.datetime(9999, 12, 31, 23, 59, 59, 999000)
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))

    def testFillval(self):
        """Test for fillval presence, type, value"""
        v = self.cdf.new('var1', data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_INT2)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(1, len(errs))
        self.assertEqual(
            'No FILLVAL attribute.',
            errs[0])
        v.attrs.new('FILLVAL', -5, type=spacepy.pycdf.const.CDF_DOUBLE)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual([
            'FILLVAL -5.0, should be -32768 for variable type CDF_INT2.',
            'FILLVAL type CDF_DOUBLE does not match variable type CDF_INT2.',
            ],
            sorted(errs))
        v.attrs['FILLVAL'] = -32768
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual([
            'FILLVAL type CDF_DOUBLE does not match variable type CDF_INT2.',
            ],
            sorted(errs))
        del v.attrs['FILLVAL']
        v.attrs.new('FILLVAL', -32768, type=spacepy.pycdf.const.CDF_INT2)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs))

    def testFillvalFloat(self):
        """Test for fillval being okay when off by float precision"""
        v = self.cdf.new('var1', data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_FLOAT)
        v.attrs.new('FILLVAL', -1e31, type=spacepy.pycdf.const.CDF_FLOAT)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs), '\n'.join(errs))

    def testFillvalString(self):
        """Test for fillval being okay in a string"""
        v = self.cdf.new('var1', data=['foo', 'bar'],
                                       type=spacepy.pycdf.const.CDF_CHAR)
        v.attrs.new('FILLVAL', ' ', type=spacepy.pycdf.const.CDF_CHAR)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs), '\n'.join(errs))

    def testFillvalEpoch(self):
        """Test for fillval being okay with epoch"""
        v = self.cdf.new('Epoch', type=spacepy.pycdf.const.CDF_EPOCH)
        v.attrs.new(
            'FILLVAL', -1e31,
            type=spacepy.pycdf.const.CDF_EPOCH)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs), '\n'.join(errs))
        del v.attrs['FILLVAL']
        v.attrs.new(
            'FILLVAL', datetime.datetime(2000, 1, 1),
            type=spacepy.pycdf.const.CDF_EPOCH)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(1, len(errs), '\n'.join(errs))
        self.assertEqual(
            'FILLVAL {} (2000-01-01 00:00:00), should be -1e+31 '
            '(9999-12-31 23:59:59.999000) for variable type CDF_EPOCH.'
            .format(6.3113904e+13), #py2k and 3k format this differently
            errs[0])

    def testMatchingRecordCount(self):
        """Same number of records for DEPEND_0"""
        e = self.cdf.new(
            'Epoch', type=spacepy.pycdf.const.CDF_EPOCH,
            data=[datetime.datetime(2010, 1, i + 1) for i in range(5)])
        v = self.cdf.new('Data', data=list(range(4)))
        v.attrs['DEPEND_0'] = 'Epoch'
        errs = spacepy.pycdf.istp.VariableChecks.recordcount(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('4 records; DEPEND_0 Epoch has 5.',
                         errs[0])
        v.append(5)
        errs = spacepy.pycdf.istp.VariableChecks.recordcount(v)
        self.assertEqual(0, len(errs))

    def testDepSize(self):
        """Depends are appropriately sized"""
        e = self.cdf.new(
            'Epoch', type=spacepy.pycdf.const.CDF_EPOCH,
            data=[datetime.datetime(2010, 1, i + 1) for i in range(5)])
        v = self.cdf.new('Data', data=[[1, 2], [3, 4], [5, 6]])
        v.attrs['DEPEND_0'] = 'Epoch'
        d = self.cdf.new('thedep', recVary=False, data=[1, 2, 3])
        v.attrs['DEPEND_1'] = 'thedep'
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Dim 1 sized 2 but DEPEND_1 thedep sized 3.', errs[0])
        del self.cdf['thedep']
        d = self.cdf.new('thedep', recVary=False, data=[1, 2])
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(0, len(errs))

    def testDepSize2(self):
        """Depends are appropriately sized, depend_2"""
        e = self.cdf.new(
            'Epoch', type=spacepy.pycdf.const.CDF_EPOCH,
            data=[datetime.datetime(2010, 1, i + 1) for i in range(5)])
        v = self.cdf.new('Data', data=[[[1, 2], [3, 4], [5, 6]]])
        v.attrs['DEPEND_0'] = 'Epoch'
        self.cdf.new('thedep', recVary=False, data=[1, 2, 3])
        v.attrs['DEPEND_1'] = 'thedep'
        self.cdf.new('dep2', recVary=False, data=[1])
        v.attrs['DEPEND_2'] = 'dep2'
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('Dim 2 sized 2 but DEPEND_2 dep2 sized 1.', errs[0])
        del self.cdf['dep2']
        self.cdf.new('dep2', recVary=False, data=[1, 2])
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(0, len(errs))

    def testMultiLayeredDep(self):
        """Check multi layered dependency."""
        # First Simple Case
        e = self.cdf.new(
            'Epoch', type=spacepy.pycdf.const.CDF_EPOCH,
            data=[datetime.datetime(2010, 1, i + 1) for i in range(5)])
        v = self.cdf.new('Data', data=[[[1,2,3],[4,5,6]],
                                       [[11,12,13],[14,15,16]],
                                       [[21,22,23],[24,25,26]],
                                       [[31,32,33],[34,35,36]],
                                       [[41,42,43],[44,45,46]]])
        v.attrs['DEPEND_0'] = 'Epoch'
        self.cdf.new('dep1', recVary=False, data=[10,20])
        v.attrs['DEPEND_1'] = 'dep1'
        self.cdf.new('dep2', recVary=False, data=[30,40,50])
        v.attrs['DEPEND_2'] = 'dep2'
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(0, len(errs))
        # Now make layered
        del self.cdf['dep2']
        vv = self.cdf.new('dep2', recVary=False, data=[[333,444,555],[666,777,888]])
        vv.attrs['DEPEND_1'] = 'dep1'
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        # Now make fail
        del self.cdf['dep1']
        vv = self.cdf.new('dep1', recVary=False, data=[10,20,30])
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual('Dim 1 sized 2 but DEPEND_1 dep1 sized 3.', errs[0])
        self.assertEqual('Dim 2 sized 3 but DEPEND_2 dep2 sized 2.', errs[1])
        
    def testDepSizeBadRV(self):
        """Depend size check with RV, non-depend 0 dep"""
        e = self.cdf.new(
            'Epoch', type=spacepy.pycdf.const.CDF_EPOCH,
            data=[datetime.datetime(2010, 1, i + 1) for i in range(5)])
        v = self.cdf.new('Data', data=[[1, 2], [3, 4], [5, 6]])
        v.attrs['DEPEND_0'] = 'Epoch'
        d = self.cdf.new('thedep', recVary=True, data=[1, 2, 3])
        v.attrs['DEPEND_1'] = 'thedep'
        errs = spacepy.pycdf.istp.VariableChecks.depsize(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('DEPEND_1 thedep is RV but has no DEPEND_0.', errs[0])

    def testValidRangeScalar(self):
        """Check validmin/max on a scalar"""
        v = self.cdf.new('var1', recVary=False, data=1)
        v.attrs['VALIDMIN'] = 0
        v.attrs['VALIDMAX'] = 2
        v.attrs['FILLVAL'] = -100
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validrange(v)))
        v.attrs['VALIDMIN'] = 2
        v.attrs['VALIDMAX'] = 3
        errs = spacepy.pycdf.istp.VariableChecks.validrange(v)
        self.assertEqual(1, len(errs))
        self.assertEqual(
            'Value 1 under VALIDMIN 2.', errs[0])

    def testValidScale(self):
        """Check scale min and max."""
        v = self.cdf.new('var1', recVary=False, data=[1, 2, 3])
        v.attrs['SCALEMIN'] = 1
        v.attrs['SCALEMAX'] = 3
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validscale(v)))
        v.attrs['SCALEMIN'] = 5
        v.attrs['SCALEMAX'] = 3
        self.assertEqual(
            1, len(spacepy.pycdf.istp.VariableChecks.validscale(v)))
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual('SCALEMIN > SCALEMAX.', errs[0])
        v.attrs['SCALEMIN'] = -200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(2, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMIN (-200) outside valid data range (-128,127).',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMIN'] = 200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(3, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMIN (200) outside valid data range (-128,127).',
             'SCALEMIN > SCALEMAX.',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMAX'] = -200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(5, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMAX (-200) outside valid data range (-128,127).',
             'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.',
             'SCALEMIN (200) outside valid data range (-128,127).',
             'SCALEMIN > SCALEMAX.',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMAX'] = 200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(4, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMAX (200) outside valid data range (-128,127).',
             'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.',
             'SCALEMIN (200) outside valid data range (-128,127).',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        
    def testValidScaleDimensioned(self):
        """Validmin/validmax with multiple elements"""
        v = self.cdf.new('var1', data=[[1, 10], [2, 20], [3, 30]])
        v.attrs['SCALEMIN'] = [2, 20]
        v.attrs['SCALEMAX'] = [300, 320]
        v.attrs['FILLVAL'] = -100
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(2, len(errs))
        errs.sort()
        self.assertEqual([
            'SCALEMAX ([300 320]) outside valid data range (-128,127).',
            'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.'
            ],
            errs)
        v.attrs['SCALEMAX'] = [30, 32]
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(1, len(errs))
        self.assertEqual([
            'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.'
            ],
            errs)

    def testValidScaleDimensionMismatch(self):
        """Validmin/validmax with something wrong in dimensionality"""
        v = self.cdf.new('var1', data=[[1, 10], [2, 20], [3, 30]])
        v.attrs['SCALEMIN'] = [1, 10, 100]
        v.attrs['SCALEMAX'] = [3, 30, 126]
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(2, len(errs))
        errs.sort()
        self.assertEqual([
            'SCALEMAX element count 3 does not match '
            'first data dimension size 2.',
            'SCALEMIN element count 3 does not match '
            'first data dimension size 2.',
            ],
            errs)

    def testValidScaleHighDimension(self):
        """scalemin/scalemax with high-dimension variables"""
        v = self.cdf.new('var1',
                         data=numpy.reshape(numpy.arange(27.), (3, 3, 3,)))
        v.attrs['SCALEMIN'] = [1, 10, 100]
        v.attrs['SCALEMAX'] = [3, 30, 300]
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(2, len(errs))
        self.assertEqual('Multi-element SCALEMIN only valid with 1D variable.',
                         errs[0])
        self.assertEqual('Multi-element SCALEMAX only valid with 1D variable.',
                         errs[1])

    def testValiddisplaytype(self):
        """Check plot type."""
        err1 = '1 dim variable with spectrogram display type.'
        err2 = 'Multi dim variable with time_series display type.'
        v = self.cdf.new('var1', recVary=True, data=[1, 2, 3])
        v.attrs['DISPLAY_TYPE'] = 'time_series'
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)))
        v.attrs['DISPLAY_TYPE'] = 'spectrogram'
        self.assertEqual(
            1, len(spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)))
        errs = spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)
        self.assertEqual(err1, errs[0])
        v = self.cdf.new('var2', recVary=True, data=[[1, 2, 3],[4,5,6]])
        v.attrs['DISPLAY_TYPE'] = 'spectrogram'
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)))
        v.attrs['DISPLAY_TYPE'] = 'time_series'
        self.assertEqual(
            1, len(spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)))
        errs = spacepy.pycdf.istp.VariableChecks.validdisplaytype(v)
        self.assertEqual(err2, errs[0])

    def testFieldnam(self):
        """Check field name matching"""
        v = self.cdf.new('var1', recVary=True, data=[1, 2, 3])
        errs = spacepy.pycdf.istp.VariableChecks.fieldnam(v)
        self.assertEqual(1, len(errs))
        self.assertEqual('No FIELDNAM attribute.', errs[0])
        v.attrs['FIELDNAM'] = 'var2'
        errs = spacepy.pycdf.istp.VariableChecks.fieldnam(v)
        self.assertEqual(1, len(errs))
        self.assertEqual(
            'FIELDNAM attribute var2 does not match var name.', errs[0])
        v.attrs['FIELDNAM'] = 'var1'
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.fieldnam(v)))


@unittest.skipIf(spacepy.pycdf.lib.version[0] < 3,
                 'Requires CDF library 3 or newer')
class VariablesTestsNew(ISTPTestsBase):
    """Tests of variable-checking functions that require CDF 3"""
    longMessage = True

    def setUp(self):
        """Disable backward-compat before making test CDF"""
        spacepy.pycdf.lib.set_backward(False)
        super(VariablesTestsNew, self).setUp()
        spacepy.pycdf.lib.set_backward(True)

    def testFillvalEpoch16(self):
        """Test for fillval being okay with epoch16"""
        v = self.cdf.new('Epoch', type=spacepy.pycdf.const.CDF_EPOCH16)
        v.attrs.new(
            'FILLVAL', (-1e31, -1e31),
            type=spacepy.pycdf.const.CDF_EPOCH16)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs), '\n'.join(errs))

    @unittest.skipIf(not spacepy.pycdf.lib.supports_int8,
                     'Requires TT2000 support in CDF library')
    def testFillvalTT2000(self):
        """Test for fillval being okay with TT2000"""
        v = self.cdf.new('Epoch', type=spacepy.pycdf.const.CDF_TIME_TT2000)
        v.attrs.new(
            'FILLVAL', -9223372036854775808,
            type=spacepy.pycdf.const.CDF_TIME_TT2000)
        errs = spacepy.pycdf.istp.VariableChecks.fillval(v)
        self.assertEqual(0, len(errs), '\n'.join(errs))

        
class FileTests(ISTPTestsBase):
    """Tests of file-checking functions"""

    def testAll(self):
        """Calling all file checks"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf['var1'].attrs['DEPEND_0'] = 'var2'
        self.cdf['var1'].attrs['FIELDNAM'] = 'var1'
        self.cdf.attrs['Logical_source'] = \
            'source_descriptor_datatype'
        self.cdf.attrs['Logical_file_id'] = \
            'source_descriptor_datatype_19990101_v00'
        errs = spacepy.pycdf.istp.FileChecks.all(self.cdf)
        self.assertEqual(2, len(errs))
        self.assertEqual([
            'var1: DEPEND_0 variable var2 missing.',
            'var1: No FILLVAL attribute.',
            ],
            sorted(errs))

    def testAllFailure(self):
        """Call file checks with a known bad one"""
        self.cdf.attrs['Logical_source'] = \
            'source_descriptor_datatype'
        self.cdf.attrs['Logical_file_id'] = \
            'source_descriptor_datatype_19990101_v00'
        class BadTestClass(spacepy.pycdf.istp.FileChecks):
            @classmethod
            def raiseserror(cls, f):
                raise RuntimeError('Bad')
        errs = BadTestClass.all(self.cdf, catch=True)
        self.assertEqual(
            ['Test raiseserror did not complete.'],
            errs)

    def testEmptyEntries(self):
        """Are there any CHAR gEntries of empty string"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf.attrs['attribute'] = ''
        errs = spacepy.pycdf.istp.FileChecks.empty_entry(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Empty CHAR entry 0 for attribute attribute.',
                         errs[0])

    def testFilename(self):
        """Compare filename to global attrs"""
        self.cdf.attrs['Logical_source'] = \
            'source_descriptor_blargh'
        self.cdf.attrs['Logical_file_id'] = \
            'source_descriptor_datatype_19990102_v00'
        errs = spacepy.pycdf.istp.FileChecks.filename(self.cdf)
        self.assertEqual(2, len(errs))
        self.assertTrue("Logical_source source_descriptor_blargh doesn't match "
                        "filename source_descriptor_datatype_19990101_v00.cdf."
                        in errs)
        self.assertTrue(
            "Logical_file_id source_descriptor_datatype_19990102_v00 doesn't "
            "match filename source_descriptor_datatype_19990101_v00.cdf."
            in errs)
        self.cdf.attrs['Logical_source'] = \
            'source_descriptor_datatype'
        self.cdf.attrs['Logical_file_id'] = \
            'source_descriptor_datatype_19990101_v00'
        errs = spacepy.pycdf.istp.FileChecks.filename(self.cdf)
        self.assertEqual(0, len(errs))
        del self.cdf.attrs['Logical_source'][0]
        errs = spacepy.pycdf.istp.FileChecks.filename(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual(['No Logical_source in global attrs.'],
                         errs)

    def testTimesMonoton(self):
        """Test monotonic time"""
        self.cdf.new('Epoch',
                     data=[datetime.datetime(1999, 1, 1, i) for i in range(3)],
                     type=spacepy.pycdf.const.CDF_EPOCH)
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 1, 5))
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 1, 4))
        errs = spacepy.pycdf.istp.FileChecks.time_monoton(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Epoch: Nonmonotonic time at record 4.', errs[0])

    def testTimes(self):
        """Compare filename to Epoch times"""
        warnings.filterwarnings(
            'ignore',
            message=r'^No type specified for time input; assuming .*$',
            category=DeprecationWarning,
            module='^spacepy.pycdf$')
        try:
            self.cdf['Epoch'] = [datetime.datetime(1999, 1, 1, i)
                                 for i in range(3)]
        finally:
            del warnings.filters[0]
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 2, 0))
        errs = spacepy.pycdf.istp.FileChecks.times(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Epoch: multiple days 19990101, 19990102.', errs[0])
        del self.cdf['Epoch'][-1]
        errs = spacepy.pycdf.istp.FileChecks.times(self.cdf)
        self.assertEqual(0, len(errs))
        warnings.filterwarnings(
            'ignore',
            message=r'^No type specified for time input; assuming .*$',
            category=DeprecationWarning,
            module='^spacepy.pycdf$')
        try:
            self.cdf['Epoch'] = [datetime.datetime(1999, 1, 2, i)
                                 for i in range(3)]
        finally:
            del warnings.filters[0]
        errs = spacepy.pycdf.istp.FileChecks.times(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Epoch: date 19990102 doesn\'t match file '
                         'source_descriptor_datatype_19990101_v00.cdf.',
                         errs[0])


class FuncTests(ISTPTestsBase):
    """Tests of simple functions"""

    def testFillval(self):
        """Set fill value"""
        expected = ((spacepy.pycdf.const.CDF_EPOCH,
                     datetime.datetime(9999, 12, 31, 23, 59, 59, 999000)),
                    (spacepy.pycdf.const.CDF_BYTE, -128),
                    (spacepy.pycdf.const.CDF_UINT1, 255),
                    (spacepy.pycdf.const.CDF_CHAR, ' '),
                    #Explicit casts to certain size so compare to same precision
                    (spacepy.pycdf.const.CDF_FLOAT, numpy.float32(-1e31)),
                    (spacepy.pycdf.const.CDF_REAL8, numpy.float64(-1e31)),
        )
        for t, e in expected:
            v = self.cdf.new('var', [1, 2,3], type=t)
            spacepy.pycdf.istp.fillval(v)
            self.assertEqual(e, v.attrs['FILLVAL'])
            del self.cdf['var']

    def testFormat(self):
        """Set the format"""
        #This is done by: type, expected, validmin, validmax (None okay)
        expected = ((spacepy.pycdf.const.CDF_EPOCH, 'A24', None, None),
                    (spacepy.pycdf.const.CDF_INT2, 'I2', -2, 9),
                    (spacepy.pycdf.const.CDF_UINT2, 'I5', None, None),
                    (spacepy.pycdf.const.CDF_INT2, 'I6', None, None),
                    (spacepy.pycdf.const.CDF_UINT2, 'I2', 0, 10),
        )
        for t, e, vmin, vmax in expected:
            v = self.cdf.new('var', type=t)
            if vmin is not None:
                v.attrs['VALIDMIN'] = vmin
            if vmin is not None:
                v.attrs['VALIDMAX'] = vmin
            spacepy.pycdf.istp.format(v)
            self.assertEqual(e, v.attrs['FORMAT'])
            del self.cdf['var']

    def testFormatChar(self):
        v = self.cdf.new('var', data=['hi', 'there'])
        spacepy.pycdf.istp.format(v)
        self.assertEqual('A5', v.attrs['FORMAT'])

    def testNanFill(self):
        """Replace fill/invalid values with nan"""
        indata = numpy.array([[5., 99., -1., 3., 4., 12.],
                              [2., 2., 2., 2., 3., -1.]])
        attrs = { 'FILLVAL': 3,
                  'VALIDMIN': 0,
                  'VALIDMAX': 12 }
        var = self.cdf.new('var', data=indata)
        var.attrs = attrs
        data = var.copy()
        expected = numpy.array(
            [[5., numpy.nan, numpy.nan, numpy.nan, 4., 12.],
             [2., 2.,        2.,        2., numpy.nan, numpy.nan]])
        spacepy.pycdf.istp.nanfill(data)
        numpy.testing.assert_almost_equal(
            data, expected, decimal=15)
        #But variable wasn't touched
        numpy.testing.assert_almost_equal(
            var[...], indata, decimal=15)
        #Impose it on the actual variable
        spacepy.pycdf.istp.nanfill(var)
        data = var[...]
        numpy.testing.assert_almost_equal(
            data, expected, decimal=15)
        #And check that integers fail
        var2 = self.cdf.new(
            'var2', data=indata, type=spacepy.pycdf.const.CDF_INT2)
        var2.attrs = attrs
        self.assertRaises(ValueError, spacepy.pycdf.istp.nanfill, var2)


class VarBundleChecksBase(unittest.TestCase):
    """Base class for VarBundle class checks"""
    testfile = 'po_l1_cam_test.cdf'

    def setUp(self):
        """Setup: make an empty, open, writeable CDF"""
        self.tempdir = tempfile.mkdtemp()
        spacepy.pycdf.lib.set_backward(False)
        self.outcdf = spacepy.pycdf.CDF(os.path.join(
            self.tempdir, 'source_descriptor_datatype_19990101_v00.cdf'),
                                     create=True)
        spacepy.pycdf.lib.set_backward(True)
        self.incdf = spacepy.pycdf.CDF(os.path.join(
            spacepy_testing.testsdir, self.testfile))

    def tearDown(self):
        """Close CDFs; delete output"""
        # Suppress did-not-compress warnings on close
        warnings.filterwarnings(
            'ignore', r'^DID_NOT_COMPRESS.*',
            spacepy.pycdf.CDFWarning, r'^spacepy\.pycdf$')
        try:
            self.incdf.close()
            self.outcdf.close()
        except:
            pass
        finally:
            del warnings.filters[0]
        shutil.rmtree(self.tempdir)


class VarBundleChecks(VarBundleChecksBase):
    """Checks for VarBundle class, CAMMICE sample file"""
    testfile = 'po_l1_cam_test.cdf'
    longMessage = True

    def testGetVarInfo(self):
        """Get dependencies, dims, etc. for a variable"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        self.assertEqual(
            [
                'ATC',
                'SectorNumbers',
                'SectorRateScalerNames',
                'SectorRateScalersCounts',
                'SectorRateScalersCountsSigma',
                'SpinNumbers',
            ],
            sorted(bundle._varinfo.keys()))
        self.assertEqual(
            [0], bundle._varinfo['ATC']['dims'])
        self.assertEqual(
            [slice(None)], bundle._varinfo['ATC']['slice'])
        self.assertEqual(
            [0, 1, 2, 3],
            bundle._varinfo['SectorRateScalersCounts']['dims'])
        self.assertEqual(
            [slice(None)] * 4,
            bundle._varinfo['SectorRateScalersCounts']['slice'])
        self.assertEqual(
            [0, 1, 2, 3],
            bundle._varinfo['SectorRateScalersCountsSigma']['dims'])
        self.assertEqual(
            [slice(None)] * 4,
            bundle._varinfo['SectorRateScalersCountsSigma']['slice'])
        self.assertEqual(
            [0, 2],
            bundle._varinfo['SectorNumbers']['dims'])
        self.assertEqual(
            'M', bundle._varinfo['SectorRateScalersCounts']['vartype'])
        self.assertEqual(
            'U', bundle._varinfo['SectorRateScalersCountsSigma']['vartype'])
        self.assertEqual(
            'D', bundle._varinfo['ATC']['vartype'])
        self.assertEqual(
            'D', bundle._varinfo['SectorNumbers']['vartype'])
        for varname, dimension in {
                'ATC': 0,
                'SectorNumbers': 2,
                'SectorRateScalerNames': 3,
                'SectorRateScalersCounts': None,
                'SectorRateScalersCountsSigma': None,
                'SpinNumbers': 1,
                }.items():
            self.assertEqual(
                dimension, bundle._varinfo[varname].get('thisdim', None),
                varname)

    def testOutputSimple(self):
        """Copy a single variable and deps with no slicing"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs,
            self.incdf['SectorRateScalersCounts'].attrs)
        self.assertEqual(
            self.incdf['ATC'].attrs['FILLVAL'],
            self.outcdf['ATC'].attrs['FILLVAL'])

    def testSimpleSlice(self):
        """Slice single element on single dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 2, single=True)
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][:, 2, ...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_2'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('SpinNumbers' in self.outcdf)
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts'].attrs)

    def testSimpleRange(self):
        """Slice a range on single dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 2, single=False)
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][:, 2:, ...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        for d in range(4):
            a = 'DEPEND_{}'.format(d)
            self.assertEqual(
                self.outcdf['SectorRateScalersCounts'].attrs[a],
                self.incdf['SectorRateScalersCounts'].attrs[a])
        numpy.testing.assert_array_equal(
            self.outcdf['SpinNumbers'][:],
            self.incdf['SpinNumbers'][2:])

    def testSliceUndo(self):
        """Slice single element on single dimension, then undo"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 2, single=True).slice(1)
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        for d in (range(4)):
            a = 'DEPEND_{}'.format(d)
            self.assertEqual(
                self.outcdf['SectorRateScalersCounts'].attrs[a],
                self.incdf['SectorRateScalersCounts'].attrs[a])
        numpy.testing.assert_array_equal(
            self.outcdf['SpinNumbers'][:],
            self.incdf['SpinNumbers'][:])

    def testSliceRecord(self):
        """Slice on the record dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(0, 2, 20)
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][2:20, ...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][2:20])
        for d in (range(4)):
            a = 'DEPEND_{}'.format(d)
            self.assertEqual(
                self.outcdf['SectorRateScalersCounts'].attrs[a],
                self.incdf['SectorRateScalersCounts'].attrs[a])
        numpy.testing.assert_array_equal(
            self.outcdf['SpinNumbers'][:],
            self.incdf['SpinNumbers'][:])

    def testSliceRecordStr(self):
        """Slice away record dimension and get str"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(0, 2, 20).mean(0)
        expected = """
        SectorRateScalersCounts: CDF_FLOAT [18, 32, 9] NRV
            SectorRateScalersCountsSigma: CDF_FLOAT [18, 32, 9] NRV
        ATC: CDF_EPOCH16 ---
        SpinNumbers: CDF_CHAR*2 [18] NRV
        SectorNumbers: CDF_CHAR*2 [32] NRV
        SectorRateScalerNames: CDF_CHAR*9 [9] NRV
        """
        expected = inspect.cleandoc(expected).split('\n')
        self.assertEqual(expected, str(bundle).split('\n'))

    def testCAMMICESortOrder(self):
        """More tests of sort order"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        for varname, sortorder in {
                'SectorRateScalersCounts': 0,
                'SectorRateScalersCountsSigma': 2,
                'ATC': 1,
                'SpinNumbers': 1,
                'SectorNumbers': 1,
                'SectorRateScalerNames': 1,
                }.items():
            self.assertEqual(
                sortorder, bundle._varinfo[varname].get('sortorder', None),
                varname)

    def testSliceMultiIDX(self):
        """Slice multiple indices"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, [2, 3, 5])
        with spacepy_testing.assertDoesntWarn(
                self, 'always',
                r'Using a non-tuple sequence for multidimensional indexing',
                FutureWarning, r'spacepy\.pycdf\.istp$'):
            bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][...][:, [2, 3, 5], ...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        for d in (range(4)):
            a = 'DEPEND_{}'.format(d)
            self.assertEqual(
                self.outcdf['SectorRateScalersCounts'].attrs[a],
                self.incdf['SectorRateScalersCounts'].attrs[a])
        numpy.testing.assert_array_equal(
            self.outcdf['SpinNumbers'][:],
            self.incdf['SpinNumbers'][:][[2, 3, 5]])

    def testSliceMultiIDXrecord(self):
        """Slice on the record dimension, multiple index"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(0, [2, 3, 5])
        bundle.output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][...][[2, 3, 5], ...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...][[2, 3, 5]])
        for d in (range(4)):
            a = 'DEPEND_{}'.format(d)
            self.assertEqual(
                self.outcdf['SectorRateScalersCounts'].attrs[a],
                self.incdf['SectorRateScalersCounts'].attrs[a])
        numpy.testing.assert_array_equal(
            self.outcdf['SpinNumbers'][:],
            self.incdf['SpinNumbers'][:])

    def testSum(self):
        """Sum over a dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.sum(2)
        self.assertEqual([False, False, True, False], bundle._summed)
        bundle.output(self.outcdf)
        counts = self.incdf['SectorRateScalersCounts'][...]
        expected = counts.sum(axis=2)
        expected[(counts < 0).max(axis=2)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCounts'][...])
        sigma = self.incdf['SectorRateScalersCountsSigma'][...]
        bad = (sigma < 0)
        sigma[bad] = 0 #avoid warning
        expected = numpy.sqrt((sigma ** 2).sum(axis=2))
        expected[bad.max(axis=2)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCountsSigma'][...])
        self.assertFalse('SectorNumbers' in self.outcdf)
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_1'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts'].attrs)
        self.assertEqual(
            self.incdf['ATC'].attrs['FILLVAL'],
            self.outcdf['ATC'].attrs['FILLVAL'])

    def testSliceSum(self):
        """Slice and sum over a dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 0, 16).sum(1)
        self.assertEqual([False, True, False, False], bundle._summed)
        bundle.output(self.outcdf)
        counts = self.incdf['SectorRateScalersCounts'][:, 0:16, ...]
        expected = counts.sum(axis=1)
        expected[(counts < 0).max(axis=1)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCounts'][...])
        sigma = self.incdf['SectorRateScalersCountsSigma'][:, 0:16, ...]
        bad = (sigma < 0)
        sigma[bad] = 0 #avoid warning
        expected = numpy.sqrt((sigma ** 2).sum(axis=1))
        expected[bad.max(axis=1)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCountsSigma'][...])
        self.assertFalse('SpinNumbers' in self.outcdf)
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_2'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts'].attrs)

    def testMean(self):
        """Average over a dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.mean(2)
        self.assertEqual([False, False, False, False], bundle._summed)
        self.assertEqual([False, False, True, False], bundle._mean)
        bundle.output(self.outcdf)
        counts = self.incdf['SectorRateScalersCounts'][...]
        counts[counts < 0] = numpy.nan
        #suppress bad value warnings
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', 'Mean of empty slice', RuntimeWarning)
            expected = numpy.nanmean(counts, axis=2)
        expected[numpy.isnan(expected)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCounts'][...])
        sigma = self.incdf['SectorRateScalersCountsSigma'][...]
        sigma[sigma < 0] = numpy.nan
        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore', r'invalid value encountered in (?:true_)divide$',
                RuntimeWarning)
            expected = numpy.sqrt(numpy.nansum(sigma ** 2, axis=2)) \
                        / (~numpy.isnan(sigma)).sum(axis=2)
        expected[numpy.isnan(expected)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCountsSigma'][...])
        self.assertFalse('SectorNumbers' in self.outcdf)
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_1'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts'].attrs)

    def testNonconflictingMultiple(self):
        """Put multiple variables without conflict in output"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.sum(1) #Sum over spin
        bundle.output(self.outcdf)
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SpinRateScalersCounts'])
        #This has some overlapping deps, but they're all the same
        bundle.output(self.outcdf)
        self.assertTrue('SpinNumbers' in self.outcdf)
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_2'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts'].attrs)
        for i in range(3):
            d = 'DEPEND_{}'.format(i)
            self.assertEqual(
                self.outcdf['SpinRateScalersCounts'].attrs[d],
                self.incdf['SpinRateScalersCounts'].attrs[d])

    def testConflictingMultiple(self):
        """Put multiple variables with conflict in output"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 0, 3)
        bundle.output(self.outcdf)
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SpinRateScalersCounts'])
        bundle.slice(1, 2, 5)
        #This is a different slice on same dim, should fail
        msg = 'Incompatible SpinNumbers already exists in output.'
        try:
            bundle.output(self.outcdf)
        except RuntimeError:
            self.assertEqual(msg, str(sys.exc_info()[1]))
        else:
            self.fail('Should have raised RuntimeError: ' + msg)

    def testNameMap(self):
        """Test name mapping"""
        #Essentially a subtest of below
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.sum(2)
        namemap = bundle._namemap(suffix="_Summed")
        expected = { n: n + '_Summed' for n in [
            'SectorRateScalersCounts', #Main var
            'SectorRateScalersCountsSigma', #its delta
            #The sum dimension. Even though it is summed completely away,
            #the check for going away is elsewhere (and if it existed, it
            #would be renamed)
            'SectorNumbers',
            #no other variables change
        ] }
        self.assertEqual(expected, namemap)

    def testSumRename(self):
        """Sum over a dimension, rename output"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.sum(2)
        bundle.output(self.outcdf, suffix='_Summed')
        counts = self.incdf['SectorRateScalersCounts'][...]
        expected = counts.sum(axis=2)
        expected[(counts < 0).max(axis=2)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCounts_Summed'][...])
        sigma = self.incdf['SectorRateScalersCountsSigma'][...]
        bad = (sigma < 0)
        sigma[bad] = 0 #avoid warning
        expected = numpy.sqrt((sigma ** 2).sum(axis=2))
        expected[bad.max(axis=2)] = -1e31
        numpy.testing.assert_allclose(
            expected, self.outcdf['SectorRateScalersCountsSigma_Summed'][...])
        self.assertFalse('SectorNumbers' in self.outcdf)
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts_Summed'].attrs['DEPEND_2'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_3'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts_Summed'].attrs['DEPEND_1'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_1'])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts_Summed'].attrs['DEPEND_0'],
            self.incdf['SectorRateScalersCounts'].attrs['DEPEND_0'])
        self.assertFalse('DEPEND_3'
                        in self.outcdf['SectorRateScalersCounts_Summed'].attrs)
        self.assertEqual('SectorRateScalersCountsSigma_Summed',
                         self.outcdf['SectorRateScalersCounts_Summed']
                         .attrs['DELTA_PLUS_VAR'])

    def testSumRenameConflict(self):
        """Sum over a dimension, rename output, with a potential conflict"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(2, 0, 2)
        bundle.output(self.outcdf, suffix='_0-1')
        bundle.slice(2, 2, 4)
        bundle.output(self.outcdf, suffix='_2-3')
        counts = self.incdf['SectorRateScalersCounts'][...]
        numpy.testing.assert_allclose(
            counts[..., 0:2, :],
            self.outcdf['SectorRateScalersCounts_0-1'][...])
        sigma = self.incdf['SectorRateScalersCountsSigma'][...]
        numpy.testing.assert_allclose(
            sigma[..., 0:2, :],
            self.outcdf['SectorRateScalersCountsSigma_0-1'][...])
        self.assertFalse('SectorNumbers' in self.outcdf)
        self.assertTrue('SectorRateScalerNames' in self.outcdf)
        self.assertFalse('SectorRateScalerNames_0-1' in self.outcdf)
        self.assertTrue('SectorNumbers_0-1' in self.outcdf)
        self.assertTrue('SectorNumbers_2-3' in self.outcdf)
        #Most depends are the same
        for which in ('0-1', '2-3'):
            for d in range(0, 4):
                if d == 2:
                    continue
                self.assertEqual(
                    self.outcdf['SectorRateScalersCounts_{}'.format(which)]
                    .attrs['DEPEND_{}'.format(d)],
                    self.incdf['SectorRateScalersCounts']
                    .attrs['DEPEND_{}'.format(d)])
        #But dim 2 is different
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts_2-3'].attrs['DEPEND_2'],
            'SectorNumbers_2-3')
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts_0-1'].attrs['DEPEND_2'],
            'SectorNumbers_0-1')
        self.assertEqual(
            'SectorRateScalersCounts_0-1',
            self.outcdf['SectorRateScalersCounts_0-1'].attrs['FIELDNAM'])
        self.assertEqual(
            'SectorNumbers_2-3',
            self.outcdf['SectorNumbers_2-3'].attrs['FIELDNAM'])


class VarBundleChecksHOPE(VarBundleChecksBase):
    """Checks for VarBundle class, HOPE sample file"""
    testfile = os.path.join('data',
                            'rbspa_rel04_ect-hope-PA-L3_20121201_v0.0.0.cdf')
    longMessage = True

    def tearDown(self):
        """Block warnings from CDF closing"""
        warnings.filterwarnings(
            'ignore', message='^DID_NOT_COMPRESS.*$',
            category=spacepy.pycdf.CDFWarning,
            module='^spacepy.pycdf')
        try:
            super(VarBundleChecksHOPE, self).tearDown()
        finally:
            del warnings.filters[0]

    def testSortOrder(self):
        """Check sort order of variables"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['Counts_P'])
        for varname, sortorder in {
                'Counts_P': 0,
                'Epoch_Ion': 1,
                'Epoch_Ion_DELTA': 2,
                'PITCH_ANGLE': 1,
                'Pitch_LABL': 3,
                'HOPE_ENERGY_Ion': 1,
                'ENERGY_Ion_DELTA': 2,
                'Energy_LABL': 3,
                }.items():
            self.assertEqual(
                sortorder, bundle._varinfo[varname].get('sortorder', None),
                varname)

    def testDepWithDelta(self):
        """Properly handle a dependency with a delta"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['Counts_P'])
        self.assertEqual('M', bundle._varinfo['Counts_P']['vartype'])
        self.assertEqual('D', bundle._varinfo['ENERGY_Ion_DELTA']['vartype'])
        bundle.slice(2, 0, 10).mean(2).output(self.outcdf)
        expected = self.incdf['Counts_P'][:, :, 0:10, ...]
        expected[expected < 0] = numpy.nan
        expected = numpy.nanmean(expected, axis=2)
        expected[numpy.isnan(expected)] = -1e31
        numpy.testing.assert_allclose(
            self.outcdf['Counts_P'], expected)

    def testOperations(self):
        """Get operations of a bundle"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['Counts_P'])
        bundle.slice(1, 1, single=True).slice(2, 0, 10).mean(2)
        ops = bundle.operations()
        self.assertEqual(
            [('slice', (1, 1), {'single': True}),
             ('slice', (2, 0, 10), {}),
             ('mean', (2,), {})],
            ops)
        #Check fancy index
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['Counts_P'])
        bundle.slice(1, [5, 6])
        ops = bundle.operations()
        self.assertEqual([('slice', (1, [5, 6]), {})], ops)

    def testSumRecord(self):
        """Sum on the record dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['Counts_P'])
        bundle.sum(0).output(self.outcdf)
        expected = numpy.sum(self.incdf['Counts_P'][...], axis=0)
        numpy.testing.assert_array_equal(
            self.outcdf['Counts_P'][...], expected)
        self.assertFalse(self.outcdf['Counts_P'].rv())
        self.assertFalse('DEPEND_0' in self.outcdf['Counts_P'].attrs)
        self.assertFalse('Epoch' in self.outcdf)
        self.assertEqual(
            'PITCH_ANGLE', self.outcdf['Counts_P'].attrs['DEPEND_1'])

    def testAvgRecord(self):
        """Average on the record dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['Counts_P'])
        bundle.mean(0).output(self.outcdf)
        expected = numpy.mean(self.incdf['Counts_P'][...], axis=0)
        numpy.testing.assert_array_equal(
            self.outcdf['Counts_P'][...], expected)
        self.assertFalse(self.outcdf['Counts_P'].rv())
        self.assertFalse('Epoch' in self.outcdf)
        self.assertFalse('DEPEND_0' in self.outcdf['Counts_P'].attrs)
        self.assertEqual(
            'PITCH_ANGLE', self.outcdf['Counts_P'].attrs['DEPEND_1'])

    def testSliceSingleRecord(self):
        """Slice single element on the record dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['Counts_P'])
        bundle.slice(0, 0, single=True).output(self.outcdf)
        expected = self.incdf['Counts_P'][0, ...]
        numpy.testing.assert_array_equal(
            self.outcdf['Counts_P'][...], expected)
        self.assertFalse(self.outcdf['Counts_P'].rv())
        self.assertFalse('DEPEND_0' in self.outcdf['Counts_P'].attrs)
        self.assertFalse('Epoch' in self.outcdf)
        self.assertEqual(
            'PITCH_ANGLE', self.outcdf['Counts_P'].attrs['DEPEND_1'])

    def testVars(self):
        """Get variables of a bundle"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['FPDU'])
        bundle.slice(1, 1, single=True).slice(2, 0, 10)
        variables = bundle.variables()
        self.assertEqual([
            [('FPDU', (100, 10))],
            [('Epoch_Ion', (100,)), ('Epoch_Ion_DELTA', (100,))],
            [('PITCH_ANGLE', None), ('Pitch_LABL', None)],
            [('HOPE_ENERGY_Ion', (100, 10)), ('ENERGY_Ion_DELTA', (100, 10)),
             ('Energy_LABL', (10,))]],
            variables)

    def testStrRepr(self):
        """Get string representation of bundle"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['FPDU'])
        bundle.slice(1, 1, single=True).slice(2, 0, 10)
        expected = """
        FPDU: CDF_FLOAT [100, 10]
        Epoch_Ion: CDF_EPOCH [100]
            Epoch_Ion_DELTA: CDF_REAL4 [100]
        PITCH_ANGLE: CDF_FLOAT ---
            Pitch_LABL: CDF_CHAR*5 ---
        HOPE_ENERGY_Ion: CDF_FLOAT [100, 10]
            ENERGY_Ion_DELTA: CDF_FLOAT [100, 10]
            Energy_LABL: CDF_CHAR*3 [10] NRV
        """
        expected = inspect.cleandoc(expected).split('\n')
        #Split on linebreak to get a better diff
        self.assertEqual(expected, str(bundle).split('\n'))
        self.assertEqual(['<VarBundle:'] + expected + ['>'],
                         repr(bundle).split('\n'))

    def testOutshape(self):
        """Get the output shape of variables"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['FPDU'])
        bundle.slice(1, 1, single=True).slice(2, 0, 10)
        expected = {
            'FPDU': (100, 10),
            'Epoch_Ion': (100,),
            'Epoch_Ion_DELTA': (100,),
            'PITCH_ANGLE': None,
            'HOPE_ENERGY_Ion': (100, 10),
            'ENERGY_Ion_DELTA': (100, 10),
            }
        for vname, shape in expected.items():
            self.assertEqual(
                shape, bundle._outshape(vname), vname)

    def testSliceNRVScalar(self):
        """Slice when the EPOCH_DELTA is NRV"""
        #Modify the input first
        self.incdf.close()
        newtest = os.path.join(self.tempdir, os.path.basename(self.testfile))
        shutil.copy2(self.testfile, newtest)
        with spacepy.pycdf.CDF(newtest, readonly=False) as cdf:
            delta = cdf['Epoch_Ion_DELTA']
            newdelta = cdf.new(
                'Epoch_Ion_DELTA_new', data=delta[0],
                type=delta.type(), recVary=False)
            newdelta.attrs.clone(delta.attrs)
            del cdf['Epoch_Ion_DELTA']
            newdelta.rename('Epoch_Ion_DELTA')
        self.incdf = spacepy.pycdf.CDF(newtest)
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['FPDU'])
        bundle.slice(0, 0, 10).output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['FPDU'][...], self.incdf['FPDU'][0:10, ...])
        numpy.testing.assert_array_equal(
            self.outcdf['Epoch_Ion'][...], self.incdf['Epoch_Ion'][0:10, ...])
        numpy.testing.assert_array_equal(
            self.outcdf['Epoch_Ion_DELTA'][...],
            self.incdf['Epoch_Ion_DELTA'][...])

    def testSliceNoRecords(self):
        """Slice when there are no records on the input"""
        #Modify the input first
        self.incdf.close()
        newtest = os.path.join(self.tempdir, os.path.basename(self.testfile))
        shutil.copy2(self.testfile, newtest)
        with spacepy.pycdf.CDF(newtest, readonly=False) as cdf:
            del cdf['FPDU'][...] #Delete data not variable
        self.incdf = spacepy.pycdf.CDF(newtest)
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['FPDU'])
        bundle.sum(1).slice(2, 0, 6).output(self.outcdf)
        self.assertEqual(
            (0, 6), self.outcdf['FPDU'].shape)


class VarBundleChecksEPILo(VarBundleChecksBase):
    """Checks for VarBundle class, EPILo sample file"""
    testfile = os.path.join('data',
                            'psp_isois-epilo_l2-ic_20190401_v0.0.0.cdf')

    def tearDown(self):
        """Block warnings from CDF closing"""
        warnings.filterwarnings(
            'ignore', message='^DID_NOT_COMPRESS.*$',
            category=spacepy.pycdf.CDFWarning,
            module='^spacepy.pycdf')
        try:
            super(VarBundleChecksEPILo, self).tearDown()
        finally:
            del warnings.filters[0]

    def testDoubleDep(self):
        """Handle a variable with a 2D depend"""
        countrate = self.incdf['H_CountRate_ChanT']
        bundle = spacepy.pycdf.istp.VarBundle(countrate)
        bundle.sum(1).slice(2, 0, 10).output(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['H_CountRate_ChanT'],
            countrate[:, :, 0:10].sum(axis=1))
        #Look direction should go away
        for v in ('Look_80_LABL', 'Look_Direction_80',
                  'Look_Direction_80_DELTAMINUS',
                  'Look_Direction_80_DELTAPLUS'):
            self.assertFalse(v in self.outcdf)

    def testDoubleDepSummed(self):
        """Handle a variable with a 2D depend, sum all dims"""
        countrate = self.incdf['H_CountRate_ChanT']
        bundle = spacepy.pycdf.istp.VarBundle(countrate)
        bundle.sum(1).slice(2, 0, 10).sum(2).output(self.outcdf, '_TS')
        numpy.testing.assert_array_equal(
            self.outcdf['H_CountRate_ChanT_TS'],
            countrate[:, :, 0:10].sum(axis=2).sum(axis=1))
        #Look direction and energy should go away
        for v in ('Look_80_LABL', 'Look_Direction_80',
                  'Look_Direction_80_DELTAMINUS',
                  'Look_Direction_80_DELTAPLUS',
                  'H_ChanT_Energy', 'H_ChanT_Energy_LABL',
                  'H_ChanT_Energy_DELTAMINUS', 'H_ChanT_Energy_DELTAPLUS'):
            self.assertFalse(v in self.outcdf)
            self.assertFalse(v + '_TS' in self.outcdf)

    def testConflictingEpoch(self):
        """Regression test for complicated name conflict"""
        bundle = spacepy.pycdf.istp.VarBundle(self.incdf['H_CountRate_ChanT'])
        bundle.sum(1).slice(2, 1).output(self.outcdf, suffix='_SP')
        #Still summed on 1!
        bundle.slice(2, 18, 32).sum(2).output(self.outcdf, suffix='_TS')
        self.assertIn('H_CountRate_ChanT_SP', self.outcdf)
        self.assertIn('H_CountRate_ChanT_TS', self.outcdf)


if __name__ == '__main__':
    unittest.main()
