#!/usr/bin/env python

"""Unit testing for istp support"""

import datetime
import os.path
import shutil
import tempfile
import unittest

import numpy
import numpy.testing
import spacepy.pycdf
import spacepy.pycdf.const
import spacepy.pycdf.istp


class ISTPTestsBase(unittest.TestCase):
    """Base class for ISTP check functions, provide a blank CDF to modify"""

    def setUp(self):
        """Setup: make an empty, open, writeable CDF"""
        self.tempdir = tempfile.mkdtemp()
        self.cdf = spacepy.pycdf.CDF(os.path.join(
            self.tempdir, 'source_descriptor_datatype_19990101_v00.cdf'),
                                     create=True)

    def tearDown(self):
        """Delete the empty cdf"""
        try:
            self.cdf.close()
        except:
            pass
        shutil.rmtree(self.tempdir)


class VariablesTests(ISTPTestsBase):
    """Tests of variable-checking functions"""

    def testDepends(self):
        """Dependencies, valid and not"""
        self.cdf['var1'] = [1, 2, 3]
        self.cdf['var1'].attrs['DEPEND_0'] = 'var2'
        errs = spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual('DEPEND_0 variable var2 missing', errs[0])
        self.cdf['var2'] = [1, 2, 3]
        self.assertEqual(
            0, len(spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])))
        self.cdf['var1'].attrs['DELTA_PLUS_VAR'] = 'foobar'
        errs = spacepy.pycdf.istp.VariableChecks.depends(self.cdf['var1'])
        self.assertEqual(1, len(errs))
        self.assertEqual('DELTA_PLUS_VAR variable foobar missing', errs[0])

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
        #This is a bit contrived to force a difference that's only
        #in terms of the precision of the float
        v = self.cdf.new('var1', recVary=False, data=[1, 2, 3],
                         type=spacepy.pycdf.const.CDF_DOUBLE)
        v.attrs['VALIDMIN'] = 0
        v.attrs['VALIDMAX'] = 10
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

    def testValidRangeFillvalDatetime(self):
        """Validmin/validmax with fillval set, Epoch var"""
        v = self.cdf.new(
            'var1', data=[datetime.datetime(2010, 1, i) for i in range(1, 6)])
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
            ['SCALEMIN (-200) outside data range (-128,127).',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMIN'] = 200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(3, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMIN (200) outside data range (-128,127).',
             'SCALEMIN > SCALEMAX.',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMAX'] = -200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(5, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMAX (-200) outside data range (-128,127).',
             'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.',
             'SCALEMIN (200) outside data range (-128,127).',
             'SCALEMIN > SCALEMAX.',
             'SCALEMIN type CDF_INT2 does not match variable type CDF_BYTE.'
             ],
             errs)
        v.attrs['SCALEMAX'] = 200
        errs = spacepy.pycdf.istp.VariableChecks.validscale(v)
        self.assertEqual(4, len(errs))
        errs.sort()
        self.assertEqual(
            ['SCALEMAX (200) outside data range (-128,127).',
             'SCALEMAX type CDF_INT2 does not match variable type CDF_BYTE.',
             'SCALEMIN (200) outside data range (-128,127).',
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
            'SCALEMAX ([300 320]) outside data range (-128,127).',
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
        self.assertEqual(1, len(errs))
        self.assertEqual('var1: DEPEND_0 variable var2 missing', errs[0])

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
        self.cdf['Epoch'] = [datetime.datetime(1999, 1, 1, i) for i in range(3)]
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 1, 5))
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 1, 4))
        errs = spacepy.pycdf.istp.FileChecks.time_monoton(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Epoch: Nonmonotonic time at record 4.', errs[0])

    def testTimes(self):
        """Compare filename to Epoch times"""
        self.cdf['Epoch'] = [datetime.datetime(1999, 1, 1, i) for i in range(3)]
        self.cdf['Epoch'].append(datetime.datetime(1999, 1, 2, 0))
        errs = spacepy.pycdf.istp.FileChecks.times(self.cdf)
        self.assertEqual(1, len(errs))
        self.assertEqual('Epoch: multiple days 19990101, 19990102.', errs[0])
        del self.cdf['Epoch'][-1]
        errs = spacepy.pycdf.istp.FileChecks.times(self.cdf)
        self.assertEqual(0, len(errs))
        self.cdf['Epoch'] = [datetime.datetime(1999, 1, 2, i) for i in range(3)]
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


class VarBundleChecks(unittest.TestCase):
    """Checks for VarBundle class"""

    def setUp(self):
        """Setup: make an empty, open, writeable CDF"""
        self.tempdir = tempfile.mkdtemp()
        spacepy.pycdf.lib.set_backward(False)
        self.outcdf = spacepy.pycdf.CDF(os.path.join(
            self.tempdir, 'source_descriptor_datatype_19990101_v00.cdf'),
                                     create=True)
        spacepy.pycdf.lib.set_backward(True)
        pth = os.path.dirname(os.path.abspath(__file__))
        self.incdf = spacepy.pycdf.CDF(os.path.join(pth, 'po_l1_cam_test.cdf'))

    def tearDown(self):
        """Close CDFs; delete output"""
        try:
            self.incdf.close()
            self.outcdf.close()
        except:
            pass
        shutil.rmtree(self.tempdir)

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
            'S', bundle._varinfo['SectorRateScalersCounts']['sumtype'])
        self.assertEqual(
            'Q', bundle._varinfo['SectorRateScalersCountsSigma']['sumtype'])
        self.assertEqual(
            'F', bundle._varinfo['ATC']['sumtype'])
        self.assertEqual(
            'F', bundle._varinfo['SectorNumbers']['sumtype'])

    def testWriteSimple(self):
        """Copy a single variable and deps with no slicing"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.write(self.outcdf)
        numpy.testing.assert_array_equal(
            self.outcdf['SectorRateScalersCounts'][...],
            self.incdf['SectorRateScalersCounts'][...])
        numpy.testing.assert_array_equal(
            self.outcdf['ATC'][...],
            self.incdf['ATC'][...])
        self.assertEqual(
            self.outcdf['SectorRateScalersCounts'].attrs,
            self.incdf['SectorRateScalersCounts'].attrs)

    def testSimpleSlice(self):
        """Slice single element on single dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 2, single=True)
        bundle.write(self.outcdf)
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
        bundle.write(self.outcdf)
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
        bundle.write(self.outcdf)
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
        bundle.write(self.outcdf)
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

    def testSliceMultiIDX(self):
        """Slice multiple indices"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, [2, 3, 5])
        bundle.write(self.outcdf)
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
        bundle.write(self.outcdf)
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
        bundle.write(self.outcdf)
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

    def testSliceSum(self):
        """Slice and sum over a dimension"""
        bundle = spacepy.pycdf.istp.VarBundle(
            self.incdf['SectorRateScalersCounts'])
        bundle.slice(1, 0, 16).sum(1)
        self.assertEqual([False, True, False, False], bundle._summed)
        bundle.write(self.outcdf)
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


if __name__ == '__main__':
    unittest.main()
