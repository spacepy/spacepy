#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's datamodel

Copyright 2010-2012 Los Alamos National Security, LLC.
"""


from __future__ import division

import datetime
import os
import tempfile
import unittest
import StringIO
import sys
import warnings

import spacepy.datamodel as dm
from spacepy import pycdf
import numpy as np

try:
    import cPickle as pickle
except:
    import pickle

__all__ = ['SpaceDataTests', 'dmarrayTests', 'converterTests', 'JSONTests']


class SpaceDataTests(unittest.TestCase):
    def setUp(self):
        super(SpaceDataTests, self).setUp()

    def tearDown(self):
        super(SpaceDataTests, self).tearDown()

    def test_dmcopy(self):
        """dmcopy should copy datamodel objects"""
        a = dm.SpaceData()
        a[1] = dm.dmarray([1,2,3], attrs={1:1})
        b = dm.dmcopy(a)
        self.assertFalse(a is b) # they are not the same memory
        np.testing.assert_allclose(a[1], b[1])
        self.assertEqual(a[1].attrs, b[1].attrs)
        b = dm.dmcopy(a[1])
        np.testing.assert_allclose(a[1], b)
        self.assertEqual(a[1].attrs, b.attrs)
        a = np.arange(10)
        b = dm.dmcopy(a)
        np.testing.assert_allclose(a, b)
        a = [1,2,3]
        b = dm.dmcopy(a)
        self.assertEqual(a, b)

    def test_SpaceData(self):
        """Spacedata dist object has certain attributes"""
        dat = dm.SpaceData()
        self.assertEqual(dat.attrs, {})
        dat = dm.SpaceData(attrs={'foo':'bar'})
        self.assertEqual(dat.attrs['foo'], 'bar')

    def test_flatten_function(self):
        """Flatten should flatted a nested dict"""
        a = dm.SpaceData()
        a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        a['4'] = dm.SpaceData(cat = 'kitty')
        a['5'] = 4
        self.assertEqual(a['1']['dog'], 5)
        b = dm.flatten(a)
        try:
            b['1']['dog']
        except KeyError:
            pass
        else:
            self.fail('KeyError not raised')
        # might be possible that list order is not preserved and this fails,
        # if so change to a bunch of self.assertTrue and in statements
        self.assertEqual(b.keys(), ['1<--pig<--fish<--a', '4<--cat', '1<--dog', '1<--pig<--fish<--b', '5'])

    def test_flatten_method(self):
        """Flatten should flatted a nested dict"""
        a = dm.SpaceData()
        a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        a['4'] = dm.SpaceData(cat = 'kitty')
        a['5'] = 4
        self.assertEqual(a['1']['dog'], 5)
        a.flatten()
        try:
            a['1']['dog']
        except KeyError:
            pass
        else:
            self.fail('KeyError not raised')
        ans =  ['4<--cat', '1<--dog', '5', '1<--pig<--fish<--a', '1<--pig<--fish<--b']
        ans.sort()
        val = a.keys()
        val.sort()
        self.assertEqual(val, ans)

    def test_numeric_key(self):
        """flatten should handle a numeric key"""
        a = dm.SpaceData()
        a[1] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        a[4] = dm.SpaceData(cat = 'kitty')
        a[5] = 4
        self.assertEqual(a[1]['dog'], 5)
        a.flatten()
        try:
            a[1]['dog']
        except KeyError:
            pass
        else:
            self.fail('KeyError not raised')
        ans = ['4<--cat', '1<--dog', 5, '1<--pig<--fish<--a', '1<--pig<--fish<--b']
        ans.sort()
        val = a.keys()
        val.sort()
        self.assertEqual(val, ans)

    def test_tree(self):
        """.tree() should call dictree"""
        a = dm.SpaceData()
        a['foo'] = dm.dmarray([1,2,3])
        realstdout = sys.stdout
        output = StringIO.StringIO()
        sys.stdout = output
        self.assertEqual(a.tree(), None)
        sys.stdout = realstdout
        result = output.getvalue()
        output.close()
        expected = "+\n|____foo\n"
        self.assertEqual(result, expected)


class dmarrayTests(unittest.TestCase):
    def setUp(self):
        super(dmarrayTests, self).setUp()
        self.dat = dm.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})

    def tearDown(self):
        super(dmarrayTests, self).tearDown()
        del self.dat

    def test_append(self):
        """append should maintain all Allowed_Attributes"""
        d2 = dm.dmarray.append(self.dat, -1)
        np.testing.assert_array_equal([1,2,3, 4, -1], d2)

    def test_count(self):
        """count should work like on a list"""
        self.assertEqual(1, self.dat.count(1))
        self.assertEqual(0, self.dat.count(10))
        self.assertEqual(3, dm.dmarray([1,1,1, 3, 4, 5, 4]).count(1))
        self.assertEqual(2, dm.dmarray([1,1,1, 3, 4, 5, 4]).count(4))

    def test_creation_dmarray(self):
        """When a dmarray is created it should have attrs empty or not"""
        self.assertTrue(hasattr(self.dat, 'attrs'))
        self.assertEqual(self.dat.attrs['a'], 'a')
        data = dm.dmarray([1,2,3])
        self.assertTrue(hasattr(data, 'attrs'))
        self.assertEqual(data.attrs, {})
        data2 = dm.dmarray([1,2,3], attrs={'coord':'GSM'})
        self.assertEqual(data.attrs, {})
        self.assertEqual(data2.attrs, {'coord':'GSM'})
        data2 = dm.dmarray([1,2,3], dtype=float, attrs={'coord':'GSM'})
        np.testing.assert_allclose([1,2,3], data2)

    def test_different_attrs(self):
        """Different instances of dmarray shouldn't share attrs"""
        a = dm.dmarray([1, 2, 3, 4])
        b = dm.dmarray([2, 3, 4, 5])
        a.attrs['hi'] = 'there'
        self.assertNotEqual(a.attrs, b.attrs)

    def test_slicing(self):
        '''Slicing a dmarray should keep the attrs'''
        dat_sl = self.dat[:-1]
        self.assertTrue(hasattr(dat_sl, 'attrs'))
        self.assertEqual(self.dat.attrs, dat_sl.attrs)
        #make sure the attrs aren't pointing at the same obj
        dat_sl.attrs = {'foo': 'bar'}
        self.assertNotEqual(self.dat.attrs, dat_sl.attrs)

    def test_pickle_dumps(self):
        """things should pickle and unpickle"""
        tmp = pickle.dumps(self.dat)
        for i, val in enumerate(self.dat):
            self.assertEqual(pickle.loads(tmp)[i], val)
        self.assertEqual(pickle.loads(tmp).attrs, self.dat.attrs)

    def test_pickle_dump(self):
        """things should pickle and unpickle to a file"""
        fname = None
        try:
            with tempfile.NamedTemporaryFile(delete=False) as fp:
                fname = fp.name
                pickle.dump(self.dat, fp)
            with open(fname, 'rb') as fp:
                dat2 = pickle.load(fp)
        finally:
            if fname != None:
                os.remove(fname)
        np.testing.assert_allclose(self.dat, dat2)
        self.assertEqual(self.dat.attrs, dat2.attrs)

    def test_attrs_only(self):
        """dmarray can only have .attrs"""
        self.assertRaises(TypeError, dm.dmarray, [1,2,3], setme = 123 )

    def test_more_attrs(self):
        """more attrs are allowed if they are predefined"""
        a = dm.dmarray([1,2,3])
        a.Allowed_Attributes = a.Allowed_Attributes + ['blabla']
        a.blabla = {}
        a.blabla['foo'] = 'you'
        self.assertEqual(a.blabla['foo'], 'you')

    def test_extra_pickle(self):
        """Extra attrs are pickled and unpicked"""
        self.dat.addAttribute('blabla', {'foo':'you'})
        val = pickle.dumps(self.dat)
        b = pickle.loads(val)
        self.assertTrue('blabla' in b.Allowed_Attributes)
        self.assertEqual(b.blabla['foo'], 'you')

    def test_extra_pickle2(self):
        """Order should not matter of Allowed_Attributes"""
        # added new one to the front
        self.dat.Allowed_Attributes = ['foo'] + self.dat.Allowed_Attributes
        self.dat.foo = 'bar'
        val = pickle.dumps(self.dat)
        b = pickle.loads(val)
        self.assertTrue('foo' in b.Allowed_Attributes)
        self.assertEqual(b.foo, 'bar')

    def test_addAttribute(self):
        """addAttribute should work"""
        a = dm.dmarray([1,2,3])
        a.addAttribute('bla')
        self.assertEqual(a.bla, None)
        a.addAttribute('bla2', {'foo': 'bar'})
        self.assertEqual(a.bla2['foo'], 'bar')
        self.assertRaises(NameError, a.addAttribute, 'bla2')

    def test_attrs(self):
        """The only attribute the can be set is attrs"""
        self.assertRaises(TypeError, dm.dmarray, [1,2,3], bbb=23)
        try:
            self.dat.bbb = 'someval'
        except TypeError:
            pass
        else:
            self.fail(
                'Assigning to arbitrary Python attribute should raise TypeError')

class converterTests(unittest.TestCase):
    def setUp(self):
        super(converterTests, self).setUp()
        self.SDobj = dm.SpaceData(attrs={'global': 'test'})
        self.SDobj['var'] = dm.dmarray([1, 2, 3], attrs={'a': 'a'})
        self.testfile = tempfile.mkstemp()
        warnings.simplefilter('error', dm.DMWarning)

    def tearDown(self):
        super(converterTests, self).tearDown()
        del self.SDobj
        os.remove(self.testfile[1])
        warnings.simplefilter('default', dm.DMWarning)

    def test_convertKeysToStr(self):
        """convertKeysToStr should give known output"""
        a = dm.SpaceData()
        a['data'] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual(a.keys(), b.keys())
        a = dm.SpaceData()
        a[50] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(a.keys()[0])], b.keys())
        a = {}
        a[50] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(a.keys()[0])], b.keys())
        a = dm.SpaceData()
        a['data'] = dm.SpaceData()
        a['data']['test'] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual(a.keys(), b.keys())
        a = dm.SpaceData()
        a[50] = dm.SpaceData()
        a[50][49] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(a.keys()[0])], b.keys())

    def test_HDF5roundtrip(self):
        """Data can go to hdf and back"""
        dm.toHDF5(self.testfile[1], self.SDobj)
        newobj = dm.fromHDF5(self.testfile[1])
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_allclose(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])
        dm.toHDF5(self.testfile[1], self.SDobj, mode='a')
        newobj = dm.fromHDF5(self.testfile[1])
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_allclose(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])

    def test_HDF5Exceptions(self):
        """HDF5 has warnings and exceptions"""
        dm.toHDF5(self.testfile[1], self.SDobj)
        self.assertRaises(IOError, dm.toHDF5, self.testfile[1], self.SDobj, overwrite=False)
        a = dm.SpaceData()
        a['foo'] = 'bar' # not an allowed type for data
        self.assertRaises(dm.DMWarning, dm.toHDF5, self.testfile[1], a)

    def test_HDF5roundtrip2(self):
        """Data can go to hdf without altering datetimes in the datamodel"""
        a = dm.SpaceData()
        a['foo'] = dm.SpaceData()
        dm.toHDF5(self.testfile[1], a)
        newobj = dm.fromHDF5(self.testfile[1])
        self.assertEqual(a['foo'], newobj['foo'])
        a['bar'] = dm.dmarray([datetime.datetime(2000, 1, 1)])
        dm.toHDF5(self.testfile[1], a)
        self.assertEqual(a['bar'], dm.dmarray([datetime.datetime(2000, 1, 1)]))

    def test_dateToISO(self):
        """dateToISO should recurse properly"""
        d1 = {'k1':datetime.datetime(2012,12,21)}
        self.assertEqual({'k1': '2012-12-21T00:00:00'}, dm._dateToISO(d1))
        d1 = {'k1':{'k2':datetime.datetime(2012,12,21)}}
        # regression, it does not traverse nested dicts
        self.assertEqual({'k1': {'k2': datetime.datetime(2012, 12, 21, 0, 0)}}, dm._dateToISO(d1))
        d1 = {'k1':[datetime.datetime(2012,12,21), datetime.datetime(2012,12,22)] }
        self.assertEqual({'k1': ['2012-12-21T00:00:00', '2012-12-22T00:00:00']}, dm._dateToISO(d1))
        d1 = datetime.datetime(2012,12,21)
        self.assertEqual('2012-12-21T00:00:00', dm._dateToISO(d1))
        d1 = [datetime.datetime(2012,12,21), datetime.datetime(2012,12,22)]
        np.testing.assert_array_equal(['2012-12-21T00:00:00', '2012-12-22T00:00:00'], dm._dateToISO(d1))


class JSONTests(unittest.TestCase):
    def setUp(self):
        super(JSONTests, self).setUp()
        self.filename = 'data/20130218_rbspa_MagEphem.txt'
        self.filename_bad = 'data/20130218_rbspa_MagEphem_bad.txt'

    def tearDown(self):
        super(JSONTests, self).tearDown()

    def test_readJSONMetadata(self):
        """readJSONMetadata should read in the file"""
        dat = dm.readJSONMetadata(self.filename)
        keys = [u'PerigeePosGeod', u'S_sc_to_pfn', u'S_pfs_to_Bmin', u'Pfs_gsm',
                u'Pfn_ED_MLAT', u'ED_R', u'Dst', u'DateTime', u'DOY', u'ED_MLON',
                u'IntModel', u'ApogeePosGeod', u'CD_MLON', u'S_sc_to_pfs',
                u'GpsTime', u'JulianDate', u'M_ref', u'ED_MLT', u'Pfs_ED_MLAT',
                u'Bfs_geo', u'Bm', u'Pfn_CD_MLON', u'CD_MLAT', u'Pfs_geo',
                u'Rsm', u'Pmin_gsm', u'Rgei', u'Rgsm', u'Pfs_CD_MLAT', u'S_total',
                u'Rgeod_Height', u'Date', u'Alpha', u'M_igrf', u'Pfs_CD_MLT',
                u'ED_MLAT', u'CD_R', u'PerigeeTimes', u'UTC', u'Pfn_ED_MLT',
                u'BoverBeq', u'Lsimple', u'Lstar', u'I', u'DipoleTiltAngle',
                u'K', u'Bmin_gsm', u'S_Bmin_to_sc', u'Bfs_gsm', u'L',
                u'ApogeeTimes', u'ExtModel', u'Kp', u'Pfs_geod_LatLon',
                u'MlatFromBoverBeq', u'Pfn_gsm', u'Loss_Cone_Alpha_n', u'Bfn_geo',
                u'Pfn_CD_MLAT', u'Rgeod_LatLon', u'Pfs_ED_MLT', u'Pfs_CD_MLON',
                u'Bsc_gsm', u'Pfn_geod_Height', u'Lm_eq', u'Rgse',
                u'Pfn_geod_LatLon', u'CD_MLT', u'FieldLineType', u'Pfn_CD_MLT',
                u'Pfs_geod_Height', u'Rgeo', u'InvLat_eq', u'M_used',
                u'Loss_Cone_Alpha_s', u'Bfn_gsm', u'Pfn_ED_MLON', u'Pfn_geo',
                u'InvLat', u'Pfs_ED_MLON']
        # make sure data has all the keys and no more or less
        for k in dat:
            self.assertTrue(k in keys)
            ind = keys.index(k)
            del keys[ind]
        self.assertEqual(len(keys), 0)

    def test_readJSONMetadata_badfile(self):
        """readJSONMetadata fails on bad files"""
        self.assertRaises(ValueError, dm.readJSONMetadata, self.filename_bad)

    def test_readJSONheadedASCII(self):
        """readJSONheadedASCII should read the test file"""
        dat = dm.readJSONheadedASCII(self.filename)
        keys = [u'PerigeePosGeod', u'S_sc_to_pfn', u'S_pfs_to_Bmin', u'Pfs_gsm',
                u'Pfn_ED_MLAT', u'ED_R', u'Dst', u'DateTime', u'DOY', u'ED_MLON',
                u'IntModel', u'ApogeePosGeod', u'CD_MLON', u'S_sc_to_pfs',
                u'GpsTime', u'JulianDate', u'M_ref', u'ED_MLT', u'Pfs_ED_MLAT',
                u'Bfs_geo', u'Bm', u'Pfn_CD_MLON', u'CD_MLAT', u'Pfs_geo',
                u'Rsm', u'Pmin_gsm', u'Rgei', u'Rgsm', u'Pfs_CD_MLAT', u'S_total',
                u'Rgeod_Height', u'Date', u'Alpha', u'M_igrf', u'Pfs_CD_MLT',
                u'ED_MLAT', u'CD_R', u'PerigeeTimes', u'UTC', u'Pfn_ED_MLT',
                u'BoverBeq', u'Lsimple', u'Lstar', u'I', u'DipoleTiltAngle',
                u'K', u'Bmin_gsm', u'S_Bmin_to_sc', u'Bfs_gsm', u'L',
                u'ApogeeTimes', u'ExtModel', u'Kp', u'Pfs_geod_LatLon',
                u'MlatFromBoverBeq', u'Pfn_gsm', u'Loss_Cone_Alpha_n', u'Bfn_geo',
                u'Pfn_CD_MLAT', u'Rgeod_LatLon', u'Pfs_ED_MLT', u'Pfs_CD_MLON',
                u'Bsc_gsm', u'Pfn_geod_Height', u'Lm_eq', u'Rgse',
                u'Pfn_geod_LatLon', u'CD_MLT', u'FieldLineType', u'Pfn_CD_MLT',
                u'Pfs_geod_Height', u'Rgeo', u'InvLat_eq', u'M_used',
                u'Loss_Cone_Alpha_s', u'Bfn_gsm', u'Pfn_ED_MLON', u'Pfn_geo',
                u'InvLat', u'Pfs_ED_MLON']
        # make sure data has all the keys and no more or less
        for k in dat:
            self.assertTrue(k in keys)
            ind = keys.index(k)
            del keys[ind]
        self.assertEqual(len(keys), 0)
        dat = dm.readJSONheadedASCII(self.filename, convert=True)
        np.testing.assert_array_equal(dat['DateTime'], [datetime.datetime(2013, 2, 18, 0, 0), datetime.datetime(2013, 2, 18, 0, 5)])

    def test_idl2html(self):
        """_idl2html should have known output"""
        self.assertEqual('R<sub>e</sub>', dm._idl2html('R!Ie'))
        self.assertEqual('R<sub>e</sub>', dm._idl2html('R!Ie!N'))
        self.assertEqual('R<sup>e</sup>', dm._idl2html('R!Ee'))
        self.assertEqual('R<sup>e</sup>', dm._idl2html('R!Ee!N'))
        self.assertEqual('Hello World', dm._idl2html('Hello World'))

    def test_toHTML(self):
        """toHTML should give known output"""
        t_file = tempfile.NamedTemporaryFile(delete=False)
        t_file.close()
        dat = dm.readJSONheadedASCII(self.filename)
        dm.toHTML(t_file.name, dat, attrs=['DESCRIPTION', 'UNITS', 'ELEMENT_LABELS'], varLinks=True)
        self.assertEqual(12834, os.path.getsize(t_file.name)) # not the best test but I am lazy
        os.remove(t_file.name)

    def test_writeJSONMetadata(self):
        """reading metadata should give same keys as original datamodel"""
        dat = dm.readJSONMetadata(self.filename)
        # make sure data has all te keys and no more or less
        t_file = tempfile.NamedTemporaryFile(delete=False)
        t_file.close()
        dm.writeJSONMetadata(t_file.name, dat)
        dat2 = dm.readJSONheadedASCII(t_file.name)
        os.remove(t_file.name)
        keylist1 = sorted(dat.keys())
        keylist2 = sorted(dat2.keys())
        self.assertTrue(keylist1==keylist2)
        #now test that values in some metadata are identical
        self.assertTrue((dat['PerigeePosGeod'] == dat2['PerigeePosGeod']).all())

    def test_toJSONheadedASCII(self):
        """Write known datamodel to JSON-headed ASCII and ensure it has right stuff added"""
        a = dm.SpaceData()
        a.attrs['Global'] = 'A global attribute'
        a['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
        a['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
        a['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
        t_file = tempfile.NamedTemporaryFile(delete=False)
        t_file.close()
        dm.toJSONheadedASCII(t_file.name, a, depend0='Var1', order=['Var1','Var2'])
        dat2 = dm.readJSONheadedASCII(t_file.name)
        #test global attr
        self.assertTrue(a.attrs==dat2.attrs)
        #test that metadata is back and all original keys are present
        for key in a['MVar'].attrs:
            self.assertTrue(key in dat2['MVar'].attrs)
        np.testing.assert_array_equal(a['MVar'], dat2['MVar'])
        #test vars are right
        np.testing.assert_allclose(a['Var1'], dat2['Var1'])
        np.testing.assert_allclose(a['Var2'], dat2['Var2'])
        #test for added dimension and start col
        self.assertTrue(dat2['Var1'].attrs['DIMENSION']==[1])
        self.assertTrue(dat2['Var2'].attrs['DIMENSION']==[2])



if __name__ == "__main__":
    unittest.main()
