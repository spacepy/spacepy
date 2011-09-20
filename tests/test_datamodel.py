#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's datamodel

Copyright ©2010 Los Alamos National Security, LLC.
"""


from __future__ import division

import datetime
import os
import tempfile
import unittest

import spacepy.toolbox as tb
import spacepy.datamodel as dm
import numpy as np

try:
    import cPickle as pickle
except:
    import pickle


class SpaceDataTests(unittest.TestCase):
    def setUp(self):
        super(SpaceDataTests, self).setUp()

    def tearDown(self):
        super(SpaceDataTests, self).tearDown()

    def test_SpaceData(self):
        """Spacedata dist object has certian attributes"""
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


class dmarrayTests(unittest.TestCase):
    def setUp(self):
        super(dmarrayTests, self).setUp()
        self.dat = dm.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})

    def tearDown(self):
        super(dmarrayTests, self).tearDown()
        del self.dat

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
        if os.path.isfile('dmh5test.h5'):
            os.remove('dmh5test.h5')

    def tearDown(self):
        super(converterTests, self).tearDown()
        del self.SDobj
        if os.path.isfile('dmh5test.h5'):
            os.remove('dmh5test.h5')

    def test_convertKeysToStr(self):
        """convertKeysToStr sjould give known output"""
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
        dm.toHDF5('dmh5test.h5', self.SDobj)
        newobj = dm.fromHDF5('dmh5test.h5')
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_allclose(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])


if __name__ == "__main__":
    unittest.main()
