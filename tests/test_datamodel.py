#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's datamodel

Copyright 2010-2012 Los Alamos National Security, LLC.
"""


from __future__ import division

import copy
import datetime
import marshal
import os
import os.path
import tempfile
import unittest

try:
    import StringIO
except ImportError:
    import io as StringIO
import sys
import warnings

import spacepy_testing
import spacepy.datamodel as dm
import spacepy.time as spt
import numpy as np

try:
    import cPickle as pickle
except:
    import pickle

# python2 python3 string wrangling
try:
    str_classes = (str, bytes, unicode)
except NameError:
    str_classes = (str, bytes)
    unicode = str

__all__ = ['SpaceDataTests', 'dmarrayTests', 'converterTests', 'JSONTests', 'converterTestsCDF']


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
        np.testing.assert_almost_equal(a[1], b[1])
        self.assertEqual(a[1].attrs, b[1].attrs)
        b = dm.dmcopy(a[1])
        np.testing.assert_almost_equal(a[1], b)
        self.assertEqual(a[1].attrs, b.attrs)
        a = np.arange(10)
        b = dm.dmcopy(a)
        np.testing.assert_almost_equal(a, b)
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
        """Flatten should flatten a nested SpaceData"""
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
        self.assertEqual(sorted(b.keys()),
                         sorted(['1<--pig<--fish<--a', '4<--cat', '1<--dog', '1<--pig<--fish<--b', '5']))

    def test_flatten_function_TypeError(self):
        """flatten, check the try:except: at the top"""
        # interesting behavior on what is created
        self.assertEqual(dm.flatten(np.arange(5)),
                         {0: 0, 1: 1, 2: 2, 3: 3, 4: 4})

    def test_unflatten_function(self):
        """Unflatten should unflatten a flattened SpaceData"""
        a = dm.SpaceData()
        a['1'] = dm.SpaceData(dog = 5, pig = dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        a['4'] = dm.SpaceData(cat = 'kitty')
        a['5'] = 4
        a[9] = dm.dmarray([1,2,3])
        b = dm.flatten(a)
        c = dm.unflatten(b)
        self.assertTrue(9 in a.keys())
        self.assertTrue(9 in c.keys())
        del a[9]
        del c[9]
        self.assertEqual(sorted(a.keys()), sorted(c.keys()))
        self.assertEqual(sorted(a['1'].keys()), sorted(c['1'].keys()))
        self.assertEqual(sorted(a['1']['pig'].keys()), sorted(c['1']['pig'].keys()))
        self.assertEqual(sorted(a['1']['pig']['fish'].keys()), sorted(c['1']['pig']['fish'].keys()))

    def test_unflatten_function_TypeError(self):
        """unflatten, check the try:except: at the top"""
        # interesting behavior on what is created
        self.assertEqual(dm.unflatten(np.arange(5)),
                         {0: 0, 1: 1, 2: 2, 3: 3, 4: 4})

    def test_unflatten_function_nonflat(self):
        """unflatten will error for a nested input"""
        a = dm.SpaceData()
        a['1'] = dm.SpaceData(dog=5, pig=dm.SpaceData(fish=dm.SpaceData(a='carp', b='perch')))
        a['4'] = dm.SpaceData(cat='kitty')
        a['5'] = 4
        with self.assertRaises(TypeError):
            dm.unflatten(a)

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
        val = sorted(a.keys())
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
        ans = ['4<--cat', '1<--dog', '1<--pig<--fish<--a', '1<--pig<--fish<--b']
        ans.sort()
        self.assertTrue(5 in a.keys())
        del a[5]
        val = sorted(a.keys())
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

    def test_fromRecArrayNames(self):
        '''given a known rec array, should get known keys in SpaceData'''
        names = ['x', 'y', 'value']
        recarr = np.zeros(3, dtype=[(names[0],'f4'),(names[1],np.float32),(names[2],'f4',(2,2))])
        sd = dm.fromRecArray(recarr)
        for kk in sd.keys():
            self.assertTrue(kk in names)

    def test_fromRecArrayShapesValues(self):
        '''given a known rec array, should get known shapes and values in SpaceData'''
        names = ['x', 'y', 'value']
        recarr = np.zeros(3, dtype=[(names[0],'f4'),(names[1],np.float32),(names[2],'f4',(2,2))])
        sd = dm.fromRecArray(recarr)
        np.testing.assert_array_equal((3,2,2), sd[names[2]].shape)
        np.testing.assert_array_equal(np.zeros(3), sd[names[0]])

    def test_multiget(self):
        '''Allow for multiple keys to be specified'''
        a = dm.SpaceData()
        a['a'] = dm.dmarray([1,2,3])
        a['b'] = dm.dmarray([1,2,3])
        a['c'] = dm.dmarray([1,2,3])
        a.attrs['foo']='bar'
        np.testing.assert_equal(a['a'],  dm.dmarray([1,2,3]))
        np.testing.assert_equal(a[['a', 'b']],  {'a': dm.dmarray([1, 2, 3]), 'b': dm.dmarray([1, 2, 3])})
        self.assertRaises(KeyError, a.__getitem__, 'NotAkey')
        self.assertRaises(KeyError, a.__getitem__, ['a', 'nokey'])

    def test_resample_input(self):
        '''resample requires SpaceData or dmarray'''
        self.assertRaises(TypeError, dm.resample, [1,2,3])

    def test_resample_shape(self):
        '''resample should give consistent results, 1d or 2d'''
        a = dm.SpaceData()
        a['a'] = dm.dmarray(range(10*3*4)).reshape(10,3,4)
        a['b'] = dm.dmarray(range(10)) + 4
        a['c'] = dm.dmarray(range(3)) + 10
        times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i) for i in range(10)]
        self.assertRaises(IndexError, dm.resample, a, times, datetime.timedelta(hours=2), datetime.timedelta(hours=0))
        
    def test_resample1(self):
        '''resample should give consistent results'''
        ans = dm.SpaceData()
        ans.attrs['foo'] = 'bar'
        ans['a'] = [ 0.5,  2.5,  4.5,  6.5, 8.5]
        ans['b'] = dm.dmarray([4.5, 6.5, 8.5, 10.5, 12.5])
        ans['b'].attrs['marco'] = 'polo'
        ans['Epoch'] = [datetime.datetime(2010, 1, 1, 1, 0),
                        datetime.datetime(2010, 1, 1, 3, 0),
                        datetime.datetime(2010, 1, 1, 5, 0),
                        datetime.datetime(2010, 1, 1, 7, 0),
                        datetime.datetime(2010, 1, 1, 9, 0)]
        
        
        a = dm.SpaceData()
        a['a'] = dm.dmarray(range(10)) #sequence 0 through 9
        a['b'] = dm.dmarray(range(10)) + 4 #sequence 4 through 13
        a['b'].attrs['marco'] = 'polo'
        a['c'] = dm.dmarray(range(3)) + 10
        a.attrs['foo'] = 'bar'
        times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i) for i in range(10)]
        out = dm.resample(a, times, winsize=datetime.timedelta(hours=2), overlap=datetime.timedelta(hours=0))
        #starting from element 0, step two hours (points) with no overlap and average
        #therefore each resulting value should be the mean of a pair of points
        #that is, ans['a'][0] should be the mean of 0 and 1
        for k, v in out.items():
            np.testing.assert_equal(v,  ans[k])
        self.assertEqual(ans.attrs, out.attrs)
        self.assertEqual(ans['b'].attrs['marco'], 'polo')
        self.assertTrue(out['b'].attrs['DEPEND_0'], 'Epoch')
        self.assertFalse('c' in out)

    def test_resample2(self):
        '''resample should give consistent results (ticktock)'''
        ans = {}
        ans['a'] = [ 0.5,  2.5,  4.5,  6.5, 8.5]
        ans['b'] = [4.5, 6.5, 8.5, 10.5, 12.5]
        ans['Epoch'] = [datetime.datetime(2010, 1, 1, 1, 0),
                        datetime.datetime(2010, 1, 1, 3, 0),
                        datetime.datetime(2010, 1, 1, 5, 0),
                        datetime.datetime(2010, 1, 1, 7, 0),
                        datetime.datetime(2010, 1, 1, 9, 0)]
        #For justification of test results see test above (test_resample1)
        a = dm.SpaceData()
        a['a'] = dm.dmarray(range(10))
        a['b'] = dm.dmarray(range(10)) + 4
        a['c'] = dm.dmarray(range(3)) + 10
        times = spt.Ticktock([datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i) for i in range(10)])
        out = dm.resample(a, times, winsize=datetime.timedelta(hours=2), overlap=datetime.timedelta(hours=0))
        for k, v in out.items():
            np.testing.assert_equal(v,  ans[k])

    def test_resample3(self):
        '''resample should give consistent results (2d)'''
        ans = {}
        ans['a'] = [[  1.,   2.],
                    [  5.,   6.],
                    [  9.,  10.],
                    [ 13.,  14.],
                    [ 17.,  18.]]
        ans['b'] = [4.5, 6.5, 8.5, 10.5, 12.5]
        ans['Epoch'] = [datetime.datetime(2010, 1, 1, 1, 0),
                        datetime.datetime(2010, 1, 1, 3, 0),
                        datetime.datetime(2010, 1, 1, 5, 0),
                        datetime.datetime(2010, 1, 1, 7, 0),
                        datetime.datetime(2010, 1, 1, 9, 0)]
        
        a = dm.SpaceData()
        a['a'] = dm.dmarray(range(10*2)).reshape(10,2)
        a['b'] = dm.dmarray(range(10)) + 4
        a['c'] = dm.dmarray(range(3)) + 10
        times = [datetime.datetime(2010, 1, 1) + datetime.timedelta(hours=i) for i in range(10)]
        out = dm.resample(a, times, winsize=datetime.timedelta(hours=2), overlap=datetime.timedelta(hours=0))
        for k, v in out.items():
            np.testing.assert_equal(v,  ans[k])

    def test_readmeta(self):
        """Check on reading from the meta property"""
        a = dm.SpaceData(attrs={'a': 1, 'b': 2})
        self.assertEqual(1, a.meta['a'])

    def test_writemeta(self):
        """Check on writing to the meta property"""
        a = dm.SpaceData(attrs={'a': 1, 'b': 2})
        a.meta['c'] = 3
        self.assertEqual(3, a.attrs['c'])
        a.meta['a'] = 99
        self.assertEqual(99, a.attrs['a'])

    def test_assignmeta(self):
        """Assign to the meta property"""
        a = dm.SpaceData(attrs={'a': 1, 'b': 2})
        a.meta = {'c': 3}
        self.assertEqual(3, a.attrs['c'])
        self.assertFalse('a' in a.attrs)

    def test_deletemeta(self):
        """Remove the meta property"""
        a = dm.SpaceData(attrs={'a': 1, 'b': 2})
        del a.meta
        self.assertFalse(hasattr(a, 'attrs'))
        self.assertFalse(hasattr(a, 'meta'))
        # Known failure, cannot delete property from instance.
#        self.assertFalse('meta' in dir(a))

    def test_pickle(self):
        """Make sure that SpaceData objects behave with pickle"""
        a = dm.SpaceData({'a': [1, 2, 3]})
        a.attrs['FILLVAL'] = 123
        p_str = pickle.dumps(a)
        ret = pickle.loads(p_str)
        assert a == ret
        assert a.attrs == ret.attrs

    @unittest.expectedFailure
    def test_marshal(self):
        """SpaceData objects don't marshal, note that here"""
        a = dm.SpaceData({'a': [1, 2, 3]})
        a.attrs['FILLVAL'] = 123
        p_str = marshal.dumps(a)
        ret = marshal.loads(p_str)
        assert a == ret
        assert a.attrs == ret.attrs


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
        self.assertEqual(d2.attrs, self.dat.attrs)
        self.assertFalse(d2.attrs is self.dat.attrs)

    def test_vstack(self):
        """vstack should maintain all Allowed_Attributes"""
        d2 = dm.dmarray.vstack(self.dat, [-1,-2,-3,-4])
        np.testing.assert_array_equal(
            np.asarray([[ 1,  2,  3,  4],[-1, -2, -3, -4]]), d2)
        self.assertEqual(d2.attrs, self.dat.attrs)
        self.assertFalse(d2.attrs is self.dat.attrs)

    def test_hstack(self):
        """hstack should maintain all Allowed_Attributes"""
        d2 = dm.dmarray.hstack(self.dat, [-1,-2,-3,-4])
        np.testing.assert_array_equal(
            np.asarray([ 1,  2,  3,  4, -1, -2, -3, -4]), d2)
        self.assertEqual(d2.attrs, self.dat.attrs)
        self.assertFalse(d2.attrs is self.dat.attrs)

    def test_dstack(self):
        """dstack should maintain all Allowed_Attributes"""
        d2 = dm.dmarray.dstack(self.dat, [-1,-2,-3,-4])
        np.testing.assert_array_equal(
            np.asarray([[[ 1, -1], [ 2, -2], [ 3, -3], [ 4, -4]]]), d2)
        self.assertEqual(d2.attrs, self.dat.attrs)
        self.assertFalse(d2.attrs is self.dat.attrs)

    def test_concatenate(self):
        """concatenate should maintain all Allowed_Attributes"""
        d2 = dm.dmarray.concatenate(self.dat, [-1,-2,-3,-4])
        np.testing.assert_array_equal(
            np.asarray([ 1,  2,  3,  4, -1, -2, -3, -4]), d2)
        self.assertEqual(d2.attrs, self.dat.attrs)
        self.assertFalse(d2.attrs is self.dat.attrs)

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
        np.testing.assert_almost_equal([1,2,3], data2)

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
        np.testing.assert_almost_equal(self.dat, dat2)
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

    def test_readmeta(self):
        """Check on reading from the meta property"""
        a = dm.dmarray([1, 2, 3], attrs={'a': 1, 'b': 2})
        self.assertEqual(1, a.meta['a'])

    def test_writemeta(self):
        """Check on writing to the meta property"""
        a = dm.dmarray([1, 2, 3], attrs={'a': 1, 'b': 2})
        a.meta['c'] = 3
        self.assertEqual(3, a.attrs['c'])
        a.meta['a'] = 99
        self.assertEqual(99, a.attrs['a'])

    def test_assignmeta(self):
        """Assign to the meta property"""
        a = dm.dmarray([1, 2, 3], attrs={'a': 1, 'b': 2})
        a.meta = {'c': 3}
        self.assertEqual(3, a.attrs['c'])
        self.assertFalse('a' in a.attrs)

    def test_deletemeta(self):
        """Remove the meta property"""
        a = dm.dmarray([1, 2, 3], attrs={'a': 1, 'b': 2})
        del a.meta
        self.assertFalse(hasattr(a, 'attrs'))
        self.assertFalse(hasattr(a, 'meta'))
        # Known failure, cannot delete property from instance.
#        self.assertFalse('meta' in dir(a))

    def test_dmfilled(self):
        """dmfilled should fill an array"""
        ans = np.asarray([1,1,1], dtype=int)
        tst = dm.dmfilled((3), fillval=1, dtype=int)
        np.testing.assert_equal(tst, ans)
        self.assertEqual(tst.dtype, ans.dtype)

        ans = np.asarray([datetime.datetime(2000, 1,1),
                          datetime.datetime(2000, 1,1),
                          datetime.datetime(2000, 1,1)], dtype=object)
        tst = dm.dmfilled((3), fillval=datetime.datetime(2000, 1,1), dtype=object)
        np.testing.assert_equal(tst, ans)
        self.assertEqual(tst.dtype, ans.dtype)

    def test_dmfilled_recarray(self):
        """dmfilled should fill a recarray"""
        dt = [('foo', 'f4'), ('bar', 'i2')]
        tst = dm.dmfilled(4, fillval=3, dtype=dt)
        a = np.empty((4, ), dt)
        a.fill(3)
        ans = dm.dmarray(a)
        np.testing.assert_equal(tst, ans)

    def test_toRecArray_contents(self):
        '''a record array can be created from a SpaceData, keys and values equal'''
        sd = dm.SpaceData()
        sd['x'] = dm.dmarray([1.0, 2.0])
        sd['y'] = dm.dmarray([2,4])
        ra = dm.toRecArray(sd)
        np.testing.assert_equal(ra['x'], [1.0, 2.0])
        np.testing.assert_equal(ra['y'], [2, 4])
        self.assertEqual(['x', 'y'], sorted(ra.dtype.fields))

    def test_toRecArray_dtypes1(self):
        '''recarray created from dmarray preserves data types (32-bit)'''
        sd = dm.SpaceData()
        sd['x'] = dm.dmarray([1.0, 2.0], dtype=np.float32)
        sd['y'] = dm.dmarray([2,4], dtype=np.int32)
        ra = dm.toRecArray(sd)
        expected = [sd[key].dtype for key in sd]
        got = [ra.dtype[name] for name in ra.dtype.names]
        self.assertEqual(expected, got)

    def test_toRecArray_dtypes2(self):
        '''recarray created from dmarray preserves data types (16-bit+str)'''
        sd = dm.SpaceData()
        sd['x'] = dm.dmarray([1.0, 2.0], dtype=np.float16)
        sd['y'] = dm.dmarray([2,4], dtype=np.int16)
        sd['str'] = dm.dmarray(['spam', 'eggs'], dtype='|S5')
        ra = dm.toRecArray(sd)
        expected = [sd[key].dtype for key in sd]
        got = [ra.dtype[name] for name in ra.dtype.names]
        self.assertEqual(expected, got)
        
class converterTests(unittest.TestCase):
    def setUp(self):
        super(converterTests, self).setUp()
        self.SDobj = dm.SpaceData(attrs={'global': 'test'})
        self.SDobj['var'] = dm.dmarray([1, 2, 3], attrs={'a': 'a'})
        self.testdir = tempfile.mkdtemp()
        self.testfile = os.path.join('test.h5')
        warnings.simplefilter('error', dm.DMWarning)

    def tearDown(self):
        super(converterTests, self).tearDown()
        del self.SDobj
        if os.path.exists(self.testfile):
            os.remove(self.testfile)
        os.rmdir(self.testdir)
        warnings.simplefilter('default', dm.DMWarning)

    def test_convertKeysToStr(self):
        """convertKeysToStr should give known output"""
        a = dm.SpaceData()
        a['data'] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual(list(a.keys()), list(b.keys()))
        a = dm.SpaceData()
        a[50] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(list(a.keys())[0])], list(b.keys()))
        a = {}
        a[50] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(list(a.keys())[0])], list(b.keys()))
        a = dm.SpaceData()
        a['data'] = dm.SpaceData()
        a['data']['test'] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual(list(a.keys()), list(b.keys()))
        a = dm.SpaceData()
        a[50] = dm.SpaceData()
        a[50][49] = dm.dmarray([1,2,3])
        b = dm.convertKeysToStr(a)
        self.assertEqual([str(list(a.keys())[0])], list(b.keys()))

    def test_toHDF5ListString(self):
        """Convert to HDF5, including a list of string in attributes"""
        a = dm.SpaceData()
        a.attrs['foo'] = ['hi']
        dm.toHDF5(self.testfile, a, mode='a')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(a.attrs['foo'], newobj.attrs['foo'])

    def test_toHDF5_method(self):
        """Convert to HDF5, using the method, catching #404"""
        a = dm.SpaceData({'dat': dm.dmarray([1, 2, 3])})
        a.toHDF5(self.testfile, mode='a')
        newobj = dm.fromHDF5(self.testfile)
        np.testing.assert_array_equal([1, 2, 3], newobj['dat'])

    def test_copy_toHDF(self):
        """Removed code in #404 for an old python bug, issue4380, make sure it works"""
        a = dm.SpaceData({'dat': dm.dmarray([1, 2, 3])})
        a2 = copy.deepcopy(a)
        a2.toHDF5(self.testfile, mode='a')
        newobj = dm.fromHDF5(self.testfile)
        np.testing.assert_array_equal([1, 2, 3], newobj['dat'])

    def test_HDF5roundtrip(self):
        """Data can go to hdf and back"""
        dm.toHDF5(self.testfile, self.SDobj)
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])
        dm.toHDF5(self.testfile, self.SDobj, mode='a')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])

    def test_HDF5roundtrip_method(self):
        """Data can go to hdf and back"""
        self.SDobj.toHDF5(self.testfile)
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])
        dm.toHDF5(self.testfile, self.SDobj, mode='a')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])

    def test_HDF5roundtripGZIP(self):
        """Data can go to hdf and back with compression"""
        dm.toHDF5(self.testfile, self.SDobj, compression='gzip')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])
        dm.toHDF5(self.testfile, self.SDobj, mode='a', compression='gzip')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(self.SDobj.attrs['global'], newobj.attrs['global'])
        np.testing.assert_almost_equal(self.SDobj['var'], newobj['var'])
        self.assertEqual(self.SDobj['var'].attrs['a'], newobj['var'].attrs['a'])

    def test_HDF5Exceptions(self):
        """HDF5 has warnings and exceptions"""
        dm.toHDF5(self.testfile, self.SDobj)
        self.assertRaises(IOError, dm.toHDF5, self.testfile, self.SDobj, overwrite=False)
        a = dm.SpaceData()
        a['foo'] = 'bar' # not an allowed type for data
        self.assertRaises(dm.DMWarning, dm.toHDF5, self.testfile, a)

    def test_HDF5roundtrip2(self):
        """Data can go to hdf without altering datetimes in the datamodel"""
        a = dm.SpaceData()
        a['foo'] = dm.SpaceData()
        dm.toHDF5(self.testfile, a)
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(a['foo'], newobj['foo'])
        a['bar'] = dm.dmarray([datetime.datetime(2000, 1, 1)])
        dm.toHDF5(self.testfile, a)
        self.assertEqual(a['bar'], dm.dmarray([datetime.datetime(2000, 1, 1)]))

    def test_HDF5roundtrip2GZIP(self):
        """Data can go to hdf without altering datetimes in the datamodel with compression"""
        a = dm.SpaceData()
        a['foo'] = dm.SpaceData()
        dm.toHDF5(self.testfile, a, compression='gzip')
        newobj = dm.fromHDF5(self.testfile)
        self.assertEqual(a['foo'], newobj['foo'])
        a['bar'] = dm.dmarray([datetime.datetime(2000, 1, 1)])
        dm.toHDF5(self.testfile, a, compression='gzip')
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
        np.testing.assert_array_equal(dm._dateToISO(d1), ['2012-12-21T00:00:00', '2012-12-22T00:00:00'])

class converterTestsCDF(unittest.TestCase):
    def setUp(self):
        super(converterTestsCDF, self).setUp()
        self.SDobj = dm.SpaceData(attrs={'global': 'test'})
        self.SDobj['var'] = dm.dmarray([1, 2, 3], attrs={'a': 'a'})
        self.testdir = tempfile.mkdtemp()
        self.testfile = os.path.join(self.testdir, 'test.cdf')
        warnings.simplefilter('error', dm.DMWarning)

    def tearDown(self):
        super(converterTestsCDF, self).tearDown()
        del self.SDobj
        if os.path.exists(self.testfile):
            os.remove(self.testfile)
        os.rmdir(self.testdir)
        warnings.simplefilter('default', dm.DMWarning)

    def test_toCDFroundtrip(self):
        """toCDF should be able to make a file and then read it in the same"""
        dm.toCDF(self.testfile, self.SDobj)
        tst = dm.fromCDF(self.testfile)
        for k in self.SDobj:
            np.testing.assert_array_equal(self.SDobj[k], tst[k])

    def test_toCDFroundtrip_method(self):
        """toCDF should be able to make a file and then read it in the same"""
        self.SDobj.toCDF(self.testfile)
        tst = dm.fromCDF(self.testfile)
        for k in self.SDobj:
            np.testing.assert_array_equal(self.SDobj[k], tst[k])

    def test_toCDF_method(self):
        """Convert to CDF, using the method, catching #404"""
        a = dm.SpaceData({'dat': dm.dmarray([1, 2, 3])})
        a.toCDF(self.testfile, mode='a')
        newobj = dm.fromCDF(self.testfile)
        np.testing.assert_array_equal([1, 2, 3], newobj['dat'])

class JSONTests(unittest.TestCase):
    def setUp(self):
        super(JSONTests, self).setUp()
        self.filename = os.path.join(
            spacepy_testing.datadir, '20130218_rbspa_MagEphem.txt')
        self.filename_bad = os.path.join(
            spacepy_testing.datadir, '20130218_rbspa_MagEphem_bad.txt')
        self.testdir = tempfile.mkdtemp()
        self.testfile = os.path.join(self.testdir, 'test.cdf')

    def tearDown(self):
        super(JSONTests, self).tearDown()
        if os.path.exists(self.testfile):
            os.remove(self.testfile)
        os.rmdir(self.testdir)

    def test_readJSONMetadata(self):
        """readJSONMetadata should read in the file"""
        dat = dm.readJSONMetadata(self.filename)
        keys = ['PerigeePosGeod', 'S_sc_to_pfn', 'S_pfs_to_Bmin', 'Pfs_gsm',
                'Pfn_ED_MLAT', 'ED_R', 'Dst', 'DateTime', 'DOY', 'ED_MLON',
                'IntModel', 'ApogeePosGeod', 'CD_MLON', 'S_sc_to_pfs',
                'GpsTime', 'JulianDate', 'M_ref', 'ED_MLT', 'Pfs_ED_MLAT',
                'Bfs_geo', 'Bm', 'Pfn_CD_MLON', 'CD_MLAT', 'Pfs_geo',
                'Rsm', 'Pmin_gsm', 'Rgei', 'Rgsm', 'Pfs_CD_MLAT', 'S_total',
                'Rgeod_Height', 'Date', 'Alpha', 'M_igrf', 'Pfs_CD_MLT',
                'ED_MLAT', 'CD_R', 'PerigeeTimes', 'UTC', 'Pfn_ED_MLT',
                'BoverBeq', 'Lsimple', 'Lstar', 'I', 'DipoleTiltAngle',
                'K', 'Bmin_gsm', 'S_Bmin_to_sc', 'Bfs_gsm', 'L',
                'ApogeeTimes', 'ExtModel', 'Kp', 'Pfs_geod_LatLon',
                'MlatFromBoverBeq', 'Pfn_gsm', 'Loss_Cone_Alpha_n', 'Bfn_geo',
                'Pfn_CD_MLAT', 'Rgeod_LatLon', 'Pfs_ED_MLT', 'Pfs_CD_MLON',
                'Bsc_gsm', 'Pfn_geod_Height', 'Lm_eq', 'Rgse',
                'Pfn_geod_LatLon', 'CD_MLT', 'FieldLineType', 'Pfn_CD_MLT',
                'Pfs_geod_Height', 'Rgeo', 'InvLat_eq', 'M_used',
                'Loss_Cone_Alpha_s', 'Bfn_gsm', 'Pfn_ED_MLON', 'Pfn_geo',
                'InvLat', 'Pfs_ED_MLON']
        if str is bytes:
            keys = [unicode(k) for k in keys]
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
        keys = ['PerigeePosGeod', 'S_sc_to_pfn', 'S_pfs_to_Bmin', 'Pfs_gsm',
                'Pfn_ED_MLAT', 'ED_R', 'Dst', 'DateTime', 'DOY', 'ED_MLON',
                'IntModel', 'ApogeePosGeod', 'CD_MLON', 'S_sc_to_pfs',
                'GpsTime', 'JulianDate', 'M_ref', 'ED_MLT', 'Pfs_ED_MLAT',
                'Bfs_geo', 'Bm', 'Pfn_CD_MLON', 'CD_MLAT', 'Pfs_geo',
                'Rsm', 'Pmin_gsm', 'Rgei', 'Rgsm', 'Pfs_CD_MLAT', 'S_total',
                'Rgeod_Height', 'Date', 'Alpha', 'M_igrf', 'Pfs_CD_MLT',
                'ED_MLAT', 'CD_R', 'PerigeeTimes', 'UTC', 'Pfn_ED_MLT',
                'BoverBeq', 'Lsimple', 'Lstar', 'I', 'DipoleTiltAngle',
                'K', 'Bmin_gsm', 'S_Bmin_to_sc', 'Bfs_gsm', 'L',
                'ApogeeTimes', 'ExtModel', 'Kp', 'Pfs_geod_LatLon',
                'MlatFromBoverBeq', 'Pfn_gsm', 'Loss_Cone_Alpha_n', 'Bfn_geo',
                'Pfn_CD_MLAT', 'Rgeod_LatLon', 'Pfs_ED_MLT', 'Pfs_CD_MLON',
                'Bsc_gsm', 'Pfn_geod_Height', 'Lm_eq', 'Rgse',
                'Pfn_geod_LatLon', 'CD_MLT', 'FieldLineType', 'Pfn_CD_MLT',
                'Pfs_geod_Height', 'Rgeo', 'InvLat_eq', 'M_used',
                'Loss_Cone_Alpha_s', 'Bfn_gsm', 'Pfn_ED_MLON', 'Pfn_geo',
                'InvLat', 'Pfs_ED_MLON']
        if str is bytes:
            keys = [unicode(k) for k in keys]
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
        if sys.platform == 'win32': #Different line endings
            expected = 12916 if str is bytes else 12892
        else:
            expected = 12834 if str is bytes else 12810 #no u on unicode strings
        self.assertEqual(expected, os.path.getsize(t_file.name)) # not the best test but I am lazy
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
        np.testing.assert_almost_equal(a['Var1'], dat2['Var1'])
        np.testing.assert_almost_equal(a['Var2'], dat2['Var2'])
        #test for added dimension and start col
        self.assertTrue(dat2['Var1'].attrs['DIMENSION']==[1])
        self.assertTrue(dat2['Var2'].attrs['DIMENSION']==[2])
        os.remove(t_file.name)

    def test_toJSONheadedASCII_method(self):
        """Write known datamodel to JSON-headed ASCII and ensure it has right stuff added"""
        a = dm.SpaceData()
        a.attrs['Global'] = 'A global attribute'
        a['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
        a['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
        a['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
        t_file = tempfile.NamedTemporaryFile(delete=False)
        t_file.close()
        a.toJSONheadedASCII(t_file.name, depend0='Var1', order=['Var1','Var2'])
        dat2 = dm.readJSONheadedASCII(t_file.name)
        #test global attr
        self.assertTrue(a.attrs==dat2.attrs)
        #test that metadata is back and all original keys are present
        for key in a['MVar'].attrs:
            self.assertTrue(key in dat2['MVar'].attrs)
        np.testing.assert_array_equal(a['MVar'], dat2['MVar'])
        #test vars are right
        np.testing.assert_almost_equal(a['Var1'], dat2['Var1'])
        np.testing.assert_almost_equal(a['Var2'], dat2['Var2'])
        #test for added dimension and start col
        self.assertTrue(dat2['Var1'].attrs['DIMENSION']==[1])
        self.assertTrue(dat2['Var2'].attrs['DIMENSION']==[2])
        os.remove(t_file.name)

    def test_toJSONheadedASCII_method_404(self):
        """Convert to toJSONheadedASCII, using the method, catching #404"""
        a = dm.SpaceData({'dat': dm.dmarray([1, 2, 3])})
        a.toJSONheadedASCII(self.testfile, mode='a')
        newobj = dm.readJSONheadedASCII(self.testfile)
        np.testing.assert_array_equal([1, 2, 3], newobj['dat'])

    def test_toJSONmetadata_globals(self):
        """Test for handling of int, float, bool, list & dict in global attrs"""
        a = dm.SpaceData()
        a.attrs['Global1'] = 'A global string attribute'
        a.attrs['Global2'] = 2
        a.attrs['Global3'] = 3.0
        a.attrs['Global4'] = [1,2,3,4]
        a.attrs['Global5'] = True
        a['Var1'] = dm.dmarray([1,2,3,4,5], attrs={'Local1': 'A local attribute'})
        a['Var2'] = dm.dmarray([[8,9],[9,1],[3,4],[8,9],[7,8]])
        a['MVar'] = dm.dmarray([7.8], attrs={'Note': 'Metadata'})
        t_file = tempfile.NamedTemporaryFile(delete=False)
        t_file.close()
        dm.writeJSONMetadata(t_file.name, a)
        dat2 = dm.readJSONMetadata(t_file.name)
        #test global attr
        self.assertTrue(a.attrs==dat2.attrs)
        #test that metadata is back and all original keys are present
        for key in a['MVar'].attrs:
            self.assertTrue(key in dat2['MVar'].attrs)
        np.testing.assert_array_equal(a['MVar'], dat2['MVar'])
        os.remove(t_file.name)

    def test_toJSON_timeunaltered(self):
        """Test to check that stored datetimes aren't changed on write"""
        data = dm.SpaceData()
        data['Epoch'] = spt.tickrange('20200101', '20200102',
                                      deltadays=datetime.timedelta(hours=1)).UTC
        exptype = type(data['Epoch'][0]) #datetime.datetime
        # save to file, then immediately clean up
        fname = None
        try:
            with tempfile.NamedTemporaryFile(delete=False) as fp:
                fname = fp.name
                data.toJSONheadedASCII(fname)
        finally:
            if fname != None:
                os.remove(fname)
        restype = type(data['Epoch'][0])
        self.assertEqual(exptype, restype)

    def test_dateToISOunaltered_SD(self):
        """Test to check that _dateToISO doesn't change datatypes, input SpaceData"""
        data = dm.SpaceData()
        data['Epoch'] = spt.tickrange('20200101', '20200102',
                                      deltadays=datetime.timedelta(hours=1)).UTC
        exptype = type(data['Epoch'][0])
        newdata = dm._dateToISO(data)
        restype = type(data['Epoch'][0])
        self.assertEqual(exptype, restype)

    def test_dateToISOunaltered_dm(self):
        """Test to check that _dateToISO doesn't change datatypes, input dmarray"""
        data = spt.tickrange('20200101', '20200102',
                             deltadays=datetime.timedelta(hours=1)).UTC
        exptype = type(data[0])
        newdata = dm._dateToISO(data)
        restype = type(data[0])
        self.assertEqual(exptype, restype)


class VariableTests(unittest.TestCase):
    def test_createISTPattrs_data(self):
        """createISTPattrs should give known results and error check for data"""
        # catch bad datatype
        with self.assertRaises(ValueError):
            dm.createISTPattrs('badval')
        # test with a datatype type and vartype
        a = dm.createISTPattrs('data', vartype='float')
        a_ans = {'CATDESC': '',
                 'DISPLAY_TYPE': 'time_series',
                 'FIELDNAM': '',
                 'FILLVAL': -1e+31,
                 'FORMAT': 'F18.6',
                 'LABLAXIS': '',
                 'SI_CONVERSION': ' > ',
                 'UNITS': ' ',
                 'VALIDMIN': '',
                 'VALIDMAX': '',
                 'VAR_TYPE': 'data',
                 'DEPEND_0': 'Epoch'}
        self.assertEqual(a, a_ans)

    def test_createISTPattrs_support_data(self):
        """createISTPattrs should give known results and error check for support_data"""
        # test with a datatype type and vartype
        a = dm.createISTPattrs('support_data', vartype='float')
        a_ans = {'CATDESC': '',
                 'FIELDNAM': '',
                 'FORMAT': 'F18.6',
                 'UNITS': ' ',
                 'VAR_TYPE': 'support_data',
                 'DEPEND_0': 'Epoch',
                 'VALIDMIN': '',
                 'VALIDMAX': '',
                 'FILLVAL': -1e+31}
        self.assertEqual(a, a_ans)
        # same for an NRV variable
        a = dm.createISTPattrs('support_data', vartype='float', NRV=True)
        a_ans = {'CATDESC': '',
                 'FIELDNAM': '',
                 'FORMAT': 'F18.6',
                 'UNITS': ' ',
                 'VAR_TYPE': 'support_data'}
        self.assertEqual(a, a_ans)

    def test_createISTPattrs_metadata(self):
        """createISTPattrs should give known results and error check for metadata"""
        # test with a datatype type and vartype
        a = dm.createISTPattrs('metadata', vartype='float')
        a_ans = {'CATDESC': '',
                 'FIELDNAM': '',
                 'FORMAT': 'F18.6',
                 'UNITS': ' ',
                 'VAR_TYPE': 'metadata',
                 'DEPEND_0': 'Epoch',
                 'FILLVAL': -1e+31}
        self.assertEqual(a, a_ans)
        # same for an NRV variable
        a = dm.createISTPattrs('metadata', vartype='float', NRV=True)
        a_ans = {'CATDESC': '',
                 'FIELDNAM': '',
                 'FORMAT': 'F18.6',
                 'UNITS': ' ',
                 'VAR_TYPE': 'metadata'}
        self.assertEqual(a, a_ans)


if __name__ == "__main__":
    unittest.main()
