#!/usr/bin/env python
from __future__ import division

import datetime
import os
import tempfile
import unittest

from spacepy import datamodel
import numpy as np

try:
    import cPickle as pickle
except:
    import pickle

__author__ = 'Brian Larsen <balarsen@lanl.gov>'

class datamodelTests(unittest.TestCase):
    def setUp(self):
        super(datamodelTests, self).setUp()
    def tearDown(self):
        super(datamodelTests, self).tearDown()

    def test_creation_dmarray(self):
        """When a dmarray is created it should have attrs empty or not"""
        data = datamodel.dmarray([1,2,3], attrs={'coord':'GSM'})
        self.assertTrue(hasattr(data, 'attrs'))
        self.assertEqual(data.attrs['coord'], 'GSM')
        data = datamodel.dmarray([1,2,3])
        self.assertTrue(hasattr(data, 'attrs'))
        self.assertEqual(data.attrs, {})
        data2 = datamodel.dmarray([1,2,3], attrs={'coord':'GSM'})
        self.assertEqual(data.attrs, {})
        self.assertEqual(data2.attrs, {'coord':'GSM'})

    def test_pickle_dumps(self):
        """things should pickle and unpickle"""
        dat = datamodel.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})
        tmp = pickle.dumps(dat)
        for i, val in enumerate(dat):
            self.assertEqual(pickle.loads(tmp)[i], val)
        self.assertEqual(pickle.loads(tmp).attrs, dat.attrs)

    def test_pickle_dump(self):
        """things should pickle and unpickle to a file"""
        dat = datamodel.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})
        fname = None
        try:
            with tempfile.NamedTemporaryFile(delete=False) as fp:
                fname = fp.name
                pickle.dump(dat, fp)
            with open(fname, 'rb') as fp:
                dat2 = pickle.load(fp)
        finally:
            if fname != None:
                os.remove(fname)
        np.testing.assert_array_almost_equal(dat, dat2)
        self.assertEqual(dat.attrs, dat2.attrs)

    def test_attrs_only(self):
        """dmarray can only have .attrs"""
        self.assertRaises(TypeError, datamodel.dmarray, [1,2,3], setme = 123 )

    def test_attrs(self):
        """The only attribute the can be set is attrs"""
        self.assertRaises(TypeError, datamodel.dmarray, [1,2,3], bbb=23)
        dat = datamodel.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})
        try:
            dat.bbb = 'someval'
        except TypeError:
            pass
        else:
            self.fail(
                'Assigning to arbitrary Python attribute should raise TypeError')


if __name__ == "__main__":
    unittest.main()
