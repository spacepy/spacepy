#!/usr/bin/env python
from __future__ import division

import unittest
import datetime
from spacepy import datamodel

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

    def test_pickle(self):
        """things should pickle and unpickle"""
        dat = datamodel.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})
        tmp = pickle.dumps(dat)
        for i, val in enumerate(dat):
            self.assertEqual(pickle.loads(tmp)[i], val)
        self.assertEqual(pickle.loads(tmp).attrs, dat.attrs)

    def test_attrs(self):
        """dmarray can only have .attrs"""
        self.assertRaises(TypeError, datamodel.dmarray, [1,2,3], setme = 123 )

    def test_attrs(self):
        """The only attribute the can be set is attrs"""
        self.assertRaises(TypeError, datamodel.dmarray, [1,2,3], bbb=23)
        dat = datamodel.dmarray([1,2,3,4], attrs={'a':'a', 'b':'b'})
        try:
            dat.bbb = 'someval'
            self.assertFail()
        except TypeError:
            self.assertTrue(True)

    #def test_SpaceData(self):
    #    """this is an abstract class that has some exceptions"""
    #    self.assertRaises(ValueError, datamodel.SpaceData)

if __name__ == "__main__":
    unittest.main()
