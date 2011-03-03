#!/usr/bin/env python
from __future__ import division

import unittest
import datetime
from spacepy import datamodel

__author__ = 'Brian Larsen <balarsen@lanl.gov>'

class datamodelTests(unittest.TestCase):
    def setUp(self):
        super(datamodelTests, self).setUp()
    def tearDown(self):
        super(datamodelTests, self).tearDown()

    def test_creation_dmarray(self):
        """When a dmarray is created it should have attrs empy or not"""
        data = datamodel.dmarray([1,2,3], attrs={'coord':'GSM'})
        self.assertTrue(hasattr(data, 'attrs'))
        self.assertEqual(data.attrs['coord'], 'GSM')
        data = datamodel.dmarray([1,2,3])
        self.assertTrue(hasattr(data, 'attrs'))

    def test_SpaceData(self):
        """this is an abstract class that has some exceptions"""
        self.assertRaises(ValueError, datamodel.SpaceData)

if __name__ == "__main__":
    unittest.main()
