#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for data_assimilation

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

import unittest

import numpy
import spacepy_testing
import spacepy.toolbox as tb
import spacepy.data_assimilation as da

__all__ = ['DataAssimilationTests']

class DataAssimilationTests(unittest.TestCase):

    def setUp(self):
        super(DataAssimilationTests, self).setUp()

    def tearDown(self):
        super(DataAssimilationTests, self).tearDown()

    def testSomething(self):
        """example test"""
        self.assertEqual(1, 1)


if __name__ == "__main__":
    unittest.main()
