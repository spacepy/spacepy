#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for qotree

Copyright ©2010 Los Alamos National Security, LLC.
"""
import unittest

import numpy

import qotree


class staticTests(unittest.TestCase):
    def setUp(self):
        super(staticTests, self).setUp()
        self.x = numpy.linspace(1, 100, 150)
        self.y = numpy.linspace(1, 100, 150)
        self.grid = numpy.vstack((self.x, self.y))

    def tearDown(self):
        super(staticTests, self).tearDown()


    def test_leftdaughter(self):
        """leftdaughter shold give known output"""
        ans = [2, 102, 202, 302, 402, 502, 602, 702]
        for i, ii in enumerate(range(1, 200, 25)):
            self.assertEqual(qotree.leftdaughter(2, ii), ans[i])


class qotreeTests(unittest.TestCase):
    def setUp(self):
        super(qotreeTests, self).setUp()
        self.x = numpy.linspace(1, 100, 150)
        self.y = numpy.linspace(1, 100, 150)
        self.grid = numpy.vstack((self.x, self.y))

    def tearDown(self):
        super(qotreeTests, self).tearDown()


    def test_QTree_imputs(self):
        """QTree does some input checking"""
        self.assertRaises(NotImplementedError, qotree.QTree, numpy.vstack((self.grid, self.x)))

    def test_creation(self):
        """create a QTree and test its contents (regression)"""
        qt = qotree.QTree(self.grid, max_depth=2)
        self.assertEqual(qt.d, 2)
        self.assertEqual(qt.keys(), [1, 2, 3, 4, 5])
        self.assertEqual(qt.getleftdaughter(1), 2)
        self.assertEqual(qt.getrightdaughter(1), 5)
        


if __name__ == "__main__":
    unittest.main()
