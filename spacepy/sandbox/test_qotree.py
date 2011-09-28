#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for qotree

Copyright ©2010 Los Alamos National Security, LLC.
"""
import unittest

import numpy

import qotree

class branchTests(unittest.TestCase):
    def setUp(self):
        super(branchTests, self).setUp()

    def tearDown(self):
        super(branchTests, self).tearDown()

    def test_inputs(self):
        """Branch does some input checking"""
        self.assertRaises(ValueError, qotree.Branch, [1])

    def test_creation(self):
        """create a Branch() and make sure it is correct"""
        br = qotree.Branch([1.2, 3.4, 1.2, 3.4])
        self.assertEqual(br.isLeaf, False)
        numpy.testing.assert_allclose(br.lim, [1.2, 3.4, 1.2, 3.4])



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

    def test_rightdaughter(self):
        """rightdaughter shold give known output"""
        ans = [5, 105, 205, 305, 405, 505, 605, 705]
        for i, ii in enumerate(range(1, 200, 25)):
            self.assertEqual(qotree.rightdaughter(2, ii), ans[i])

    def test_boxes_through_level(self):
        """boxes_through_level should give known results"""
        ans = [1, 5, 21, 85, 341, 1365, 5461]
        for i, ii in enumerate(range(1, 7)):
            self.assertEqual(qotree.boxes_through_level(2, ii), ans[i])

    def test_mother(self):
        """mother should give known results"""
        ans = [0, 7, 13, 19, 25, 32, 38, 44]
        for i, ii in enumerate(range(1, 200, 25)):
            self.assertEqual(qotree.mother(2, ii), ans[i])


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
        self.assertRaises(NotImplementedError, qotree.QTree, self.grid, grid_edges=True)

    def test_creation(self):
        """create a QTree and test its contents (regression)"""
        qt = qotree.QTree(self.grid, max_depth=2)
        self.assertEqual(qt.d, 2)
        self.assertEqual(qt.keys(), [1, 2, 3, 4, 5])
        self.assertEqual(qt.getleftdaughter(1), 2)
        self.assertEqual(qt.getrightdaughter(1), 5)
        self.assertEqual(qt.getboxes_through_level(2), 5)
        self.assertEqual(qt.getmother(5), 1)
        numpy.testing.assert_allclose(qt.grid, self.grid)
        numpy.testing.assert_array_equal(qt.locs, range(150))
        self.assertEqual(qt.max_depth, 2)
        self.assertEqual(qt.npoints, 150)

    def test_creation2(self):
        """create a QTree and test its contents (regression)"""
        qt = qotree.QTree(self.grid, max_depth=3)
        self.assertEqual(qt.d, 2)
        self.assertEqual(qt.keys(), range(1, 22))
        self.assertEqual(qt.max_depth, 3)
        self.assertEqual(qt.npoints, 150)

if __name__ == "__main__":
    unittest.main()
