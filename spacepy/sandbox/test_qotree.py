#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for qotree

Copyright ©2010 Los Alamos National Security, LLC.
"""
import unittest

import numpy

import qotree


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
        qt = qotree.QTree(self.grid)


if __name__ == "__main__":
    unittest.main()
