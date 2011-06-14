#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for borg

Copyright ©2010 Los Alamos National Security, LLC.
"""

import unittest

import numpy
import spacepy.toolbox as tb

class BorgTests(unittest.TestCase):

    def setUp(self):
        super(BorgTests, self).setUp()

    def tearDown(self):
        super(BorgTests, self).tearDown()

    def testSomething(self):
        """example test"""
        self.assertEqual(1, 1)


if __name__ == "__main__":
    unittest.main()
