#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for libspacepy

Copyright 2014 Los Alamos National Security, LLC.
"""

import unittest
import warnings

import spacepy.lib


class SpacepyLibTests(unittest.TestCase):
    """Tests for libspacepy"""

    def setUp(self):
        warnings.simplefilter('always')

    def testExists(self):
        """Make sure we're finding libspacepy"""
        self.assertTrue(spacepy.lib.have_libspacepy)


if __name__ == '__main__':
    unittest.main()
