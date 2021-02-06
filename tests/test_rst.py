#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for toolbox

Copyright 2010-2012 Los Alamos National Security, LLC.
"""
import unittest

import spacepy_testing
import spacepy.rst as rst

__all__ = ['RSTTests', ]

class RSTTests(unittest.TestCase):

    def setUp(self):
        super(RSTTests, self).setUp()
        self.list1 = ['list1', 'list2', 'list3']
        self.list2 = [ ['list11', 'list21', 'list31'],
                       ['list12', 'list22', 'list32'],
                       ['list13', 'list23', 'list33'] ]

    def tearDown(self):
        super(RSTTests, self).tearDown()

    def test_listToEnumerate(self):
        """listToEnumerate"""
        self.assertEqual('1. list1\n2. list2\n3. list3\n\n', rst.listToEnumerate(self.list1))

    def test_listToList(self):
        """listToList"""
        self.assertEqual('- list1\n- list2\n- list3\n\n', rst.listToList(self.list1))

    def test_listToTable(self):
        """listToTable"""
        self.assertEqual('.. csv-table:: \n\t:header: \n\t:widths: 6, 6, 6\n\n\tlist11, list21, list31\n\tlist12, list22, list32\n\tlist13, list23, list33\n',
                         rst.listToTable(self.list2))

    def test_strToHeading(self):
        """strToHeading"""
        ans = ['hello world\n===========\n', 'hello world\n-----------\n',
               'hello world\n~~~~~~~~~~~\n', '-----------\nhello world\n-----------\n',
               '===========\nhello world\n===========\n']
        levels = [0,1,2,-1,-2]
        for ii, level in enumerate(levels):
            self.assertEqual(ans[ii], rst.strToHeading('hello world', level) )
        self.assertRaises(ValueError, rst.strToHeading, 'hello world', 999)


if __name__ == "__main__":
    unittest.main()
