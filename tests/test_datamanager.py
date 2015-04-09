#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's data manager

Copyright 2015 University System of New Hampshire
"""

import datetime
import unittest

import numpy.testing

import spacepy.datamanager


__all__ = ["RePathTests"]


class RePathTests(unittest.TestCase):
    def test_path_split(self):
        """Verify the simple splitting works"""
        path = "foo/bar/baz"
        self.assertEqual(["foo", "bar", "baz"],
                         spacepy.datamanager.RePath.path_split(path))

    def test_path_slice(self):
        """Verify slicing on a path"""
        path = "foo/bar/baz"
        self.assertEqual("bar",
                         spacepy.datamanager.RePath.path_slice(path, 1))
        self.assertEqual(["foo", "baz"],
                         spacepy.datamanager.RePath.path_slice(
                             path, 0, step=2))
        self.assertEqual(["bar", "baz"],
                         spacepy.datamanager.RePath.path_slice(
                             path, 1, 3))


class DataManagerFunctionTests(unittest.TestCase):
    def test_insert_fill(self):
        e = [1, 2, 3, 5, 6, 7]
        d = [1, 2, 3, 3, 2, 1]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5, 6, 7]
        expected_df = [1, 2, 3, -1, 3, 2, 1]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1, 2], [3, 4], [4, 3], [2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2], [3, 4], [-1, -1], [4, 3], [2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2, -1, 3, 4], [4, 3, -1, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 4, 5]
        d = [[1., 2, 3, 4], [4, 3, 2, 1]] #has to be float for nan output
        ef, df = spacepy.datamanager.insert_fill(e, d)
        expected_ef = [1, 2, 3, 4, 5]
        expected_df = [[1, 2, numpy.nan, 3, 4], [4, 3, numpy.nan, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [1, 2, 3, 4]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [1, 2, 3, 4]
        expected_df = [[1, 2, 3, 4], [4, 3, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, 1, i, 0) for i in (1, 2, 4, 5)]
        d = [[1, 2, 3, 4], [4, 3, 2, 1]]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, 1, i, 0)
                       for i in (1, 2, 3, 4, 5)]
        expected_df = [[1, 2, -1, 3, 4], [4, 3, -1, 2, 1]]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, 1, 0, 0, 0, i)
             for i in range(1, 4) + range(5, 8) + range(9, 12)]
        d = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, 1, 0, 0, 0, i)
                       for i in range(1, 12)]
        expected_df = [1, 2, 3, -1, 5, 6, 7, -1, 9, 10, 11]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)

        e = [datetime.datetime(2000, 1, i, 0, 0, 0, 0)
             for i in range(1, 4) + range(5, 8) + range(9, 12)]
        d = [1, 2, 3, 5, 6, 7, 9, 10, 11]
        ef, df = spacepy.datamanager.insert_fill(e, d, -1)
        expected_ef = [datetime.datetime(2000, 1, i, 0, 0, 0, 0)
                       for i in range(1, 12)]
        expected_df = [1, 2, 3, -1, 5, 6, 7, -1, 9, 10, 11]
        numpy.testing.assert_array_equal(ef, expected_ef)
        numpy.testing.assert_array_equal(df, expected_df)


if __name__ == "__main__":
    unittest.main()
