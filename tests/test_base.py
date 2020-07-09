#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for base spacepy

Copyright 2012 Los Alamos National Security, LLC.
"""

import unittest
import warnings

import spacepy


class SpacepyFuncTests(unittest.TestCase):
    """Tests for functions in core spacepy library"""

    def testDeprecation(self):
        """Test the deprecation decorator"""
        @spacepy.deprecated(0.1, 'pithy message')
        def testfunc(x):
            """
            test function

            this will test things
            """
            return x + 1
        self.assertEqual(
            "\n"
            "            test function\n"
            "\n"
            "            .. deprecated:: 0.1\n"
            "               pithy message\n"
            "\n"
            "            this will test things\n"
            "            ",
            testfunc.__doc__)
        with warnings.catch_warnings(record=True) as w:
            #make sure to catch expected warnings
            warnings.filterwarnings('always', 'pithy message',
                                    DeprecationWarning, '^spacepy')
            self.assertEqual(2, testfunc(1))
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual('pithy message', str(w[0].message))

    def testDeprecationNone(self):
        """Test the deprecation decorator with no docstring"""
        @spacepy.deprecated(0.1, 'pithy message')
        def testfunc(x):
            return x + 1
        self.assertEqual(
            "\n"
            "    .. deprecated:: 0.1\n"
            "       pithy message",
            testfunc.__doc__)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'pithy message',
                                    DeprecationWarning, '^spacepy')
            self.assertEqual(2, testfunc(1))
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual('pithy message', str(w[0].message))

    def testDeprecationDifferentIndent(self):
        """Test the deprecation decorator, first line indented differently"""
        @spacepy.deprecated(0.1, 'pithy message')
        def testfunc(x):
            """test function

            this will test things
            """
            return x + 1
        self.assertEqual(
            "test function\n"
            "\n"
            "            .. deprecated:: 0.1\n"
            "               pithy message\n"
            "\n"
            "            this will test things\n"
            "            ",
            testfunc.__doc__)

    def testDeprecationDifferentMessage(self):
        """Test the deprecation decorator, docstring different from message"""
        @spacepy.deprecated(0.1, 'pithy message', docstring='foo\nbar')
        def testfunc(x):
            """test function

            this will test things
            """
            return x + 1
        self.assertEqual(
            "test function\n"
            "\n"
            "            .. deprecated:: 0.1\n"
            "               foo\n"
            "               bar\n"
            "\n"
            "            this will test things\n"
            "            ",
            testfunc.__doc__)
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', 'pithy message',
                                    DeprecationWarning, '^spacepy')
            self.assertEqual(2, testfunc(1))
        self.assertEqual(1, len(w))
        self.assertEqual(DeprecationWarning, w[0].category)
        self.assertEqual('pithy message', str(w[0].message))


if __name__ == '__main__':
    unittest.main()
