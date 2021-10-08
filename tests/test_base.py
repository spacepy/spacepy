#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for base spacepy

Copyright 2012 Los Alamos National Security, LLC.
"""

import os
import shutil
import tempfile
import unittest
import warnings

import spacepy_testing
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
        with spacepy_testing.assertWarns(self, 'always', r'pithy message$',
                                         DeprecationWarning, r'spacepy$'):
            self.assertEqual(2, testfunc(1))

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
        with spacepy_testing.assertWarns(self, 'always', r'pithy message$',
                                         DeprecationWarning, r'spacepy$'):
            self.assertEqual(2, testfunc(1))

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
        with spacepy_testing.assertWarns(self, 'always', r'pithy message$',
                                         DeprecationWarning, r'spacepy$'):
            self.assertEqual(2, testfunc(1))

    def testDotfln(self):
        """Checks DOT_FLN calculations"""
        old_env = { k: os.environ.get(k, None) for k in ('SPACEPY', 'HOME') }
        td = None
        try:
            td = tempfile.mkdtemp()
            os.environ['SPACEPY'] = os.path.join(td, 'spacepy')
            os.environ['HOME'] = os.path.join(td, 'notspacepy')
            self.assertEqual(os.path.join(td, 'spacepy', '.spacepy'),
                             spacepy._find_spacepy_dir())
            del os.environ['SPACEPY']
            self.assertEqual(os.path.join(td, 'notspacepy', '.spacepy'),
                             spacepy._find_spacepy_dir())
        finally:
            if td:
                shutil.rmtree(td)
            for k, v in old_env.items():
                if v is None:
                    if k in os.environ:
                        del os.environ[k]
                else:
                    os.environ[k] = v


if __name__ == '__main__':
    unittest.main()
