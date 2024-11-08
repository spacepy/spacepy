#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit test suite for base spacepy

Copyright 2012 Los Alamos National Security, LLC.
"""

import inspect
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
            "test function\n"
            "\n"
            ".. deprecated:: 0.1\n"
            "   pithy message\n"
            "\n"
            "this will test things",
            inspect.getdoc(testfunc))
        with spacepy_testing.assertWarns(self, 'always', r'.*pithy message$',
                                         DeprecationWarning):
            self.assertEqual(2, testfunc(1))

    def testDeprecationNone(self):
        """Test the deprecation decorator with no docstring"""
        @spacepy.deprecated(0.1, 'pithy message')
        def testfunc(x):
            return x + 1
        self.assertEqual(
            ".. deprecated:: 0.1\n"
            "   pithy message",
            inspect.getdoc(testfunc))
        with spacepy_testing.assertWarns(self, 'always', r'pithy message$',
                                         DeprecationWarning):
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
            ".. deprecated:: 0.1\n"
            "   pithy message\n"
            "\n"
            "this will test things",
            inspect.getdoc(testfunc))

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
            ".. deprecated:: 0.1\n"
            "   foo\n"
            "   bar\n"
            "\n"
            "this will test things",
            inspect.getdoc(testfunc))
        with spacepy_testing.assertWarns(self, 'always', r'pithy message$',
                                         DeprecationWarning):
            self.assertEqual(2, testfunc(1))

    def testDeprecationOneParagraph(self):
        """Test the deprecation decorator, docstring is one paragraph"""
        @spacepy.deprecated(0.1, 'pithy message')
        def testfunc(x):
            """test function"""
            return x + 1
        self.assertEqual(
            "test function\n"
            "\n"
            ".. deprecated:: 0.1\n"
            "   pithy message",
            inspect.getdoc(testfunc))


class SpacepyDirTests(unittest.TestCase):
    """Tests on the .spacepy directory and related."""

    def setUp(self):
        self.td = None
        self.old_env = {}
        super(SpacepyDirTests, self).setUp()
        self.old_env = { k: os.environ.get(k, None)
                         for k in ('SPACEPY', 'HOME') }
        self.td = tempfile.mkdtemp()

    def tearDown(self):
        if self.td:
            shutil.rmtree(self.td)
        for k, v in self.old_env.items():
            if v is None:
                if k in os.environ:
                    del os.environ[k]
            else:
                os.environ[k] = v
        super(SpacepyDirTests, self).tearDown()

    def testDotfln(self):
        """Checks DOT_FLN calculations"""
        os.environ['SPACEPY'] = os.path.join(self.td, 'spacepy')
        os.environ['HOME'] = os.path.join(self.td, 'notspacepy')
        self.assertEqual(os.path.join(self.td, 'spacepy', '.spacepy'),
                         spacepy._find_spacepy_dir())
        self.assertTrue(os.path.isdir(os.path.join(self.td, 'spacepy')))
        self.assertEqual(os.path.join(self.td, 'spacepy', '.spacepy'),
                         spacepy._find_spacepy_dir())
        del os.environ['SPACEPY']
        self.assertEqual(os.path.join(self.td, 'notspacepy', '.spacepy'),
                         spacepy._find_spacepy_dir())

    def testDotflnRelative(self):
        """Checks DOT_FLN with a relative path"""
        wd = os.getcwd()
        try:
            os.chdir(self.td)
            cwd = os.path.abspath('.')  # Mac puts you in a different dir
            os.environ['SPACEPY'] = ''
            self.assertEqual(os.path.join(cwd, '.spacepy'),
                             spacepy._find_spacepy_dir())
            os.environ['SPACEPY'] = 'spacepy'
            self.assertEqual(os.path.join(cwd, 'spacepy', '.spacepy'),
                             spacepy._find_spacepy_dir())
            self.assertTrue(os.path.isdir(os.path.join(cwd, 'spacepy')))
        finally:
            os.chdir(wd)

    def testDotflnTilde(self):
        """Checks DOT_FLN with no available environment variables"""
        keys = ['SPACEPY', 'HOME', 'HOMEDRIVE']
        backup = {key: os.environ[key] for key in keys if key in os.environ}
        try:
            for key in keys:
                if key in os.environ:
                    del os.environ[key]
            self.assertEqual(os.path.join(os.path.expanduser('~'), '.spacepy'),
                             spacepy._find_spacepy_dir())
        finally:
            os.environ.update(backup)

    def testNoDotfln(self):
        """Check creating .spacepy"""
        spdir = os.path.join(self.td, 'spacepy')
        os.mkdir(spdir)
        spacepy._populate_spacepy_dir(os.path.join(spdir, '.spacepy'))
        self.assertTrue(os.path.isdir(os.path.join(
            self.td, 'spacepy', '.spacepy', 'data')))

    def testNoDataDir(self):
        """Check creating data directory only"""
        spdir = os.path.join(self.td, 'spacepy')
        os.mkdir(spdir)
        spacepy._populate_spacepy_dir(os.path.join(spdir, '.spacepy'))
        self.assertTrue(os.path.isdir(os.path.join(
            self.td, 'spacepy', '.spacepy', 'data')))

    def testEmptyConfig(self):
        """Treat an empty config file as corrupt"""
        configfile = os.path.join(self.td, 'spacepy.rc')
        open(configfile, 'w').close()
        spacepy._read_config(configfile)
        self.assertIn('enable_old_data_warning', spacepy.config)
        self.assertTrue(os.stat(configfile).st_size > 100)
        spacepy._read_config(spacepy.rcfile)  # Restore the previous config


class SpacepyConfigTests(unittest.TestCase):
    """Tests for _read_config and _write_defaults"""
    def testWriteDefaultsMultipleSections(self):
        """Test _write_defaults where the config file has multiple sections"""
        td = tempfile.mkdtemp()
        configfile = os.path.join(td, 'spacepy.rc')
        expected = ['[spacepy]\n',
                    'some sample text\n',
                    f'#SpacePy {spacepy.__version__} default test entry: value\n',
                    '#test entry: value\n',
                    '[test section]\n']
        try:
            with open(configfile, 'w') as cf:
                cf.write('[spacepy]\n')
                cf.write('some sample text\n')
                cf.write('[test section]\n')
            spacepy._write_defaults(configfile, {"test entry": "value"})
            with open(configfile, 'r') as cf:
                actual_lines = cf.readlines()
            self.assertEqual(expected, actual_lines)
        finally:
            shutil.rmtree(td)

    def testReadConfigError(self):
        """Test generic error cp for _read_config"""
        td = tempfile.mkdtemp()
        configfile = os.path.join(td, 'spacepy.rc')
        try:
            with open(configfile, 'w') as cf:
                cf.write("this text is not in a section\n")
            spacepy._read_config(configfile)
            self.assertIn('enable_old_data_warning', spacepy.config)
        finally:
            shutil.rmtree(td)
            spacepy._read_config(spacepy.rcfile)  # Restore the previous config



if __name__ == '__main__':
    unittest.main()
