#!/usr/bin/env python

"""Test of test-suite helper functions

Tests to ensure that the more involved support in spacepy_testing works.
"""

import unittest
import warnings

import spacepy_testing


class SPTestWarning(Warning):
    """Used to have a unique type of warning"""
    pass

class SpacepyTestingTests(unittest.TestCase):
    """Simple tests for the testing support"""

    def testassertWarnsSimple(self):
        """Simple test of assertWarns"""
        # This is most similar to how it would normally be used
        with spacepy_testing.assertWarns(
                self, 'always', r'f.*', SPTestWarning):
            warnings.warn('foo', SPTestWarning)

    def testassertWarnsfilter(self):
        """assertWarns only acts on desired warnings"""
        with warnings.catch_warnings(record=True) as cm:
            with spacepy_testing.assertWarns(
                    self, 'always', r'f.*', SPTestWarning):
                warnings.warn('foo', SPTestWarning)
                warnings.warn('foo')
                warnings.warn('bar', SPTestWarning)
        # Verify that the other two warnings were issued as expected
        self.assertEqual(2, len(cm))
        for w in cm:
            self.assertIn(w.category, (UserWarning, SPTestWarning))
            if w.category is UserWarning:
                self.assertEqual('foo', str(w.message))
            else:
                self.assertEqual('bar', str(w.message))

    def testassertWarnsTwice(self):
        """Issue warning twice, test fails"""
        with self.assertRaises(AssertionError):
            with spacepy_testing.assertWarns(
                self, 'always', r'f.*', SPTestWarning):
                warnings.warn('foo', SPTestWarning)
                warnings.warn('foo2', SPTestWarning)

    def testassertDoesntWarnSimple(self):
        """Simple test of doesn't warn"""
        with spacepy_testing.assertDoesntWarn(
                self, 'always', r'f.*', SPTestWarning):
            pass

    def testassertDoesntWarnIssued(self):
        """Issue warning, test fails"""
        with self.assertRaises(AssertionError):
            with spacepy_testing.assertDoesntWarn(
                self, 'always', r'f.*', SPTestWarning):
                warnings.warn('foo', SPTestWarning)

    def testassertDoesntWarnfilter(self):
        """assertDoesntWarn only acts on desired warnings"""
        with warnings.catch_warnings(record=True) as cm:
            with spacepy_testing.assertDoesntWarn(
                    self, 'always', r'f.*', SPTestWarning):
                warnings.warn('foo')
                warnings.warn('bar', SPTestWarning)
        # Verify that the other two warnings were issued as expected
        self.assertEqual(2, len(cm))
        for w in cm:
            self.assertIn(w.category, (UserWarning, SPTestWarning))
            if w.category is UserWarning:
                self.assertEqual('foo', str(w.message))
            else:
                self.assertEqual('bar', str(w.message))


if __name__ == '__main__':
    unittest.main()
