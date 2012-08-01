#!/usr/bin/env python
"""functions useful for finding multiple instances of a condition"""

#Jon Niehof, 2012 Aug 1

import unittest

import numpy
import numpy.testing

def find_ascending(x, val):
    """In a list x, find indexes AFTER or AT point where x
    transitions from below val to above val
    
    i.e., will return index OF val if present, otherwise of
    first number greater than val.
    """
    return numpy.nonzero(numpy.diff((x >= val).astype(numpy.int)) > 0)[0] + 1

def find_descending(x, val):
    """In a list x, find indexes AFTER or AT point where x
    transitions from above x to below x.
    
    i.e., will return index OF val if present, otherwise of
    first number less than val.
    """
    return numpy.nonzero(numpy.diff((x <= val).astype(numpy.int)) > 0)[0] + 1

class TestStuff(unittest.TestCase):
    def testAscending(self):
        """Test ascending"""
        x = numpy.array([1, 2, 3, 4, 5, 4, 3, 2] * 10)
        numpy.testing.assert_array_equal(list(range(3, 83,8)),
                                         find_ascending(x, 3.5))
        numpy.testing.assert_array_equal(list(range(3, 83,8)),
                                         find_ascending(x, 4))

    def testDescending(self):
        """Test ascending"""
        x = numpy.array([1, 2, 3, 4, 5, 4, 3, 2] * 10)
        numpy.testing.assert_array_equal(list(range(6, 83,8)),
                                         find_descending(x, 3.5))
        numpy.testing.assert_array_equal(list(range(5, 83,8)),
                                         find_descending(x, 4))

#These NAN do not work properly
#compares to nan always false, can we use that?
    def testAscendingNAN(self):
        """Test nan handling"""
        x = numpy.array([1, 2, 3, 4, 5, 4, 3, 2] * 10, dtype=numpy.float)
        x[5] = numpy.nan
        numpy.testing.assert_array_equal(list(range(3, 83,8)),
                                         find_ascending(x, 3.5))
        numpy.testing.assert_array_equal(list(range(3, 83,8)),
                                         find_ascending(x, 4))  

    def testDescendingNAN(self):
        """Test nan handling on descending"""
        x = numpy.array([1, 2, 3, 4, 5, 4, 3, 2] * 10, dtype=numpy.float)
        x[5] = numpy.nan
        numpy.testing.assert_array_equal(list(range(6, 83,8)),
                                         find_descending(x, 3.5))
        numpy.testing.assert_array_equal(list(range(5, 83,8)),
                                         find_descending(x, 4))


if __name__ == '__main__':
    unittest.main()
