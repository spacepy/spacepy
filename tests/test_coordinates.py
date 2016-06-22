
import unittest
import spacepy.coordinates as spc
import glob
import os
import datetime
from numpy import array
import numpy as np
from spacepy.time import Ticktock
try:
    import spacepy.irbempy as ib
except ImportError:
    pass #tests will fail, but won't bring down the entire suite

__all__ = ['coordsTest']


class coordsTest(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        try:
            self.cvals = spc.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        except ImportError:
            pass #tests will fail, but won't bring down the entire suite

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        np.testing.assert_equal([1,1], self.cvals.x)
        np.testing.assert_equal([2,2], self.cvals.y)
        np.testing.assert_equal([4,2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')

    def test_append(self):
        c2 = spc.Coords([[6,7,8],[9,10,11]], 'GEO', 'car')
        actual = self.cvals.append(c2)
        expected = [[1,2,4],[1,2,2],[6,7,8],[9,10,11]]
        np.testing.assert_equal(expected, actual.data.tolist())


if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(coordsTest)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()
