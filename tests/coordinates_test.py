
import unittest
import spacepy.coordinates as c
import glob
import os
import datetime
from numpy import array
import numpy as n
from spacepy.time import Ticktock
        

class coordsTest(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        self.cvals = c.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""        
        n.testing.assert_equal([1,1], self.cvals.x)
        n.testing.assert_equal([2,2], self.cvals.y)
        n.testing.assert_equal([4,2], self.cvals.z)
        self.cvals.ticks = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = self.cvals.convert('GSM', 'sph')
        
    def test_append(self):
        c2 = c.Coords([[6,7,8],[9,10,11]], 'GEO', 'car')
        actual = self.cvals.append(c2)
        expected = [[1,2,4],[1,2,2],[6,7,8],[9,10,11]]
        n.testing.assert_equal(expected, actual.data.tolist()) 




if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    ## suite = unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)


    unittest.main()








