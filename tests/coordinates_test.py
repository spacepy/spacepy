
import unittest
import spacepy.coordinates as c
import glob
import os
import datetime
from numpy import array
import numpy
from spacepy.time import Ticktock
        

class coordsTest(unittest.TestCase):
    def setUp(self):
        #super(tFunctionTests, self).setUp()
        pass

    def tearDown(self):
        #super(tFunctionTests, self).tearDown()
        pass

    def test_coords(self):
        """Coords should create and do simple conversions"""
        cvals = c.Coords([[1,2,4],[1,2,2]], 'GEO', 'car')
        self.assertEqual(np.array([1,1]), cvals.x)
        self.assertEqual(np.array([2,2]), cvals.y)
        self.assertEqual(np.array([4,2]), cvals.y)
        cvals.ticktock = Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:00:00'], 'ISO') # add ticktock
        newcoord = cvals.convert('GSM', 'sph')
        







if __name__ == "__main__":
    ## suite = unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)

    ## suite = unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ## unittest.TextTestRunner(verbosity=2).run(suite)


    unittest.main()








