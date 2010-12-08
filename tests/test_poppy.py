#!/usr/bin/env python

"""Unit test suite for PoPPy"""

__version__ = "0.0"
__author__ = "Jonathan Niehof <jniehof@lanl.gov>"


import math
import unittest

import numpy.random
import scipy.special
try:
    import poppy
except ImportError:
    from spacepy import poppy
from spacepy import toolbox


class BootstrapTests(unittest.TestCase):
    """Tests of the boots_ci function

    @ivar long_test: Run a very long test
    @type long_test: bool
    """
    
    def __init__(self, *args, **kwargs):
        super(BootstrapTests, self).__init__(*args, **kwargs)
        self.long_test = False

    def testGaussian(self):
        """Check output on a Gaussian"""
        gauss = lambda x: math.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                          (5 * math.sqrt(2 * math.pi))
        mean = lambda x: sum(x) / len(x)
        if self.long_test:
            gauss_values = toolbox.dist_to_list(gauss, 5000)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(gauss_values, 50000, 95.0, mean)
            #standard error on the mean for 5000 samples, standard dev 5 is
            #0.0707, 95% confidence interval is 1.96 standard dev or
            #0.13859292911256332 (+/-)
            self.assertAlmostEqual(-0.138373060749, ci_low, places=10)
            self.assertAlmostEqual(0.137630863579, ci_high, places=10)
        else:
            gauss_values = toolbox.dist_to_list(gauss, 100)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(gauss_values, 100, 95.0, mean)
            #for 100 samples, stderr is 0.5, 95% is +/- 0.98
            self.assertAlmostEqual(-1.03977357727, ci_low, places=10)
            self.assertAlmostEqual(0.914472603387, ci_high, places=10)

    def testPoisson(self):
        """Check output on a continuous Poisson-like distribution"""
        l = 4 #Expected number of counts
        poiss = lambda k: l **k * math.exp(-l) / scipy.special.gamma(k + 1)
        mean = lambda x: sum(x) / len(x)
        if self.long_test:
            poiss_values = toolbox.dist_to_list(poiss, 5000, 0, 100)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(poiss_values, 50000, 95.0, mean)
            #IF this were normal (which it isn't), expected confidence interval
            #3.94456 - 4.05544, in reality should be skewed right
            self.assertAlmostEqual(3.97171801164, ci_low, places=10)
            self.assertAlmostEqual(4.08051711605, ci_high, places=10)
        else:
            poiss_values = toolbox.dist_to_list(poiss, 100, 0, 100)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(poiss_values, 100, 95.0, mean)
            #'Expected' 3.608 - 4.392
            self.assertAlmostEqual(3.57505503117, ci_low, places=10)
            self.assertAlmostEqual(4.4100232162, ci_high, places=10)


if __name__ == '__main__':
    unittest.main()
