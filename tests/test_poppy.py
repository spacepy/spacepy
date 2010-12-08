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
            self.assertAlmostEqual(-0.13950902739218726, ci_low, places=12)
            self.assertAlmostEqual(0.13799083192825787, ci_high, places=12)
        else:
            gauss_values = toolbox.dist_to_list(gauss, 100)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(gauss_values, 100, 95.0, mean)
            #for 100 samples, stderr is 0.5, 95% is +/- 0.98
            self.assertAlmostEqual(-1.0499918089210114, ci_low, places=12)
            self.assertAlmostEqual(0.88439433514210331, ci_high, places=12)

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
            self.assertAlmostEqual(3.9711311149820565, ci_low, places=12)
            self.assertAlmostEqual(4.0806294481588035, ci_high, places=12)
        else:
            poiss_values = toolbox.dist_to_list(poiss, 100, 0, 100)
            numpy.random.seed(0)
            (ci_low, ci_high) = poppy.boots_ci(poiss_values, 100, 95.0, mean)
            #'Expected' 3.608 - 4.392
            self.assertAlmostEqual(3.6309861974417017, ci_low, places=12)
            self.assertAlmostEqual(4.345026226622057, ci_high, places=12)


if __name__ == '__main__':
    unittest.main()
