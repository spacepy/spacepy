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


class AssocTests(unittest.TestCase):
    """Tests of association analysis"""

    def testAssocRegress(self):
        """Test entire association analysis

        This is a regression test: values are taken from existing
        implementation and assumed correct.
        """
        #Created from numpy.random.randint(0, 864000, [100])
        t1 = [655963,  70837, 681914, 671054, 120856,  52827, 674833, 727548,
              569743, 168496, 417843, 846365, 124030, 149118, 331591, 499616,
              401455, 484069, 470137, 476552, 105777, 711075, 688800, 544110,
              730248,  90896, 249113, 132259, 261126, 723237, 197802, 236209,
              239956, 231539, 537703, 528239, 354486, 135203, 357335, 683024,
              258661, 316019, 214481, 672555, 212034, 103213, 746137, 111065,
              436613, 286559, 511857, 489026,  32363, 275374, 342386, 192878,
              266724, 116734, 176916, 818199, 842720, 863230,  44519, 439228,
              419808, 179134, 639349, 763944, 605673, 507410, 618672, 828078,
              85643, 438461, 737909,  39394, 708392, 652198, 588335, 670095,
              675992, 580726, 455917, 267198, 708865, 423934, 126758, 537249,
              228845, 364992, 843165, 687482, 162219, 107074, 263547, 272363,
              838316, 574703, 421124, 484203]
        t2 = [757071, 820424, 522664, 399644, 558150, 342836, 834185,   9550,
              278019, 110750, 514961, 473701, 673490, 830021, 762821, 231852,
              662373, 729429, 545901, 409830, 432443, 649515, 446492, 286360,
              364346, 794408, 270865, 291723, 160511, 376046, 483439, 522531,
              92429, 642585, 803893,  61784, 116165, 405721, 565018,   3538,
              815125, 311551, 850973, 629556, 701310, 490674, 183441, 116949,
              388805, 457620, 302912,  75785, 717289, 424186, 370460,  93986,
              194428, 125804,  95628, 382477, 234520,  34429, 568429, 110523,
              519464, 530399, 244645, 345020, 690005, 750812, 237726, 549233,
              297069,  60590, 779392, 120764, 298320, 587738, 141891, 114935,
              585671, 138104, 752052, 585814, 670661, 281514, 148099, 682492,
              660800, 429724, 832390, 536037, 618901, 363413, 257753, 858464,
              674609, 191279, 337199, 193586]
        pop = poppy.PPro(t1, t2)
        lags = [(i - 144) * 120 for i in range(289)]
        pop.assoc(lags, 120)
        expected = [3.0, 2.0, 1.0, 1.0, 3.0, 2.0, 2.0, 4.0, 3.0, 3.0, 3.0, 1.0,
                    0.0, 2.0, 2.0, 2.0, 5.0, 5.0, 5.0, 3.0, 0.0, 2.0, 5.0, 5.0,
                    4.0, 4.0, 5.0, 7.0, 4.0, 1.0, 2.0, 2.0, 3.0, 4.0, 6.0, 4.0,
                    3.0, 3.0, 2.0, 3.0, 4.0, 6.0, 5.0, 2.0, 1.0, 3.0, 4.0, 3.0,
                    3.0, 4.0, 3.0, 1.0, 1.0, 1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 7.0,
                    5.0, 2.0, 4.0, 5.0, 4.0, 7.0, 6.0, 2.0, 3.0, 5.0, 4.0, 3.0,
                    4.0, 5.0, 4.0, 1.0, 2.0, 4.0, 3.0, 2.0, 3.0, 3.0, 1.0, 3.0,
                    5.0, 2.0, 3.0, 4.0, 1.0, 1.0, 3.0, 4.0, 3.0, 2.0, 6.0, 6.0,
                    2.0, 3.0, 3.0, 2.0, 2.0, 2.0, 3.0, 2.0, 3.0, 3.0, 3.0, 4.0,
                    2.0, 4.0, 3.0, 1.0, 1.0, 0.0, 0.0, 2.0, 3.0, 1.0, 1.0, 1.0,
                    1.0, 3.0, 3.0, 2.0, 1.0, 1.0, 1.0, 0.0, 2.0, 5.0, 6.0, 4.0,
                    3.0, 4.0, 3.0, 2.0, 3.0, 4.0, 4.0, 6.0, 5.0, 2.0, 3.0, 3.0,
                    1.0, 2.0, 4.0, 3.0, 3.0, 3.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
                    3.0, 4.0, 4.0, 4.0, 2.0, 4.0, 4.0, 2.0, 1.0, 3.0, 3.0, 1.0,
                    2.0, 6.0, 5.0, 1.0, 2.0, 4.0, 4.0, 2.0, 1.0, 2.0, 2.0, 4.0,
                    4.0, 2.0, 3.0, 5.0, 4.0, 5.0, 7.0, 4.0, 1.0, 1.0, 2.0, 3.0,
                    4.0, 4.0, 4.0, 6.0, 5.0, 3.0, 4.0, 4.0, 3.0, 2.0, 4.0, 3.0,
                    3.0, 3.0, 2.0, 4.0, 4.0, 4.0, 3.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                    4.0, 2.0, 4.0, 6.0, 4.0, 2.0, 0.0, 0.0, 2.0, 3.0, 6.0, 6.0,
                    3.0, 2.0, 1.0, 2.0, 2.0, 1.0, 3.0, 6.0, 5.0, 4.0, 6.0, 5.0,
                    2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 2.0, 1.0, 2.0, 3.0, 2.0, 3.0,
                    4.0, 3.0, 4.0, 3.0, 1.0, 0.0, 2.0, 3.0, 4.0, 6.0, 5.0, 3.0,
                    4.0, 3.0, 3.0, 3.0, 3.0, 2.0, 3.0, 7.0, 5.0, 1.0, 0.0, 3.0,
                    4.0, 6.0, 4.0, 0.0, 2.0, 2.0, 1.0, 2.0, 2.0, 1.0, 0.0, 2.0,
                    3.0]
        self.assertEqual(expected, list(pop.assoc_total))


if __name__ == '__main__':
    unittest.main()
