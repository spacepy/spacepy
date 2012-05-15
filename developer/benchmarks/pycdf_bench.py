

from __future__ import division

import matplotlib.pyplot as plt 
import os
import datetime

import numpy as np
import spacepy.toolbox as tb

# this is a git access module useful in here.
import dulwich
from benchmarker import Benchmarker

import spacepy.pycdf as cdf

#==============================================================================
# utility stuff from test_pycdf
#==============================================================================

class est_tz(datetime.tzinfo):
    """Eastern Standard timezone (no daylight time"""

    def utcoffset(self, dt):
        """Offset from UTC"""
        return datetime.timedelta(hours=-5)

    def dst(self, dt):
        """Minute offset for DST"""
        return datetime.timedelta(0)

    def tzname(self, dt):
        """Name of this time zone"""
        return 'EST'



#==============================================================================
# test a function vs numpy
#==============================================================================

three_d = [[[90, 90, 96, 71, 90], [29, 18, 90, 78, 51],
           [14, 29, 41, 25, 50], [73, 59, 83, 92, 24],
           [10, 1, 4, 61, 54]],
          [[40, 8, 0, 28, 47], [3, 98, 28, 9, 38],
           [34, 95, 7, 87, 9], [11, 73, 71, 54, 69],
           [42, 75, 82, 16, 73]],
          [[88, 40, 5, 69, 41], [35, 15, 32, 68, 8],
           [68, 74, 6, 30, 9], [86, 48, 52, 49, 100],
           [8, 35, 26, 16, 61]],
          [[49, 81, 57, 37, 98], [54, 64, 28, 21, 17],
           [73, 100, 90, 8, 25], [40, 75, 52, 41, 40],
           [42, 72, 55, 16, 39]],
          [[24, 38, 26, 85, 25], [5, 98, 63, 29, 33],
           [91, 100, 17, 85, 9], [59, 50, 50, 41, 82],
           [21, 45, 65, 51, 90]]]

four_d = [[[[14, 84, 79, 74, 45], [39, 47, 93, 32, 59],
            [15, 47, 1, 84, 44], [13, 43, 13, 88, 3]],
           [[65, 75, 36, 90, 93], [64, 36, 59, 39, 42],
            [59, 85, 21, 88, 61], [64, 29, 62, 33, 35]],
           [[46, 69, 3, 50, 44], [86, 15, 32, 17, 51],
            [79, 20, 29, 10, 55], [29, 10, 79, 7, 58]]],
          [[[20, 76, 81, 40, 85], [44, 56, 5, 83, 32],
            [34, 88, 23, 57, 74], [24, 55, 83, 39, 60]],
           [[79, 56, 5, 98, 29], [28, 50, 77, 33, 45],
            [38, 82, 82, 28, 97], [42, 14, 56, 48, 38]],
           [[58, 27, 38, 43, 25], [72, 91, 85, 44, 43],
            [17, 57, 91, 19, 35], [98, 62, 61, 14, 60]]]]

one_d = [1, 2, 3, 4]

two_d = [[6, 7, 48, 81], [61, 67, 90, 99], [71, 96, 58, 85],
         [35, 31, 71, 73], [77, 41, 71, 92], [74, 89, 94, 64],
         [64, 30, 66, 94]]

zero_d = 1

zero_d_n = np.require(zero_d, requirements='C')

three_d_n = np.require(three_d, requirements='C')
four_d_n = np.require(four_d, requirements='C')
one_d_n = np.require(one_d, requirements='C')
two_d_n = np.require(two_d, requirements='C')


def testFlipMajority(zero_d, one_d, two_d, three_d, four_d):
    """Python"""
    #Code to generate this 5x5x5:
    #[[[random.randint(0,100) for i in range(5)]
    #  for j in range(5)] for k in range(5)]

    flipped = cdf._pycdf._Hyperslice.flip_majority(three_d)

    #[[[[random.randint(0,100) for i in range(5)]
    #  for j in range(4)] for k in range(3)] for l in range(2)]

    flipped = cdf._pycdf._Hyperslice.flip_majority(four_d)

    zero_d = 1
    flipped = cdf._pycdf._Hyperslice.flip_majority(zero_d)

    flipped = cdf._pycdf._Hyperslice.flip_majority(one_d)


    flipped = cdf._pycdf._Hyperslice.flip_majority(two_d)

def testFlipMajority_numpy(zero_d, one_d, two_d, three_d, four_d):
    """Numpy"""
    #Code to generate this 5x5x5:
    #[[[random.randint(0,100) for i in range(5)]
    #  for j in range(5)] for k in range(5)]

    flipped = np.require(three_d, requirements='F')

    #[[[[random.randint(0,100) for i in range(5)]
    #  for j in range(4)] for k in range(3)] for l in range(2)]

    flipped = np.require(four_d, requirements='F')

    flipped = np.require(zero_d, requirements='F')

    flipped = np.require(one_d, requirements='F')

    flipped = np.require(two_d, requirements='F')


loop = 1 # not used in ths test
for bm in Benchmarker(width=25, cycle=50, extra=2):
    bm.run(testFlipMajority, zero_d, one_d, two_d, three_d, four_d)
    bm.run(testFlipMajority_numpy, zero_d_n, one_d_n, two_d_n, three_d_n, four_d_n)


#f6bdbf9106ec3bacdf0917b6ecf0f3d8b168a0dc   (rbsp2)
### Average of 50 (=54-2*2)     user       sys     total      real
#Python                      0.0004    0.0000    0.0004    0.0002
#Numpy                       0.0002    0.0000    0.0002    0.0001
#
### Ranking                    real
#Numpy                       0.0001 (100.0%) *************************
#Python                      0.0002 ( 38.4%) **********
#
### Ratio Matrix               real    [01]    [02]
#[01] Numpy                  0.0001  100.0%  260.1%
#[02] Python                 0.0002   38.4%  100.0%



#==============================================================================
# Just bench some of the functions, taken from the tests
#==============================================================================
