#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
test the LANLstar module

Copyright 2011-2012 Los Alamos National Security, LLC.
"""

import os.path
import sys
import unittest

import spacepy_testing
import spacepy
import spacepy.LANLstar as sl
import numpy
from numpy import array, hstack

__all__ = ['LANLStarFunctionsTest', 'lanlstarTest', 'lanlstarTestMulti']


class LANLStarFunctionsTest(unittest.TestCase):
    """Tests of simple support functions for LANLStar"""

    def testNetPath(self):
        """Get the path to a .net file"""
        #okay, this is stupid since it's really the same code,
        #but at least it checks for syntax errors....
        self.assertEqual(
            os.path.join(
            os.path.split(spacepy.__file__)[0], 'data', 'LANLstar', 'Lmax_T89.txt'),
            sl._get_net_path('Lmax_T89.txt'))
        try:
            q = sl._get_net_path('sample.txt')
        except RuntimeError:
            self.assertEqual(
                'Could not find neural network file sample.txt',
                str(sys.exc_info()[1]))


class lanlstarTest(unittest.TestCase):
    nel = 1
    def setUp(self):
        self.dat = {
            'Kp'     : array([2.7    ]*self.nel),
            'Dst'    : array([7.7777 ]*self.nel),
            'dens'   : array([4.1011 ]*self.nel),
            'velo'   : array([400.101]*self.nel),
            'Pdyn'   : array([4.1011 ]*self.nel),
            'ByIMF'  : array([3.7244 ]*self.nel),
            'BzIMF'  : array([-0.1266]*self.nel),
            'G1'     : array([1.02966]*self.nel),
            'G2'     : array([0.54933]*self.nel),
            'G3'     : array([0.81399]*self.nel),
            'W1'     : array([0.12244]*self.nel),
            'W2'     : array([0.2514 ]*self.nel),
            'W3'     : array([0.0892 ]*self.nel),
            'W4'     : array([0.0478 ]*self.nel),
            'W5'     : array([0.2258 ]*self.nel),
            'W6'     : array([1.0461 ]*self.nel),
            'Year'   : array([1996   ]*self.nel),
            'DOY'    : array([6      ]*self.nel),
            'Hr'     : array([1.2444 ]*self.nel),
            'Lm'     : array([4.9360 ]*self.nel),
            'Bmirr'  : array([315.620]*self.nel),
            'rGSM'   : array([4.8341 ]*self.nel),
            'lonGSM' : array([-40.266]*self.nel),
            'latGSM' : array([36.4469]*self.nel),
            'PA'     : array([57.3875]*self.nel),
            'SMx'    : array([3.9783 ]*self.nel),
            'SMy'    : array([-2.51335]*self.nel),
            'SMz'    : array([1.106617]*self.nel)}
                                                                                
    def test_get_lanlstar(self):
        expected_lstar = {'OPDYN'   : array([4.7171]*self.nel),
                          'OPQUIET' : array([4.6673]*self.nel),
                          'T01QUIET': array([4.8427]*self.nel),
                          'T01STORM': array([4.8669]*self.nel),
                          'T89'     : array([4.5187]*self.nel),
                          'T96'     : array([4.6439]*self.nel),
                          'T05'     : array([4.7174]*self.nel),
                          'RAMSCB'  : array([5.9610]*self.nel)}

        Bmodels = ['OPDYN', 'OPQUIET', 'T01QUIET', 'T01STORM', 'T89', 'T96', 'T05', 'RAMSCB']
        actual = sl.LANLstar(self.dat, Bmodels)
        for key in Bmodels:
            numpy.testing.assert_almost_equal(expected_lstar[key], actual[key], decimal=4)

    def test_get_lanlmax(self):
        expected_lmax = {'OPDYN'   : array([10.6278]*self.nel),
                          'OPQUIET' : array([9.3352]*self.nel),
                          'T01QUIET': array([10.0538]*self.nel),
                          'T01STORM': array([9.9300]*self.nel),
                          'T89'     : array([8.2888]*self.nel),
                          'T96'     : array([9.2410]*self.nel),
                          'T05'    : array([9.9295]*self.nel)}

        Bmodels = ['OPDYN', 'OPQUIET', 'T01QUIET', 'T01STORM', 'T89', 'T96', 'T05']
        actual = sl.LANLmax(self.dat, Bmodels)
        for key in Bmodels:
            numpy.testing.assert_almost_equal(expected_lmax[key], actual[key], decimal=4)
            
    def test_get_lanlmax_G(self):
        if self.nel > 1:
            self.dat['G'] = numpy.repeat([hstack([self.dat['G1'][0], self.dat['G2'][0], self.dat['G3'][0]])], self.nel, axis=0)
            self.dat['W'] = numpy.repeat([hstack([self.dat['W1'][0], self.dat['W2'][0], self.dat['W3'][0],
                                          self.dat['W4'][0], self.dat['W5'][0], self.dat['W6'][0]])], self.nel, axis=0)
        else:
            self.dat['G'] = hstack([self.dat['G1'], self.dat['G2'], self.dat['G3']])
            self.dat['W'] = hstack([self.dat['W1'], self.dat['W2'], self.dat['W3'], self.dat['W4'], self.dat['W5'], self.dat['W6']])
        for i in range(1,4): del self.dat['G{0}'.format(i)]
        for i in range(1,7): del self.dat['W{0}'.format(i)]
        expected_lmax = {'T01QUIET' : array([10.0538]*self.nel),
                          'T01STORM': array([9.9300]*self.nel),
                          'T05'     : array([9.9295]*self.nel)}

        Bmodels = ['T01QUIET','T01STORM','T05']
        actual = sl.LANLmax(self.dat, Bmodels)
        for key in Bmodels:
            numpy.testing.assert_almost_equal(expected_lmax[key], actual[key], decimal=4)


class lanlstarTestMulti(lanlstarTest):
    nel = 3


if __name__=="__main__":
    unittest.main()
