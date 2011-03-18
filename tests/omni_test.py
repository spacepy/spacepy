#!/usr/bin/env	python
#	-*-	coding:	utf-8	-*-
"""
testing	the	omni	module
"""
__version__	=	"$Revision: 1.2 $,	$Date: 2011/03/18 16:39:02 $"
__author__	=	'Josef	Koller,	Los	Alamos	National	Lab	(jkoller@lanl.gov)'

# -----------------------------------------------------------------------

import unittest
import spacepy
import spacepy.omni as om
import datetime
import glob
import os
import numpy as n
import numpy.testing
from numpy import array

class BigTests(unittest.TestCase):

    def setUp(self):
        self.ticks = spacepy.time.Ticktock(['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO')
        #self.loci = spacepy.coordinates.Coords([[3,0,0],[2,0,0]], 'GEO', 'car')

    def tearDown(self):
        pass

    def test_get_omni(self):
        expected = {'ByIMF': array([-3.        , -3.41666667]),
 'Bz1': array([-5.        , -5.18333333]),
 'Bz2': array([-5.        , -5.18333333]),
 'Bz3': array([-5.        , -5.18333333]),
 'Bz4': array([-5.        , -5.18333333]),
 'Bz5': array([-5.        , -5.18333333]),
 'Bz6': array([-5.        , -5.18333333]),
 'BzIMF': array([-5.        , -5.18333333]),
 'DOY': array([ 33.,  33.]),
 'Dst': array([-75.        , -74.83333333]),
 'G1': array([ 7.69      ,  7.82666667]),
 'G2': array([ 10.59      ,  10.56666667]),
 'G3': array([ 10.87      ,  10.59666666]),
 'Hr': array([ 12.        ,  12.16666667]),
 'Kp': array([ 3.7,  3.7]),
 'Pdyn': array([ 2.2       ,  2.23666667]),
 'Qbits': {'ByIMF': array([2, 2]),
           'BzIMF': array([2, 2]),
           'G1': array([2, 2]),
           'G2': array([2, 2]),
           'G3': array([2, 2]),
           'Pdyn': array([2, 2]),
           'W1': array([0, 0]),
           'W2': array([2, 2]),
           'W3': array([0, 0]),
           'W4': array([0, 0]),
           'W5': array([2, 2]),
           'W6': array([2, 2]),
           'dens': array([2, 2]),
           'velo': array([2, 2])},
 'RDT': array([ 730883.5       ,  730883.50694444]),
 'UTC': [datetime.datetime(2002, 2, 2, 12, 0),
         datetime.datetime(2002, 2, 2, 12, 10)],
 'W1': array([ 1.835     ,  1.81683333]),
 'W2': array([ 1.828     ,  1.80666667]),
 'W3': array([ 1.17      ,  1.17216667]),
 'W4': array([ 2.052 ,  2.0235]),
 'W5': array([ 1.45      ,  1.44783333]),
 'W6': array([ 2.944 ,  2.8965]),
 'Year': array([ 2002.,  2002.]),
 'akp3': array([ 3.95      ,  3.94333333]),
 'dens': array([ 9.  ,  8.95]),
 'ticks': spacepy.time.Ticktock( ['2002-02-02T12:00:00', '2002-02-02T12:10:00'], 'ISO'),
 'velo': array([ 358.        ,  362.66666671])}
        actual = om.get_omni(self.ticks)
        keylist = expected.keys()
        keylist.remove('ticks')
        keylist.remove('UTC')
        keylist.remove('Qbits')
        for key in keylist:
            numpy.testing.assert_almost_equal(expected[key], actual[key], decimal=6)
        for key in actual['Qbits'].keys():
            numpy.testing.assert_almost_equal(expected['Qbits'][key], actual['Qbits'][key])
        self.assertEqual(expected['ticks'].data, actual['ticks'].data)
            
        
# -----------------------------------------------------------------------


if	__name__	==	"__main__":
    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(SimpleFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    ##	suite	=	unittest.TestLoader().loadTestsFromTestCase(tFunctionTests)
    ##	unittest.TextTestRunner(verbosity=2).run(suite)

    unittest.main()





