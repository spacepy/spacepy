#!/usr/bin/env python

"""Integration tests for omni

Copyright 2019 SpacePy contributors
"""

import unittest

import numpy.testing
import spacepy_testing
import spacepy
import spacepy.omni
import spacepy.time


class OMNIIntegration(unittest.TestCase):
    """Tests of OMNI data assuming up-to-date database"""

    def testFillValues(self):
        """Test that expected fill values are nan"""
        #This is a check that, after reading the CDF from SPDF, fill is
        #replaced with nan to match usage in ViRBO
        #In order to be accurate, it should be run after an update()
        tt = spacepy.time.tickrange('2004-12-02T12:00:00', '2004-12-03', 3./24)
        omni = spacepy.omni.get_omni(tt, dbase='OMNI2hourly')
        expected = spacepy.dmarray([0.305, 0.295, numpy.nan, numpy.nan, 0.295])
        #There is a very small difference resulting from float/double mismatch
        numpy.testing.assert_array_almost_equal(
            omni['Proton_flux_gt_10_MeV'], expected, decimal=8)


if __name__ == "__main__":
    unittest.main()
