#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's PyBats subpackage.

Copyright 2015 University of Michigan
"""

import unittest

import spacepy.pybats as pb
import spacepy.pybats.bats as pbs

__all__ = ['TestIdlFile']


class TestIdlFile(unittest.TestCase):
    '''
    Test the class :class:`spacepy.pybats.IdlFile` for different output
    file types and formats (ascii and binary).
    '''

    varnames = 'x z rho ux uy uz bx by bz p b1x b1y b1z e jx jy jz'
    units='R R Mp/cc km/s km/s km/s nT nT nT nPa nT nT nT J/m3 uA/m2 uA/m2 uA/m2'
    knownMhdUnits = dict(zip(varnames.split(), units.split()))
    knownMhdXmax = 31.0
    knownMhdXmin = -220.0
    knownMhdZlim = 124.0
    
    def testBinary(self):
        # Open file:
        mhd = pb.IdlFile('data/pybats_test/y0_binary.out')

        # Test units are loaded correctly:
        for v in mhd:
            if v not in self.varnames: continue
            self.assertEqual(self.knownMhdUnits[v], mhd[v].attrs['units'])

        # Test values inside of mhd:
        self.assertEqual(self.knownMhdXmax, mhd['x'].max())
        self.assertEqual(self.knownMhdXmin, mhd['x'].min())
        self.assertEqual(self.knownMhdZlim, mhd['z'].max())
        self.assertEqual(self.knownMhdZlim*-1, mhd['z'].min())
            
    def testAscii(self):
        # Open file:
        mhd = pb.IdlFile('data/pybats_test/y0_ascii.out', format='ascii')

        # Test units are loaded correctly:
        for v in mhd:
            if v not in self.varnames: continue
            self.assertEqual(self.knownMhdUnits[v], mhd[v].attrs['units'])

        # Test values inside of mhd:
        self.assertEqual(self.knownMhdXmax, mhd['x'].max())
        self.assertEqual(self.knownMhdXmin, mhd['x'].min())
        self.assertEqual(self.knownMhdZlim, mhd['z'].max())
        self.assertEqual(self.knownMhdZlim*-1, mhd['z'].min())


class TestMagGrid(unittest.TestCase):
    '''
    Test the class :class:`spacepy.pybats.bats.MagGridFile` to ensure opening,
    handling, and calculations are working correctly.
    '''

    knownDbnMax = 8.0770346781
    knownDbhMax = 8.07703468843
    knownPedMax = 2.783385368
    
    def testOpen(self):
        # Open both binary and ascii versions of same data.
        # Ensure expected values are loaded.
        m1 = pbs.MagGridFile('data/pybats_test/mag_grid_ascii.out',
                              format='ascii')
        m2 = pbs.MagGridFile('data/pybats_test/mag_grid_binary.out')

        self.assertAlmostEqual(self.knownDbnMax, m1['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m1['dBnPed'].max())
        self.assertAlmostEqual(self.knownDbnMax, m2['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m2['dBnPed'].max())

    def testCalc(self):
        # Open both binary and ascii versions of same data.
        # Ensure calculations give expected values.
        m1 = pbs.MagGridFile('data/pybats_test/mag_grid_ascii.out',
                              format='ascii')
        m2 = pbs.MagGridFile('data/pybats_test/mag_grid_binary.out')

        # Calculation of H-component:
        m1.calc_h()
        m2.calc_h()
        self.assertAlmostEqual(self.knownDbhMax, m1['dBh'].max())
        self.assertAlmostEqual(self.knownDbhMax, m2['dBh'].max())
        
class TestVirtSat(unittest.TestCase):
    '''
    Test the class :class:`spacepy.pybats.VirtSat` to ensure opening, handling,
    and calculations are working correctly.
    '''

    knownSatXmax = -1.6397
    knownSatPmax = 0.00509526
    knownSatOmax = 1.72270E-04
    knownSatHmax = 1.72235
    
    def testOpen(self):
        # Open file, ensure values are read properly.
        sat = pbs.VirtSat('data/pybats_test/sat_multispecies.sat')
        self.assertEqual(self.knownSatXmax, sat['x'].max())
        self.assertEqual(self.knownSatPmax, sat['p'].max())
        self.assertEqual(self.knownSatHmax, sat['rhoh'].max())
        self.assertEqual(self.knownSatOmax, sat['rhoo'].max())

    def testCalc(self):
        # Test various unit calculations
        sat = pbs.VirtSat('data/pybats_test/sat_multispecies.sat')

        # Test calculation of species number density:
        sat.calc_ndens()
        self.assertTrue('N' in sat)
        self.assertEqual(100, sat['oFrac'][0]+sat['hFrac'][0]+sat['heFrac'][0])
        
class TestImfInput(unittest.TestCase):
    '''
    Test reading, writing, and plotting from ImfInput files.
    '''

    # Files to open:
    file_single = 'data/pybats_test/imf_single.dat'
    file_multi  = 'data/pybats_test/imf_multi.dat'

    # Known values:
    knownImfBz   = [3, -10.]
    knownImfRho  = [5., 15.]
    knownImfTemp = [.80E+05, 1.20E+05]
    knownImfIono = [4.99, 0.01]

    def setUp(self):
        self.sing = pb.ImfInput(self.file_single)
        self.mult = pb.ImfInput(self.file_multi)
    
    def testOpen(self):
        # Test single fluid/default variable names:
        self.assertEqual(self.knownImfBz[0],   self.sing['bz'][0])
        self.assertEqual(self.knownImfBz[-1],  self.sing['bz'][-1])
        self.assertEqual(self.knownImfRho[0],  self.sing['rho'][0])
        self.assertEqual(self.knownImfRho[-1], self.sing['rho'][-1])
        self.assertEqual(self.knownImfTemp[0], self.sing['temp'][0])
        self.assertEqual(self.knownImfTemp[-1],self.sing['temp'][-1])

        # Open and test multi-fluid/non-standard variable names:
        self.assertEqual(self.knownImfBz[0],   self.mult['bz'][0])
        self.assertEqual(self.knownImfBz[-1],  self.mult['bz'][-1])
        self.assertEqual(self.knownImfRho[0],  self.mult['n'][0])
        self.assertEqual(self.knownImfTemp[0], self.mult['t'][0])
        self.assertEqual(self.knownImfTemp[-1],self.mult['t'][-1])
        self.assertEqual(self.knownImfIono[0], self.mult['IonoRho'][0])
        self.assertEqual(self.knownImfIono[-1],self.mult['IonoRho'][-1])

    def testPlot(self):
        import matplotlib.pyplot as plt

        # Create plots:
        f1 = self.sing.quicklook()
        f2 = self.sing.quicklook()

        # Test figures:
        self.assertTrue(isinstance(f1, plt.Figure))
        self.assertTrue(isinstance(f2, plt.Figure))

        plt.close(f1)
        plt.close(f2)
        
class TestExtraction(unittest.TestCase):
    '''
    Test Extraction class by opening a file with known solution.
    '''
    mhd = pbs.Bats2d('data/pybats_test/z0_sine.out')
    
    def testExtract(self):
        from numpy import pi, cos
        
        analytic = lambda(x): 1.+.5*cos(x*pi/10.)
        extr = self.mhd.extract(range(-5, 6),[-8]*11)
        for x, rho in zip(extr['x'], extr['rho']):
            self.assertAlmostEqual(rho, analytic(x), 2)
        
        
if __name__=='__main__':
    unittest.main()

