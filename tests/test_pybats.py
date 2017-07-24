#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's PyBats subpackage.

Copyright 2015 University of Michigan
"""


import datetime as dt
import unittest

import numpy as np

import spacepy.pybats as pb
import spacepy.pybats.bats as pbs

__all__ = ['TestIdlFile']

class TestParseFileTime(unittest.TestCase):
    '''
    Test the parse_filename_time function, which attempts to extract
    the datetime/runtime/run iteration from a standard SWMF file name.
    '''

    from datetime import datetime as dt
    
    files = ['mag_grid_e20130924-232600.out',
             'y=0_mhd_1_e20130924-220500-054.out',
             'y=0_mhd_2_t00001430_n00031073.out',
             'z=0_mhd_2_t00050000_n00249620.out']
    dates = [dt(2013,9,24,23,26,0), dt(2013,9,24,22, 5,0),
             None, None]
    times = [None, None, 870, 18000]
    iters = [None, None, 31073, 249620]

    def testParse(self):
        from spacepy.pybats import parse_filename_time
        for f, d, t, i in zip(self.files, self.dates, self.times, self.iters):
            self.assertEqual( parse_filename_time(f), (i,t,d) )
        
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

class TestRim(unittest.TestCase):

    def testRead(self):
        from spacepy.pybats import rim

        # Open file:
        iono=rim.Iono('data/pybats_test/it000321_104510_000.idl.gz')

class TestBats2d(unittest.TestCase):
    '''
    Test functionality of Bats2d objects.
    '''

    mhd = pbs.Bats2d('data/pybats_test/y0_binary.out')
    
    def testCalc(self):
        # Test all calculations:
        self.mhd.calc_all()
        
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
        
class TestSatOrbit(unittest.TestCase):
    '''
    Test reading and writing of satellite orbit files.
    '''
    def setUp(self):
        # Create 5 minutes of fake data:
        self.start = dt.datetime(2000,1,1)
        self.time = np.array([self.start+dt.timedelta(minutes=i)
                              for i in range(5)])

        self.pos = np.zeros( (3,5) )
        for i in range(5):
            self.pos[:,i]=[i, 10.+i, 100.+i]

    def testWrite(self):
        '''
        Create a file from scratch.
        '''
        from os import remove
        from spacepy.pybats import SatOrbit
        from numpy.testing import assert_array_equal as assert_array

        # New empty object:
        sat = SatOrbit()

        # Load in fake data:
        sat['time'] = self.time
        sat['xyz']  = self.pos

        # Add some info:
        sat.attrs['coor'] = 'SMG'
        sat.attrs['file'] = 'testsat.dat'
        sat.attrs['head'] = ['test','header','values']

        sat.write()
        import os
        os.listdir('.')

        #### Now, test the integrity of the file:
        sat = SatOrbit('./testsat.dat')

        # Check header:
        self.assertEqual(sorted(sat.attrs['head']),
                         sorted(['test','header','values']))
        self.assertEqual(sat.attrs['coor'], 'SMG')

        # Check time and position:
        assert_array(sat['time'], self.time)
        assert_array(sat['xyz'],  self.pos)

        # Get rid of file:
        remove('./testsat.dat')

    def testRead(self):
        from spacepy.pybats import SatOrbit
        from numpy.testing import assert_array_equal as assert_array

        sat = SatOrbit('data/pybats_test/testsat.dat')

        # Check header:
        self.assertEqual(sorted(sat.attrs['head']),
                         sorted(['test','header','values']))
        self.assertEqual(sat.attrs['coor'], 'SMG')

        # Check time and position:
        assert_array(sat['time'], self.time)
        assert_array(sat['xyz'],  self.pos)
    
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

    def tearDown(self):
        import os
        import glob

        # Remove temporary files.
        for f in glob.glob('*.tmp'):
            os.remove(f)

    def testWrite(self):
        # Test that files are correctly written to file.
        
        from numpy.testing import assert_array_equal as assert_array
        
        # Save original file names:
        old_file_1 = self.sing.attrs['file']
        old_file_2 = self.mult.attrs['file']

        # Rename files, write:
        self.sing.attrs['file'] = './imf_sing.tmp'
        self.mult.attrs['file'] = './imf_mult.tmp'
        self.sing.write()
        self.mult.write()

        # Reopen files:
        sing = pb.ImfInput('./imf_sing.tmp')
        mult = pb.ImfInput('./imf_mult.tmp')

        # Ensure files were written correctly:
        for v in sing:
            assert_array(self.sing[v], sing[v])
        for v in mult:
            assert_array(self.mult[v], mult[v])
        
        
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

    def testAppend(self):
        '''
        Test combining two files via timeseries_append
        '''

        # Combine two imf files, get size arrays before append:
        npts = self.mult['time'].size
        self.mult.timeseries_append(self.sing)

        # Check that the arrays were combined appropriately.
        for v in self.mult:
            if v in self.sing:
                self.assertEqual(self.mult['time'].size, self.mult[v].size)
            else:
                self.assertEqual(self.mult[v].size, npts)
        
class TestExtraction(unittest.TestCase):
    '''
    Test Extraction class by opening a file with known solution.
    '''
    mhd = pbs.Bats2d('data/pybats_test/z0_sine.out')
    
    def testExtract(self):
        from numpy import pi, cos
        
        analytic = lambda x: 1.+.5*cos(x*pi/10.)
        extr = self.mhd.extract(range(-5, 6),[-8]*11)
        for x, rho in zip(extr['x'], extr['rho']):
            self.assertAlmostEqual(rho, analytic(x), 2)
        
        
if __name__=='__main__':
    unittest.main()

