#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Test suite for SpacePy's PyBats subpackage.

Copyright 2015 University of Michigan
"""

import matplotlib
matplotlib.use('Agg')
from matplotlib import dates

import datetime as dt
import glob
import os
import unittest

import numpy as np
import numpy.testing

import spacepy_testing
import spacepy.pybats      as pb
import spacepy.pybats.bats as pbs
import spacepy.pybats.ram  as ram
import spacepy.pybats.gitm as gitm

__all__ = ['TestParseFileTime', 'TestIdlFile', 'TestRim', 'TestBats2d',
           'TestMagGrid', 'TestSatOrbit', 'TestVirtSat', 'TestImfInput',
           'TestExtraction', 'TestProbeIdlFile']

class TestParseFileTime(unittest.TestCase):
    '''
    Test the parse_filename_time function, which attempts to extract
    the datetime/runtime/run iteration from a standard SWMF file name.
    '''

    from datetime import datetime as dt

    files = ['mag_grid_e20130924-232600.out',
             'y=0_mhd_1_e20130924-220500-054.out',
             'y=0_mhd_2_t00001430_n00031073.out',
             'z=0_mhd_2_t00050000_n00249620.out',
             os.path.join(spacepy_testing.datadir,
                          'pybats_test', 'mag_grid_ascii.out'),
             'y=0_mhd_1_t20140410000000_n00001500.out',
             'z=0_mhd_2_e20140410-000000-000_20140410-000300-000.outs',
             'z=0_mhd_2_n00001500_00001889.outs'
             ]
    iters = [None, None, 31073, 249620, None, 1500, None, [1500, 1889]]
    times = [None, None, 870, 18000, None, None, None, None]
    dates = [dt(2013,9,24,23,26,0), dt(2013,9,24,22, 5,0),
             None, None, None, dt(2014, 4, 10, 0, 0),
             [dt(2014,4,10,0,0,0), dt(2014,4,10,0,3,0)], None]

    def testParse(self):
        from spacepy.pybats import parse_filename_time
        for f, d, t, i in zip(self.files, self.dates, self.times, self.iters):
            self.assertEqual( parse_filename_time(f), (i,t,d) )

class TestProbeIdlFile(unittest.TestCase):
    '''
    Test the function :func:`spacepy.pybats._probe_idlfile` across many
    different compatible files.
    '''

    filelist = [os.path.join(spacepy_testing.datadir, 'pybats_test', 'y=0_mhd_1_e20140410-000050.out'),
                os.path.join(spacepy_testing.datadir, 'pybats_test', 'y0_ascii.out'),
                os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_binary.out'),
                os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_ascii.out')]

    knownResponses = [('bin', '<', np.dtype('int32'), np.dtype('float32')),
                      ('asc', False, False, False),
                      ('bin', '<', np.dtype('int32'), np.dtype('float64')),
                      ('asc', False, False, False)]

    def testProbe(self):
        # Loop over various IdlFile files:
        for f, known in zip(self.filelist, self.knownResponses):
            # Probe to get format, binary information:
            response = pb._probe_idlfile(f)
            # Test for correct responses:
            for r, k in zip(response, known):
                self.assertEqual(r, k)

class TestScanBinHeader(unittest.TestCase):
    '''
    Test to ensure _scan_bin_header is properly working.
    '''

    # Files that have a single frame:
    files_single = [os.path.join(spacepy_testing.datadir, 'pybats_test', 'y=0_mhd_1_e20140410-000050.out'),
                    os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_binary.out')]
    # File that has multiple frames:
    file_multi = os.path.join(spacepy_testing.datadir, 'pybats_test',
                              'y=0_mhd_1_e20140410-000000-000_20140410-000200-000.outs')

    # Responses from single-frame files:
    knownSingle = [{'iter':68, 'runtime':10., 'ndim':2, 'nvar':15, 'end':175288},
                   {'iter':0,  'runtime':10., 'ndim':2, 'nvar':15, 'end':2416}]

    # Responses from multi-frame file:
    knownMulti = [{'start': 0, 'iter':2500, 'runtime':0.0, 'ndim':2, 'nvar':11, 'end':4512},
                  {'start': 4512, 'iter':2512, 'runtime': 120.0, 'ndim':2, 'nvar':11, 'end':9024}]

    def testOneFrame(self):
        '''Test files that only have one epoch frame.'''

        # Open files, get file properties, then probe header and test
        # against known values:
        for f, known in zip(self.files_single, self.knownSingle):
            props = pb._probe_idlfile(f)
            with open(f, 'rb') as data:
                info = pb._scan_bin_header(data, *props[1:])
                for key in known:
                    self.assertEqual(info[key], known[key])

    def testMultiFrame(self):
        '''Test files that have more than one epoch frame.'''

        # Test the only two frames in the file.
        # The ability to scan the number of entries in a file is
        # located within the IdlBin class.
        props = pb._probe_idlfile(self.file_multi)
        with open(self.file_multi, 'rb') as data:
            for known in self.knownMulti:
                info = pb._scan_bin_header(data, *props[1:])
                for key in known:
                    self.assertEqual(info[key], known[key])

class TestIdlFile(unittest.TestCase):
    '''
    Test the class :class:`spacepy.pybats.IdlFile` for different output
    file types and formats (ascii and binary).
    '''

    # Known values for single-frame *.out files:
    varnames = 'x z rho ux uy uz bx by bz p b1x b1y b1z e jx jy jz'
    units='R R Mp/cc km/s km/s km/s nT nT nT nPa nT nT nT J/m3 uA/m2 uA/m2 uA/m2'
    knownMhdUnits = dict(zip(varnames.split(), units.split()))
    knownMhdXmax = 31.0
    knownMhdXmin = -220.0
    knownMhdZlim = 124.0
    knownMhdTime = dt.datetime(2014, 4, 10, 0, 0, 50)
    
    # Known values for multi-frame *.outs files:
    # Time/iteration range covered by files:
    knownIterRng1  = [2500, 2512]
    knownIterRng2  = [2500, 2512]
    knownRtimeRng1 = [0.0, 120.0]
    knownRtimeRng2 = [0.0, 120.0]
    knownTimeRng1  = [dt.datetime(2014, 4, 10, 0, 0), dt.datetime(2014, 4, 10, 0, 2)]
    knownTimeRng2  = [dt.datetime(2014, 4, 10, 0, 0), dt.datetime(2014, 4, 10, 0, 2)]

    knownMax1 = {"Rho":14.871581077575684, "Ux":-1093.626953125,
                 "Bz":4.795100212097168, "P":2.5764927864074707,
                 "jy":0.0006453984533436596}
    knownMax2 = {"Rho":14.71767520904541,"Ux":-1107.3873291015625,
                 "Bz":4.878035068511963,"P":2.2604243755340576,
                 "jy":0.0007115363841876388}

    def testBinary(self):
        # Open file:
        mhd = pb.IdlFile(os.path.join(spacepy_testing.datadir,
                                      'pybats_test', 'y=0_mhd_1_e20140410-000050.out'))

        # Test time attribute:
        self.assertEqual(self.knownMhdTime, mhd.attrs['time'])
        
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
        mhd = pb.IdlFile(os.path.join(
            spacepy_testing.datadir, 'pybats_test', 'y0_ascii.out'))

        # Test units are loaded correctly:
        for v in mhd:
            if v not in self.varnames: continue
            self.assertEqual(self.knownMhdUnits[v], mhd[v].attrs['units'])

        # Test values inside of mhd:
        self.assertEqual(self.knownMhdXmax, mhd['x'].max())
        self.assertEqual(self.knownMhdXmin, mhd['x'].min())
        self.assertEqual(self.knownMhdZlim, mhd['z'].max())
        self.assertEqual(self.knownMhdZlim*-1, mhd['z'].min())

    def testReadOuts(self):
        # Start with y=0 slice MHD outs file:
        mhd = pb.IdlFile(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                      'y=0_mhd_1_e20140410-000000-000_20140410-000200-000.outs'))

        for i in range(len(self.knownIterRng1)):
            self.assertEqual(self.knownIterRng1[i],  mhd.attrs['iter_range'][i])
            self.assertEqual(self.knownRtimeRng1[i], mhd.attrs['runtime_range'][i])
            self.assertEqual(self.knownTimeRng1[i],  mhd.attrs['time_range'][i])

    def testSwitchFrame(self):
        '''Test our ability to open on arbitrary frame and change frame'''

        # Start with y=0 slice MHD outs file, but start on 2nd frame:
        f = os.path.join(spacepy_testing.datadir, 'pybats_test',
                         'y=0_mhd_1_e20140410-000000-000_20140410-000200-000.outs')
        mhd = pb.IdlFile(f, iframe=1)

        # Check against 2nd frame info:
        self.assertEqual(self.knownIterRng2[1],  mhd.attrs['iter'])
        self.assertEqual(self.knownRtimeRng2[1], mhd.attrs['runtime'])
        self.assertEqual(self.knownTimeRng2[1],  mhd.attrs['time'])

        # Check against 2nd frame data:
        for v in ['Rho', 'Ux', 'Bz', 'P', 'jy']:
            self.assertAlmostEqual(self.knownMax2[v], mhd[v].max(), places=14)

        # Switch frames, check for successful update of attributes:
        mhd.switch_frame(0)
        self.assertEqual(self.knownIterRng1[0],  mhd.attrs['iter'])
        self.assertEqual(self.knownRtimeRng1[0], mhd.attrs['runtime'])
        self.assertEqual(self.knownTimeRng1[0],  mhd.attrs['time'])

        # Check against 1st frame data:
        for v in ['Rho', 'Ux', 'Bz', 'P', 'jy']:
            self.assertAlmostEqual(self.knownMax1[v], mhd[v].max(), places=14)

class TestRim(unittest.TestCase):

    # Solutions for calc_I:
    knownI = {'n_I'    :2.4567986751877576e-10,
              'n_Iup'  :0.25176036603984825,
              'n_Idown':-0.25176036579416844,
              's_I'    :-8.211648588249157e-10,
              's_Iup'  :0.2517603660687805,
              's_Idown':-0.2517603668899454}

    def testReadZip(self):
        from spacepy.pybats import rim

        # Open file:
        iono=rim.Iono(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                   'it000321_104510_000.idl.gz'))


    def testReadAscii(self):
        import gzip
        from os import remove
        from shutil import copyfileobj
        from spacepy.pybats import rim

        try:
            # Unzip file and create a copy of it:
            name_in = os.path.join(spacepy_testing.datadir, 'pybats_test',
                                   'it000321_104510_000.idl.gz')
            name_out= name_in[:-3]
            with gzip.open(name_in, 'rb') as f_in, open(name_out, 'wb') as f_out:
                copyfileobj(f_in, f_out)

            # Test file:
            iono = rim.Iono(name_out)
        finally:
            # Remove temp file:
            remove(name_out)

    def testReadWrapped(self):
        '''Test reading files where entries are wrapped over multiple lines.'''
        from spacepy.pybats import rim
        iono = rim.Iono(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                     'it_wrapped.idl.gz'))

    def testIonoCalc(self):
        '''Test calculations made by rim.Iono objects.'''
        from spacepy.pybats import rim

        iono = rim.Iono(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                     'it000321_104510_000.idl.gz'))
        iono.calc_I()
        for key in self.knownI:
            self.assertAlmostEqual(self.knownI[key], iono[key])

    def testAddCont(self):
        from spacepy.pybats import rim
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        iono = rim.Iono(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                     'it000321_104510_000.idl.gz'))
        out = iono.add_cont('n_jr', add_cbar=True)

        self.assertTrue(isinstance(out[0], plt.Figure))
        self.assertTrue(isinstance(out[1], plt.Axes))
        self.assertTrue(isinstance(out[2], mpl.contour.QuadContourSet))
        self.assertTrue(isinstance(out[3], mpl.colorbar.Colorbar))

class TestBats2d(unittest.TestCase):
    '''
    Test functionality of Bats2d objects.
    '''

    knownMax1 = {'jx':1.4496836229227483e-05, 'jbz':7.309692051649108e-08,
                 'wy':0.0, 'u':1285.6114501953125}
    knownMax2 = {'jx': 1.680669083725661e-05, 'jbz': 8.276679608343329e-08,
                 'wy': 0.0, 'u': 1285.6114501953125}

    def setUp(self):
        self.mhd = pbs.Bats2d(os.path.join(spacepy_testing.datadir, 'pybats_test', 'y=0_mhd_1_e20140410-000050.out'))
        self.outs = pbs.Bats2d(os.path.join(spacepy_testing.datadir, 'pybats_test',
                        'y=0_mhd_1_e20140410-000000-000_20140410-000200-000.outs'))

    def testCalc(self):
        # Test all calculations:
        self.mhd.calc_all()

    def testSwitchFrame(self):
        '''Test switching frames and associated calculations'''
        # Perform some calculations on the first data frame:
        self.outs.calc_utotal()
        self.outs.calc_jxb()
        self.outs.calc_vort(conv=0) # Test maintaing kwargs between frames

        # Check initial calculated max values against reference:
        for k in self.knownMax1:
            self.assertAlmostEqual(self.outs[k].max(), self.knownMax1[k])

        # Switch frames and ensure that values update:
        self.outs.switch_frame(1)
        for k in self.knownMax2:
            self.assertAlmostEqual(self.outs[k].max(), self.knownMax2[k])

    def testMultispecies(self):
        # Open file:
        mhd = pbs.Bats2d(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                      'cut_multispecies.out'))
        mspec_varnames='x Rho Ux Uy Uz Bx By Bz P OpRho OpUx OpUy OpUz OpP jx jy jz g rbody cuty cutz'.split()
        mspec_units='km Mp/cc km/s km/s km/s nT nT nT nPa Mp/cc km/s km/s km/s nPa uA/m2 uA/m2 uA/m2'.split()
        knownMultispecUnits = dict(zip(mspec_varnames,
                                       mspec_units))

        # Test units are loaded correctly:
        for v in mhd:
            if v not in mspec_varnames: continue
            self.assertEqual(knownMultispecUnits[v], mhd[v].attrs['units'])

        mhd.calc_all(exclude=['calc_gradP', 'calc_vort'])

    def testPlotting(self):
        '''
        Create a contour plot, add stream traces with arrows,
        close plot, pass if no Exceptions.  This is a basic test that
        ensures that all methods and functions underlying contours and
        field line tracing are at least operating to completion.
        '''

        import matplotlib.pyplot as plt

        # Test adding a basic contour:
        fig, ax, cnt, cbar = self.mhd.add_contour('x', 'z', 'p', add_cbar=True)

        # Test adding field lines via "add_b_magsphere":
        self.mhd.add_b_magsphere(target=ax, narrow=5)

        # Test adding streams via "add_stream_scatter":
        self.mhd.add_stream_scatter('ux','uz',target=ax,narrow=1)

    def testGetStreamBad(self):
        """Get a streamline with a start point outside the input grid"""
        with self.assertRaises(ValueError) as cm:
            # X range is -220 to 31 (Z==0 is in range)
            self.mhd.get_stream(50, 0, 'ux', 'uz', method='rk4')
        self.assertEqual(
            'Start value 50 out of range for variable x.',
            str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            # Z range is -124 to 124 (X==0 is in range)
            self.mhd.get_stream(0, 150, 'ux', 'uz', method='rk4')
        self.assertEqual(
            'Start value 150 out of range for variable z.',
            str(cm.exception))

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

        m1 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_ascii.out'))
        m2 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_binary.out'))

        self.assertAlmostEqual(self.knownDbnMax, m1['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m1['dBnPed'].max())
        self.assertAlmostEqual(self.knownDbnMax, m2['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m2['dBnPed'].max())

    def testOpenTypeGuess(self):
        # Open both binary and ascii versions of same data
        # without specifying the type.
        # Ensure expected values are loaded.
        m1 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir,
                                          'pybats_test', 'mag_grid_ascii.out'))
        m2 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir,
                                          'pybats_test', 'mag_grid_binary.out'))

        self.assertAlmostEqual(self.knownDbnMax, m1['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m1['dBnPed'].max())
        self.assertAlmostEqual(self.knownDbnMax, m2['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, m2['dBnPed'].max())

    def testCalc(self):
        # Open both binary and ascii versions of same data.
        # Ensure calculations give expected values.
        m1 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_ascii.out'))
        m2 = pbs.MagGridFile(os.path.join(spacepy_testing.datadir, 'pybats_test', 'mag_grid_binary.out'))

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
        try:
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
        finally:
            # Get rid of file:
            remove('./testsat.dat')

    def testRead(self):
        from spacepy.pybats import SatOrbit
        from numpy.testing import assert_array_equal as assert_array

        sat = SatOrbit(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                    'testsat.dat'))

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
        sat = pbs.VirtSat(os.path.join(spacepy_testing.datadir,
                                       'pybats_test', 'sat_multispecies.sat'))
        self.assertEqual(self.knownSatXmax, sat['x'].max())
        self.assertEqual(self.knownSatPmax, sat['p'].max())
        self.assertEqual(self.knownSatHmax, sat['rhoh'].max())
        self.assertEqual(self.knownSatOmax, sat['rhoo'].max())

    def testCalc(self):
        # Test various unit calculations
        sat = pbs.VirtSat(os.path.join(spacepy_testing.datadir,
                                       'pybats_test', 'sat_multispecies.sat'))

        # Test calculation of species number density:
        sat.calc_ndens()
        self.assertTrue('N' in sat)
        self.assertEqual(100, sat['oFrac'][0]+sat['hFrac'][0]+sat['heFrac'][0])

class TestImfInput(unittest.TestCase):
    '''
    Test reading, writing, and plotting from ImfInput files.
    '''
    def setUp(self):
        # Files to open:
        self.file_single = os.path.join(spacepy_testing.datadir,
                                        'pybats_test', 'imf_single.dat')
        self.file_multi  = os.path.join(spacepy_testing.datadir,
                                        'pybats_test', 'imf_multi.dat')
        self.sing = pb.ImfInput(self.file_single)
        self.mult = pb.ImfInput(self.file_multi)

    # Known values:
    knownImfBz   = [3, -10.]
    knownImfRho  = [5., 15.]
    knownImfTemp = [.80E+05, 1.20E+05]
    knownImfIono = [4.99, 0.01]
    knownSubMilli= dt.datetime(2017, 9, 6, 16, 42, 37, 0)

    def tearDown(self):
        # Remove temporary files.
        for f in glob.glob('*.tmp'):
            os.remove(f)

    def testSubMillisec(self):
        '''
        Test case where sub-millisecond time values can create problems on
        read/write.
        '''
        # Create an IMF object from scratch, fill with zeros.
        imf = pb.ImfInput()
        for key in imf: imf[key]=[0]

        # Add a sub-millisecond time:
        imf['time'] = [dt.datetime(2017,9,6,16,42,36,999600)]

        # Write and test for non-failure on re-read:
        imf.write('testme.tmp')
        imf2 = pb.ImfInput('testme.tmp')

        # Test for floor of sub-millisecond times:
        self.assertEqual(self.knownSubMilli, imf2['time'][0])

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
    def setUp(self):
        self.mhd = pbs.Bats2d(os.path.join(spacepy_testing.datadir,
                                           'pybats_test', 'z0_sine.out'))

    def testExtract(self):
        analytic = lambda x: 1.+.5*np.cos(x*np.pi/10.)
        extr = self.mhd.extract(range(-5, 6),[-8]*11)
        for x, rho in zip(extr['x'], extr['rho']):
            self.assertAlmostEqual(rho, analytic(x), 2)


class TestGitm(unittest.TestCase):
    '''
    Test opening GITM binary files, handling files.
    '''
    nVars = 13
    nLat  = 18
    vers  = 4.03
    time  = dt.datetime(2015, 3, 16, 20, 1, 8)
    shape = (18,18)
    lat1  = 1.48352986

    def testBinary(self):
        '''
        This tests the ability to open a file and correctly read the attributes and
        variables as well as properly reshape the arrays and remove unused dimensions.
        '''
        # Open 2D file:
        f = gitm.GitmBin(os.path.join(spacepy_testing.datadir, 'pybats_test',
                                      'gitm_2D.bin'))
        # Check some critical attributes/values:
        self.assertEqual(self.nVars, f.attrs['nVars'])
        self.assertEqual(self.nLat,  f.attrs['nLat'])
        self.assertEqual(self.vers,  f.attrs['version'])
        self.assertEqual(self.time,  f['time'])
        self.assertEqual(self.shape, f['Longitude'].shape)
        self.assertAlmostEqual(   self.lat1, f['Latitude'][0,-1], 6)
        self.assertAlmostEqual(-1*self.lat1, f['Latitude'][0, 0], 6)

class RampyTests(unittest.TestCase):
    '''
    Tests for pybats.rampy
    '''
    def setUp(self):
        super(RampyTests, self).setUp()
        self.testfile = os.path.join(spacepy_testing.datadir, 'pybats_test',
                                     'ramsat_test.nc')
        self.p_test = os.path.join(spacepy_testing.datadir, 'pybats_test',
                                   'rampress_test.dat')
        self.p_test_noe = os.path.join(spacepy_testing.datadir, 'pybats_test',
                                   'rampress_test_noe.dat')

    def tearDown(self):
        super(RampyTests, self).tearDown()

    def test_RamSat_load(self):
        '''Test that RAM satellite netCDF will load'''
        data = ram.RamSat(self.testfile)

    def test_RamSat_contents_time(self):
        '''Test that start time attribute and Time variable are as expected'''
        data = ram.RamSat(self.testfile)
        self.assertEqual(data.starttime, dt.datetime(2012, 10, 29))
        tst = [dt.datetime(2012, 10, 29) + dt.timedelta(seconds=sc) for sc in [60, 120, 180]]
        numpy.testing.assert_array_equal(data['Time'], tst)

    def test_RamSat_contents_dims(self):
        '''Test for expected shape of loaded data array'''
        data = ram.RamSat(self.testfile)
        numpy.testing.assert_array_equal(data['FluxH+'].shape, [3, 72, 35])

    def test_RamSat_preserveFlux(self):
        '''Ensure that original flux is unchanged on calculating omni-flux'''
        data = ram.RamSat(self.testfile)
        flux_h = data['FluxH+'].copy()
        data.create_omniflux(check=False)
        numpy.testing.assert_array_equal(flux_h, data['FluxH+'])

    def test_RamSat_omni_pabin(self):
        '''Check that internal PA bin calc is consistent'''
        data = ram.RamSat(self.testfile)
        data.create_omniflux(check=False)
        omni1 = np.asarray(data['omniHe'].copy())
        origwid = np.asarray(data['pa_width'].copy())
        del data['pa_width']
        data.create_omniflux(check=False)
        omni2 = np.asarray(data['omniHe'])
        # test that calculated omni is close -- bin widths are
        # not fully recoverable from grid, so this is approximate
        numpy.testing.assert_allclose(omni1, omni2, rtol=1e-1)
        # and test that "pa_width" is consistent with simulation
        # (again, calculation from pa_grid isn't fully recovering
        # actual widths - so test is approximate)
        newwid = np.asarray(data['pa_width'])
        numpy.testing.assert_array_almost_equal(origwid, newwid, decimal=2)

    def test_RamSat_omnicalc_regress(self):
        '''Regression test for omni flux calculation'''
        data = ram.RamSat(self.testfile)
        #remove any precalculated omniflux
        rmkeys = [key for key in data if key.lower().startswith('omni')]
        for rmk in rmkeys:
            del data[rmk]
        regrH = np.array([0.0000000e+00, 1.4857327e+07, 1.4904620e+07, 1.4523555e+07,
                          1.3417377e+07, 1.1787100e+07, 9.8903560e+06, 8.0948140e+06,
                          6.6876425e+06, 5.7112065e+06, 5.0520345e+06, 4.5785050e+06,
                          4.2270450e+06, 3.9780850e+06, 3.8351362e+06, 3.7873392e+06,
                          3.9363230e+06, 4.0524422e+06, 2.9526408e+06, 9.0399219e+05,
                          5.2025672e+05, 2.4965306e+05, 9.1892180e+04, 3.0133383e+04,
                          1.6105718e+04, 1.4831358e+04, 1.7809768e+04, 2.2456926e+04,
                          2.5075262e+04, 2.3687787e+04, 1.8142951e+04, 1.0017233e+04,
                          3.8184448e+03, 1.1123656e+03, 2.2883945e+02], dtype=np.float32)
        data.create_omniflux(check=False)
        testarr = np.array(data['omniH'][0].data)
        numpy.testing.assert_array_almost_equal(regrH, testarr)

    def test_RamSat_orbit_formatter(self):
        '''Test label variables are as expected'''
        data = ram.RamSat(self.testfile)
        tst_time = matplotlib.dates.date2num(dt.datetime(2012, 10, 29, 0, 2))
        fmtstr = data._orbit_formatter(tst_time, None)
        expected = '00:02 UT\n04:35 MLT\n-13.4$^{\circ}$ MLat\nR=5.05 $R_{E}$'
        self.assertEqual(expected, fmtstr)

    def test_PressureFile_load(self):
        '''Make sure pressure file loads with no error'''
        data = ram.PressureFile(self.p_test)
        self.assertIn('aniHe', data)
        # single point checks
        self.assertAlmostEqual(data['pere'][0], 3.5930268833807135)
        self.assertAlmostEqual(data['perH'][0], 0.17925744886042114)
        self.assertAlmostEqual(data['total'][0], 2.615855029414899)

    def test_PressureFile_no_elec(self):
        '''Make sure pressure file loads with electrons turned off'''
        data = ram.PressureFile(self.p_test_noe)
        self.assertIn('aniH', data)
        self.assertIn('total', data)
        self.assertNotIn('pere', data)
        self.assertNotIn('pare', data)
        self.assertAlmostEqual(data['perH'][0], 0.17925744886042114)
        self.assertAlmostEqual(data['total'][0], 2.615855029414899)

    def test_pcolor_plot(self):
        '''Pcolor plot should not fail and should return known types'''
        import matplotlib.pyplot as plt
        data = ram.PressureFile(self.p_test)
        out = data.add_pcol_press('totH')
        self.assertTrue(isinstance(out[0], plt.Figure))
        self.assertTrue(isinstance(out[1], plt.Axes))
        self.assertTrue(isinstance(out[2], matplotlib.collections.QuadMesh))
        self.assertTrue(out[3] is None)

    def test_contour_plot(self):
        '''Contour plot should not fail and should return known types'''
        import matplotlib.pyplot as plt
        data = ram.PressureFile(self.p_test)
        out = data.add_cont_press('totH', add_cbar=True)
        self.assertTrue(isinstance(out[0], plt.Figure))
        self.assertTrue(isinstance(out[1], plt.Axes))
        self.assertTrue(isinstance(out[2], matplotlib.tri.TriContourSet))
        self.assertTrue(isinstance(out[3], matplotlib.colorbar.Colorbar))


if __name__=='__main__':
    unittest.main()

