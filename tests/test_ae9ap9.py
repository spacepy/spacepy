#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for toolbox

Copyright 2010-2012 Los Alamos National Security, LLC.
"""
import datetime
import glob
import gzip
import io
import os
import sys
import shutil
import tempfile
import unittest

import matplotlib.collections
import matplotlib.patches
import numpy.testing

import spacepy_testing
from spacepy import ae9ap9

__all__ = ['ae9ap9Tests', 'ae9ap9PlotTests',]

class ae9ap9Tests(unittest.TestCase):
    """
    This class mostly provides regression tests as there are few first principles things that
    can be tested in the file reader and parser
    """
    
    def setUp(self):
        super(ae9ap9Tests, self).setUp()
        self.datafiles = sorted(glob.glob(os.path.join(
            spacepy_testing.datadir,
            'Run1.AE9.CLoutput_mc_fluence_agg_pctile_??.txt')))
        self.orbitfile = os.path.join(
            spacepy_testing.datadir, 'ephem_AE9test.dat')

    def tearDown(self):
        super(ae9ap9Tests, self).tearDown()

    def test_setUnits_error(self):
        """Invalid units raise the correct error and message"""
        ans = ae9ap9.readFile(self.datafiles[0])
        self.assertRaisesRegex(ValueError, '^(Units of FeV)', ans.setUnits, 'FeV')

    def test_setUnits_print(self):
        """Print out current units"""
        ans = ae9ap9.readFile(self.datafiles[0])
        output = io.StringIO()
        realstdout = sys.stdout
        try:
            sys.stdout = output
            ans.setUnits()
            result = output.getvalue()
        finally:
            output.close()
            sys.stdout = realstdout
        self.assertEqual(
            "Energy units: MeV\nFluence units: #/cm^2/MeV\n", result)

    def test_setUnits_convert(self):
        """Conversion correctly changes flux/fluence values, energy values and units"""
        ans = ae9ap9.readFile(self.datafiles[0])
        curr = ans['Energy'].attrs['UNITS']
        self.assertEqual(curr, 'MeV') #test files should have MeV units
        E0 = ans['Energy'][0]
        F0 = ans['Fluence'][0,0]
        ans.setUnits('keV')
        self.assertEqual(ans['Energy'].attrs['UNITS'], 'keV')
        self.assertEqual(E0*1000.0, ans['Energy'][0])
        self.assertTrue('/keV' in ans['Fluence'].attrs['UNITS'])
        self.assertEqual(F0/1000.0, ans['Fluence'][0,0])

    def test_readFile(self):
        """Can read a file in and get the same answer"""
        ans = ae9ap9.readFile(self.datafiles[0])
        self.assertEqual((21, ), ans['Energy'].shape)
        self.assertEqual((121, ), ans['Epoch'].shape)
        self.assertEqual((121, 21 ), ans['Fluence'].shape)
        self.assertEqual((121, 3), ans['Coords'].shape)
        self.assertEqual((121, ), ans['MJD'].shape)
        self.assertEqual((3, ), ans['posComp'].shape)

    def test_readFileGzip(self):
        """can read a gzipped file as well"""
        # take the input file and gzip it to a temp file
        tmpf = tempfile.NamedTemporaryFile(suffix='.gz')
        tmpf.close()
        tmpname = tmpf.name
        try:
            with open(self.datafiles[0], 'rb') as f_in:
                with gzip.open(tmpname, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            ans = ae9ap9.readFile(tmpname)
            self.assertEqual((21, ), ans['Energy'].shape)
            self.assertEqual((121, ), ans['Epoch'].shape)
            self.assertEqual((121, 21 ), ans['Fluence'].shape)
            self.assertEqual((121, 3), ans['Coords'].shape)
            self.assertEqual((121, ), ans['MJD'].shape)
            self.assertEqual((3, ), ans['posComp'].shape)
        finally:
            os.remove(tmpname)

    def test_combinePercentiles(self):
        """Can read and combine percentile files"""
        realstdout = sys.stdout
        output = io.StringIO()
        sys.stdout = output
        ans = ae9ap9.combinePercentiles(self.datafiles)
        output.close()
        sys.stdout = realstdout
        self.assertEqual((21, ), ans['Energy'].shape)
        self.assertEqual((21, 2), ans['Fluence'].shape)
        self.assertEqual((2, ), ans['Percentile'].shape)

    def test_parseInfo_orbit(self):
        """Read in an orbit data file"""
        with open(self.orbitfile, 'r') as fp:
            dat = fp.readlines()
        dat = [v.strip() for v in dat]
        out = ae9ap9._parseInfo(dat)
        ans = {'propagator': 'Kepler with J2',
               'time_format': 'MJD',
               'coord_system': ('GEI', 'Re'),
               'delimiter': ','}
        for k in out:
            self.assertEqual(out[k], ans[k])

    def test_Lm(self):
        """Calculate McIlwain L"""
        ans = ae9ap9.readFile(self.datafiles[0])
        ans.getLm(model='OPQUIET')
        self.assertEqual('OPQUIET', ans['Lm'].attrs['MODEL'])
        self.assertEqual((121,), ans['Lm'].shape)
        numpy.testing.assert_allclose(6.6, ans['Lm'], atol=.2)

    def test_readDOYFrac(self):
        """Read times in year, DOY, fraction"""
        infile = glob.glob(os.path.join(
            spacepy_testing.datadir, 'Run2.AE9*'))[0]
        ans = ae9ap9.readFile(infile)
        numpy.testing.assert_array_equal(
            [datetime.datetime(2015, 1, 1, 0, *ms)
             for ms in [(0,), (4, 59, 998080), (9, 59, 996160)]],
            ans['Epoch'])

    def test_readUTC(self):
        """Read times in ymdhms"""
        infile = glob.glob(os.path.join(
            spacepy_testing.datadir, 'Run3.AE9*'))[0]
        ans = ae9ap9.readFile(infile)
        numpy.testing.assert_array_equal(
            [datetime.datetime(2015, 1, 1, 0, minute)
             for minute in [0, 5, 10]],
            ans['Epoch'])


class ae9ap9PlotTests(spacepy_testing.TestPlot):
    """Test plotting functions of ae9ap9"""

    def setUp(self):
        super().setUp()
        self.datafile = sorted(glob.glob(os.path.join(
            spacepy_testing.datadir,
            'Run1.AE9.CLoutput_mc_fluence_agg_pctile_??.txt')))[0]

    def test_plotOrbit(self):
        """Simple orbit plot, all the defaults"""
        ans = ae9ap9.readFile(self.datafile)
        fig = ans.plotOrbit()
        ax = fig.get_axes()
        self.assertEqual(2, len(ax))
        ax_xy, ax_xz = ax
        gse_xy = ax_xy.get_lines()[0].get_xydata()
        numpy.testing.assert_array_equal(gse_xy, ans['Coords'][:, :2])
        wedges = [c for c in ax_xy.get_children()
                  if isinstance(c, matplotlib.patches.Wedge)]
        self.assertEqual(2, len(wedges))  # half-lit Earth
        self.assertEqual((0, 0), wedges[0].center)

    def test_plotSpectrogram(self):
        """Plot simple spectrogram, all defaults"""
        ans = ae9ap9.readFile(self.datafile)
        # if not calculated before plot, uses T89 and needs omni
        ans.getLm(model='OPQUIET')
        ax = ans.plotSpectrogram()
        qm = [c for c in ax.get_children()
              if isinstance(c, matplotlib.collections.QuadMesh)]
        self.assertEqual(1, len(qm))
        data = qm[0].get_array()
        # default 10 time bins, L 2 to 8 (exclusive) in 1/4 Re bins
        # exact shape depends on mpl version
        self.assertEqual(10 * 23, data.size)
        self.assertTrue((data > 0).any())

    def test_plotSummary(self):
        """Summary plot, all defaults"""
        ans = ae9ap9.readFile(self.datafile)
        fig = ans.plotSummary()
        ax = fig.get_axes()
        self.assertEqual(3, len(ax))
        ax_xy, ax_xz, ax_fluence = ax
        # 21 energies
        self.assertEqual(21, len(ax_fluence.get_lines()))

    def test_plotSummaryTimerange(self):
        """Summary plot with a timerange"""
        ans = ae9ap9.readFile(self.datafile)
        fig = ans.plotSummary(timerange=[datetime.datetime(2010, 1, 19, 8, 0),
                                         datetime.datetime(2010, 1, 19, 8, 30)])
        ax = fig.get_axes()
        self.assertEqual(3, len(ax))
        ax_xy, ax_xz, ax_fluence = ax
        l = ax_xy.get_lines()[0]
        x, y = l.get_data()
        numpy.testing.assert_array_equal(x, ans['Coords'][60:90, 0])
        # Is the line on-screen?
        xlim = ax_xy.get_xlim()
        self.assertLess(xlim[1], x.min())  # axis inverted
        self.assertGreater(xlim[0], x.max())
        ylim = ax_xy.get_ylim()
        self.assertLess(ylim[1], y.min())
        self.assertGreater(ylim[0], y.max())


if __name__ == "__main__":
    unittest.main()
