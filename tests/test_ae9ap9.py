#!/usr/bin/env python2.6
# -*- coding: utf-8 -*-

"""
Test suite for toolbox

Copyright 2010-2012 Los Alamos National Security, LLC.
"""
import glob
import gzip
import os
import shutil
import tempfile
import unittest

from spacepy import ae9ap9

__all__ = ['ae9ap9Tests', ]

class ae9ap9Tests(unittest.TestCase):
    """
    This class mostly provides regression tests as there are few first principle things that
    can be tested in the file reader and parser
    """
    
    def setUp(self):
        super(ae9ap9Tests, self).setUp()
        self.datafiles = glob.glob('data/Run1.AE9.CLoutput_mc_fluence_agg_pctile_??.txt')
        
    def tearDown(self):
        super(ae9ap9Tests, self).tearDown()

    def test_readFile(self):
        """Can read a file in and get the same answer"""
        ans = ae9ap9.readFile(self.datafiles[0])
        self.assertEqual((21, ), ans['Energy'].shape)
        self.assertEqual((121, ), ans['Epoch'].shape)
        self.assertEqual((121, 21 ), ans['Fluence'].shape)
        self.assertEqual((121, 3), ans['GSE'].shape)
        self.assertEqual((121, ), ans['MJD'].shape)
        self.assertEqual((3, ), ans['posComp'].shape)

    def test_readFileGzip(self):
        """can read a gziped file as well"""
        # take the input file and gzip it to a temp file
        tmpf = tempfile.NamedTemporaryFile(suffix='.gz')
        tmpf.close()
        tmpname = tmpf.name
        try:
            with open(self.datafiles[0], 'rb') as f_in, gzip.open(tmpname, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            ans = ae9ap9.readFile(tmpname)
            self.assertEqual((21, ), ans['Energy'].shape)
            self.assertEqual((121, ), ans['Epoch'].shape)
            self.assertEqual((121, 21 ), ans['Fluence'].shape)
            self.assertEqual((121, 3), ans['GSE'].shape)
            self.assertEqual((121, ), ans['MJD'].shape)
            self.assertEqual((3, ), ans['posComp'].shape)
        finally:
            os.remove(tmpname)

    def test_combinePercentiles(self):
        """Can read and combine percentile files"""
        ans = ae9ap9.combinePercentiles(self.datafiles)
        import spacepy.toolbox as tb
        #tb.dictree(ans, verbose=1)
        self.assertEqual((21, ), ans['Energy'].shape)
        self.assertEqual((21, 2), ans['Fluence'].shape)
        self.assertEqual((2, ), ans['Percentile'].shape)
        
if __name__ == "__main__":
    unittest.main()
