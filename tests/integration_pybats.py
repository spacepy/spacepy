#!/usr/bin/env python3

"""Integration tests for Pybats.

Copyright 2025 SpacePy contributors
"""

import os
import unittest

import spacepy_testing
import spacepy.pybats as pb


class PybatsFileIntegration(unittest.TestCase):
    '''
    Large and cumbersome tests related to :class:`spacepy.pybats.IdlFile`
    file types and formats (ascii and binary).
    '''

    # This single file is 2.4K in size.
    filebase = os.path.join(spacepy_testing.datadir, 'pybats_test',
                            'mag_grid_binary.out')

    # Reference values:
    knownDbnMax = 8.0770346781
    knownPedMax = 2.783385368

    def testBinaryLarge(self):
        '''Test ability to open large (>2Gb) file.'''

        # Set location for temporary large file:
        tempfile = 'testfile_large.outs'

        # Concat small file into a large one (>2gb)
        with open(tempfile, 'wb') as f_out:
            with open(self.filebase, 'rb') as f_in:
                rawbits = f_in.read()
            for i in range(1000000):
                f_out.write(rawbits)

        # Read the large file to test for failures:
        data = pb.IdlFile(tempfile)

        # Switch to final frame and check values:
        data.switch_frame(-1)
        self.assertAlmostEqual(self.knownDbnMax, data['dBn'].max())
        self.assertAlmostEqual(self.knownPedMax, data['dBnPed'].max())

        # Remove temporary test file.
        os.remove(tempfile)


if __name__ == "__main__":
    unittest.main()
