#!/usr/bin/env python

"""Helper functions for SpacePy unit tests

This module is importable from test scripts that are in this directory.
"""

import os.path


testsdir = os.path.dirname(os.path.abspath(__file__))
"""Directory containing the unit test scripts"""

datadir = os.path.join(testsdir, 'data')
"""Directory containing unit test data"""
