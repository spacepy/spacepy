#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Master test suite for SpacePy

Copyright 2010-2014 Triad National Security, LLC.
Copyright 2015-2020 the contributors
"""

import glob
import importlib
import os.path
import sys
try:
    import unittest_pretty as ut
except ImportError:
    import unittest as ut
import warnings

# Warnings attributed to SpacePy code should be errors
warnings.filterwarnings("error", module=r"spacepy\.")
warnings.filterwarnings("error", module="spacepy_testing$")
# omni is imported from spacepy and may not have this data in testing
warnings.filterwarnings("default", "Qin-Denton/OMNI2 data not found",
                        UserWarning, r"spacepy\.")

import spacepy_testing

# Duplicative of test discovery, but avoids pulling the integration tests
for testfile in glob.glob(os.path.join(os.path.dirname(__file__), 'test_*.py')):
    modname = os.path.basename(testfile)[:-3]
    if modname == 'test_all':
        continue
    mod = importlib.import_module(modname)
    for name in dir(mod):
        cls = getattr(mod, name)
        if isinstance(cls, type) and issubclass(cls, ut.TestCase):
            globals()[name] = cls


if __name__ == '__main__':
    ut.main(warnings=False)  # do not reset warning filter
