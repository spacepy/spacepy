#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Master test suite for SpacePy

author: Steve Morley
organization: LANL
contact: smorley@lanl.gov

version: V1: 28-May-2011 
"""

import sys

try:
    import unittest_pretty as ut
except ImportError:
    import unittest as ut

from test_time import *
from test_empiricals import *
from test_toolbox import *
from test_omni import *
from test_coordinates import *
from test_seapy import *
from test_poppy import *
from test_pycdf import *
# add others here as they are written 

if __name__ == '__main__':
    ut.main()
