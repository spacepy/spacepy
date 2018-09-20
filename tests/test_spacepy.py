#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Master test suite for SpacePy

version: V1: 28-May-2011


Copyright 2010-2014 Los Alamos National Security, LLC.
"""

import sys

try:
    import unittest_pretty as ut
except ImportError:
    import unittest as ut

from test_pybats import *
from test_time import *
from test_empiricals import *
from test_toolbox import *
from test_omni import *
from test_coordinates import *
from test_seapy import *
from test_poppy import *
from test_pycdf import *
from test_data_assimilation import *
from test_spectrogram import *
from test_irbempy import *
from test_datamanager import *
from test_datamodel import *
from test_base import *
from test_plot_utils import *
from test_rst import *
from test_lib import *
from test_ae9ap9 import *
# add others here as they are written

if __name__ == '__main__':
    ut.main()
