#!/usr/bin/env python

"""Integration tests for SpacePy

These are heavier tests than the unit tests that may e.g. temporarily
affect the state of the machine, rather than the fully isolated unit
tests.

Copyright 2019 SpacePy contributors
"""


import unittest

from integration_omni import *
from integration_toolbox import *


if __name__ == '__main__':
    unittest.main()
