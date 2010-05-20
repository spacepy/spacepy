#!/usr/bin/python
# -*- coding: utf-8 -*-

"""SpacePy: Space Science Tools for Python

SpacePy is a package of tools primarily aimed at the space science community.
This __init__.py file sets the parameters for import statements.

If running the ipython shell, simply type '?' after any command for help.
ipython also offers tab completion, so hitting tab after '<module name>.'
will list all available functions, classes and variables.

Detailed HTML documentation is available locally in the spacepy/doc directory
and can be launched by typing:
>>> spacepy.help()
"""

def help():
    """Launches web browser with local HTML help"""
    
    import webbrowser
    path=__path__[0]+'/doc/'
    webbrowser.open(path+'index.html')

# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", "onerapy"]

# Expose definitions from modules in this package.
from toolbox import loadpickle, savepickle, dictree, printfig
from time import tickrange, doy2date

#package info
__version__ = '0.1dev'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley (smorley@lanl.com/morley_steve@hotmail.com)',
'Josef Koller (jkoller@lanl.gov )']
__license__ = """SpacePy: Space Science Tools for Python

SpacePy is a package of tools primarily aimed at the space science community.
Copyright (C) 2010  Steven Morley, Josef Koller

SpacePy is released under GPL v3.0:
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__notice__ = """SpacePy: Space Science Tools for Python
SpacePy is released under GPL v3.0. See __licence__ for details, and help() for HTML help."""

__licence__ = __license__ #for those who speak English, rather than an odd dialect

try: #if in iPython interactive shell, print licence notice
    assert __IPYTHON__active
    print __notice__
except: #otherwise print single line notice
    print "SpacePy is released under GPL v3.0. See __licence__ for details, and help() for HTML help."


# -----------------------------------------------
# some settings
NCPUS = 1

# -----------------------------------------------
def test_all():
    """
    test all spacepy routines

    Returns:
    ========
        - nFAIL (int) : number of failures

    Example:
    ========
    >>> spacepy.test()
    test_ticktock: PASSED TEST 1
    test_ticktock: PASSED TEST 2
    0

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)

    Version:
    ========
    V1: 24-Jan-2010
    """
    
    import spacepy.time

    nFAIL = spacepy.time.test()
    nFAIL =+ spacepy.toolbox.test()

    
    
    return nFAIL
