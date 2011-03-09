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
    path = __path__[0]+'/doc/'
    webbrowser.open(path+'index.html')

# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
    "irbempy", "empiricals", "radbelt"]

# Expose definitions from modules in this package.
from .toolbox import loadpickle, savepickle, dictree, printfig
import spacepy.time, spacepy.coordinates

#package info
__version__ = '0.1dev'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley', 'Josef Koller', 'Dan Welling', 'Brian Larsen', 'Jon Niehof', 'Mike Henderson']
__contact__ = 'spacepy@lanl.gov'
__license__ = """SpacePy: Space Science Tools for Python

SpacePy is a package of tools primarily aimed at the space science community.
Copyright (C) 2010 The SpacePy Team.

For license details, please see the license supplied with this package
"""

__notice__ = """SpacePy: Space Science Tools for Python
SpacePy is released under license. See __licence__ for details, and help() for HTML help."""

__licence__ = __license__ #for those who speak English, rather than an odd dialect

try: #if in iPython interactive shell, print licence notice
    assert __IPYTHON__active
    print(__notice__)
except NameError: #otherwise print single line notice
    print("SpacePy is released under license. See __licence__ for details, and help() for HTML help.")

# import some settings
from os import environ as ENVIRON
if 'SPACEPY' in ENVIRON:
    exec(compile(open(ENVIRON['SPACEPY']+'/.spacepy/spacepy.rc').read(), ENVIRON['SPACEPY']+'/.spacepy/spacepy.rc', 'exec'))
    DOT_FLN = ENVIRON['SPACEPY']+'/.spacepy'
else:
    exec(compile(open(ENVIRON['HOME']+'/.spacepy/spacepy.rc').read(), ENVIRON['HOME']+'/.spacepy/spacepy.rc', 'exec'))
    DOT_FLN = ENVIRON['HOME']+'/.spacepy'
