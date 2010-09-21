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
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
    "irbempy", "constants", "empiricals", "radbelt"]

# Expose definitions from modules in this package.
from toolbox import loadpickle, savepickle, dictree, printfig
import time, coordinates, constants

#package info
__version__ = '0.1dev'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley (smorley@lanl.com/morley_steve@hotmail.com)',
'Josef Koller (jkoller@lanl.gov )']
__license__ = """SpacePy: Space Science Tools for Python

SpacePy is a package of tools primarily aimed at the space science community.
Copyright (C) 2010  Steven Morley, Josef Koller

LANL LICENSE AGREEMENT FOR SPACEPY 1.0

This LICENSE AGREEMENT is between the Los Alamos National Laboratory ("LANL"), and the Individual or Organization (“Licensee”) accessing and otherwise using SpacePy 1.0 software in source or binary form and its associated documentation.
Subject to the terms and conditions of this License Agreement, LANL hereby grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce, analyze, test, perform and/or display publicly, prepare derivative works, distribute, and otherwise use SpacePy 1.0 alone or in any derivative version, provided, however, that LANL’s License Agreement and LANL’s notice of copyright, i.e., “Copyright © 2001-2010 Los Alamos National Laboratory; All Rights Reserved” are retained in SpacePy 1.0 alone or in any derivative version prepared by Licensee.
In the event Licensee prepares a derivative work that is based on or incorporates SpacePy 1.0 or any part thereof, and wants to make the derivative work available to others as provided herein, then Licensee hereby agrees to include in any such work a brief summary of the changes made to SpacePy.
LANL is making SpacePy 1.0 available to Licensee on an “AS IS” basis. LANL MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE, BUT NOT LIMITATION, LANL MAKES NO AND DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF SPACEPY 1.0 WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.
LANL SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF SPACEPY 1.0 FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING SPACEPY 1.0, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
This License Agreement will automatically terminate upon a material breach of its terms and conditions.
Nothing in this License Agreement shall be deemed to create any relationship of agency, partnership, or joint venture between LANL and Licensee. This License Agreement does not grant permission to use LANL trademarks or trade name in a trademark sense to endorse or promote products or services of Licensee, or any third party.
By copying, installing or otherwise using SpacePy 1.0, Licensee agrees to be bound by the terms and conditions of this License Agreement.

"""

__notice__ = """SpacePy: Space Science Tools for Python
SpacePy is released under license. See __licence__ for details, and help() for HTML help."""

__licence__ = __license__ #for those who speak English, rather than an odd dialect

try: #if in iPython interactive shell, print licence notice
    assert __IPYTHON__active
    print __notice__
except: #otherwise print single line notice
    print "SpacePy is released under license. See __licence__ for details, and help() for HTML help."

# import some settings
from os import environ as ENVIRON
if ENVIRON.has_key('SPACEPY'):
	execfile(ENVIRON['SPACEPY']+'/.spacepy/spacepy.rc')
	DOT_FLN = ENVIRON['SPACEPY']+'/.spacepy'
else:
	execfile(ENVIRON['HOME']+'/.spacepy/spacepy.rc')
	DOT_FLN = ENVIRON['HOME']+'/.spacepy'


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

