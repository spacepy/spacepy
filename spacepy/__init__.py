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

Copyright ©2010 Los Alamos National Security, LLC.
"""

def help():
    """Launches web browser with local HTML help"""
    
    import webbrowser
    path = __path__[0]+'/doc/'
    webbrowser.open(path+'index.html')

# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
    "irbempy", "empiricals", "radbelt", "borg"]

# Expose definitions from modules in this package.
from .toolbox import loadpickle, savepickle, dictree, printfig
import spacepy.time, spacepy.coordinates

#package info
__version__ = '0.1'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley', 'Josef Koller', 'Dan Welling', 'Brian Larsen', 'Jon Niehof', 'Mike Henderson']
__contact__ = 'spacepy@lanl.gov'
__license__ = """SpacePy: Space Science Tools for Python


Copyright ©2010 Los Alamos National Security, LLC.
All Rights Reserved.

 This material was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE


 1. This LICENSE AGREEMENT is between the Los Alamos National Security, LLC ("LANS"), and the Individual or Organization ("Licensee") accessing and otherwise using SpacePy 0.1 software in source or binary form and its associated documentation.

 2. Subject to the terms and conditions of this License Agreement, LANS hereby grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce, analyze, test, perform and/or display publicly, prepare derivative works, distribute, and otherwise use SpacePy 0.1 alone or in any derivative version, provided, however, that LANS’ License Agreement and LANS’ notice of copyright, i.e., "Copyright (c) 2010 Los Alamos National Security, LLC; All Rights Reserved" are retained in SpacePy 0.1 alone or in any derivative version prepared by Licensee.

 3. In the event Licensee prepares a derivative work that is based on or incorporates SpacePy 0.1 or any part thereof, and wants to make the derivative work available to others as provided herein, then Licensee hereby agrees to include in any such work a brief summary of the changes made to SpacePy 0.1.

 4. LANS is making SpacePy 0.1 available to Licensee on an "AS IS" basis. LANS MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE, BUT NOT LIMITATION, LANS MAKES NO AND DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF SPACEPY 0.1 WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.

 5. LANS SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF SPACEPY 0.1 FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING SPACEPY 0.1, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

 6. This License Agreement will automatically terminate upon a material breach of its terms and conditions.

 7. Nothing in this License Agreement shall be deemed to create any relationship of agency, partnership, or joint venture between LANS and Licensee. This License Agreement does not grant permission to use LANS trademarks or trade name in a trademark sense to endorse or promote products or services of Licensee, or any third party.

 8. By copying, installing or otherwise using SpacePy 0.1, Licensee agrees to be bound by the terms and conditions of this License Agreement.
"""

__citation__ = """When publishing research which used SpacePy, please provide appropriate
credit to the SpacePy team via citation or acknowledgement.

To cite SpacePy in publications, use (BibTeX code):
@INPROCEEDINGS{spacepy11,
author = {{Morley}, S.~K. and {Koller}, J. and {Welling}, D.~T. and {Larsen}, B.~A. and {Henderson}, M.~G. and {Niehof}, J.~T.},
title = "{Spacepy - A Python-based library of tools for the space sciences}",
booktitle = "{Proceedings of the 9th Python in science conference (SciPy 2010)}",
year = 2011,
address = {Austin, TX}
}

Certain modules may provide additional citations in the __citation__
attribute. Contact a module's author before publication or public
presentation of analysis performed by that module. This allows the author
to validate the analysis and receive appropriate credit for his or her
work.
"""

__notice__ = """SpacePy: Space Science Tools for Python
SpacePy is released under license.
See __licence__ for details, __citation__ for citation information,
and help() for HTML help."""

__licence__ = __license__ #for those who speak English, rather than an odd dialect

try: #if in iPython interactive shell, print licence notice
    if __IPYTHON__active == 1: print(__notice__)
except NameError: #otherwise print single line notice
    print("SpacePy is released under license. See __licence__ and __citation__ for details, and help() for HTML help.")

# import some settings
from os import environ as ENVIRON
if 'SPACEPY' in ENVIRON:
    exec(compile(open(ENVIRON['SPACEPY']+'/.spacepy/spacepy.rc').read(), ENVIRON['SPACEPY']+'/.spacepy/spacepy.rc', 'exec'))
    DOT_FLN = ENVIRON['SPACEPY']+'/.spacepy'
else:
    exec(compile(open(ENVIRON['HOME']+'/.spacepy/spacepy.rc').read(), ENVIRON['HOME']+'/.spacepy/spacepy.rc', 'exec'))
    DOT_FLN = ENVIRON['HOME']+'/.spacepy'
