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

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import multiprocessing
import os
import os.path
import warnings

def help():
    """Launches web browser with local HTML help"""
    
    import webbrowser
    fspec = os.path.join(__path__[0], 'Doc', 'index.html')
    if not os.path.exists(fspec):
        fspec = os.path.join(__path__[0], '..',
                             'Doc', 'build', 'html', 'index.html')
    if os.path.exists(fspec):
        webbrowser.open(fspec)
    else:
        print("Can't find help files in {0}".format(__path__[0]))

# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
    "irbempy", "empiricals", "radbelt", "data_assimilation"]

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


 1. This LICENSE AGREEMENT is between the Los Alamos National Security, LLC ("LANS"), and the Individual or Organization ("Licensee") accessing and otherwise using SpacePy 0.1.0 software in source or binary form and its associated documentation.

 2. Subject to the terms and conditions of this License Agreement, LANS hereby grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce, analyze, test, perform and/or display publicly, prepare derivative works, distribute, and otherwise use SpacePy 0.1.0 alone or in any derivative version, provided, however, that LANS’ License Agreement and LANS’ notice of copyright, i.e., "Copyright (c) 2010 Los Alamos National Security, LLC; All Rights Reserved" are retained in SpacePy 0.1.0 alone or in any derivative version prepared by Licensee.

 3. In the event Licensee prepares a derivative work that is based on or incorporates SpacePy 0.1.0 or any part thereof, and wants to make the derivative work available to others as provided herein, then Licensee hereby agrees to include in any such work a brief summary of the changes made to SpacePy 0.1.0.

 4. LANS is making SpacePy 0.1.0 available to Licensee on an "AS IS" basis. LANS MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE, BUT NOT LIMITATION, LANS MAKES NO AND DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF SPACEPY 0.1.0 WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.

 5. LANS SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF SPACEPY 0.1.0 FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING SPACEPY 0.1.0, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

 6. This License Agreement will automatically terminate upon a material breach of its terms and conditions.

 7. Nothing in this License Agreement shall be deemed to create any relationship of agency, partnership, or joint venture between LANS and Licensee. This License Agreement does not grant permission to use LANS trademarks or trade name in a trademark sense to endorse or promote products or services of Licensee, or any third party.

 8. By copying, installing or otherwise using SpacePy 0.1.0, Licensee agrees to be bound by the terms and conditions of this License Agreement.
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
def _read_config(rcfile):
    """Read configuration information from a file"""
    global ENABLE_DEPRECATION_WARNING, NCPUS, OMNI_URL, LEAPSEC_URL, PSDDATA_URL
    defaults = {'enable_deprecation_warning': True,
                'ncpus': multiprocessing.cpu_count(),
                'omni_url': 'ftp://virbo.org/QinDenton/hour/merged/latest/WGhour-latest.d.zip',
                'leapsec_url': 'ftp://maia.usno.navy.mil/ser7/tai-utc.dat',
                'psddata_url': 'http://spacepy.lanl.gov/repository/psd_dat.sqlite',
                }
    cp = ConfigParser.SafeConfigParser(defaults)
    try:
        successful = cp.read([rcfile])
    except ConfigParser.Error:
        successful = []
    if successful: #New file structure
        ENABLE_DEPRECATION_WARNING = cp.getboolean('spacepy', 'enable_deprecation_warning')
        NCPUS = cp.getint('spacepy', 'ncpus')
        OMNI_URL = cp.get('spacepy', 'omni_url')
        LEAPSEC_URL = cp.get('spacepy', 'leapsec_url')
        PSDDATA_URL = cp.get('spacepy', 'psddata_url')
    else:
        output_vals = {} #default to writing nothing
        if os.path.exists(rcfile):
            with open(rcfile, 'r') as f: #old file structure, read and convert
                for l in f:
                    if l.count('=') != 1:
                        continue
                    name, val = l.split('=')
                    name = name.strip(' "\'\n')
                    val = val.strip(' "\'\n')
                    if not name in ('ENABLE_DEPRECATION_WARNING', 'NCPUS',
                                    'OMNI_URL', 'LEAPSEC_URL', 'PSDDATA_URL'):
                        continue
                    output_vals[name.lower()] = val
                    if name == 'ENABLE_DEPRECATION_WARNING':
                        ENABLE_DEPRECATION_WARNING = val
                    elif name == 'NCPUS':
                        NCPUS = val
                    elif name == 'OMNI_URL':
                        OMNI_URL = val
                    elif name == 'LEAPSEC_URL':
                        LEAPSEC_URL = val
                    elif name == 'PSDDATA_URL':
                        PSDDATA_URL = val
        cp = ConfigParser.SafeConfigParser()
        cp.add_section('spacepy')
        for k in output_vals:
            cp.set('spacepy', k, str(output_vals[k]))
        with open(rcfile, 'wb') as cf:
            cp.write(cf)

from os import environ as ENVIRON
if 'SPACEPY' in ENVIRON:
    DOT_FLN = os.path.join(ENVIRON['SPACEPY'], '.spacepy')
else:
    if 'HOME' in ENVIRON:
        DOT_FLN = os.path.join(ENVIRON['HOME'], '.spacepy')
    else:
        DOT_FLN = os.path.join(ENVIRON['HOMEDRIVE'],
                               ENVIRON['HOMEPATH'],
                               '.spacepy')
rcfile = os.path.join(DOT_FLN, 'spacepy.rc')
if not os.path.exists(DOT_FLN):
    import shutil, sys
    from . import toolbox
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data')
    dataout = os.path.join(DOT_FLN, 'data')
    os.mkdir(DOT_FLN)
    os.chmod(DOT_FLN, 0o777)
    shutil.copy(os.path.join(datadir, 'spacepy.rc'), DOT_FLN)
    os.mkdir(dataout)
    os.chmod(dataout, 0o777)
    shutil.copy(os.path.join(datadir, 'tai-utc.dat'), dataout)
    print('SpacePy data installed to ' + DOT_FLN)
    print('If you wish to start fresh in the future, delete this directory.')
    _read_config(rcfile)
    if sys.version_info[0] < 3 and \
           toolbox.query_yes_no("\nDo you want to update OMNI database and leap seconds table? (Internet connection required)", default = "no") == 'yes':
        toolbox.update()
    print('Regular OMNI updates are recommended: spacepy.toolbox.update()')
    print('Thanks for using SpacePy!')
else:
    _read_config(rcfile)

#Set up a filter to always warn on deprecation
try:
    ENABLE_DEPRECATION_WARNING
except:
    ENABLE_DEPRECATION_WARNING = True
if ENABLE_DEPRECATION_WARNING:
    warnings.filterwarnings('default', '', DeprecationWarning,
                            '^spacepy', 0, False)
