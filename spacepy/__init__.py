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

Most functionality is in spacepy's submodules. Each module has specific
help available:

    coordinates
    data_assimilation
    datamodel
    empiricals
    irbempy
    LANLstar
    omni
    poppy
    pybats
    pycdf
    radbelt
    seapy
    time
    toolbox

Copyright 2010-2012 Los Alamos National Security, LLC.
"""

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import functools
import multiprocessing
import os
import os.path
import sys
import warnings

def help():
    """Launches web browser with local HTML help"""
    
    import webbrowser
    fspec = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'Doc', 'index.html')
    if not os.path.exists(fspec):
        fspec = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                             'Doc', 'build', 'html', 'index.html')
    if os.path.exists(fspec):
        webbrowser.open('file://' + fspec)
    else:
        print("Can't find help files in {0}".format(__path__[0]))

# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
    "irbempy", "empiricals", "radbelt", "data_assimilation"]

# on windows, make sure the Fortran libs are findable
if sys.platform == 'win32':
    fortlibs = os.path.join(os.path.dirname(__file__), 'mingw')
    if 'PATH' in os.environ:
        if not fortlibs in os.environ['PATH']:
            if os.environ['PATH']:
                os.environ['PATH'] += (';' + fortlibs)
            else: #empth PATH
                os.environ['PATH'] = fortlibs
    else:
        os.environ['PATH'] = fortlibs

def deprecated(version, message):
    """Decorator to deprecate a function/method

    Parameters
    ==========
    version : str
        What is the first version where this was deprecated?

    message : str
        Message to include in the documentation and the warning message.
    """
    message = str(message)
    version = str(version)
    #actual decorator, with version and message curried in
    def _deprecator(func):
        #this is the actual, deprecated function
        @functools.wraps(func)
        def _deprecated(*args, **kwargs):
            warnings.warn(message, DeprecationWarning)
            return func(*args, **kwargs)
        if func.__doc__ is None:
            doclines = []
        else:
            doclines = func.__doc__.split('\n')
        #first non-blank line
        idx = next((i for i in range(len(doclines)) if doclines[i].strip()),
                   None)
        if idx is None: #no non-blank
            leading = '    '
            idx = len(doclines) #insert at end
        else:
            first = doclines[idx]
            #copy whitespace
            leading = first[0:len(first) - len(first.lstrip())]
            idx += 1 #insert just after first non-blank
        #REVERSE order since insert goes before
        doclines.insert(idx, leading)
        doclines.insert(idx, leading + '   ' + message)
        doclines.insert(idx, leading + '.. deprecated:: ' + version)
        doclines.insert(idx, leading)
        _deprecated.__doc__ = '\n'.join(doclines) + '\n'
        return _deprecated
    return _deprecator

# Expose definitions from modules in this package.
# since datamodel depends on top level, delay the variable binding
from . import datamodel
dmarray = datamodel.dmarray
SpaceData = datamodel.SpaceData

#package info
__version__ = '0.1.4'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley', 'Josef Koller', 'Dan Welling', 'Brian Larsen', 'Jon Niehof', 'Mike Henderson']
__contact__ = 'spacepy@lanl.gov'
__license__ = """SpacePy: Space Science Tools for Python


Copyright 2010-2013 Los Alamos National Security, LLC.
All Rights Reserved.

 This material was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE


 1. This LICENSE AGREEMENT is between the Los Alamos National Security, LLC ("LANS"), and the Individual or Organization ("Licensee") accessing and otherwise using SpacePy 0.1.4 software in source or binary form and its associated documentation.

 2. Subject to the terms and conditions of this License Agreement, LANS hereby grants Licensee a nonexclusive, royalty-free, world-wide license to reproduce, analyze, test, perform and/or display publicly, prepare derivative works, distribute, and otherwise use SpacePy 0.1.4 alone or in any derivative version, provided, however, that LANS' License Agreement and LANS' notice of copyright, i.e., "Copyright (c) 2010 Los Alamos National Security, LLC; All Rights Reserved" are retained in SpacePy 0.1.4 alone or in any derivative version prepared by Licensee.

 3. In the event Licensee prepares a derivative work that is based on or incorporates SpacePy 0.1.4 or any part thereof, and wants to make the derivative work available to others as provided herein, then Licensee hereby agrees to include in any such work a brief summary of the changes made to SpacePy 0.1.4.

 4. LANS is making SpacePy 0.1.4 available to Licensee on an "AS IS" basis. LANS MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED. BY WAY OF EXAMPLE, BUT NOT LIMITATION, LANS MAKES NO AND DISCLAIMS ANY REPRESENTATION OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF SPACEPY 0.1.4 WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.

 5. LANS SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF SPACEPY 0.1.4 FOR ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF MODIFYING, DISTRIBUTING, OR OTHERWISE USING SPACEPY 0.1.4, OR ANY DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.

 6. This License Agreement will automatically terminate upon a material breach of its terms and conditions.

 7. Nothing in this License Agreement shall be deemed to create any relationship of agency, partnership, or joint venture between LANS and Licensee. This License Agreement does not grant permission to use LANS trademarks or trade name in a trademark sense to endorse or promote products or services of Licensee, or any third party.

 8. By copying, installing or otherwise using SpacePy 0.1.4, Licensee agrees to be bound by the terms and conditions of this License Agreement.

 The modified version of IRBEMlib distributed with SpacePy is covered by the Lesser GNU Public License (LGPL).

 The LANLstar module depends on the ffnet package which is distributed under the GNU Public License (GPL). The use of LANLstar is therefore covered by the GPL.
"""

if sys.platform == 'win32':
    __license__ += \
        """
Fortran library support provided by MinGW. The MinGW base runtime package has been placed in the public domain, and is not governed by copyright.
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

#Global spacepy configuration information
config = {}

# import some settings
def _read_config(rcfile):
    """Read configuration information from a file"""
    global config
    defaults = {'enable_deprecation_warning': str(True),
                'ncpus': str(multiprocessing.cpu_count()),
                'qindenton_url': 'ftp://virbo.org/QinDenton/hour/merged/latest/WGhour-latest.d.zip',
                'omni2_url': 'ftp://virbo.org/OMNI/OMNI2/merged/latest/OMNI_OMNI2-latest.cdf.zip',
                'leapsec_url': 'ftp://maia.usno.navy.mil/ser7/tai-utc.dat',
                'psddata_url': 'http://spacepy.lanl.gov/repository/psd_dat.sqlite',
                }
    #Functions to cast a config value; if not specified, value is a string
    str2bool = lambda x: x.lower() in ('1', 'yes', 'true', 'on')
    caster = {'enable_deprecation_warning': str2bool,
              'ncpus': int,
              }
    cp = ConfigParser.SafeConfigParser(defaults)
    try:
        successful = cp.read([rcfile])
    except ConfigParser.Error:
        successful = []
    if successful: #New file structure
        config = dict(cp.items('spacepy'))
    else: #old file structure, wipe it out
        cp = ConfigParser.SafeConfigParser()
        cp.add_section('spacepy')
        with open(rcfile, 'wb') as cf:
            cp.write(cf)
        for k in defaults:
            if not k in config:
                config[k] = defaults[k]
    for k in caster:
        config[k] = caster[k](config[k])

if 'SPACEPY' in os.environ:
    DOT_FLN = os.path.join(os.environ['SPACEPY'], '.spacepy')
else:
    if 'HOME' in os.environ:
        DOT_FLN = os.path.join(os.environ['HOME'], '.spacepy')
    elif 'HOMEDRIVE' in os.environ and 'HOMEPATH' in os.environ:
        DOT_FLN = os.path.join(os.environ['HOMEDRIVE'],
                               os.environ['HOMEPATH'],
                               '.spacepy')
    else:
        DOT_FLN = os.path.expanduser(os.path.join('~', '.spacepy'))
rcfile = os.path.join(DOT_FLN, 'spacepy.rc')
if not os.path.exists(DOT_FLN):
    print("""SpacePy: Space Science Tools for Python
  See __licence__ and __citation__ for licensing, and help() for HTML help.""")
    import shutil, sys
    datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'data')
    dataout = os.path.join(DOT_FLN, 'data')
    os.mkdir(DOT_FLN)
    os.mkdir(dataout)
    shutil.copy(os.path.join(datadir, 'tai-utc.dat'), dataout)
    print('Data and configuration installed to ' + DOT_FLN)
    _read_config(rcfile)
    print('Downloading OMNI database and leap seconds table is recommended:'
          '\n\timport spacepy.toolbox; spacepy.toolbox.update()')
    print('Thanks for using SpacePy!')
else:
    _read_config(rcfile)

#Set up a filter to always warn on deprecation
if config['enable_deprecation_warning']:
    warnings.filterwarnings('default', '', DeprecationWarning,
                            '^spacepy', 0, False)
