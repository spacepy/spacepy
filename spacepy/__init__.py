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

.. autosummary::

    ~spacepy.coordinates
    ~spacepy.data_assimilation
    ~spacepy.datamodel
    ~spacepy.empiricals
    ~spacepy.irbempy
    ~spacepy.LANLstar
    ~spacepy.omni
    ~spacepy.poppy
    ~spacepy.pybats
    ~spacepy.pycdf
    ~spacepy.radbelt
    ~spacepy.seapy
    ~spacepy.time
    ~spacepy.toolbox
    ~spacepy.ae9ap9

Copyright 2010-2016 Los Alamos National Security, LLC.
"""

try:
    import ConfigParser
except ImportError:
    import configparser as ConfigParser
import functools
import multiprocessing
import os
import os.path
import re
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
           "irbempy", "empiricals", "radbelt", "data_assimilation", "pycdf",
           "datamanager", "datamodel", "ae9ap9"]

# on windows, make sure the mingw runtime libs are findable
if sys.platform == 'win32':
    minglibs = os.path.join(os.path.dirname(__file__), 'mingw')
    if 'PATH' in os.environ:
        if not minglibs in os.environ['PATH']:
            if os.environ['PATH']:
                os.environ['PATH'] += (';' + minglibs)
            else: #empth PATH
                os.environ['PATH'] = minglibs
    else:
        os.environ['PATH'] = minglibs

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
__version__ = '0.2.0'
__author__ = 'The SpacePy Team'
__team__ = ['Steve Morley', 'Josef Koller', 'Dan Welling', 'Brian Larsen', 'Jon Niehof', 'Mike Henderson']
__contact__ = 'spacepy@lanl.gov'
__license__ = """SpacePy: Space Science Tools for Python


Copyright 2010 Triad National Security, LLC.
All Rights Reserved.

This is open source software; you can redistribute it and/or modify it under the 
terms of the Python Software Foundation License. If software is modified to
produce derivative works, such modified software should be clearly marked, so
as not to confuse it with the version available from LANL. Full text of the 
Python Software Foundation License can be found in the LICENSE.md file in the
main development branch of the repository (https://github.com/spacepy/spacepy).
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

def _write_defaults(rcfile, defaults, section='spacepy'):
    """Write configuration defaults out to a file if they aren't there"""
    f = open(rcfile, 'r+t') #Avoid race condition, open for read and write
    try:
        startpos = f.tell()
        rclines = f.readlines()
        writeme = [k for k in defaults.keys()]
        #Places where sections start
        secidx = [i for i in range(len(rclines))
                  if re.match(r"^\[[^\[\]]*\]\n$", rclines[i])]
        #Line containing the start of this section
        thissec = next((i for i in secidx
                        if rclines[i] == '[{section}]\n'.format(
                                section=section)), -1)
        if thissec != -1: #this section exists, read it
            #Find the line that starts the next section
            nextsecidx = secidx.index(thissec) + 1
            if nextsecidx < len(secidx):
                nextsec = secidx[nextsecidx]
            else:
                nextsec = None
            #Read only this section for lines matching default values
            present = []
            for l in rclines[thissec:nextsec]:
                #For each key, does the line match, commented or not?
                for k in writeme:
                    if re.match(r'(#\s?)?{key}\s?[:=]'.format(key=k),
                                l):
                        writeme.remove(k)
                        break
        f.seek(startpos) #Back to start of file
        #Scan to just after section header. Text mode, can only seek() to tell()
        if thissec != -1:
            while f.readline() != '[{section}]\n'.format(section=section):
                pass
        if sys.platform == 'win32':
            #Tickle the file for read/write switch
            #https://stackoverflow.com/questions/11176724/python-file-operations
            f.seek(0, 2)
        #Write default values for anything not read
        for k in sorted(writeme):
            f.write("#{key}: {value} #default in SpacePy {ver}\n".format(
                key=k, value=defaults[k], ver=__version__))
        #And write all the remaining lines from the section header to end
        if writeme:
            f.writelines(rclines[thissec+1:])
    finally:
        f.close()


# import some settings
def _read_config(rcfile):
    """Read configuration information from a file"""
    global config
    defaults = {'enable_deprecation_warning': str(True),
                'ncpus': str(multiprocessing.cpu_count()),
                'qindenton_url': 'http://virbo.org/ftp/QinDenton/hour/merged/latest/WGhour-latest.d.zip',
                'omni2_url': 'http://virbo.org/ftp/OMNI/OMNI2/merged/latest/OMNI_OMNI2-latest.cdf.zip',
                'leapsec_url': 'http://maia.usno.navy.mil/ser7/tai-utc.dat',
                'psddata_url': 'http://spacepy.lanl.gov/repository/psd_dat.sqlite',
                'support_notice': str(True),
                'apply_plot_styles': str(True),
                }
    #Functions to cast a config value; if not specified, value is a string
    str2bool = lambda x: x.lower() in ('1', 'yes', 'true', 'on')
    caster = {'enable_deprecation_warning': str2bool,
              'ncpus': int,
              'support_notice': str2bool,
              'apply_plot_styles': str2bool,
              }
    #SafeConfigParser deprecated in 3.2. And this is hideous, but...
    if hasattr(ConfigParser, 'SafeConfigParser'):
        cp_class = ConfigParser.SafeConfigParser
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings(
                'always', 'The SafeConfigParser class has been renamed.*',
                DeprecationWarning,
                '^spacepy$') #configparser lies about source of warnings
            ConfigParser.SafeConfigParser()
        for this_w in w:
            if isinstance(this_w.message, DeprecationWarning):
                cp_class = ConfigParser.ConfigParser
            else:
                warnings.showwarning(this_w.message, this_w.category,
                                     this_w.filename, this_w.lineno,
                                     this_w.file, this_w.line)
    else:
        cp_class = ConfigParser.ConfigParser
    cp = cp_class(defaults)
    try:
        successful = cp.read([rcfile])
    except ConfigParser.Error:
        successful = []
    if successful: #New file structure
        config = dict(cp.items('spacepy'))
    else: #old file structure, wipe it out
        cp = cp_class()
        cp.add_section('spacepy')
        with open(rcfile, 'w') as cf:
            cp.write(cf)
        for k in defaults:
            if not k in config:
                config[k] = defaults[k]
    _write_defaults(rcfile, defaults)
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

if __version__ == 'UNRELEASED' and config['support_notice']:
    print('This unreleased version of SpacePy is not supported '
          'by the SpacePy team.')


#Set up a filter to always warn on deprecation
if config['enable_deprecation_warning']:
    warnings.filterwarnings('default', '', DeprecationWarning,
                            '^spacepy', 0, False)
