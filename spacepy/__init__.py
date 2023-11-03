#!/usr/bin/python
# -*- coding: utf-8 -*-

"""SpacePy: Space Science Tools for Python

SpacePy is a package of tools primarily aimed at the space science community.
This __init__.py file sets the parameters for import statements.

If running the ipython shell, simply type '?' after any command for help.
ipython also offers tab completion, so hitting tab after '<module name>.'
will list all available functions, classes and variables.

Detailed HTML documentation is available online by typing:

    >>> spacepy.help()

Most functionality is in spacepy's submodules. Each module has specific
help available:

.. rubric:: Submodules

.. autosummary::
    :template: clean_module.rst

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

.. rubric:: Functions

.. autosummary::

    deprecated
    help

.. autofunction:: deprecated
.. autofunction:: help

Copyright 2010-2016 Los Alamos National Security, LLC.
"""

import configparser
import errno
import functools
import multiprocessing
import os
import os.path
import re
import shutil
import sys
import warnings
import webbrowser


def help():
    """Launches web browser with SpacePy documentation

    Online help is always for the latest release of SpacePy.
    """
    print('Opening docs for latest release. Installed SpacePy is {}.'.format(
        __version__))
    webbrowser.open('https://spacepy.github.io/')


# put modules here that you want to be accessible through 'from spacepy import *'
__all__ = ["seapy", "toolbox", "poppy", "coordinates", "time", "omni", 
           "irbempy", "empiricals", "radbelt", "data_assimilation", "pycdf",
           "datamanager", "datamodel", "ae9ap9"]

# Make sure the mingw runtime libs from our binary wheel are findable
minglibs = os.path.join(os.path.dirname(__file__), 'mingw')
if sys.platform == 'win32' and os.path.isdir(minglibs):
    if os.environ.get('PATH'):
        if not minglibs in os.environ['PATH']:
            os.environ['PATH'] += (';' + minglibs)
    else:  # empty or nonexistent PATH
        os.environ['PATH'] = minglibs
    try:
        os.add_dll_directory(minglibs)
    except AttributeError:  # Python 3.8+ only
        pass

#actual deprecation decorator
def _deprecator(version, message, docstring, func):
    #this is the actual, deprecated function
    @functools.wraps(func)
    def _deprecated(*args, **kwargs):
        warnings.warn(message, DeprecationWarning)
        return func(*args, **kwargs)
    if func.__doc__ is None:
        doclines = []
    else:
        doclines = func.__doc__.split('\n')
    # Docstring SHOULD be a single non-blank line with summary, blank line,
    # and then the rest of the content. Want to put the deprecation
    # information just before the first line of "rest of content"
    isblank = [not bool(l.strip()) for l in doclines]
    # All places with a non-blank line following a blank
    paragraphs =  [i for i in range(len(isblank))
                   if not isblank[i] and
                   (i == 0 or isblank[i - 1])]
    if not paragraphs: # No non-blank, guess indentation, insert at end
        leading = '    '
        insert_at = len(doclines) #insert at end
    elif len(paragraphs) == 1: # Get indent from only para, insert at end
        l = doclines[paragraphs[0]]
        leading = l[:len(l) - len(l.lstrip())]
        insert_at = len(doclines)
    else: # Get indent from 2nd paragraph, insert just before it
        l = doclines[paragraphs[1]]
        leading = l[:len(l) - len(l.lstrip())]
        # Insert before blank line before the paragraph.
        insert_at = paragraphs[1] - 1
    to_insert = [
        '',
        leading + '.. deprecated:: ' + version,
    ] \
    + [leading + '   ' + d for d in docstring.split('\n')]
    doclines[insert_at:insert_at] = to_insert
    _deprecated.__doc__ = '\n'.join(doclines)
    return _deprecated

def deprecated(version, message, docstring=None):
    """Decorator to deprecate a function/method

    Modifies a function so that calls to it raise
    ``DeprecationWarning`` and the docstring has a deprecation
    note added in accordance with `numpydoc format
    <https://numpydoc.readthedocs.io/en/latest/format.html#sections>`_

    Parameters
    ==========
    version : str
        What is the first version where this was deprecated?

    message : str
        Message to include in the deprecation warning and in the
        docstring.

    Other Parameters
    ================
    docstring : str

        .. versionadded:: 0.2.2

        If specified, ``docstring`` will be added to the modified function's
        docstring instead of ``message`` (which will only be used in the
        deprecation warning.) It can be a multi-line string (separated with
        ``\\n``). It will be indented to match the existing docstring.

    Notes
    =====
    On Python 2, the deprecated function's signature won't be preserved.
    The function will work but will not have proper argument names listed
    in e.g. ``help``.

    This warning will show as coming from SpacePy, not the deprecated
    function.

    Examples
    ========
        >>> import spacepy
        >>> @spacepy.deprecated('0.2.1', 'Use a different function instead',
        ...                     docstring='A different function is better\\n'
        ...                               'because of reasons xyz')
        ... def foo(x):
        ...     '''This is a test function
        ...
        ...     It may do many useful things.
        ...     '''
        ...     return x + 1
        >>> help(foo)
        Help on function foo in module __main__:
        foo(x)
            This is a test function
            .. deprecated:: 0.2.1
               A different function is better
               because of reasons xyz
            It may do many useful things.
        >>> foo(2)
        DeprecationWarning: Use a different function instead
        3
    """
    message = str(message)
    version = str(version)
    if docstring is None:
        docstring = message
    return functools.partial(_deprecator, version, message, docstring)

# Expose definitions from modules in this package.
# since datamodel depends on top level, delay the variable binding
from . import datamodel
dmarray = datamodel.dmarray
SpaceData = datamodel.SpaceData

#package info
__version__ = 'UNRELEASED'
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
            f.write(("#SpacePy {ver} default {key}: {value}\n"
                     "#{key}: {value}\n").format(
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
                'keepalive': str(True),
                'ncpus': str(multiprocessing.cpu_count()),
                'qindenton_url': 'http://virbo.org/ftp/QinDenton/hour/merged/latest/WGhour-latest.d.zip',
                'qd_daily_url': 'https://rbsp-ect.newmexicoconsortium.org/data_pub/QinDenton/',
                'omni2_url': 'https://spdf.gsfc.nasa.gov/pub/data/omni/omni_cdaweb/hourly/',
                'leapsec_url': 'https://maia.usno.navy.mil/ser7/tai-utc.dat',
                'psddata_url': 'http://spacepy.lanl.gov/repository/psd_dat.sqlite',
                'support_notice': str(True),
                'enable_old_data_warning': str(True),
                }
    #Functions to cast a config value; if not specified, value is a string
    str2bool = lambda x: x.lower() in ('1', 'yes', 'true', 'on')
    caster = {'enable_deprecation_warning': str2bool,
              'keepalive': str2bool,
              'ncpus': int,
              'support_notice': str2bool,
              'enable_old_data_warning': str2bool,
              }
    cp = configparser.ConfigParser(defaults)
    try:
        successful = cp.read([rcfile])
    except configparser.Error:
        successful = []
    if successful:  # New file structure
        try:
            config = dict(cp.items('spacepy'))
        except configparser.NoSectionError:
            successful = []
            config = {}
    if not successful:  # Old or bad file structure, wipe it out
        cp = configparser.ConfigParser()
        cp.add_section('spacepy')
        with open(rcfile, 'w') as cf:
            cp.write(cf)
        for k in defaults:
            if not k in config:
                config[k] = defaults[k]
    _write_defaults(rcfile, defaults)
    for k in caster:
        config[k] = caster[k](config[k])


def _find_spacepy_dir():
    """Determine the .spacepy directory (DOT_FLN)

    This does not create the directory (although will create the parent
    directory if ``SPACEPY`` environment variable is set).

    Returns
    -------
    DOT_FLN : str
        Full path to the .spacepy directory.
    """
    if 'SPACEPY' in os.environ:
        parentdir = os.path.abspath(os.path.expanduser(os.environ['SPACEPY']))
        if not os.path.exists(parentdir):
            try:
                os.makedirs(parentdir)
            except OSError as e:
                if e.errno != errno.EEXIST:  # FileExistsError in 3.3+
                    raise
    elif 'HOME' in os.environ:
        parentdir = os.environ['HOME']
    elif 'HOMEDRIVE' in os.environ and 'HOMEPATH' in os.environ:
        parentdir = os.path.join(
            os.environ['HOMEDRIVE'], os.environ['HOMEPATH'])
    else:
        parentdir = os.path.expanduser('~')
    return os.path.join(parentdir, '.spacepy')


def _populate_spacepy_dir(DOT_FLN):
    """Create .spacepy directory and populate.

    Makes sure created and data directory exists.

    Parameters
    ----------
    DOT_FLN : str
        Full path to the .spacepy directory.
    """
    if not os.path.exists(DOT_FLN):
        try:
            os.mkdir(DOT_FLN)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    dataout = os.path.join(DOT_FLN, 'data')
    if not os.path.exists(dataout):
        datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'data')
        try:
            os.mkdir(dataout)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        shutil.copy2(os.path.join(datadir, 'tai-utc.dat'), dataout)
        # Set permissions based on umask, not perms in SpacePy install
        u = os.umask(0)
        os.umask(u)
        os.chmod(os.path.join(dataout, 'tai-utc.dat'), 0o666 & ~u)


DOT_FLN = _find_spacepy_dir()
_populate_spacepy_dir(DOT_FLN)
rcfile = os.path.join(DOT_FLN, 'spacepy.rc')
_read_config(rcfile)

if __version__ == 'UNRELEASED' and config['support_notice']:
    print('This unreleased version of SpacePy is not supported '
          'by the SpacePy team.')
# Set up a filter to always warn on deprecation
if config['enable_deprecation_warning']:
    warnings.filterwarnings('default', '', DeprecationWarning,
                            '^spacepy', 0, False)
