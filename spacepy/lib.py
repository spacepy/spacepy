#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Support for fast C-based versions of spacepy routines.

To add functions, put a .c file containing the desired functions in
the libspacepy directory. Any functions for export to python should
be added to the functions dictionary below. These functions will
be callable as attributes of spacepy.lib. If the C library is not
found, spacepy.lib.have_libspacepy will be False. See Poppy
for an example.

Authors: Jon Niehof
Institution: Los Alamos National Laboratory
Contact: jniehof@lanl.gov

Copyright 2010-2014 Los Alamos National Security, LLC.
"""

import ctypes
import os.path
import sys
import sysconfig

import numpy

__contact__ = 'Jon Niehof, jniehof@lanl.gov'

#Handy shortcut types
dptr = ctypes.POINTER(ctypes.c_double)
"""alias of LP_c_double"""
ulptr = ctypes.POINTER(ctypes.c_ulong)
"""alias of LP_c_ulong"""
lptr = ctypes.POINTER(ctypes.c_long)
"""alias of LP_c_long"""

#Dictionary of function signatures
#Key is name of function, value is a list of types
#First element in list is return type; rest are parameter types.
functions = {
    'boots': [None, dptr, dptr, ctypes.c_ulong,
              ctypes.c_ulong, ctypes.c_ulong,
              ctypes.c_ulong, ctypes.c_int],
    'assoc': [None, dptr, dptr, dptr, lptr,
              ctypes.c_double, ctypes.c_long, ctypes.c_long,
              ctypes.c_long],
    'aa_ci': [None, ulptr, ulptr, ctypes.c_ulong, ctypes.c_ulong,
              ctypes.c_ulong, ulptr, ctypes.c_int],
    'solve_cn': [None, dptr, dptr, dptr, dptr, dptr, dptr,
                 ctypes.c_double, ctypes.c_int, dptr],
    'hypot_tb': [ctypes.c_double,
                 numpy.ctypeslib.ndpointer(
                     dtype=ctypes.c_double, flags='C_CONTIGUOUS'),
                 ctypes.c_long],
    'cEuler': [ctypes.c_int, #return value, number of points used
               ctypes.c_int, ctypes.c_int, #grid size
               ctypes.c_int, ctypes.c_double, #maxsteps, step size
               ctypes.c_double, ctypes.c_double, #x,y of start
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x grid
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y grid
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x field
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y field
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x of stream
               numpy.ctypeslib.ndpointer(
                   dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y of stream
               ],
    'cRk4': [ctypes.c_int, #return value, number of points used
             ctypes.c_int, ctypes.c_int, #grid size
             ctypes.c_int, ctypes.c_double, #maxsteps, step size
             ctypes.c_double, ctypes.c_double, #x,y of start
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x grid
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y grid
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x field
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y field
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # x of stream
             numpy.ctypeslib.ndpointer(
                 dtype=ctypes.c_double, flags='C_CONTIGUOUS'), # y of stream
             ],
    }


def load_lib():
    """Find and load libspacepy

    Normally this will be in the directory where spacepy is installed,
    under libspacepy.

    Returns
    =======
    library : ctypes.CDLL or ctypes.winDLL
        the open library
    """
    libdir = os.path.dirname(os.path.abspath(__file__))
    libnames = {
        'win32': ['libspacepy.dll.a', 'spacepy.dll'],
        'darwin':  ['libspacepy.dylib', 'libspacepy.so',
                    'spacepy.dylib', 'spacepy.so'],
        }.get(sys.platform, ['libspacepy.so'])
    ext = sysconfig.get_config_var('EXT_SUFFIX')
    if ext is None:
        ext = sysconfig.get_config_var('SO')
    if ext:
        libnames.append('libspacepy' + ext)
        libnames.append('spacepy' + ext)
    libpaths = [os.path.join(libdir, n) for n in libnames]
    libpaths = [p for p in libpaths if os.path.exists(p)]
    for p in libpaths:
        try:
            return ctypes.CDLL(p)
        except OSError:
            pass
    return None


def load_call_dict(call_dict, lib):
    """Loads argument/return types from the call dictionary

    Parameters
    ==========
    call_dict : dict
        call dictionary. Keyed by function name;
        values are [return type, argtype0, argtype 1...]
    
    lib : ctypes.CDLL
        library where functions specified in ``call_dict`` live.

    """
    for funcname in call_dict:
        func = getattr(lib, funcname)
        args = call_dict[funcname]
        func.restype = args[0]
        if len(args) <= 1:
            func.argtypes = None
        else:
            func.argtypes = args[1:]


library = load_lib()
if library == None:
    have_libspacepy = False
else:
    have_libspacepy = True
    load_call_dict(functions, library)
    for f in functions:
        globals()[f] = getattr(library, f)
    #Clean up no-longer-useful bits
    del f, load_call_dict, load_lib, functions, ctypes, os, sys
