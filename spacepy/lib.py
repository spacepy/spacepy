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

Copyright 2010 Los Alamos National Security, LLC.
"""

import ctypes
import os.path
import sys

import numpy as np

__contact__ = 'Jon Niehof, jniehof@lanl.gov'

#Handy shortcut types
dptr = ctypes.POINTER(ctypes.c_double)
ulptr = ctypes.POINTER(ctypes.c_ulong)
lptr = ctypes.POINTER(ctypes.c_long)

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
                 np.ctypeslib.ndpointer(dtype=ctypes.c_double, flags='C_CONTIGUOUS'),
                 ctypes.c_long]
    }

def load_lib():
    """Find and load libspacepy

    Normally this will be in the directory where spacepy is installed,
    under libspacepy.

    @return: the open library
    @rtype: ctypes.CDLL or ctypes.WinDLL
    """
    libdir = os.path.dirname(os.path.abspath(__file__))
    if sys.platform == 'win32':
        libpath = os.path.join(libdir, 'spacepy.dll')
    elif sys.platform == 'darwin':
        for n in ('libspacepy.dylib', 'libspacepy.so',
                  'spacepy.dylib', 'libspacepy.so'):
            libpath = os.path.join(libdir, n)
            if os.path.exists(libpath):
                break
    else:
        libpath = os.path.join(libdir, 'libspacepy.so')
    if not os.path.exists(libpath):
        return None
    return ctypes.CDLL(libpath)


def load_call_dict(call_dict, lib):
    """Loads argument/return types from the call dictionary

    @param call_dict: call dictionary. Keyed by function name;
                      values are [return type, argtype0, argtype 1...]
    @type call_dict: dict
    @param lib: library where functions specified in L{call_dict} live.
    @type lib: ctypes.WinDLL or ctypes.CDLL
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
