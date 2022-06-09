#!/usr/bin/env python

"""Test the loading of the CDF library and calls to each function"""

import ctypes
import os
import subprocess
import sys


# 32-bit arm needs special treatment for variadic floats
typepun = os.uname()[4].startswith('arm') and sys.maxsize <= 2 ** 32

call_dict = {
    'breakdownTT2000': [None, ctypes.c_longlong]\
        + [ctypes.POINTER(ctypes.c_double)] * 3,
    'CDF_TT2000_from_UTC_EPOCH': [ctypes.c_longlong, ctypes.c_double],
    'CDF_TT2000_from_UTC_EPOCH16': [ctypes.c_longlong,
                                    ctypes.POINTER(ctypes.c_double * 2)],
    'CDF_TT2000_to_UTC_EPOCH': [ctypes.c_double, ctypes.c_longlong],
    'CDF_TT2000_to_UTC_EPOCH16': [ctypes.c_double, ctypes.c_longlong,
                                  ctypes.POINTER(ctypes.c_double * 2)],
    'CDFlib': [ctypes.c_long, ctypes.c_long],
    'CDFsetFileBackward': [None, ctypes.c_long],
    'computeEPOCH': [ctypes.c_double] + [ctypes.c_long] * 7,
    'computeEPOCH16': [ctypes.c_double] + [ctypes.c_long] * 10\
        + [ctypes.POINTER(ctypes.c_double * 2)],
    'computeTT2000': [ctypes.c_longlong]\
        + [ctypes.c_longlong if typepun else ctypes.c_double] * 3,
    'EPOCH16breakdown': [None, ctypes.c_double * 2]\
        + [ctypes.POINTER(ctypes.c_long)] * 10,
    'EPOCHbreakdown': [ctypes.c_long, ctypes.c_double]\
        + [ctypes.POINTER(ctypes.c_long)] * 7,
}


def load_lib():
    # Full path to the dylib. Must be 3.7.1
    libpath = "/usr/local/lib/libcdf.so"
    lib = ctypes.CDLL(libpath)

    for funcname in call_dict:
        func = getattr(lib, funcname)
        args = call_dict[funcname]
        func.restype = args[0]
        func.argtypes = None if len(args) <= 1 else args[1:]
    return lib


def cast_ll(x):
    """Reinterpret-cast a float to a long long"""
    return ctypes.cast(ctypes.pointer(ctypes.c_double(x)),
                       ctypes.POINTER(ctypes.c_longlong)).contents


def test_breakdownTT2000(lib):
    yyyy = ctypes.c_double(0)
    mm = ctypes.c_double(0)
    dd = ctypes.c_double(0)
    hh = ctypes.c_double(0)
    mn = ctypes.c_double(0)
    sec = ctypes.c_double(0)
    msec = ctypes.c_double(0)
    usec = ctypes.c_double(0)
    nsec = ctypes.c_double(0)
    tt2000 = 315673512171654000
    lib.breakdownTT2000(tt2000, yyyy, mm, dd, ctypes.byref(hh), ctypes.byref(mn), ctypes.byref(sec), ctypes.byref(msec), ctypes.byref(usec), ctypes.byref(nsec))
    print('Expect 2010-1-2 3:4:5.987654000')
    print(
        'Actual {:.0f}-{:.0f}-{:.0f} {:.0f}:{:.0f}:{:.0f}'\
        '.{:03.0f}{:03.0f}{:03.0f}'
        .format(
            yyyy.value, mm.value, dd.value, hh.value, mn.value, sec.value,
            msec.value, usec.value, nsec.value))

    
def test_CDF_TT2000_from_UTC_EPOCH(lib):
    epoch = 63429523200000.0
    res = lib.CDF_TT2000_from_UTC_EPOCH(epoch)
    print('Expect 315576066184000000')
    print('Actual {}'.format(res))


def test_CDF_TT2000_from_UTC_EPOCH16(lib):
    epoch = (63429523200.0, 1000.0)
    res = lib.CDF_TT2000_from_UTC_EPOCH16((ctypes.c_double * 2)(*epoch))
    print('Expect 315576066184000001')
    print('Actual {}'.format(res))


def test_CDF_TT2000_to_UTC_EPOCH(lib):
    tt = 315576066184000000
    res = lib.CDF_TT2000_to_UTC_EPOCH(tt)
    print('Expect 63429523200000.0')
    print('Actual {}'.format(res))


def test_CDF_TT2000_to_UTC_EPOCH16(lib):
    tt = 315576066184000001    
    epoch16 = (ctypes.c_double * 2)(-1., -1.)
    res = lib.CDF_TT2000_to_UTC_EPOCH16(tt, epoch16)
    print('Expect (63429523200.0, 1000.0)')
    print('Actual ({}, {})'.format(epoch16[0], epoch16[1]))


def test_CDFlib(lib):
    ver = ctypes.c_long(0)
    rel = ctypes.c_long(0)
    GET_ = ctypes.c_long(1007)
    LIB_VERSION_ = ctypes.c_long(21)
    LIB_RELEASE_ = ctypes.c_long(22)
    NULL_ = ctypes.c_long(1000)
    lib.CDFlib(GET_, LIB_VERSION_, ctypes.byref(ver),
               GET_, LIB_RELEASE_, ctypes.byref(rel),
               NULL_)
    print('Version {}.{}'.format(ver.value, rel.value))


def test_CDFsetFileBackward(lib):
    lib.CDFsetFileBackward(1)


def test_computeEPOCH(lib):
    epoch = lib.computeEPOCH(2010, 1, 2,
                             3, 4, 5,
                             987)
    print('Expect 63429620645987.0')
    print('Actual {}'.format(epoch))


def test_computeEPOCH16(lib):
    epoch16 = (ctypes.c_double * 2)(-1., -1.)
    epoch = lib.computeEPOCH16(2010, 1, 2,
                               3, 4, 5,
                               987, 654, 0, 0,
                               epoch16)
    print('Expect (63429620645.0, 987654000000.0)')
    print('Actual ({}, {})'.format(*epoch16))


def test_computeTT2000(lib):
    args = (2010, 1, 2, 3, 4, 5, 987, 654, 0)
    tt = lib.computeTT2000(*[(cast_ll if typepun else ctypes.c_double)(a)
                             for a in args])
    print('Expect 315673512171654000')
    print('Actual {}'.format(tt))


def test_EPOCH16breakdown(lib):
    epoch16 = (ctypes.c_double * 2)(63429620645.0, 987654321000.0)
    yyyy = ctypes.c_long(0)
    mm = ctypes.c_long(0)
    dd = ctypes.c_long(0)
    hh = ctypes.c_long(0)
    mn = ctypes.c_long(0)
    sec = ctypes.c_long(0)
    msec = ctypes.c_long(0)
    usec = ctypes.c_long(0)
    nsec = ctypes.c_long(0)
    psec = ctypes.c_long(0)
    lib.EPOCH16breakdown(epoch16, yyyy, mm, dd, hh, mn, sec, msec, usec, nsec, psec)
    print('Expect 2010-1-2 3:4:5.987654321000')
    print(
        'Actual {:.0f}-{:.0f}-{:.0f} {:.0f}:{:.0f}:{:.0f}'\
        '.{:03.0f}{:03.0f}{:03.0f}{:03.0f}'
        .format(
            yyyy.value, mm.value, dd.value, hh.value, mn.value, sec.value,
            msec.value, usec.value, nsec.value, psec.value))


def test_EPOCHbreakdown(lib):
    epoch = 63429620645987.
    yyyy = ctypes.c_long(0)
    mm = ctypes.c_long(0)
    dd = ctypes.c_long(0)
    hh = ctypes.c_long(0)
    mn = ctypes.c_long(0)
    sec = ctypes.c_long(0)
    msec = ctypes.c_long(0)
    lib.EPOCHbreakdown(epoch, yyyy, mm, dd, hh, mn, sec, msec)
    print('Expect 2010-1-2 3:4:5.987')
    print(
        'Actual {:.0f}-{:.0f}-{:.0f} {:.0f}:{:.0f}:{:.0f}'\
        '.{:03.0f}'
        .format(
            yyyy.value, mm.value, dd.value, hh.value, mn.value, sec.value,
            msec.value))


def main():
    if len(sys.argv) == 1:
        for k in sorted(call_dict, key=str.lower):
            print(k)
            res = subprocess.run([sys.executable, __file__, k])
            print('{} {}\n'.format(
                k, 'failed' if res.returncode else 'completed'))
    else:
        k = sys.argv[1]
        f = globals().get('test_{}'.format(k), None)
        if f is None:
            print('skipping {}, no test'.format(k))
        else:
            lib = load_lib()
            f(lib)
        

if __name__ == '__main__':
    main()
