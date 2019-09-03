#!/usr/bin/env python

import sys

import spacepy.pycdf
import spacepy.pycdf.istp


def main(f):
    """Do all checks on a file

    :param str f: filename to check
    """
    with spacepy.pycdf.CDF(f) as cdffile:
        errs = spacepy.pycdf.istp.FileChecks.all(cdffile)
    if errs:
        print('\n'.join(errs))


if __name__ == '__main__':
    main(*sys.argv[1:])
