#!/usr/bin/env python

import argparse

import spacepy.pycdf
import spacepy.pycdf.istp


def main(f):
    """Do all checks on a file

    :param str f: filename to check
    """
    with spacepy.pycdf.CDF(f) as cdffile:
        errs = spacepy.pycdf.istp.FileChecks.all(cdffile, catch=True)
    if errs:
        print('\n'.join(errs))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Check a CDF file for ISTP compliance.")
    parser.add_argument('cdffile', action='store',
                        help='Path to CDF file to check')
    args = parser.parse_args()
    main(args.cdffile)
