import argparse
import filecmp
import glob
import os
import subprocess


def not_ascii(string):
    """Check if the characters in string s are not ASCII, U+0-U+7F."""
    return len(string) != len(string.encode())


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--irbem_dir', dest='irbem_dir',
                    default='../../../IRBEM/')
opts = parser.parse_args()

# Compare fortran sources
print('Comparing source directory ({})\n'.format(os.path.join(opts.irbem_dir,
                                                              'source')))
ffiles = glob.glob(os.path.join(opts.irbem_dir,
                                'source', '*.f'))
ifiles = glob.glob(os.path.join(opts.irbem_dir,
                                'source', '*.inc'))
cfiles = glob.glob(os.path.join(opts.irbem_dir,
                                'source', '*.c'))
combined_files = ffiles + ifiles + cfiles

for ff in combined_files:
    local = os.path.split(ff)[-1]
    wrapdir = '../../spacepy/irbempy'
    cand = os.listdir(wrapdir)
    spaceirb = [cc for cc in cand if 'irbem-lib-' in cc][0]
    spacefn = os.path.join(wrapdir, spaceirb, 'source/', local)
    status = filecmp.cmp(ff, spacefn, shallow=False)
    if not status:
        print('Compared {}: Files are different'.format(local))
        # Test for unicode and flag lines that need characters changing
        # only test if the file differs from the spacepy version as
        # anything in our repo is assumed to be fine with f2py
        with open(ff, 'r', errors='ignore') as fh:
            cont = fh.read()
        cont = cont.split('\n')
        for lineno, line in enumerate(cont, start=1):
            unic_present = not_ascii(line)
            if unic_present:
                print('Line {} of {} has unicode characters'.format(lineno, ff))

# suggest new name for irbem source dir and remind me where I need to change it
shorthash = subprocess.getoutput('cd {}; git rev-parse --short HEAD'.format(opts.irbem_dir))
comdate = subprocess.getoutput(r'cd {}; git log -1 --format=%cd --date=short'.format(opts.irbem_dir))
comdate = comdate.replace('-', '')
newirb = 'irbem-lib-{}-{}'.format(comdate, shorthash)
print('Change {}\n    to {}'.format(os.path.join(wrapdir, spaceirb),
                               os.path.join(wrapdir, newirb)))