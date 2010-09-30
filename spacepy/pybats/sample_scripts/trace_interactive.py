#!/usr/bin/env python
# trace_interactive.py

'''
Open a BATS-R-US output file into an interactive plot window to 
examine log(pressure) and magnetic field.

Syntax:
trace_interactive [options] filename

Options:
-h or -help:
     Print this help.

-r=[resolution]:
     If input file has not been resampled to a regular grid, regrid using
     the given resolution.  Default is 0.5 Earth Radii.  This may be given
     as a fraction, e.g. 1.0/8.0.
'''

# Set defaults.
res = 0.5
infile = False

# Handle options.
import sys
import glob
for option in sys.argv[1:]:
    # Handle options:
    if option[0] == '-':
        if option[0:2] == '-r':
            try:
                res = float(option[3:])
            except ValueError:
                res = eval(option[3:])
        if option == '-h' or option == '-help':
            print(__doc__)
            exit()
    # Search for files that match option.
    else:
        infile = option

import math

import spacepy.pybats.bats as pb
import spacepy.pybats.interact as int
from spacepy.pybats import load

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import (LogNorm, Normalize)
from matplotlib.ticker import (LogLocator, LogFormatter, LogFormatterMathtext, 
                               MultipleLocator)

# Read file.
i = infile.rfind('.')
if infile[i:] == '.pb':     #PyBats Save file
    mhd = load(infile)
elif infile[i:] == '.out':  #IDL binary file.
    mhd = pb.Bats2d(infile)
    print('Regridding to %.3f Earth Radii' % res)
    mhd.regrid(res, [-35, 15], [-30,30])

# Set up color scales.
# Pressure log:
minlog = 0.01
maxlog = 50.0
lev_exp = np.arange(np.log10(minlog), np.log10(maxlog), 0.1)
levs = np.power(10, lev_exp)

# Initialize plot.
fig = plt.figure(figsize=(10,7))
#fig.subplots_adjust(left=0.08, right=0.99, wspace=0.10, hspace=0.18)
hour = math.floor(mhd.time / 3600.0)
mint = (mhd.time / 3600.0 - hour)*60.0
time = '%04ih %05.2fm' % (hour, mint)

fig.text(0.5, 0.01, 'File = '+mhd.filename, ha='center')

# Pressure plot.
ax1 = fig.add_subplot(111)
mhd.data['p_log'] = np.where(mhd.data['p']>minlog, 
                             mhd.data['p'], 1.01*minlog)
cont = mhd.contourf(ax1, 'x','z','p_log',levs, norm=LogNorm() )
cbar = plt.colorbar(cont, pad=0.01, ticks=LogLocator(), 
                    format=LogFormatterMathtext())
cbar.set_label('Pressure ($nPa$)')
mhd.add_body(ax1)
ax1.set_title('Results at iter=%i, t=%s' % (mhd.iter, time))

ax1.set_xlabel('GSM X ($R_{E}$)')
ax1.set_ylabel('GSM Z ($R_{E}$)')
ax1.set_ylim([-20,20])
ax1.set_xlim([-30,15])
ax1.invert_xaxis()

tracer = int.ClickTracer(ax1, mhd, 'bx', 'bz')

plt.show()
