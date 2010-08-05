#!/usr/bin/env python
# plot-mhdquicklook.py

'''
This plot creates quick-look, 4x4 plots of BATS output.
Usage: 
     plot_mhdquicklook [options] [files]
     
Files can be whole names or use wildcards such as '*' and '?' for
examining many files at once.

This script can accept raw idl files (*.out, long load times but small) or
saved pybats files (*.pb, faster load times but very large.)  Be sure to
include the extension in your wildcard!

Options:
-h or -help:
      Show this help text.
-b or -batch:
      Peform batch processing; do not print to screen.
      All plots are automatically saved to file.
-f[format code]:
      Set output format (batch mode only.)  Default is png file.
      Options include pdf, ps, or png.
-s or -save:
      If reading raw binary files (*.out), save the regridded Bats2d objects
      using Pybats' save function.  This uses cPickle to save the files.
      Note that this requires much disk space and will increase the run
      time, but saving these files will allow for the use of the load 
      option in the future.
-iono:
      Set the fourth frame (lower right) to "iono" mode, where matching IE
      files are searched for and used to plot ionospheric potential and 
      field-aligned currents.  This will slow down the plotting, but adds
      more information.  The default is to either plot oxygen composition
      (for multispecies runs) or leave this frame blank.

Examples:
     Plot all files in PWD as pdfs in batch mode, 
     pickle (save) the regridded files:
     >>> plot_mhdquicklook -b -s -fpdf y*.out

     Plot all saved files in PWD as PNGs to screen:
     >>> plot_mhdquicklook 
'''

# Default values and settings:
batchmode = False     # Run in batch mode?
usesave   = False     # Use save (pickled) files?
DoIono    = False     # Plot ionospheric stuff.
format = 'png'
files = []

# Handle options.
import sys
import glob
for option in sys.argv[1:]:
    # Handle options:
    if option[0] == '-':
        if option[0:2] == '-b':
            batchmode = True
            import matplotlib
            matplotlib.use('Agg')  # Turn off Tk for speed.
        if option[0:2] == '-f':
            format = option[2:].lower()
        if option == '-h' or option == '-help':
            print __doc__
            exit()
    # Search for files that match option.
    else:
        files = files + glob.glob(option)

# Imports:
import math

import spacepy.pybats.bats as pb
from spacepy.pybats import save, load

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import (LogNorm, Normalize)
from matplotlib.ticker import (LogLocator, LogFormatter, LogFormatterMathtext, 
                               MultipleLocator)

def ax_details(ax):
    '''
    Quickly set plot details common to each frame.
    '''
    loc = ax.get_geometry()
    if loc[2] == 1:
        ax.set_ylabel('GSM Z ($R_{E}$)')
    elif loc[2] == 3:
        ax.set_xlabel('GSM X ($R_{E}$)')
        ax.set_ylabel('GSM Z ($R_{E}$)')
    elif loc[2] == 4:
        ax.set_xlabel('GSM X ($R_{E}$)')
    ax.set_ylim([-25,25])
    ax.set_xlim([-70,15])
    ax.invert_xaxis()

# Set up color scales.
# Pressure log:
minlog = 0.01
maxlog = 50.0
lev_exp = np.arange(np.log10(minlog), np.log10(maxlog), 0.1)
levs = np.power(10, lev_exp)
# Number density or Temperature
minden = 0.0
maxden = 50.0
denrng = Normalize(vmin=minden, vmax=maxden)
denlev = np.linspace(minden, maxden, 50)
denlct = MultipleLocator(10)
# Percents:
perrng = Normalize(vmin=0.0, vmax=100.0)
perlev = np.linspace(0.0, 100.0, 50)

# Cycle through all files and plot them.
for infile in files:
    print 'Reading file %s' % infile
    # Read file.
    i = infile.rfind('.')
    if infile[i:] == '.pb':     #PyBats Save file
        mhd = load(infile)
    elif infile[i:] == '.out':  #IDL binary file.
        mhd = pb.Bats2d(infile)
        mhd.regrid(0.25, [-70, 15], [-30,30])

    # Initialize plot.
    fig = plt.figure(figsize=(10.6666666,8)) #3/4 Aspect Ratio
    fig.subplots_adjust(left=0.08, right=0.99, wspace=0.10, hspace=0.18)
    hour = math.floor(mhd.time / 3600.0)
    mint = (mhd.time / 3600.0 - hour)*60.0
    time = '%04ih %05.2fm' % (hour, mint)
    fig.suptitle('Results at iter=%i, t=%s' % (mhd.iter, time))
    fig.text(0.5, 0.01, 'File = '+mhd.filename, ha='center')

    # Pressure plot.
    ax1 = fig.add_subplot(221)
    mhd.data['p_log'] = np.where(mhd.data['p']>minlog, 
                                 mhd.data['p'], 1.01*minlog)
    mhd.planet_add_closed(ax1, DoImf=True, DoOpen=True)
    cont = mhd.contourf(ax1, 'x','z','p_log',levs, norm=LogNorm() )
    cbar = plt.colorbar(cont, pad=0.01, ticks=LogLocator(), 
                        format=LogFormatterMathtext())
    ax1.set_title('Pressure ($nPa$)')

    # Density plot
    ax2 = fig.add_subplot(222)
    cont2 = mhd.contourf(ax2, 'x', 'z', 'rho', denlev, crange=denrng,
                          cmap=plt.get_cmap('RdYlBu_r'), extend='max')
    cbar = plt.colorbar(cont2, pad=0.01, ticks=denlct)
    ax2.set_title('Density ($cm^{-3}$)')

    # Temperature Plot
    ax3 = fig.add_subplot(223)
    mhd.calc_temp(units='kev')
    cont3 = mhd.contourf(ax3, 'x', 'z', 'temp', denlev, crange=denrng,
                         cmap=plt.get_cmap('gist_ncar'), extend='max')
    cbar = plt.colorbar(cont3, pad=0.01, ticks=denlct)
    ax3.set_title('Temperature ($KeV$)')
    
    # Ox Plot
    if mhd.data.has_key('rhoo'):
        mhd.data['comp'] = mhd.data['rhoo']/mhd.data['rho'] * 100.0
        ax4 = fig.add_subplot(224)
        cont4 = mhd.contourf(ax4, 'x', 'z', 'comp', perlev, crange=perrng,
                             cmap=plt.get_cmap('Greens'))
        cbar = plt.colorbar(cont4, pad=0.01, ticks=denlct)
        ax4.set_title('Percent O$^{+}$ by Mass')
        ax_details(ax4)
        mhd.add_body(ax4)

    # Adjust all plots.
    ax_details(ax1)
    ax_details(ax2)
    ax_details(ax3)
    mhd.add_body(ax1)
    mhd.add_body(ax2)
    mhd.add_body(ax3)

    # Save if in batchmode.
    if batchmode:
        # Build output file name by removing file type.
        loc = infile.find('.',-5)
        if loc != -1:
            outname = infile[0:loc] + '.' + format
        else:
            outname = infile + '.' + format
        print 'Saving %s' % outname
        fig.savefig(outname, format=format)
    
    # Close figure to save memory.
    if batchmode:
        plt.close(fig)
    

# Show if not in batch mode.    
if not batchmode:
    plt.show()
