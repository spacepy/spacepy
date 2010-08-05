#!/usr/bin/env python
'''
This script demonstrates how Python and PyBats combine to make a powerful
script that searches an SWMF run directory OR an SWMF results directory
compiled by the PostProc.pl script, finds GM and IE results, and plots
the information on a single plot.
'''
import sys
import glob

from rampy import apply_smart_timeticks
import pybats.bats as pb
import pybats.rim as rim
from pybats import ImfInput

import datetime as dt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import (LogNorm, Normalize)
from matplotlib.ticker import (LogLocator, LogFormatter, LogFormatterMathtext, 
                               MultipleLocator)

# Files to open:
swfile = 'ace20050831.dat'
bfiles = glob.glob('GM/IO2/y*t*.out')
imf = ImfInput(swfile)

# CONSTANTS AND CRAP.
bats_rect = (.1, 1./3.+.07, 2./3.*.91, 2./3.*.8)
minlog = 0.01
maxlog = 50.0
lev_exp = np.arange(np.log10(minlog), np.log10(maxlog), 0.1)
levs = np.power(10, lev_exp)
timerange = [dt.datetime(2005,8,31,9),dt.datetime(2005,9,2,1)]

print "Found %d files to plot." % (len(bfiles))

# MAIN LOOP:
for bfile in bfiles:
    print 'Working on file ', bfile
    # Open Batsfile
    mhd = pb.Bats2d(bfile)
    mhd.regrid(0.25, [-40, 15], [-30,30])

    # What time is it now?
    timenow = timerange[0] + dt.timedelta(seconds=mhd.time)

    # Get corresponding Iono file
    ifile = 'IE/ionosphere/it05%02d%02d_%02d%02d00_000.idl' % (
        timenow.month, timenow.day, timenow.hour, timenow.minute)
    print '\t...with iono file ', ifile
    iono = rim.Iono(ifile)

    # Start by creating a figure:
    fig = plt.figure(figsize=(10.6666666,8))
    fig.subplots_adjust(bottom=0.07, top=0.93, hspace=0.3, left=0.1, right=0.94)
    # BATS plot:
    ax1 = fig.add_axes(bats_rect)
    mhd.data['p_log'] = np.where(mhd.data['p']>minlog, 
                                 mhd.data['p'], 1.01*minlog)
    mhd.planet_add_closed(ax1, DoOpen=False, DoImf=False)
    mhd.add_body(ax1)
    cont = mhd.contourf(ax1, 'x','z','p_log',levs, norm=LogNorm() )
    cbar = plt.colorbar(cont, pad=0.01, ticks=LogLocator(), 
                        format=LogFormatterMathtext())
    ax1.set_title('Pressure ($nPa$)')
    ax1.set_ylim([-20,20])
    ax1.set_xlim([15,-32])
    ax1.set_xlabel('GSM X ($R_{E}$)')
    ax1.set_ylabel('GSM Z ($R_{E}$)')

    # Plot 2 and 3: IE stuff.
    ax2 = fig.add_subplot(333, polar=True)
    cnt = iono.add_cont(ax2,'phi',add_cbar=True)
    ax3 = fig.add_subplot(336, polar=True)
    cnt = iono.add_cont(ax3, 'jr', add_cbar=True)

    # Plot 4: SW Input
    ax4b = fig.add_subplot(313)
    ax4b.plot(imf.time, imf.data['bz'], 'k-')
    ax4b.plot([imf.time[0], imf.time[-1]], [0,0], 'k:')
    ax4b.set_ylabel('IMF $B_{Z} (nT)$')
    ax4b.set_title('Solar Wind Conditions')
    
    # Nudge just a bit.
    pos = ax4b.get_position()
    pos.x1=0.92
    ax4b.set_position(pos)

    ax4v = ax4b.twinx()
    ax4v.plot(imf.time, -1.0*imf.data['vx'], 'b-')
    ax4v.set_ylabel('$V_{x}$ ($km/s$)', color='b')
    labels = ax4v.get_yticklabels()
    for t1 in labels:
        t1.set_color('b')
        
    ax4v.plot([timenow,timenow],[100,600], 'r--', lw=3)
    ax4v.set_ylim([350,550])
    apply_smart_timeticks(ax4v, timerange)
    ax4b.set_xlabel('Universal Time from %s' % timerange[0].isoformat())


    # Save file
    loc = bfile.rfind('.')
    outname = bfile[0:loc] + '_multi.png'
    fig.savefig(outname)
    plt.close(fig)
