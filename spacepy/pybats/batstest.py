#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# Let's test some pybats goodness.
################################
# TEMPORARY USE ONLY.
# CANDIDATE FOR REMOVAL.
# USE AS SANDBOX BUT DON'T SAVE.
################################

Copyright ©2010 Los Alamos National Security, LLC.
"""

import spacepy.pybats.bats as pb
import pylab as plt
import numpy as np
from matplotlib.colors import LogNorm

infile = 'y=0_mhd_1_t00070000_n0872925.out'
mhd1 = pb.Bats2d(infile)

# Regrid.
mhd1.regrid(0.25, [-30,20], [-30,30])

# Plot.
fig1 = plt.figure(figsize=(12,8))
plt.subplots_adjust(left=0.05, right=0.95)
ax1 = fig1.add_subplot(221)
lev_exp = np.arange(np.log10(0.001), np.log10(50.0), 0.10)
levs = np.power(10, lev_exp)
cont = mhd1.contourf(ax1, 'x','z','p',levs, 
                     norm=LogNorm() )
plt.colorbar(cont)

ax2 = fig1.add_subplot(222)
cont2 = mhd1.contourf(ax2, 'x', 'z', 'rho', 35, 
                      cmap=plt.get_cmap('RdYlBu_r'))
plt.colorbar(cont2)

cmap = plt.get_cmap('gist_heat_r')
ax3 = fig1.add_subplot(223)
mhd1.calc_temp(units='kev')
cont3 = mhd1.contourf(ax3, 'x', 'z', 'temp', 35,
                      cmap=plt.get_cmap('gist_ncar'))#'Spectral_r'))
plt.colorbar(cont3)

# Adjust the plots so they look better.
ax1.invert_xaxis()
ax1.autoscale_view(tight=True)
ax2.invert_xaxis()
ax2.autoscale_view(tight=True)
ax3.invert_xaxis()
ax3.autoscale_view(tight=True)
mhd1.add_body(ax1)
mhd1.add_body(ax2)
mhd1.add_body(ax3)

plt.show()
