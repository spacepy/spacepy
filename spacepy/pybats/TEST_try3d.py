#!/usr/bin/env python

import matplotlib.pyplot as plt
#import matplotlib.tri as tri
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

import pybats.bats as pb

ydat = pb.Bats2d('y=0_mhd_1_t00070000_n0872925.out')
zdat = pb.Bats2d('z=0_mhd_2_t00070000_n0872925.out')

fig = plt.figure()
ax = Axes3D(fig)
ax.tricontour(ydat.grid['x'], ydat.data['rho'], ydat.grid['z'], 
              zdir='y', N=100, offset=-125)
ax.tricontour(zdat.grid['x'], zdat.grid['y'], 
              zdat.data['rho'], N=1, zdir='z',offset=-125)

