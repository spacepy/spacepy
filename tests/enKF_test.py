"""
test python script for EnKF on spacepy

Copyright 2010 Los Alamos National Security, LLC.
"""

# load modules
import numpy
import pylab
import spacepy
import datetime
import sys
from spacepy import radbelt


# radiation belt model
rmod = radbelt.RBmodel()

# Next, set the start time, end time, and the size of the timestep:
start = datetime.datetime(2002,10,23)
end = datetime.datetime(2002,11,4)
#start = datetime.datetime(2002,11,03)
#end = datetime.datetime(2002,11,03,1,0)
delta = datetime.timedelta(hours=0.5)
rmod.setup_ticks(start, end, delta, dtype='UTC')
rmod.SRC_model = False

# add PSD data
rmod.add_PSD()

# Now, run the model over the entire time range using the evolve method:
rmod.evolve()

# visualize the results of rmod
rmod.plot(values=rmod.PSD,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

# visualize data
rmod.plot_obs(clims=[-10,-6],Lmax=False,Kp=False,Dst=False,title='Observations Plot')

print('==================================================')
print('                   ASSIMILATING')
print('==================================================')

# ASSIMILATE DATA
# there are three different inflation methodologies within this data
# assimilation scheme:
#
#       inflation == 0: Add model error (perturbation for the ensemble)
#       around model state values only where observations are available.
#
#       inflation == 1: Add model error (perturbation for the ensemble)
#       around observation values only where observations are available.
#
#       inflation == 2: Inflate around ensemble average for EnKF.
rmod.assimilate(inflation=0)

# visualize the results
rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

pylab.show()
