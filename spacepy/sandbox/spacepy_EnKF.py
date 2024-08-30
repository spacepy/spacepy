# test python script for EnKF on spacepy

# load modules
import numpy
import pylab
import spacepy
import datetime
import sys
from spacepy import radbelt
import fp_equality
import spacepy.time as st
import pdb

__contact__ = 'Josef Koller, jkoller@lanl.gov'

# flags for script
# -----------------------------------------

# assimilation flag
#       iassim = 0 : no assimilation
#       iassim = 1 : EnKF assimilation
#       iassim = 2 : data insertion
iassim = 0

# gaussian source term
#       igauss = 0 : no artificial 
#                    gaussian source term
#       igauss = 1 : add artificial 
#                    gaussian source term
igauss = 1
# -----------------------------------------

# radiation belt model
rmod = radbelt.RBmodel()

# Next, set the start time, end time, and the size of the timestep:
start = datetime.datetime(2002,10,23)
end = datetime.datetime(2002,11,4)
#start = datetime.datetime(2002,11,03)
#end = datetime.datetime(2002,11,03,1,0)
delta = datetime.timedelta(hours=0.5)
rmod.setup_ticks(start, end, delta, dtype='UTC')


if (igauss == 1):
    # add Gaussian source term
    rmod.add_source()
else:
    # add PSD data
    rmod.add_PSD_obs()

# Now, run the model over the enitre time range using the evolve method:
rmod.evolve()

if (igauss == 1):
    # observation time delta
    dt  = datetime.timedelta(hours=6.0)
    time = st.tickrange(start, end, dt, dtype='UTC')

    # observation space delta
    Lt = 2

    Tgrid = rmod.ticks
    nTAI = len(Tgrid)

    # compute time delta
    delta = Tgrid[1].UTC[0]-Tgrid[0].UTC[0]
    delta = delta.seconds

    # initialize arrays
    Lstar = ['']*(len(time))
    PSD = ['']*(len(time))

    icount=0
    # loop for data
    for i in numpy.arange(nTAI):

        # get observations every dt
        if (numpy.mod(i*delta,dt.seconds) == 0):
            # determine position of observations
            Lstar[icount] = rmod.Lgrid[0:len(rmod.Lgrid):Lt]
            # determine observations PSD
            PSD[icount] = rmod.PSD[0:len(rmod.Lgrid):Lt,i]

            # add to counter index
            icount = icount + 1

    # add PSD data
    rmod.add_PSD_obs(time=time, Lstar=Lstar, PSD=PSD)

    # observation time delta
    #dt = 6*60*60
    # add PSD data from previous experiment
    #rmod.add_PSD_twin(dt=dt,Lt=2)

# visualize the results of rmod
rmod.plot(values=rmod.PSD,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

# visualize data
rmod.plot_obs(clims=[-10,-6],Lmax=False,Kp=False,Dst=False,title='Observations Plot')


if (iassim > 0):
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

    # determine inflation technique by argument given
    if (len(sys.argv) > 1):
        inflation_opt = int(sys.argv[1])
    else:
        print('default inflation technique')
        inflation_opt = 1
    print('inflation option',inflation_opt)

    if (igauss == 1):
        # neglect Gaussian source term
        rmod.add_source(source=False)

    # assimilate
    if (iassim == 1):
        rmod.assimilate(method='EnKF',inflation=inflation_opt)
    elif (iassim == 2):
        rmod.assimilate(method='insert')

    # visualize the results

    rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

pylab.show()
