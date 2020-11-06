#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions supporting radiation belt diffusion codes

Authors: Josef Koller
Institution: Los Alamos National Laboratory
Contact: jkoller@lanl.gov

Copyright 2010 Los Alamos National Security, LLC.
"""

import ctypes
import datetime

from spacepy import help
import numpy as np
import spacepy.time as st
import pdb

__contact__ = 'Josef Koller, jkoller@lanl.gov'


# -----------------------------------------------
# RBmodel class
# -----------------------------------------------
class RBmodel(object):
    """1-D Radial diffusion class

    This module contains a class for performing and visualizing 1-D radial
    diffusion simulations of the radiation belts.

    Here is an example using the default settings of the model.
    Each instance must be initialized with (assuming import radbelt as rb):

    >>> rmod = rb.RBmodel()

    Next, set the start time, end time, and the size of the timestep:

    >>> import datetime
    >>> start = datetime.datetime(2003,10,14)
    >>> end = datetime.datetime(2003,12,26)
    >>> delta = datetime.timedelta(hours=1)
    >>> rmod.setup_ticks(start, end, delta, dtype='UTC')

    Now, run the model over the entire time range using the evolve method:

    >>> rmod.evolve()

    Finally, visualize the results:

    >>> rmod.plot_summary()

    .. currentmodule:: spacepy.radbelt
    .. autosummary::
        ~RBmodel.Gaussian_source
        ~RBmodel.add_Lmax
        ~RBmodel.add_Lpp
        ~RBmodel.add_PSD_obs
        ~RBmodel.add_PSD_twin
        ~RBmodel.add_omni
        ~RBmodel.add_source
        ~RBmodel.assimilate
        ~RBmodel.evolve
        ~RBmodel.get_DLL
        ~RBmodel.plot
        ~RBmodel.plot_obs
        ~RBmodel.set_lgrid
        ~RBmodel.setup_ticks
    .. automethod:: Gaussian_source
    .. automethod:: add_Lmax
    .. automethod:: add_Lpp
    .. automethod:: add_PSD_obs
    .. automethod:: add_PSD_twin
    .. automethod:: add_omni
    .. automethod:: add_source
    .. automethod:: assimilate
    .. automethod:: evolve
    .. automethod:: get_DLL
    .. automethod:: plot
    .. automethod:: plot_obs
    .. automethod:: set_lgrid
    .. automethod:: setup_ticks
    """

    def __init__(self, grid='L', NL=91, const_kp=False):
        """
        format for grid e.g., L-PA-E
        """
        import spacepy.lib

        self.const_kp=const_kp

        # Initialize the code to use to advance the solutionp.
        if spacepy.lib.have_libspacepy and spacepy.lib.solve_cn:
            self.advance=spacepy.lib.solve_cn
        else:
            self.advance=get_modelop_L

        grid = grid.upper()
        for type in grid.split():
            if type == 'L':
                self.MU = 2083
                self.K = 0.03
                self.DLL_model = 'BA2000'
                self.Lmax_model = 'JKemp'
                self.Lpp_model = 'CA1992'
                #self.SRC_model = 'JK1'
                #self.SRCmagn = datetime.timedelta(days=1e-1) # relative acceleration per day
                self.MPloss = datetime.timedelta(minutes=0.1) # minutes time scale
                self.PPloss = datetime.timedelta(days=10.) # days time scale
                self.set_lgrid(NL)
                self.MODerr = 5. # relative factor * PSD
                #self.PSDerr = 0.3 # relative observational error
                self.PSDerr = 0.25 # relative observational error
                self.MIN_PSD = 1e-99

        # source term flag set to false
        self.source = False

    # -----------------------------------------------
    def __str__(self):
        return '<RB Model; mu=%f, k=%f, DLL_model=%s >' % \
            (self.MU, self.K, self.DLL_model)

    __repr__ = __str__

    # -----------------------------------------------
    def __getitem__(self, idx):
        """
        """

        return self.PSD[:,idx]

    # -----------------------------------------------
    def set_lgrid(self, NL=91):
        '''
        Using NL grid points, create grid in L.
        Default number of points is 91 (dL=0.1).
        '''
        from numpy import linspace, zeros
        self.NL = NL
        self.Lgrid = linspace(1,10,self.NL)
        self.PSDinit = zeros(self.NL, dtype=ctypes.c_double)

    # -----------------------------------------------
    def setup_ticks(self, start, end, delta, dtype='ISO'):
        """
        Add time information to the simulation by specifying
        a start and end time, timestep, and time type (optional).

        Examples
        ========
        >>> start = datetime.datetime(2003,10,14)
        >>> end = datetime.datetime(2003,12,26)
        >>> delta = datetime.timedelta(hours=1)
        >>> rmod.setup_ticks(start, end, delta, dtype='UTC')
        """
        self.ticks = st.tickrange(start, end, delta, dtype)

    # -----------------------------------------------
    def add_omni(self, keylist=None):

        """
        add omni data to instance according to the tickrange in ticks
        """
        import spacepy.omni as om

        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        omni = om.get_omni(self.ticks)

        if 'params' not in self.__dict__:
            self.params = {}

        if keylist is None:
            # add all keys
            self.params = omni
        else:
            # just add requested keys
            for key in keylist:
                self.params[key] = omni[key]

    # -----------------------------------------------
    def add_Lmax(self, Lmax_model):

        """
        add last closed drift shell Lmax
        """
        import spacepy.empiricals as em

        if 'params' not in self.__dict__:
            self.params = {}

        assert self.ticks, "Provide tick range with 'setup_ticks'"
        self.params['Lmax'] = em.get_Lmax(self.ticks, Lmax_model)

    # -----------------------------------------------
    def add_Lpp(self, Lpp_model):

        """
        add last closed drift shell Lmax
        """
        import spacepy.empiricals as em

        if 'params' not in self.__dict__:
            self.params = {}

        assert self.ticks, "Provide tick range with 'setup_ticks'"
        self.params['Lpp'] = em.get_plasma_pause(self.ticks, Lpp_model)

    # -----------------------------------------------
    def add_PSD_obs(self, time=None, PSD=None, Lstar=None, satlist=None):

        """
        add PSD observations

        Parameters
        ----------

        time : Ticktock datetime array
            array of observation times

        PSD : list of numpy arrays
            PSD observational data for each time. Each entry in the list is a
            numpy array with the observations for the corresponding time

        Lstar : list of numpy arrays
            Lstar location of each PSD observations. Each entry in the list is
            a numpy array with the location of the observations for the
            corresponding time

        satlist : list of satellite names

        Returns
        -------
        out : list of dicts
            Information of the observational data, where each entry
            contains the observations and locations of observations for each
            time specified in the time array. Each list entry is a dictionary
            with the following information:

        Ticks : Ticktock array
            time of observations
        Lstar : numpy array
            location of observations
        PSD   : numpy array
            PSD observation values
        sat   : list of strings
            satellite names
        MU    : scalar value
            Mu value for the observations
        K     : scalar value
            K value for the observations
        """

        import pdb

        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        Tgrid = self.ticks
        nTAI = len(Tgrid)

        # initialize PSDdata list
        self.PSDdata = ['']*(nTAI-1)

        if (PSD == None):
        # PSD data not provided,
        # extract from database
            import spacepy.sandbox.PSDdata as PD

            for i, Tnow, Tfut in zip(np.arange(nTAI-1), Tgrid[:-1], Tgrid[1:]):
                start_end = spacepy.time.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
                self.PSDdata[i] = PD.get_PSD(start_end, self.MU, self.K, satlist)

        else:
        # PSD data arrays provided

            # model grid
            Lgrid = self.Lgrid

            itime = 0
            # loop over time defined for the model integration
            for i, Tnow, Tfut in zip(np.arange(nTAI-1), Tgrid[:-1], Tgrid[1:]):

                # empty array
                lstar = np.array([],dtype=float)
                time_idx = np.array([],dtype=int)

                # loop over observation time
                for itime in np.arange(len(time)):

                    if (Tnow <= time[itime] and time[itime] <= Tfut):
                        #print 'match!! '
                        #print itime
                        #print i
                        #print time[itime]
                        #print Tnow
                        #print Tfut

                        # concatenate to lstar
                        lstar = np.concatenate((lstar,Lstar[itime]))
                        lstar = np.unique(lstar)
                        #idx = lstar.argsort()

                        #pdb.set_trace()
                        # add time index
                        time_idx = np.append(time_idx,itime)
                #end loop

                #pdb.set_trace()
                if (time_idx.shape[0] > 0):
                    # initialize PSD array
                    psd = np.zeros_like(lstar)

                    # initialize number of obs array
                    num_obs = np.zeros_like(lstar)

                    # sort time index
                    time_idx = np.unique(time_idx)

                    # loop over time index
                    for itime in time_idx:
                        # sort observations
                        idx = Lstar[itime].argsort()
                        tmplstar = Lstar[itime][idx]
                        tmppsd = PSD[itime][idx]

                        # run through all unique grid-points and compute
                        # average observation
                        for j, iL in enumerate(lstar):
                            # identify idex of grid-point
                            idx = np.where(np.isclose(iL,tmplstar))
                            # assign observation for grid-point
                            psd[j] = psd[j] + tmppsd[idx]
                            # add for number of observations
                            num_obs[j] = num_obs[j] + 1.0

                    psd = psd/num_obs

                    # assign time for observations
                    #Ticks = time[itime]
                    Ticks = Tfut
                    # determine position of observations
                    #lstar = Lstar[itime]
                    # determine observations PSD
                    #psd = PSD[itime]
                    # provide MU
                    MU = self.MU*np.ones_like(lstar)
                    # provide K
                    K = self.K*np.ones_like(lstar)
                    # empy satellite
                    sat = ['']

                    # add to dictionary
                    self.PSDdata[i] = {'Ticks':Ticks, 'Lstar':lstar, \
                                       'PSD':psd, 'sat':sat, \
                                       'MU':MU, 'K':K}

        # adjust initial conditions to these PSD values
        mval = np.mean(self.PSDdata[0]['PSD'])
        self.PSDinit = mval*np.exp(-(self.Lgrid - 5.5)**2/0.8)

        return
    # -----------------------------------------------
###    def add_PSD(self, satlist=None):
###
###        """
###        add observations from PSD database using the ticks list
###        """
###
###        import spacepy.sandbox.PSDdata as PD
###        import spacepy.time
###
###        assert 'ticks' in self.__dict__ , \
###            "Provide tick range with 'setup_ticks'"
###        Tgrid = self.ticks
###        nTAI = len(Tgrid)
###
###        self.PSDdata = ['']*(nTAI-1)
###        for i, Tnow, Tfut in zip(np.arange(nTAI-1), Tgrid[:-1], Tgrid[1:]):
###            start_end = spacepy.time.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
###            self.PSDdata[i] = PD.get_PSD(start_end, self.MU, self.K, satlist)
###
###        # adjust initial conditions to these PSD values
###        mval = np.mean(self.PSDdata[0]['PSD'])
###        #self.PSDinit = mval*np.exp(-(self.Lgrid - 5.5)**2/0.2)
###        self.PSDinit = mval*np.exp(-(self.Lgrid - 5.5)**2/0.8)
###
###        return
###    # -----------------------------------------------
    def add_PSD_twin(self,dt=0,Lt=1):

        """
        add observations from PSD database using the ticks list
        the arguments are the following:

            dt = observation time delta in seconds
            Lt = observation space delta
        """
        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        Tgrid = self.ticks
        nTAI = len(Tgrid)

        # compute time delta
        delta = Tgrid[1].UTC[0]-Tgrid[0].UTC[0]
        delta = delta.seconds

        # compute observations time delta
        #dt = 2*60*60

        self.PSDdata = ['']*(nTAI-1)
        for i, Tnow, Tfut in zip(np.arange(nTAI-1), Tgrid[:-1], Tgrid[1:]):
            # get observations every dt
            if (np.mod(i*delta,dt) == 0):
                # assign time for observations
                Ticks = Tnow
                # determine position of observations
                lstar = self.Lgrid[0:len(self.Lgrid):Lt]
                # determine observations PSD
                psd = self.PSD[0:len(self.Lgrid):Lt,i]
                # provide MU
                MU = self.MU*np.ones_like(lstar)
                # provide K
                K = self.K*np.ones_like(lstar)
                # empy satellite
                sat = np.zeros_like(lstar)
                # add to dictionary
                self.PSDdata[i] = {'Ticks':Ticks, 'Lstar':lstar, \
                                   'PSD':psd, 'sat':sat, \
                                   'MU':MU, 'K':K}
            else:
                self.PSDdata[i] = {}

        # adjust initial conditions to these PSD values
        mval = np.mean(self.PSDdata[0]['PSD'])
        self.PSDinit = mval*np.exp(-(self.Lgrid - 5.5)**2/0.2)

        return
    # -----------------------------------------------
    def add_source(self,source=True,A=1.0e-8,mu=5.0,sigma=0.5):

        """
        add source parameters A, mu, and sigma for the Gaussian source function
        """
        # set source term flag
        self.source = source

        # define height of Gaussian function
        self.source_A = A

        # define center of Gaussian function
        self.source_mu = mu

        # define amplitude
        self.source_sigma = sigma

        return
    # -----------------------------------------------
    def Gaussian_source(self):
        """
        Gaussian source term added to radiation belt model. The source term is
        given by the equation:

        S = A exp{-(L-mu)^2/(2*sigma^2)}

        with A=10^(-8), mu=5.0, and sigma=0.5 as default values
        """
        Lgrid = self.Lgrid

        # determine whether to include source or not
        if self.source:
            # define height of Gaussian function
            A = self.source_A

            # define center of Gaussian function
            mu = self.source_mu

            # define amplitude
            sigma = self.source_sigma
        else:
            A = 0.0
            mu = 1.0
            sigma = 1.0

        # Gaussian function
        S = A*np.exp(-(Lgrid-mu)**2/(2.0*sigma**2))

        return S

    # -----------------------------------------------
    def evolve(self):
        """
        calculate the diffusion in L at constant mu,K coordinates
        """
        from . import radbelt as rb

        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"

        f = self.PSDinit
        Tgrid = self.ticks.TAI
        nTAI = len(Tgrid)
        Lgrid = self.Lgrid
        self.PSD  = np.zeros( (len(f),len(Tgrid)), dtype=ctypes.c_double)
        self.PSD[:,0] = f.copy()

        # add omni if not already given
        if 'omni' not in self.__dict__:
            self.add_omni(keylist=['Kp', 'Dst'])

        # run Lmax model
        self.add_Lmax(self.Lmax_model)
        # run plasma pause model
        self.add_Lpp(self.Lpp_model)

        # setup params dictionary
        params = {}
        params['DLL_model'] = self.DLL_model
        #params['SRC_model'] = self.SRC_model
        #params['SRCmagn'] = self.SRCmagn
        params['MPloss'] = self.MPloss
        params['PPloss'] = self.PPloss
        if 'SRCartif' in self.__dict__:
            params['SRCartif'] = self.SRCartif
        keylist = ['Kp', 'Dst', 'Lmax', 'Lpp']

        # start with the first one since 0 is the initial condition
        for i, Tnow, Tfut in zip(np.arange(nTAI-1)+1, Tgrid[:-1], Tgrid[1:]):
            Tdelta = Tfut - Tnow
            Telapse= Tnow - Tgrid[0]
            # copy over parameter list
            for key in keylist:
                params[key] = self.params[key][i]
            # now integrate from Tnow to Tfut
            f = diff_LL(self, Lgrid, f, Tdelta, Telapse, params=params)

            # add Gaussian source term
            f = f + self.Gaussian_source()

            # copy to PSD
            self.PSD[:,i] = f.copy()

    # -----------------------------------------------
    def assimilate(self, method='EnKF',inflation=0):
        """
        Assimilates data for the radiation belt model using the Ensemble
        Kalman Filter. The algorithm used is the SVD method presented by
        Evensen in 2003 (Evensen, G., Ocean dynamics, 53, pp.343--367, 2003).
        To compensate for model errors, three inflation algorithms are
        implemented. The inflation methodology is specified by the
        'inflation' argument, and the options are the following:

            inflation == 0: Add model error (perturbation for the ensemble)
            around model state values only where observations are available
            (DEFAULT).

            inflation == 1: Add model error (perturbation for the ensemble)
            around observation values only where observations are available.

            inflation == 2: Inflate around ensemble average for EnKF.

        Prior to assimilation, a set of data values has to be speficied by
        setting the start and end dates, and time step, using the setup_ticks
        funcion of the radiation belt model:

        >>> import spacepy
        >>> import datetime
        >>> from spacepy import radbelt

        >>> start = datetime.datetime(2002,10,23)
        >>> end = datetime.datetime(2002,11,4)
        >>> delta = datetime.timedelta(hours=0.5)
        >>> rmod.setup_ticks(start, end, delta, dtype='UTC')

        Once the dates and time step are specified, the data is added using the
        add_PSD function:

        >>> rmod.add_PSD()

        The observations are averaged over the time windows, whose interval is
        give by the time step.

        Once the dates and data are set, the assimiation is performed using the
        'assimilate' function:

        >>> rmod.assimilate(inflation=1)

        This function will add the PSDa values, which are the analysis state of
        the radiation belt using the observations within the dates. To plot the
        analysis simply use the plot funtion:

        >>> rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

        """
        import spacepy.data_assimilation
        import spacepy.sandbox.PSDdata as PD
        import copy as c

        # add PSD observations with add_PSD,
        # this has to be done to the class
        # module when running the RBmodel.
        # add_PSD will add the PSDdata to the class
        # which is a dictionary for the data that has been added

        # debugging command
        #pdb.set_trace()

        # setup method
        assert method in ['EnKF','insert'], 'data assimilation method='+method+' not implemented'

        nTAI = len(self.ticks)

        # enKF method
        if method == 'EnKF':
            da = spacepy.data_assimilation.ensemble()

            # initialize A with initial condition
            # the initial condition is ones, have to change this to initialize
            # the ensemble with perturbed reference state
            A = np.ones( (self.NL, da.Nens) )*self.PSDinit[:,np.newaxis]

            self.PSDf = np.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa = np.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa[:,0] = self.PSDinit
            self.PSDf[:,0] = self.PSDinit

            # diagnostic tools:
            # observations-minus-background
            self.PSD_omb = ['']*(nTAI)
            # observations-minus-analysis
            self.PSD_oma = ['']*(nTAI)
            # analysis-minus-background
            self.PSD_amb = np.zeros( (self.PSDinit.shape[0], nTAI) )

            # add model error (perturbation) in the ensemble initial conditionp.
            #A = da.add_model_error(self, A, self.PSDdata[0])
            std = 0.35
            normal = np.random.randn( da.Nens )
            for iens in np.arange(da.Nens):
                A[:,iens] = A[:,iens] + std*normal[iens]*A[:,iens]
            A[np.where(A < self.MIN_PSD)] = self.MIN_PSD

            # ==========================================
            # DEBUG
            #np.savetxt('ensemble_IC.dat',A)
            #np.savetxt('model_IC.dat',self.PSDinit)
            #np.savetxt('model_grid_IC.dat',self.Lgrid)
            # ==========================================

            # create temporary RB class instance
            rbtemp = c.copy(self)

            # time loop
            for i, Tnow, Tfut in zip(np.arange(nTAI-1)+1, self.ticks[:-1], self.ticks[1:]):

                # make forcast and add model error
                # make forecast using all ensembles in A
                iens = 0
                for f in A.T:
                    rbtemp.ticks = st.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
                    rbtemp.PSDinit = f.copy()
                    rbtemp.evolve()
                    #rbtemp.PSDdata = self.PSDdata[i-1]
                    A[:,iens] = rbtemp.PSD[:,1].copy()
                    iens += 1

                # save result in ff
                Tnow = Tfut
                self.PSDf[:,i] = np.mean(A, axis=1)

                # verify that there are data points within the interval, if
                # there are data points then extract average observations
                # within the window, if not return empty observation array y
                # and Lobs.

                if len(self.PSDdata[i-1]) > 0:
                    # get observations for time window ]Tnow-Twindow,Tnow]
                    Lobs, y = spacepy.data_assimilation.average_window(self.PSDdata[i-1], self.Lgrid)
                else:
                    y = np.array([])
                    Lobs = np.array([])

                print(Lobs)
                print(y)

                # ==========================================
                # DEBUG
                #np.savetxt('obs_location_IC.dat',Lobs)
                #np.savetxt('obs_IC.dat',y)
                #pdb.set_trace()
                # ==========================================

                # then assimilate otherwise do another forcast
                if len(y) > 0:

                    ### check for minimum PSD values
                    ### A[np.where(A<self.MIN_PSD)] = self.MIN_PSD
                    ### # insert observations directly
                    ### A = da.add_model_error_obs(self, A, Lobs, y)
                    ### # dictionary
                    ### self.PSDa[:,i] = np.mean(A, axis=1)

                    # INFLATION SCHEMES
                    if inflation == 0:
                        print('inflation around model state values at observation locations')
                        # Add model error (perturbation for the ensemble) around model
                        # state values.  This acts as an inflation scheme for EnKF
                        A = da.add_model_error(self, A, self.PSDdata[i-1])
                    elif inflation == 1:
                        print('inflation around observation values at observation locations')
                        # Add model error (perturbation for the ensemble) around
                        # observation values. This acts as an inflation scheme for EnKF
                        A = da.add_model_error_obs(self, A, Lobs, y)
                    elif inflation == 2:
                        print('inflation around ensemble average')
                        # Inflate around ensemble average for EnKF
                        # ensemble average
                        ens_avg = np.mean(A, axis=1)

                        # inflation factor
                        inflation_factor = 1.8

                        # loop over ensemble members and inflate
                        iens = 0
                        for ens in A.T:
                            # inflate ensemble
                            A[:,iens] = inflation_factor*(ens - ens_avg) + ens_avg
                            iens += 1
                    A[np.where(A < self.MIN_PSD)] = self.MIN_PSD

                    # prepare assimilation analysis
                    # project ensemble states to obs. grid
                    HA = da.getHA(self, Lobs, A)

                    # measurement perturbations ensemble
                    Psi = da.getperturb(self, y)

                    # ensemble of innovation vectors
                    Inn = da.getInnovation(y, Psi, HA)

                    # calculate ensemble perturbation HA' = HA-HA_mean
                    HAp = da.getHAprime(HA)

                    # calculate prior diagnostics
                    # observation minus background
                    omb = y-np.average(HA,axis=1)
                    self.PSD_omb[i] = {
                            'Lobs':Lobs,
                            'y':y,
                            'omb':omb,
                            }
                    xf_avg = np.average(A,axis=1)

                    # now call the main analysis routine
                    if len(y) == 1:
                        A = da.EnKF_oneobs(A, Psi, Inn, HAp)
                    else:
                        A = da.EnKF(A, Psi, Inn, HAp)

                    # check for minimum PSD values
                    A[np.where(A<self.MIN_PSD)] = self.MIN_PSD

                    # average A from analysis step and save in results
                    # dictionary
                    self.PSDa[:,i] = np.mean(A, axis=1)

                    # calculate posterior diagnostics
                    # observation minus analysis
                    Hanalysis = da.getHA(self, Lobs, A)
                    oma = y-np.average(Hanalysis,axis=1)
                    self.PSD_oma[i] = {
                            'Lobs':Lobs,
                            'y':y,
                            'omb':oma,
                            }
                    # analysis minus background
                    xa_avg = np.average(A,axis=1)
                    self.PSD_amb[:,i] = xa_avg - xf_avg
                    #pdb.set_trace()

                    # print assimilated result
                    Hx = np.zeros_like(y)
                    for iL,Lstar in enumerate(Lobs):
                        idx = np.where(np.isclose(self.Lgrid,Lstar))
                        Hx[iL] = self.PSDa[idx,i]
                    print(Hx)
                    print(HA)

                elif len(y) == 0:
                    print('no observations within this window')
                    self.PSDa[:,i] = self.PSDf[:,i]
                    continue #

                # print message
                print('Tnow: ', self.ticks[i].ISO)
        # insert obsrvations for data assimilation
        elif method == 'insert':
            da = spacepy.data_assimilation.ensemble()

            A = np.ones( (self.NL, 1) )*self.PSDinit[:,np.newaxis]

            self.PSDf = np.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa = np.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa[:,0] = self.PSDinit
            self.PSDf[:,0] = self.PSDinit

            # add model error (perturbation) in the ensemble initial condition.
            std = 0.15
            normal = np.random.randn( self.NL )
            A[:,0] = (1.0 + std*normal)*A[:,0]
            A[np.where(A < self.MIN_PSD)] = self.MIN_PSD

            rbtemp = c.copy(self)

            # time loop
            for i, Tnow, Tfut in zip(np.arange(nTAI-1)+1, self.ticks[:-1], self.ticks[1:]):
                # evolve solution
                rbtemp.ticks = st.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
                rbtemp.PSDinit = A[:,0]
                rbtemp.evolve()
                A[:,0] = rbtemp.PSD[:,1]

                # save result in ff
                Tnow = Tfut
                self.PSDf[:,i] = A[:,0]

                # verify that there are data points within the interval, if
                # there are data points then extract average observations
                # within the window, if not return empty observation array y
                # and Lobs.
                if len(self.PSDdata[i-1]) > 0:
                    Lobs, y = spacepy.data_assimilation.average_window(self.PSDdata[i-1], self.Lgrid)
                else:
                    y = np.array([])
                    Lobs = np.array([])

                print(Lobs)
                print(y)

                #pdb.set_trace()
                # then assimilate otherwise do another forcast
                if len(y) > 0:
                    # check for minimum PSD values
                    A[np.where(A<self.MIN_PSD)] = self.MIN_PSD

                    # create index of obs location
                    idx = np.array([],dtype=int)
                    for Lval in Lobs:
                        Lidx = np.where( fp_equality.eq(Lval,self.Lgrid) )[0]
                        idx = np.append(idx,np.array([int(Lidx)]))

                    # insert observations directly
                    A[idx,0] = y

                    #A = da.add_model_error_obs(self, A, Lobs, y)
                    # dictionary
                    self.PSDa[:,i] = A[:,0]

                elif len(y) == 0:
                    print('no observations within this window')
                    self.PSDa[:,i] = self.PSDf[:,i]
                    continue #

                # print message
                print('Tnow: ', self.ticks[i].ISO)

    # -----------------------------------------------
    def plot(self, Lmax=True, Lpp=False, Kp=True, Dst=True,
             clims=[0,10], title=None, values=None):
        """
        Create a summary plot of the RadBelt object distribution function.
        For reference, the last closed drift shell, Dst, and Kp are all
        included.  These can be disabled individually using the corresponding
        Boolean kwargs.

        The clims kwarg can be used to manually set the color bar range.
        To use, set it equal to a two-element list containing minimum and
        maximum Log_10 value to plot.  Default action is to use [0,10] as
        the log_10 of the color range.  This is good enough for most
        applications.

        The title of the top most plot defaults to 'Summary Plot' but can be
        customized using the title kwarg.

        The figure object and all three axis objects (PSD axis, Dst axis,
        and Kp axis) are all returned to allow the user to further customize
        the plots as necessary.  If any of the plots are excluded, None is
        returned in their stead.

        Examples
        ========

        >>> rb.plot(Lmax=False, Kp=False, clims=[2,10], title='Good work!')

        This command would create the summary plot with a color bar range
        of 100 to 10^10.  The Lmax line and Kp values would be excluded.
        The title of the topmost plot (phase space density) would be set to
        'Good work!'.
        """
        import matplotlib.pyplot as p
        from matplotlib.colors  import LogNorm
        from matplotlib.ticker  import (LogLocator, LogFormatter,
                                        LogFormatterMathtext)
        import pdb

        # debugging command
        #pdb.set_trace()

        # test for default values
        if values is None:
            values = self.PSD

        # Initialize axis variables so that they can be returned,
        # even if not used.
        ax1 = None
        ax2 = None
        ax3 = None

        fig = p.figure()
        fig.subplots_adjust(left=0.10, right=0.999, top=0.92)

        # PLOT PSD
        if Kp or Dst:
            ax1 = p.subplot(2,1,1)
        else:
            ax1 = p.subplot(1,1,1)
        # Plot phase space density, masking out values of 0.
        map = ax1.pcolorfast(self.ticks.eDOY, self.Lgrid,
                             np.where(values > 0.0, self.PSD, 10.0**-39),
                             vmin=10.0**clims[0], vmax=10.0**clims[1],
                             norm=LogNorm())
        ax1.set_ylabel('L*')

        if title is not None:
            ax1.set_title(title)
        # Add color bar.
        cbar = p.colorbar(map, pad=0.01, shrink=.85, ticks=LogLocator(),
                          format=LogFormatterMathtext())
        cbar.set_label('Phase Space Density')
        # add Lmax line
        if Lmax:
            p.plot(self.ticks.eDOY, self.params['Lmax'], 'w')
        # Minimize time range.
        ax1.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])
        # Finally, save the position of the plot to make next plot match.
        pos = ax1.get_position()

        if Dst is True:
            ax2 = p.subplot(2,1,2)
            pos2 = ax2.get_position()
            pos2.x1 = pos.x1
            ax2.set_position(pos2)
            ax2.plot(self.ticks.eDOY, self.params['Dst'], color='r')
            ax2.set_xlabel('DOY in '+str(self.ticks.UTC[0].year))
            ax2.set_ylabel('Dst', color='r')
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
            ax2.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])

        if Kp is True:
            if Dst is True:
                p.subplot(2,1,2)
                ax3 = p.twinx()
            else:
                ax3 = p.subplot(2,1,2)
                pos3 = ax3.get_position()
                pos3.x1 = pos.x1
                ax3.set_position(pos3)
            ax3.plot(self.ticks.eDOY, self.params['Kp'], 'k:')
            if Dst is True:
                ax3.yaxis.tick_right()
            ax3.set_xlabel('DOY in '+str(self.ticks.UTC[0].year))
            ax3.set_ylabel('Kp')
            ax3.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])

        #p.show()

        return fig, ax1, ax2, ax3

    # -----------------------------------------------
    def plot_obs(self, Lmax=True, Lpp=False, Kp=True, Dst=True,
             clims=[0,10], title=None, values=None):
        """
        Create a summary plot of the observations.  For reference, the last
        closed drift shell, Dst, and Kp are all included.  These can be
        disabled individually using the corresponding boolean kwargs.

        The clims kwarg can be used to manually set the color bar range.
        To use, set it equal to a two-element list containing minimum and
        maximum Log_10 value to plot.  Default action is to use [0,10] as
        the log_10 of the color range.  This is good enough for most
        applications.

        The title of the top most plot defaults to 'Summary Plot' but can be
        customized using the title kwarg.

        The figure object and all three axis objects (PSD axis, Dst axis,
        and Kp axis) are all returned to allow the user to further customize
        the plots as necessary.  If any of the plots are excluded, None is
        returned in their stead.

        Examples
        ========

        >>> rb.plot_obs(Lmax=False, Kp=False, clims=[2,10], title='Observations Plot')

        This command would create the summary plot with a color bar range
        of 100 to 10^10.  The Lmax line and Kp values would be excluded.
        The title of the topmost plot (phase space density) would be set to
        'Good work!'.
        """
        import spacepy.data_assimilation
        import matplotlib.pyplot as p
        from matplotlib.colors  import LogNorm
        from matplotlib.ticker  import (LogLocator, LogFormatter,
                                        LogFormatterMathtext)
        import pdb

        # debugging command
        #pdb.set_trace()

        # test for default values
        if values is None:
            values = self.PSDdata

        # Initialize axis variables so that they can be returned,
        # even if not used.
        ax1 = None
        ax2 = None
        ax3 = None

        fig = p.figure()
        fig.subplots_adjust(left=0.10, right=0.999, top=0.92)

        # compute time-window average observation value
        y = np.array([],dtype=float)
        Lobs = np.array([],dtype=float)
        eDOYobs = np.array([],dtype=float)
        nTAI = len(self.ticks)
        # time loop
        for i, Tnow, Tfut in zip(np.arange(nTAI-1)+1, self.ticks[:-1], self.ticks[1:]):
            if len(values[i-1]) > 0:
                # get observations for time window ]Tnow-Twindow,Tnow]
                Lobs_tmp, y_tmp = spacepy.data_assimilation.average_window(values[i-1], self.Lgrid)
                y = np.append(y,y_tmp)
                Lobs = np.append(Lobs,Lobs_tmp)
                eDOYobs = np.append(eDOYobs,np.ones(len(y_tmp))*self.ticks.eDOY[i-1])

        # PLOT PSDdata
        if Kp or Dst:
            ax1 = p.subplot(2,1,1)
        else:
            ax1 = p.subplot(1,1,1)
        # Plot phase space density observations.
        map = ax1.scatter(eDOYobs,Lobs,c=y,norm=LogNorm(),
                          vmin=10.0**clims[0], vmax=10.0**clims[1],
                          edgecolor='none')
        ax1.set_ylabel('L*')

        if title is not None:
            ax1.set_title(title)

        ax1.set_ylim(self.Lgrid[0],self.Lgrid[len(self.Lgrid)-1])
        # Add color bar.
        cbar = p.colorbar(map, pad=0.01, shrink=.85, ticks=LogLocator(),
                          format=LogFormatterMathtext())
        cbar.set_label('Phase Space Density')
        # add Lmax line
        if Lmax:
            p.plot(self.ticks.eDOY, self.params['Lmax'], 'w')
        # Minimize time range.
        ax1.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])
        # Finally, save the position of the plot to make next plot match.
        pos = ax1.get_position()

        if Dst is True:
            ax2 = p.subplot(2,1,2)
            pos2 = ax2.get_position()
            pos2.x1 = pos.x1
            ax2.set_position(pos2)
            ax2.plot(self.ticks.eDOY, self.params['Dst'], color='r')
            ax2.set_xlabel('DOY in '+str(self.ticks.UTC[0].year))
            ax2.set_ylabel('Dst', color='r')
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
            ax2.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])

        if Kp is True:
            if Dst is True:
                p.subplot(2,1,2)
                ax3 = p.twinx()
            else:
                ax3 = p.subplot(2,1,2)
                pos3 = ax3.get_position()
                pos3.x1 = pos.x1
                ax3.set_position(pos3)
            ax3.plot(self.ticks.eDOY, self.params['Kp'], 'k:')
            if Dst is True:
                ax3.yaxis.tick_right()
            ax3.set_xlabel('DOY in '+str(self.ticks.UTC[0].year))
            ax3.set_ylabel('Kp')
            ax3.set_xlim([self.ticks.eDOY[0], self.ticks.eDOY[-1]])

        #p.show()

        return fig, ax1, ax2, ax3


    def get_DLL(self, Lgrid, params, DLL_model='BA2000'):
        """
        Calculate DLL as a simple power law function (alpha*L**Bbta)
        using alpha/beta values from popular models found in the
        literature and chosen with the kwarg "DLL_model".

        The calculated DLL is returned, as is the alpha and beta
        values used in the calculationp.

        The output DLL is in units of units/day.
        """


        if DLL_model == 'BA2000': # Brautigam and Albert (2000)
            if type(self.const_kp) == type(0.0):
                Kp=self.const_kp
            else:
                Kp = params['Kp']
            alpha = 10.0**(0.506*Kp-9.325)
            beta = 10.0

        elif DLL_model == 'FC2006': # Fei and Chan (2006)
            alpha = 1.5e-6
            beta  = 8.5

        elif DLL_model == 'U2008': # Ukhorskiy (2008)
            alpha = 7.7e-6
            beta  = 6.0

        elif DLL_model == 'S1997': # Selesnick (1997)
            alpha = 1.9e-10
            beta  = 11.7
        elif DLL_model == 'const': # Constant DLL.
            alpha= 1.0
            beta = 1.0
            DLL  = np.zeros(len(Lgrid), dtype=ctypes.c_double)+10.
            # approximately BA2000 for Kp=1, L=1.
        else:
            raise ValueError("Radial diffusion model %s not implemented" % DLL_model)


        if (DLL_model!='const'): DLL = alpha * Lgrid ** beta

        return DLL, alpha, beta
# -----------------------------------------------
def get_modelop_L(f, L, Dm_old, Dm_new, Dp_old, Dp_new, Tdelta, NL):

    """
    Advance the distribution function, f, discretized into the Lgrid, L, forward
    in time by a timestep, Tdelta.  The off-grid current and next diffusion
    coefficients, D[m,p]_[old,new] will be used.  The number of grid points is set
    by NL.

    This function performs the same calculation as the C-based code,
    spacepy.lib.solve_cnp.  This code is very slow and should only be used when
    the C code fails to compile.
    """

    import numpy.linalg as nlin

    # get grid and setup centered grid Lm=L_i-1/2, Lp=L_i+1/2
    Lmin = L[0]
    Lmax = L[-1]
    dL = L[1] - L[0]
    Lp = L + 0.5*dL
    Lm = L - 0.5*dL

    # setup the diffusion coefficient on each grid and center grid
    #Dllm = np.zeros(NL); Dllp = np.zeros(NL)
    betam = np.zeros(NL); betap = np.zeros(NL)
    for i in range(1,int(NL)-1,1):
        Dllp[i] = 0.5*(DLL[i]+DLL[i+1])
        Dllm[i] = 0.5*(DLL[i]+DLL[i-1])
        betam[i] = Dllm[i]*dt / (2*dL*dL*Lm[i]*Lm[i])
        betap[i] = Dllp[i]*dt / (2*dL*dL*Lp[i]*Lp[i])

    # initialize some arrays
    A = np.zeros( (NL,NL) )
    B = np.zeros( (NL,NL) )
    C = np.zeros( (NL,NL) )

    # off diagonal elements
    for i in np.arange(1,int(NL)-1,1):
        # diagonal elements
        A[i,i] = 1 + betam[i]*L[i]*L[i] + betap[i]*L[i]*L[i]
        B[i,i] = 1 - betam[i]*L[i]*L[i] - betap[i]*L[i]*L[i]
        # off diagonal elements
        A[i,i+1] = -betap[i]*L[i]*L[i]
        A[i,i-1] = -betam[i]*L[i]*L[i]
        B[i,i+1] =  betap[i]*L[i]*L[i]
        B[i,i-1] =  betam[i]*L[i]*L[i]


    # boundary condition
    A[0,0] = 1.
    A[-1,-1] = 1.
    B[0,0] = 1.
    B[-1,-1] = 1.

    # get inverse and multipy
    Ai = nlinp.inv(A)
    C = np.dot(Ai,B)

    return np.dot(C, f)

# -----------------------------------------------
def diff_LL(r, grid, f, Tdelta, Telapsed, params=None):
    """
    calculate the diffusion in L at constant mu,K coordinates
    time units
    """

    import ctypes as ct

    # prepare some parameter variables
    if params:
        DLL_model = params['DLL_model']
        Kp = params['Kp']
        Dst = params['Dst']
        Lmax = params['Lmax']
        Lpp = params['Lpp']
    #else: # use standard config
        # but which one?

    Lgrid = grid
    NL = len(Lgrid)
    dL = Lgrid[1]-Lgrid[0]

    # get DLL and model operator
    (DLL,  alpha, beta) = r.get_DLL(Lgrid,          params, DLL_model)
    (DLLp, alpha, beta) = r.get_DLL(Lgrid + dL/2.0, params, DLL_model)
    (DLLm, alpha, beta) = r.get_DLL(Lgrid - dL/2.0, params, DLL_model)
    DLL = DLL/86400.; DLLp = DLLp/86400.; DLLm = DLLm/86400.

    # Set default of NO sources:
    src = np.zeros(NL)

    # Create source operators (source splitting) using implicit
    # trapezoidal method to solve source ODE.
    # Add artificial sources
    if 'SRCartif' in params:
        # Call the artificial source function, sending info as
        # key word arguments.  Note that such functions should be
        # able to handle extra kwargs through the use of **kwargs!
        sfunc = params['SRCartif']

        # Apply using correct CN-method.
        src=0.5*Tdelta*(
            sfunc(Telapsed, Lgrid, alpha, DLL, beta) +
            sfunc(Telapsed+Tdelta, Lgrid, alpha, DLL, beta))
        src[0], src[-1] = 0.,0.

    else:
        src1 = np.zeros(NL)
        src2 = np.zeros(NL)

    # Apply solution operators to f.
    dptr = ctypes.POINTER(ctypes.c_double)
    r.advance(f.ctypes.data_as(dptr),
              Lgrid.ctypes.data_as(dptr),
              DLLm.ctypes.data_as(dptr),
              DLLm.ctypes.data_as(dptr),
              DLLp.ctypes.data_as(dptr),
              DLLp.ctypes.data_as(dptr),
              ct.c_double(Tdelta), NL,
              src.ctypes.data_as(dptr))

   # add source according to values in SRC...
#    if params['SRC_model']:
#        # setup source vector
#        S = get_local_accel(Lgrid, params, SRC_model='JK1')
#        f = f + S*Tdelta

    # add losses through magnetopause shadowing, time scale taken from MPloss
    mp_sec = params['MPloss'].days * 86400.0 + params['MPloss'].seconds + \
             params['MPloss'].microseconds / 1000000.0
    if mp_sec > 0.0:
        # setup loss vector
        LSS = np.zeros(NL)
        LSS[np.where(Lgrid>params['Lmax'])] = -1./mp_sec
        f = f*np.exp(LSS*Tdelta)

    # add losses inside plasma pause, time scale taken from PPloss
    pp_sec = params['PPloss'].days * 86400.0 + params['PPloss'].seconds + \
             params['PPloss'].microseconds / 1000000.0
    if pp_sec > 0.0:
        # calculate plasma pause location after Carpenter & Anderson 1992
        LSS = np.zeros(NL)
        LSS[np.where(Lgrid<params['Lpp'])] = -1./pp_sec
        f = f*np.exp(LSS*Tdelta)

    return f

# -----------------------------------------------
def get_local_accel(Lgrid, params, SRC_model='JK1'):
    """
    calculate the diffusion coefficient D_LL
    """

    if SRC_model == 'JK1':
        magn = params['SRCmagn'].seconds
        Lcenter = 5.6
        Lwidth = 0.3
        Kp = params['Kp']
        S = magn*np.exp(-(Lgrid-Lcenter)**2/(2*Lwidth**2))*Kp*Kp

    return S
