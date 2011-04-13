#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions supporting radiation belt diffusion codes
"""

from spacepy import help
import ctypes
__version__ = "$Revision: 1.19 $, $Date: 2011/04/13 22:10:21 $"
__author__ = 'J. Koller, Los Alamos National Lab (jkoller@lanl.gov)'


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

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 17-Mar-2010 (JK)
    v1.01 12-May-2010 (dtw)
    
    """
    
    def __init__(self, grid='L', NL=91, const_kp=False):
        """
        format for grid e.g., L-PA-E
        """
        import numpy as n
        import spacepy.time as st
        import spacepy.lib

        self.const_kp=const_kp

        # Initialize the code to use to advance the solution.
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
                self.SRC_model = 'JK1'
                self.SRCmagn = st.Tickdelta(days=1e-1) # relative acceleration per day
                self.MPloss = st.Tickdelta(minutes=0.1) # minutes time scale
                self.PPloss = st.Tickdelta(days=10.) # days time scale
                self.set_lgrid(NL)
                self.MODerr = 5. # relative factor * PSD
                self.PSDerr = 0.3 # relative observational error
                self.MIN_PSD = 1e-99
                        
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

        Example:

        >>> start = datetime.datetime(2003,10,14)
        >>> end = datetime.datetime(2003,12,26)
        >>> delta = datetime.timedelta(hours=1)
        >>> rmod.setup_ticks(start, end, delta, dtype='UTC')
        """

        import spacepy.time as st
        
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
    def add_PSD(self, satlist=None):
        
        """
        add observations from PSD database using the ticks list        
        """
        
        import numpy as n
        import spacepy.sandbox.PSDdata as PD
        import spacepy.time
        
        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        Tgrid = self.ticks
        nTAI = len(Tgrid)
        
        self.PSDdata = ['']*(nTAI-1)
        for i, Tnow, Tfut in zip(n.arange(nTAI-1), Tgrid[:-1], Tgrid[1:]):
            start_end = spacepy.time.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
            self.PSDdata[i] = PD.get_PSD(start_end, self.MU, self.K, satlist)
            
        # adjust initial conditions to these PSD values
        mval = n.mean(self.PSDdata[0]['PSD'])
        self.PSDinit = mval*n.exp(-(self.Lgrid - 5.5)**2/0.2)
        
        return
        
    # -----------------------------------------------
    def evolve(self):
        """
        calculate the diffusion in L at constant mu,K coordinates    
        """
        from . import radbelt as rb
        import numpy as n
        
        assert 'ticks' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        
        f = self.PSDinit
        Tgrid = self.ticks.TAI
        nTAI = len(Tgrid)
        Lgrid = self.Lgrid
        self.PSD  = n.zeros( (len(f),len(Tgrid)), dtype=ctypes.c_double)
        self.PSD[:,0] = f

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
        params['SRC_model'] = self.SRC_model
        params['SRCmagn'] = self.SRCmagn
        params['MPloss'] = self.MPloss
        params['PPloss'] = self.PPloss
        if 'SRCartif' in self.__dict__:
            params['SRCartif'] = self.SRCartif
        keylist = ['Kp', 'Dst', 'Lmax', 'Lpp']
                
        # start with the first one since 0 is the initial condition
        for i, Tnow, Tfut in zip(n.arange(nTAI-1)+1, Tgrid[:-1], Tgrid[1:]):
            Tdelta = Tfut - Tnow
            Telapse= Tnow - Tgrid[0]
            # copy over parameter list
            for key in keylist:
                params[key] = self.params[key][i]
            # now integrate from Tnow to Tfut
            f = diff_LL(self, Lgrid, f, Tdelta, Telapse, params=params)
            self.PSD[:,i] = f

    # -----------------------------------------------
    def assimilate(self, method='enKF',inflation=0):
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

        >> import spacepy
        >> import datetime
        >> from spacepy import radbelt

        >> start = datetime.datetime(2002,10,23)
        >> end = datetime.datetime(2002,11,4)
        >> delta = datetime.timedelta(hours=0.5)
        >> rmod.setup_ticks(start, end, delta, dtype='UTC')

        Once the dates and time step are specified, the data is added using the
        add_PSD function:

        >> rmod.add_PSD()

        The observations are averaged over the time windows, whose interval is
        give by the time step. 

        Once the dates and data are set, the assimiation is performed using the
        'assimilate' function:

        >> rmod.assimilate(inflation=1)

        This function will add the PSDa values, which are the analysis state of
        the radiation belt using the observations within the dates. To plot the
        analysis simply use the plot funtion:

        >> rmod.plot(values=rmod.PSDa,clims=[-10,-6],Lmax=False,Kp=False,Dst=False)

        """
        import spacepy.borg
        import spacepy.sandbox.PSDdata as PD
        import spacepy.time as st
        import numpy as n
        import copy as c
        import pdb
        import spacepy.toolbox as tb

        # add PSD observations with add_PSD, 
        # this has to be done to the class
        # module when running the RBmodel.
        # add_PSD will add the PSDdata to the class
        # which is a dictionary for the data that has been added
        
        # debugging command 
        #pdb.set_trace()

        # setup method
        assert method in ['enKF'], 'DA method='+method+' not implemented'

        nTAI = len(self.ticks)

        # enKF method
        if method == 'enKF':
            da = spacepy.borg.enKF()

            # initialize A with initial condition
            # the initial condition is ones, have to change this to initialize
            # the ensemble with perturbed reference state
            A = n.ones( (self.NL, da.Nens) )*self.PSDinit[:,n.newaxis]

            self.PSDf = n.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa = n.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa[:,0] = self.PSDinit
            self.PSDf[:,0] = self.PSDinit

            # add model error (perturbation) in the ensemble initial condition.
            A = da.add_model_error(self, A, self.PSDdata[0])
            A[n.where(A < self.MIN_PSD)] = self.MIN_PSD

            # time loop
            for i, Tnow, Tfut in zip(n.arange(nTAI-1)+1, self.ticks[:-1], self.ticks[1:]):
                    
                # make forcast and add model error
                # make forecast using all ensembles in A
                iens = 0
                for f in A.T:
                    # create temporary RB class instance
                    rbtemp = c.copy(self)
                    rbtemp.ticks = st.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
                    rbtemp.PSDinit = f
                    rbtemp.evolve()
                    #rbtemp.PSDdata = self.PSDdata[i-1]
                    A[:,iens] = rbtemp.PSD[:,1]
                    iens += 1
        
                # save result in ff
                Tnow = Tfut
                self.PSDf[:,i] = n.mean(A, axis=1)

                # verify that there are data points within the interval, if
                # there are data points then extract average observations
                # within the window, if not return empty observation array y
                # and Lobs.
        
                if len(self.PSDdata[i-1]) > 0:
                    # get observations for time window ]Tnow-Twindow,Tnow]
                    Lobs, y = spacepy.borg.average_window(self.PSDdata[i-1], self.Lgrid)
                else:
                    y = n.array([])
                    Lobs = n.array([])

                print Lobs
                print y
                
                # then assimilate otherwise do another forcast
                if len(y) > 0:
        
                    # INFLATION SCHEMES
                    if inflation == 0:

                        print 'inflation around model state values at observation locations'
                        # Add model error (perturbation for the ensemble) around model
                        # state values.  This acts as an inflation scheme for EnKF
                        A = da.add_model_error(self, A, self.PSDdata[i-1])

                    elif inflation == 1:

                        print 'inflation around observation values at observation locations'
                        # Add model error (perturbation for the ensemble) around
                        # observation values. This acts as an inflation scheme for EnKF
                        A = da.add_model_error_obs(self, A, Lobs, y)

                    elif inflation == 2:


                        print 'inflation around ensemble average'
                        # Inflate around ensemble average for EnKF

                        # ensemble average
                        ens_avg = n.mean(A, axis=1)

                        # inflation factor
                        inflation_factor = 1.8

                        # loop over ensemble members and inflate
                        iens = 0
                        for ens in A.T:
                            # inflate ensemble
                            A[:,iens] = inflation_factor*(ens - ens_avg) + ens_avg
                            iens += 1

                        A[n.where(A < self.MIN_PSD)] = self.MIN_PSD

                
                    # prepare assimilation analysis

                    # project ensemble states to obs. grid
                    HA = da.getHA(self, Lobs, A)
                    
                    # measurement perturbations ensemble
                    Psi = da.getperturb(self, y)

                    # ensemble of innovation vectors
                    Inn = da.getInnovation(y, Psi, HA)

                    # calculate ensemble perturbation HA' = HA-HA_mean 
                    HAp = da.getHAprime(HA)

                    # now call the main analysis routine
                    if len(y) == 1:
                        A = da.analysis_oneobs(A, Psi, Inn, HAp)
                    else:
                        #A = da.analysis(A, Psi, Inn, HAp)
                        A = da.analysis_Evensen(A, Psi, Inn, HAp)
                
                    # check for minimum PSD values
                    A[n.where(A<self.MIN_PSD)] = self.MIN_PSD
                
                    # average A from analysis step and save in results
                    # dictionary
                    self.PSDa[:,i] = n.mean(A, axis=1)

                    # print assimilated result
                    Hx = n.zeros_like(y)
                    for iL,Lstar in enumerate(Lobs):
                        idx = n.where(tb.feq(self.Lgrid,Lstar))
                        Hx[iL] = self.PSDa[idx,i]
                    print Hx
                    print HA
                    
                elif len(y) == 0:

                    print 'no observations within this window'

                    self.PSDa[:,i] = self.PSDf[:,i]

                    continue # 
                    
                # print message
                print('Tnow: ', self.ticks[i].ISO)

    # -----------------------------------------------
    def plot(self, Lmax=True, Lpp=False, Kp=True, Dst=True, 
             clims=[0,10], title='Summary Plot', values=None):
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

        Example:
        
        >>> rb.plot(Lmax=False, Kp=False, clims=[2,10], title='Good work!')

        This command would create the summary plot with a color bar range
        of 100 to 10^10.  The Lmax line and Kp values would be excluded.
        The title of the topmost plot (phase space density) would be set to
        'Good work!'.
        """
        import numpy as n
        import matplotlib.pyplot as p
        from matplotlib.colors  import LogNorm
        from matplotlib.ticker  import (LogLocator, LogFormatter,
                                        LogFormatterMathtext)
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
                             n.where(values > 0.0, self.PSD, 10.0**-39),
                             vmin=10.0**clims[0], vmax=10.0**clims[1], 
                             norm=LogNorm())
        ax1.set_ylabel('L*')
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

        p.show()
        
        return fig, ax1, ax2, ax3

    # -----------------------------------------------
    def plot_obs(self, Lmax=True, Lpp=False, Kp=True, Dst=True, 
             clims=[0,10], title='Summary Plot', values=None):
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

        Example:
        
        >>> rb.plot_obs(Lmax=False, Kp=False, clims=[2,10], title='Observations Plot')

        This command would create the summary plot with a color bar range
        of 100 to 10^10.  The Lmax line and Kp values would be excluded.
        The title of the topmost plot (phase space density) would be set to
        'Good work!'.
        """
        import spacepy.borg
        import numpy as n
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
        y = n.array([],dtype=float)
        Lobs = n.array([],dtype=float)
        eDOYobs = n.array([],dtype=float)
        nTAI = len(self.ticks)
        # time loop
        for i, Tnow, Tfut in zip(n.arange(nTAI-1)+1, self.ticks[:-1], self.ticks[1:]):
            if len(values[i-1]) > 0:
                # get observations for time window ]Tnow-Twindow,Tnow]
                Lobs_tmp, y_tmp = spacepy.borg.average_window(values[i-1], self.Lgrid)
                y = n.append(y,y_tmp)
                Lobs = n.append(Lobs,Lobs_tmp)
                eDOYobs = n.append(eDOYobs,n.ones(len(y_tmp))*self.ticks.eDOY[i-1])

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
        values used in the calculation. 
        
        The output DLL is in units of units/day.
        """

        import numpy as n

        if DLL_model is 'BA2000': # Brautigam and Albert (2000)
            if type(self.const_kp) == type(0.0):
                Kp=self.const_kp
            else:
                Kp = params['Kp']
            alpha = 10.0**(0.506*Kp-9.325)
            beta = 10.0

        elif DLL_model is 'FC2006': # Fei and Chan (2006)
            alpha = 1.5e-6
            beta  = 8.5
        
        elif DLL_model is 'U2008': # Ukhorskiy (2008)
            alpha = 7.7e-6
            beta  = 6.0
        
        elif DLL_model is 'S1997': # Selesnick (1997)
            alpha = 1.9e-10
            beta  = 11.7
        elif DLL_model is 'const': # Constant DLL.
            alpha= 1.0
            beta = 1.0
            DLL  = n.zeros(len(Lgrid), dtype=ctypes.c_double)+10.
            # approximately BA2000 for Kp=1, L=1.
        else:
            raise ValueError, \
                "Radial diffusion model %s not implemented" % DLL_model
        

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
    spacepy.lib.solve_cn.  This code is very slow and should only be used when 
    the C code fails to compile.
    """
    
    import numpy as n
    import numpy.linalg as nlin

    # get grid and setup centered grid Lm=L_i-1/2, Lp=L_i+1/2
    Lmin = L[0]
    Lmax = L[-1]
    dL = L[1] - L[0]
    Lp = L + 0.5*dL
    Lm = L - 0.5*dL

    # setup the diffusion coefficient on each grid and center grid
    #Dllm = n.zeros(NL); Dllp = n.zeros(NL)
    betam = n.zeros(NL); betap = n.zeros(NL)
    for i in range(1,int(NL)-1,1):
        Dllp[i] = 0.5*(DLL[i]+DLL[i+1])
        Dllm[i] = 0.5*(DLL[i]+DLL[i-1])
        betam[i] = Dllm[i]*dt / (2*dL*dL*Lm[i]*Lm[i])
        betap[i] = Dllp[i]*dt / (2*dL*dL*Lp[i]*Lp[i])

    # initialize some arrays
    A = n.zeros( (NL,NL) )
    B = n.zeros( (NL,NL) )
    C = n.zeros( (NL,NL) )
    
    # off diagonal elements
    for i in n.arange(1,int(NL)-1,1):
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
    Ai = nlin.inv(A)
    C = n.dot(Ai,B)

    return n.dot(C, f)

# -----------------------------------------------
def diff_LL(r, grid, f, Tdelta, Telapsed, params=None):
    """
    calculate the diffusion in L at constant mu,K coordinates
    time units
    """

    import numpy as n
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
    src = n.zeros(NL)

    # Create source operators (source splitting) using implicit 
    # trapezoidal method to solve source ODE.
    # Add artificial sources
    if params.has_key('SRCartif'):
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
        src1 = n.zeros(NL)
        src2 = n.zeros(NL)

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
    if params['SRC_model']:
        # setup source vector
        S = get_local_accel(Lgrid, params, SRC_model='JK1')
        f = f + S*Tdelta
     
    # add losses through magnetopause shadowing, time scale taken from MPloss
    if params['MPloss'].seconds > 0.0:
        # setup loss vector
        LSS = n.zeros(NL)
        LSS[n.where(Lgrid>params['Lmax'])] = -1./params['MPloss'].seconds
        f = f*n.exp(LSS*Tdelta)
        
    # add losses inside plasma pause, time scale taken from PPloss
    if params['PPloss'].seconds > 0.0:
        # calculate plasma pause location after Carpenter & Anderson 1992
        LSS = n.zeros(NL)
        LSS[n.where(Lgrid<params['Lpp'])] = -1./params['PPloss'].seconds
        f = f*n.exp(LSS*Tdelta)

    return f
    
# -----------------------------------------------
def get_local_accel(Lgrid, params, SRC_model='JK1'):
    """
    calculate the diffusion coefficient D_LL
    """
    import numpy as n
    
    if SRC_model is 'JK1':        
        magn = params['SRCmagn'].seconds
        Lcenter = 5.6
        Lwidth = 0.3
        Kp = params['Kp']
        S = magn*n.exp(-(Lgrid-Lcenter)**2/(2*Lwidth**2))*Kp*Kp        

    
    return S
