#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions supporting radiation belt diffusion codes
"""

from spacepy import help
__version__ = "$Revision: 1.12 $, $Date: 2010/10/26 22:03:15 $"
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
    
    >>> start = datetime.datetime(2003,10,14)
    >>> end = datetime.datetime(2003,12,26)
    >>> delta = datetime.timedelta(hours=1)
    >>> rmod.setup_ticks(start, end, delta, dtype='UTC')

    Now, run the model over the enitre time range using the evolve method:

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
    
    def __init__(self, grid='L', NL=91):
        """
        format for grid e.g., L-PA-E
        """
        import numpy as n
        import spacepy.time as st
        
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
        self.PSDinit = zeros(self.NL)

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
        
        self.ticktock = st.tickrange(start, end, delta, dtype)
        
    # -----------------------------------------------    
    def add_omni(self, keylist=None):
        
        """
        add omni data to instance according to the tickrange in ticktock
        """
        import spacepy.omni as om

        assert 'ticktock' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        omni = om.get_omni(self.ticktock)
        
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
            
        assert self.ticktock, "Provide tick range with 'setup_ticks'"
        self.params['Lmax'] = em.get_Lmax(self.ticktock, Lmax_model)
                
    # -----------------------------------------------    
    def add_Lpp(self, Lpp_model):
        
        """
        add last closed drift shell Lmax
        """
        import spacepy.empiricals as em
 
        if 'params' not in self.__dict__:           
            self.params = {}
        
        assert self.ticktock, "Provide tick range with 'setup_ticks'"
        self.params['Lpp'] = em.get_plasma_pause(self.ticktock, Lpp_model)

    # -----------------------------------------------    
    def add_PSD(self, satlist=None):
        
        """
        add observations from PSD database using the ticktock list        
        """
        
        import numpy as n
        import spacepy.sandbox.PSDdata as PD
        import spacepy.time
        
        assert 'ticktock' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        Tgrid = self.ticktock
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
        import copy as c
        
        assert 'ticktock' in self.__dict__ , \
            "Provide tick range with 'setup_ticks'"
        
        f = self.PSDinit
        Tgrid = self.ticktock.TAI
        nTAI = len(Tgrid)
        Lgrid = self.Lgrid
        self.PSD  = n.zeros( (len(f),len(Tgrid)) )
        self.PSD[:,0] = c.copy(f)
        
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
            f = rb.diff_LL(Lgrid, f, Tdelta, Telapse, params)
            self.PSD[:,i] = c.copy(f)

    # -----------------------------------------------
    def assimilate(self, method='enKF'):
        """
        call data assimilation function in assimilate.py
        """
        import spacepy.borg
        import spacepy.sandbox.PSDdata as PD
        import spacepy.time as st
        import numpy as n
        import copy as c
        
        # setup method
        assert method in ['enKF'], 'DA method='+method+' not implemented'

        nTAI = len(self.ticktock)
        

        # enKF method
        if method == 'enKF':
            da = spacepy.borg.enKF()
            # initialize A with initial condition
            A = n.ones( (self.NL, da.Nens) )*self.PSDinit[:,n.newaxis]
            self.PSDf = n.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa = n.zeros( (self.PSDinit.shape[0], nTAI) )
            self.PSDa[:,0] = self.PSDinit
            self.PSDf[:,0] = self.PSDinit

            # time loop
            for i, Tnow, Tfut in zip(n.arange(nTAI-1)+1, self.ticktock[:-1], self.ticktock[1:]):
                    
                # make forcast and add model error
                # make forecast using all ensembles in A
                iens = 0
                for f in A.T:
                    # create temporary RB class instance
                    rbtemp = c.copy(self)
                    rbtemp.ticktock = st.Ticktock([Tnow.UTC[0], Tfut.UTC[0]], 'UTC')
                    rbtemp.PSDinit = f
                    rbtemp.evolve()
                    #rbtemp.PSDdata = self.PSDdata[i-1]
                    A[:,iens] = rbtemp.PSD[:,1]
                    iens += 1
        
                # save result in ff
                Tnow = Tfut
                self.PSDf[:,i] = n.mean(A, axis=1)
        
                # add model error 
                A = da.add_model_error(self, A, self.PSDdata[i-1])
                A[n.where(A < self.MIN_PSD)] = self.MIN_PSD
                
                # get observations for time window ]Tnow-Twindow,Tnow]
                Lobs, y = spacepy.borg.average_window(self.PSDdata[i-1], self.Lgrid)
                
                #Lobs = n.array(self.PSDdata[i-1]['Lstar'])
                #y = n.array(self.PSDdata[i-1]['PSD'])
                print(Lobs, y)
                
                if len(y) > 1:  # then assimilate otherwise do another forcast
        
                    # prepare assimilation analysis
                    HA = da.getHA(self, Lobs, A) # project ensemble states to obs. grid
                    Psi = da.getperturb(self, y) # measurement perturbations ensemble
                    Inn = da.getInnovation(y, Psi, HA) # ensemble of innovation vectors
                    HAp = da.getHAprime(HA) # calculate ensemble perturbation HA' = HA-HA_mean 
                    # now call the main analysis routine
                    A = da.analysis(A, Psi, Inn, HAp)
                
                    # check for minimum PSD values
                    A[n.where(A<self.MIN_PSD)] = self.MIN_PSD
                
                    # average A from analysis step and save in results dictionary
                    self.PSDa[:,i] = n.mean(A, axis=1)            
                    
                elif len(y) == 0:
                
                    self.PSDa[:,i] = n.mean(A, axis=1)
                    # save average innovation vector
                    self.Inn[i] = 0.0
                    continue # 
                    
                # print message
                print('Tnow: ', self.ticktock[i].ISO)


    # -----------------------------------------------
    def plot(self, Lmax=True, Lpp=False, Kp=True, Dst=True, 
             clims=[0,10], title='Summary Plot'):
        """
        Create a summary plot of the RadBelt object distribution function.
        For reference, the last closed drift shell, Dst, and Kp are all
        included.  These can be disabled individually using the corresponding 
        boolean kwargs.

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
        map = ax1.pcolorfast(self.ticktock.eDOY, self.Lgrid, 
                             n.where(self.PSD > 0.0, self.PSD, 10.0**-39),
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
            p.plot(self.ticktock.eDOY, self.params['Lmax'], 'w')
        # Minimize time range.
        ax1.set_xlim([self.ticktock.eDOY[0], self.ticktock.eDOY[-1]])
        # Finally, save the position of the plot to make next plot match.
        pos = ax1.get_position()
        
        if Dst is True:
            ax2 = p.subplot(2,1,2)
            pos2 = ax2.get_position()
            pos2.x1 = pos.x1
            ax2.set_position(pos2)
            ax2.plot(self.ticktock.eDOY, self.params['Dst'], color='r')
            ax2.set_xlabel('DOY in '+str(self.ticktock.UTC[0].year))
            ax2.set_ylabel('Dst', color='r')
            for tl in ax2.get_yticklabels():
                tl.set_color('r')
            ax2.set_xlim([self.ticktock.eDOY[0], self.ticktock.eDOY[-1]])

        if Kp is True:
            if Dst is True:
                p.subplot(2,1,2)
                ax3 = p.twinx()
            else:
                ax3 = p.subplot(2,1,2)
                pos3 = ax3.get_position()
                pos3.x1 = pos.x1
                ax3.set_position(pos3)
            ax3.plot(self.ticktock.eDOY, self.params['Kp'], 'k:')            
            if Dst is True:
                ax3.yaxis.tick_right()
            ax3.set_xlabel('DOY in '+str(self.ticktock.UTC[0].year))
            ax3.set_ylabel('Kp')                    
            ax3.set_xlim([self.ticktock.eDOY[0], self.ticktock.eDOY[-1]])

        p.show()
        
        return fig, ax1, ax2, ax3
        
# -----------------------------------------------
def diff_LL(grid, f, Tdelta, Telapsed, params=None):
    """
    calculate the diffusion in L at constant mu,K coordinates
    time units
    """

    import numpy as n
    
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
    (DLL, alpha, beta) = get_DLL(Lgrid, params, DLL_model)
    DLL = DLL/86400.
    C = get_modelop_L(Lgrid, Tdelta, DLL)

    # set boundaries
    f[-1] = f[-2]
    f[0] = f[1]

    # apply diffusion to f
    f = n.dot(C,f)
    
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

    # Add artificial sources
    if 'SRCartif' in params:
        # Call the artificial source function, sending info as 
        # key word arguments.  Note that such functions should be
        # able to handle extra kwargs through the use of **kwargs!
        f = f + params['SRCartif'](Telapsed, Lgrid, dll=DLL, 
                                   alpha=alpha, beta=beta) * Tdelta

    return f

# -----------------------------------------------
def get_modelop_L(Lgrid, Tdelta, DLL):
    """
    calculate the model oparator for a single small timestep
    """
    
    import numpy as n
    import numpy.linalg as nlin

    # get grid and setup centered grid Lm=L_i-1/2, Lp=L_i+1/2
    L = Lgrid
    NL = len(L)
    Lmin = L[0]
    Lmax = L[-1]
    dL = L[1] - L[0]
    Lp = L + 0.5*dL
    Lm = L - 0.5*dL

    # setup smaller timesteps
    k = 100
    dt = Tdelta/float(k)
        
    # setup the diffusion coefficient on each grid and center grid
    Dllm = n.zeros(NL); Dllp = n.zeros(NL)
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

    # now complete full timestep by ktimes
    C = n.array(  n.mat(C)**int(k) )
    
    return C

# -----------------------------------------------
def get_DLL(Lgrid, params, DLL_model='BA2000'):
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
        
    else:
        raise ValueError("Radial diffusion model %s not implemented" % DLL_model)
        

    DLL = alpha * Lgrid ** beta
    return DLL, alpha, beta
    
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
