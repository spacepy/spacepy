'''
A PyBats module for handling input, output, and visualization of 
binary SWMF output files taylored to BATS-R-US-type data.

.. currentmodule:: spacepy.pybats.bats

.. rubric:: Classes

.. autosummary::
    :template: clean_class.rst
    :toctree: autosummary

    BatsLog
    Stream
    Bats2d
    Mag
    MagFile
    GeoIndexFile
    VirtSat
'''

import numpy as np
from spacepy.pybats import PbData, IdlBin, LogFile, set_target
from spacepy.datamodel import dmarray

class BatsLog(LogFile):
    '''
    A specialized version of :class:`~spacepy.pybats.LogFile` that includes
    special methods for plotting common BATS-R-US log file values, such as
    D$_{ST}$.

    .. autosummary::

        ~BatsLog.add_dst_quicklook

    .. automethod:: add_dst_quicklook
    '''

    def add_dst_quicklook(self, target=None, loc=111, showObs=True, **kwargs):
        '''
        Create a quick-look plot of Dst (if variable present in file) 
        and compare against observations.
        
        Like all *add_\* * methods in Pybats, the *target* kwarg determines
        where to place the plot.
        If kwarg *target* is **None** (default), a new figure is 
        generated from scratch.  If *target* is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot location
        *loc* (defaults to 111).  If target is a matplotlib Axes object, 
        the plot is placed into that axis at subplot location *loc*.

        Observed Dst is automatically fetched from the Kyoto World Data Center
        via the :mod:`spacepy.pybats.kyoto` module.  The associated 
        :class:`spacepy.pybats.kyoto.KyotoDst` object, which holds the observed
        Dst, is stored as *self.obs_dst* for future use.

        The figure and axes objects are returned to the user.
        '''
        
        import matplotlib.pyplot as plt
        from spacepy.pybats import apply_smart_timeticks
        
        if 'dst' not in self:
            return None, None

        fig, ax = set_target(target, figsize=(10,4), loc=loc)

        ax.plot(self['time'], self['dst'], 
                label='BATS-R-US $D_{ST}$ (Biot-Savart)', **kwargs)
        ax.hlines(0.0, self['time'][0], self['time'][-1], 
                  'k', ':', label='_nolegend_')
        apply_smart_timeticks(ax, self['time'])
        ax.set_ylabel('D_${ST}$ ($nT$)')
        ax.set_xlabel('Time from '+ self['time'][0].isoformat()+' UTC')

        if(showObs):
            try:
                import spacepy.pybats.kyoto as kt
            except ImportError:
                return fig, ax
        
            try:
                stime = self['time'][0]; etime = self['time'][-1]
                if not hasattr(self, 'obs_dst'):
                    self.obs_dst = kt.fetch('dst', stime, etime)

            except BaseException as args:
                print('WARNING! Failed to fetch Kyoto Dst: ', args)
            else:
                ax.plot(self.obs_dst['time'], self.obs_dst['dst'], 
                        'k--', label='Obs. Dst')
                ax.legend(loc='best')
                apply_smart_timeticks(ax, self['time'])
        else:
            ax.legend(loc='best')
            

        return fig, ax

    
class Stream(object):
    '''
    A class for streamlines.  Contains all of the information about
    the streamline, including extracted variables.

    Upon instantiation, the object will trace through the vector
    field determined by the "[x/y]field" values and the Bats object
    "bats".

    Parameters
    ----------
    bats : Bats
        Bats object to trace through. (Should this be Bats2D?)
    xstart : float
        X value of location to start the trace.
    ystart : float
        Y value of location to start the trace.
    xfield : str
        Name of variable in ``bats`` which contains X values of the field
    yfield : str
        Name of variable in ``bats`` which contains Y values of the field

    Other Parameters
    ----------------
    style : str
        Sets line style, including colors. See :meth:`set_style` for details.
        (Default 'mag')
    type : str
        (Default 'streamline')
    method : str
        Integration method. The default is Runge-Kutta 4 ('rk4') which gives
        a good blend of speed and accuracy. See the test functions in
        :mod:`~spacepy.pybats.trace2d` for more info.  The other option is
        a simple Euler's method approach ('eul'). (Default 'rk4')
    extract : bool
        (Default: False)
    maxPoints : int
        (Default : 20000)

    Notes
    -----
    .. Not really "notes" but need to keep this section from being parsed as parameters

    .. rubric:: Methods

    .. autosummary::

        ~Stream.set_style
        ~Stream.treetrace
        ~Stream.trace
        ~Stream.plot

    .. automethod:: set_style
    .. automethod:: treetrace
    .. automethod:: trace
    .. automethod:: plot
    '''
    
    def __init__(self, bats, xstart, ystart, xfield, yfield,
                 style = 'mag', type='streamline', method='rk4',
                 extract=False, maxPoints=20000):
        # Key values:
        self.xstart = xstart #X and Y starting
        self.ystart = ystart #points in the field.
        self.xvar = xfield #Strings that list the variables
        self.yvar = yfield #that will be used for tracing.
        # Descriptors:
        self.type   = type
        self.open   = True
        self.method = method

        # Do tracing:
        self.treetrace(bats, extract, maxPoints=maxPoints)

        # Set style
        self.set_style(style)

    #def __repr__(self):
    #    pass
    #
    #def __str__(self):
    #    pass

    def set_style(self, style):
        '''
        Set the line style either using a simple matplotlib-type style
        string or using a preset style type.  Current types include:
        
        'mag' : treat line as a magnetic field line.  Closed lines are
                white, other lines are black.
        '''
        
        if style == 'mag':
            if self.open:
                self.style = 'k-'
            else:
                self.style = 'w-'
        else:
            self.style = style

    def treetrace(self, bats, extract=False, maxPoints=20000):
        '''
        Trace through the vector field using the quad tree.
        '''
        from numpy import array, sqrt, append
        if self.method == 'euler' or self.method == 'eul':
            from spacepy.pybats.trace2d import trace2d_eul as trc
        elif self.method == 'rk4':
            from spacepy.pybats.trace2d import trace2d_rk4 as trc
        else:
            raise ValueError('Tracing method {} not recognized!'.format(
                self.method))

        # Get name of dimensions in order.
        grid = bats['grid'].attrs['dims']

        # Find starting block and set starting locations.
        block = bats.find_block(self.xstart, self.ystart)
        xfwd=[self.xstart]; yfwd=[self.ystart]
        xnow, ynow = self.xstart, self.ystart

        # Trace forwards.
        while(block):
            # Grab indices of all values inside current block.
            # Ghost cell check is for experimental testing:
            if hasattr(bats.qtree[block], 'ghost'):
                loc = bats.qtree[block].ghost
            # Typical approach is no ghost cell info:
            else:
                loc = bats.qtree[block].locs
            # Trace through this block.
            x, y = trc(bats[self.xvar][loc], bats[self.yvar][loc], 
                       xnow, ynow, bats[grid[0]][loc][0,:], 
                       bats[grid[1]][loc][:,0], ds=0.01)
            # Update location and block:
            xnow, ynow = x[-1], y[-1]
            newblock = bats.find_block(xnow,ynow)
            # If we didn't leave the block, stop tracing.
            # Additionally, if inside rBody, stop.
            if(block==newblock) or (xnow**2+ynow**2)<bats.attrs['rbody']*.8 :
                block=False
            elif newblock:
                block=newblock
            else:
                block=False
            # Append to full trace vectors.
            xfwd = np.append(xfwd, x[1:])
            yfwd = np.append(yfwd, y[1:])
            del(x); del(y)
            # It's possible to get stuck swirling around across
            # a few blocks.  If we spend a lot of time tracing,
            # call it quits.
            if xfwd.size>maxPoints: block=False

        # Trace backwards.  Same Procedure as above.
        block = bats.find_block(self.xstart, self.ystart)
        xbwd=[self.xstart]; ybwd=[self.ystart]
        xnow, ynow = self.xstart, self.ystart
        while(block):
            if hasattr(bats.qtree[block], 'ghost'):
                loc = bats.qtree[block].ghost
            else:
                loc = bats.qtree[block].locs
            x, y = trc(bats[self.xvar][loc], bats[self.yvar][loc], 
                       xnow, ynow, bats[grid[0]][loc][0,:], 
                       bats[grid[1]][loc][:,0], ds=-0.01)
            xnow, ynow = x[-1], y[-1]
            newblock = bats.find_block(xnow,ynow)
            if(block==newblock) or (xnow**2+ynow**2)<bats.attrs['rbody']*.8 :
                block=False
            elif newblock:
                block=newblock
            else:
                block=False
            # Append to full trace vectors.
            xbwd = np.append(x[::-1], xbwd)
            ybwd = np.append(y[::-1], ybwd)
            if xbwd.size>maxPoints: 
                block=False

        # Combine foward and backward traces.
        self.x = np.append(xbwd[:-1],xfwd)
        self.y = np.append(ybwd[:-1],yfwd)

        # If planetary run w/ body:
        # 1) Check if line is closed to body.
        # 2) Trim points within body.
        if 'rbody' in bats.attrs:
            # Radial distance:
            r = sqrt(self.x**2.0  + self.y**2.0)
            # Closed field line?
            if (r[0] < bats.attrs['rbody']) and (r[-1] < bats.attrs['rbody']):
                self.open = False
            # Trim the fat!
            limit = bats.attrs['rbody']*.8
            self.x, self.y = self.x[r>limit], self.y[r>limit]

    def trace(self, bats, extract=False):
        '''
        Trace through the vector field.
        '''
        from numpy import array, sqrt
        if self.method == 'euler':
            from spacepy.pybats.trace2d import trace2d_eul as trc
        elif self.method == 'rk4':
            from spacepy.pybats.trace2d import trace2d_rk4 as trc
        
        # Get name of dimensions in order.
        grid = bats['grid'].attrs['dims']

        # Trace forward
        x1, y1 = trc(bats[self.xvar], bats[self.yvar], 
                     self.xstart, self.ystart,
                     bats[grid[0]], bats[grid[1]])
        # Trace backward
        x2, y2 = trc(bats[self.xvar], bats[self.yvar], 
                     self.xstart, self.ystart,
                     bats[grid[0]], bats[grid[1]], ds=-0.1)
        # Join lines together such that point 0 is beginning of line
        # and point -1 is the end (when moving parallel to line.)
        self.x = array(x2[::-1].tolist() + x1[1:].tolist())
        self.y = array(y2[::-1].tolist() + y1[1:].tolist())

        # Check if line is closed to body.
        if 'rbody' in bats.attrs:
            r1 = sqrt(self.x[0]**2.0  + self.y[0]**2.0)
            r2 = sqrt(self.x[-1]**2.0 + self.y[-1]**2.0)
            if (r1 < bats.attrs['rbody']) and (r2 < bats.attrs['rbody']):
                self.open = False

    def plot(self, ax, *args, **kwargs):
        '''
        Add streamline to axes object "ax". 
        '''
        ax.plot(self.x, self.y, self.style, *args, **kwargs)

class Bats2d(IdlBin):
    '''
    An object class of parent pybats.idlbin taylored to BATS-R-US output.
    This function requires the Matplotlib griddata function.
    
    '''
    # Init by calling IdlBin init and then building qotree, etc.
    def __init__(self, filename):
        import spacepy.pybats.qotree as qo
        reload(qo)
        from numpy import array
        # Read file.
        IdlBin.__init__(self, filename)

        # Parse grid into quad tree.
        if self['grid'].attrs['gtype'] != 'Regular':
            xdim, ydim = self['grid'].attrs['dims'][0:2]
            try:
                self.qtree=qo.QTree(array([self[xdim],self[ydim]]))
                self.find_block = self.qtree.find_leaf
            except:
                self.qtree=False
                self.find_block=lambda: False
            

    ####################
    # CALCULATIONS
    ####################

    def calc_temp(self, units='eV'):
        '''
        Calculate plasma temperature for each fluid.  Number density is
        calculated using *calc_ndens* if it hasn't been done so already.
        
        Temperature is obtained via density and pressure through the simple
        relationship P=nkT.

        Use the units kwarg to set output units.  Current choices are
        KeV, eV, and K.  Default is eV.
        '''

        from spacepy.datamodel import dmarray

        units = units.lower()

        # Create dictionary of unit conversions.
        conv  = {'ev' : 6241.50935,  # nPa/cm^3 --> eV.
                 'kev': 6.24150935,  # nPa/cm^3 --> KeV.
                 'k'  : 72429626.47} # nPa/cm^3 --> K.

        # Calculate number density if not done already.
        if not 'N' in self:
            self.calc_ndens()
        
        # Find all number density variables.
        for key in list(self.keys()):
            # Next variable if not number density:
            if key[-1] != 'N':
                continue
            # Next variable if no matching pressure:
            if not key[:-1]+'p' in self:
                continue
            self[key[:-1]+'t'] = dmarray(
                conv[units] * self[key[:-1]+'p']/self[key],
                attrs = {'units':units})

    def calc_b(self):
        '''
        Calculates total B-field strength using all three B components.
        Retains units of components.  Additionally, the unit vector
        b-hat is calculated as well.
        '''
        from numpy import sqrt

        self['b'] = sqrt(self['bx']**2.0 + self['by']**2.0 + self['bz']**2.0)
        self['b'].attrs={'units':self['bx'].attrs['units']}

        self['bx_hat'] = self['bx'] / self['b']
        self['by_hat'] = self['by'] / self['b']
        self['bz_hat'] = self['bz'] / self['b']

        self['bx_hat'].attrs={'units':'unitless'}
        self['by_hat'].attrs={'units':'unitless'}
        self['bz_hat'].attrs={'units':'unitless'}

    def calc_Econv(self):
        '''
        Calculates the convective electric field, -UxB.  Works for default
        MHD units of nT and km/s; if these units are not correct, an 
        exception will be raised.  Returns E_convective in mV/m.
        '''
        from copy import copy

        # Some quick declarations for more readable code.
        ux = self['ux']; uy=self['uy']; uz=self['uz']
        bx = self['bx']; by=self['by']; bz=self['bz']

        # Check units.  Should be nT(=Volt*s/m^2) and km/s.
        if (bx.attrs['units']!='nT') or (ux.attrs['units']!='km/s'):
            raise Exception('Incorrect units!  Should be km/s and nT.')

        # Calculate; return in millivolts per meter
        self['ex'] = -1.0*(uy*bz - uz*by) / 1000.0
        self['ey'] = -1.0*(uz*bx - ux*bz) / 1000.0
        self['ez'] = -1.0*(ux*by - uy*bx) / 1000.0
        self['ex'].attrs={'units':'mV/m'}
        self['ey'].attrs={'units':'mV/m'}
        self['ez'].attrs={'units':'mV/m'}
        

    def calc_ndens(self):
        '''
        Calculate number densities for each fluid.  Species mass is 
        ascertained either by searching file parameters or by assuming
        mass density from the species name (e.g. OpRho is clearly oxygen).
        If neither can be found, an exception is raised.

        New values are saved using the keys *SpeciesN* (e.g. *OpN*).
        '''

        from spacepy.datamodel import dmarray

        mass={'hp':1.0, 'op':16.0, 'he':4.0, 
              'sw':1.0, 'o':16.0, 'h':1.0, 'iono':1.0}
        species = []
        names   = []

        # Find all species, the variable names end in "Rho".
        for k in self:
            if (k[-3:] == 'rho')   \
                    and (k!='rho') \
                    and (k[:-3]+'N' not in self):
                species.append(k)
                names.append(k[:-3])
            if (k[:3] == 'rho')   \
                    and (k!='rho') \
                    and (k[3:]+'N' not in self):
                species.append(k)
                names.append(k[3:])
        # Individual number density
        for s, n in zip(species, names):
            self[n+'N'] = dmarray(self[s] / mass[n], 
                                       attrs={'units':'$cm^{-3}$'})

        # Total N is sum of individual  number densities.
        self['N'] = dmarray(np.zeros(self['rho'].shape), 
                            attrs={'units':'$cm^{-3}$'}) 
        if species:
            for n in names: self['N']+=self[n+'N']
        else:
            self['N'] += self['rho']
                                
    def calc_comp(self):
        '''
        Calculate composition by percent number for "minor" species
        (e.g. oxygen, helium).  If the required values are not present
        (i.e. 'rhoo', 'rhohe'), then the calculation is skipped without
        raising exceptions.
        '''

        from spacepy.datamodel import dmarray

        # Recognized mass densities.
        mass={'h':1.0,'he':4.0,'o':16.0}

        # Find all mass density variables that are species-specific.
        species = []
        for k in self:
            if ('rho' in k) and (len(k)>3):
                species.append(k.lower())
        if not k: return # No composition?  No calculation!

        # Create total number density.
        self['n'] = dmarray(np.zeros(self['rho'].shape), attrs={'units':'cm-3'})
        for s in species:
            self['n']+=self[s]/mass[s[3:]]

        # Loop through species calculating composition for each.
        # Create new variables with new attrs.
        for s in species:
            self['comp_'+s[3:]] = 100.0*self[s]/mass[s[3:]]/self['n']
            self['comp_'+s[3:]].attrs = {'units':'Percent'}

    def calc_beta(self):
        '''
        Calculates plasma beta (ratio of plasma to magnetic pressure, 
        indicative of who - plasma or B-field - is "in charge" with regards
        to flow.
        Assumes:
        -pressure in units of nPa
        -B in units of nT.
        Resulting value is unitless.
        Values where b_total = zero are set to -1.0 in the final array.
        '''
        from numpy import pi

        if not self.has_key('b'):
            self.calc_b()
        mu_naught = 4.0E2 * pi # Mu_0 x unit conversion (nPa->Pa, nT->T)
        temp_b = self['b']**2.0
        temp_b[temp_b<1E-8] =  -1.0*mu_naught*self['p'][temp_b==0.0]
        temp_beta=mu_naught*self['p']/temp_b
        #temp_beta[self['b']<1E-9] = -1.0
        self['beta']=temp_beta
        self['beta'].attrs={'units':'unitless'}

    def calc_jxb(self):
        '''
        Calculates the JxB force assuming:
        -Units of J are uA/m2, units of B are nT.
        Under these assumptions, the value returned is force density (nN/m^3).
        '''
        from numpy import sqrt
        from spacepy.datamodel import dmarray

        # Unit conversion (nT, uA, cm^-3 -> nT, A, m^-3) to nN/m^3.
        conv = 1E-6
        # Calculate cross product, convert units.
        self['jbx']=dmarray( (self['jy']*self['bz']-self['jz']*self['by'])*conv,
                             {'units':'nN/m^3'})
        self['jby']=dmarray( (self['jz']*self['bx']-self['jx']*self['bz'])*conv,
                             {'units':'nN/m^3'})
        self['jbz']=dmarray( (self['jx']*self['by']-self['jy']*self['bx'])*conv,
                             {'units':'nN/m^3'})
        self['jb'] =dmarray( sqrt(self['jbx']**2 +
                                  self['jby']**2 +
                                  self['jbz']**2), {'units':'nN/m^3'})

    def calc_alfven(self):
        '''
        Calculate the Alfven speed, B/(mu*rho)^(1/2) in km/s.  This is performed
        for each species and the total fluid.
        The variable is saved under key "alfven" in self.data.
        '''
        from numpy import sqrt, pi
        from spacepy.datamodel import dmarray
        
        if not self.has_key('b'):
            self.calc_b()
        #M_naught * conversion from #/cm^3 to kg/m^3
        mu_naught = 4.0E-7 * pi * 1.6726E-27 * 1.0E6

        for k in self:
            if (k[-3:]) == 'rho':
                # Alfven speed in km/s:
                self[k[:-3]+'alfven'] = dmarray(self['b']*1E-12 / 
                                                 sqrt(mu_naught*self[k]),
                                                 attrs={'units':'km/s'})

    def calc_divmomen(self):
        '''
        Calculate the divergence of momentum, i.e. 
        $\rho(u \dot \nabla)u$.
        '''

        from spacepy.datamodel import dmarray
        from spacepy.pybats.batsmath import d_dx, d_dy
        
        print("SUPER WARNING!  This is very, very exploratory.")

        # Create empty arrays to hold new values.
        size = self['ux'].shape
        self['divmomx'] = dmarray(np.zeros(size), {'units':'nN/m3'})
        self['divmomz'] = dmarray(np.zeros(size), {'units':'nN/m3'})

        # Units!
        c1 = 1000./6371.0 # km2/Re/s2 -> m/s2
        c2 = 1.6726E-21   # AMU/cm3 -> kg/m3
        c3 = 1E9          # N/m3 -> nN/m3.

        for k in self.qtree:
            # Calculate only on leafs of quadtree.
            if not self.qtree[k].isLeaf: continue
            
            # Extract values from current leaf.
            leaf = self.qtree[k]
            ux   = self['ux'][leaf.locs]
            uz   = self['uz'][leaf.locs]

            self['divmomx'][leaf.locs]=ux*d_dx(ux, leaf.dx)+uz*d_dy(ux, leaf.dx)
            self['divmomz'][leaf.locs]=ux*d_dx(uz, leaf.dx)+uz*d_dy(uz, leaf.dx)

        # Unit conversion.
        self['divmomx']*=self['rho']*c1*c2*c3
        self['divmomz']*=self['rho']*c1*c2*c3


    def calc_gradP(self):
        '''
        Calculate the pressure gradient force.
        '''

        from spacepy.datamodel import dmarray
        from spacepy.pybats.batsmath import d_dx, d_dy

        if 'p' not in self:
            raise KeyError('Pressure not found in object!')

        # Create new arrays to hold pressure.
        dims = self['grid'].attrs['dims']
        size = self['p'].shape
        self['gradP'] = dmarray(np.zeros(size), {'units':'nN/cm^3'})
        for d in dims:
            self['gradP_'+d] = dmarray(np.zeros(size), {'units':'nN/m^3'})

        for k in self.qtree:
            # Plot only leafs of the tree.
            if not self.qtree[k].isLeaf: continue

            # Extract leaf; place pressure into 2D array.
            leaf=self.qtree[k]
            z=self['p'][leaf.locs]
            
            # Calculate derivatives; place into new dmarrays.
            # Unit conversion: P in nPa => gradP in nN/m3, dx=Re
            # Convert to nN/m3 by multiplying by 6378000m**-1.
            # Clever indexing maps values back to "unordered" array so that
            # values can be plotted like all others.
            conv = 1.0/-6378000.0
            self['gradP_'+dims[0]][leaf.locs] = d_dx(z, leaf.dx)*conv
            self['gradP_'+dims[1]][leaf.locs] = d_dy(z, leaf.dx)*conv

        
        # Scalar magnitude:
        for d in dims:
            self['gradP'] += self['gradP_'+d]**2
        self['gradP'] = np.sqrt(self['gradP'])

    def calc_utotal(self):
        '''
        Calculate bulk velocity magnitude: $u^2 = u_X^2 + u_Y^2 + u_Z^2$.
        This is done on a per-fluid basis.
        '''

        from numpy import sqrt
        from spacepy.datamodel import dmarray

        species = []

        # Find all species, the variable names end in "ux".
        for k in self:
            if (k[-2:] == 'ux')   \
                    and (k[:-2]+'u' not in self):
                species.append(k[:-2])

        units = self['ux'].attrs['units']

        for s in species:
            self[s+'u'] = dmarray(sqrt( self[s+'ux']**2+
                                        self[s+'uy']**2+
                                        self[s+'uz']**2), 
                                  attrs={'units':units})

    def calc_Ekin(self, units='eV'):
        '''
        Calculate average kinetic energy per particle using 
        $E=\frac{1}{2}mv^2$.  Note that this is not the same as energy
        density.  Units are $eV$.
        '''
        from numpy import sqrt
        from spacepy.datamodel import dmarray

        print("WARNING - self.calc_Ekin: I think this is wrong!")
        
        conv =  0.5 * 0.0103783625 # km^2-->m^2, amu-->kg, J-->eV.
        if units.lower == 'kev':
            conv=conv/1000.0

        mass={'hp':1.0, 'op':16.0, 'he':4.0, 'sw':1.0, '':1.0}
        species = []

        # Find all species, the variable names end in "Rho".
        for k in self:
            #and (k!='rho') \
            if (k[-3:] == 'rho') and (k[:-3]+'Ekin' not in self):
                species.append(k[:-3])

        for s in species:
            #THIS IS WRONG HERE: 1/2mV**2?  Notsomuch.
            self[s+'Ekin'] = dmarray(sqrt( self[s+'ux']**2+
                                           self[s+'uy']**2+
                                           self[s+'uz']**2)
                                     * conv * mass[s],
                                     attrs={'units':units})


    def calc_all(self):
        '''
        Perform all variable calculations (e.g. calculations that
        begin with 'calc_').  Any exceptions raised by functions that
        could not be peformed (typicaly from missing variables) are
        discarded.
        '''
        for command in dir(self):
            if (command[0:5] == 'calc_') and (command != 'calc_all'):
                try:
                    eval('self.'+command+'()')
                except AttributeError, Error:
                    print('WARNING: Did not perform %s: %s' % (command, Error))

    #####################
    # Other calculations
    #####################
    def gradP_regular(self, cellsize=None, dim1range=-1, dim2range=-1):
        '''
        Calculate pressure gradient on a regular grid.
        Note that if the Bats2d object is not on a regular grid, one of 
        two things will happen.  If kwarg cellsize is set, the value of 
        cellsize will be used to call self.regrid and the object will
        be regridded using a cellsize of cellsize.  Kwargs dim1range and
        dim2range can be used in the same way they are used in self.regrid
        to restrict the regridding to a smaller domain.  If cellsize is
        not set and the object is on an irregular grid, an exception is
        raised.

        The gradient is calculated using numpy.gradient.  The output units
        are force density (N/m^3).  Three variables are added to self.data:
        gradp(dim1), gradp(dim2), gradp.  For example, if the object is an
        equatorial cut, the variables gradpx, gradpy, and gradp would be
        added representing the gradient in each direction and then the 
        magnitude of the vector.
        '''
        from numpy import gradient, sqrt

        if self.gridtype != 'Regular':
            if not cellsize:
                raise ValueError,('Grid must be regular or ' +
                                  'cellsize must be given.')
            self.regrid(cellsize, dim1range=dim1range, dim2range=dim2range)

        # Order our dimensions alphabetically.
        newvars = []
        for key in sorted(self.grid.keys()):
            newvars.append('gradp'+key)

        dx=self.resolution*6378000.0 # RE to meters
        p =self['p']*10E-9        # nPa to Pa
        self[newvars[0]], self[newvars[1]]=gradient(p, dx, dx)
        self['gradp']=sqrt(self[newvars[0]]**2.0 + self[newvars[1]]**2.0)

    def regrid(self, cellsize=1.0, dim1range=-1, dim2range=-1, debug=False):
        '''
        Re-bin data to regular grid of spacing cellsize.  Action is 
        performed on all data entries in the bats2d object.

        '''
        from matplotlib.mlab import griddata

        if self['grid'].attrs['gtype'] == 'Regular': return
        
        # Order our dimensions alphabetically.
        dims = self['grid'].attrs['dims']
        if debug: print("Ordered dimensions: ", dims)

        # Check to see if dimranges are 2-element lists.
        # If not, either set defaults or raise exceptions.
        if dim1range == -1:
            dim1range = [self[dims[0]].min(),self[dims[0]].max()]
        else:
            if isinstance(dim1range, ( type(()), type([]) ) ):
                if len(dim1range) != 2:
                    raise ValueError('dim1range must have two elements!')
            else: raise TypeError('dim1range must be a tuple or list!')
        if dim2range == -1:
            dim2range = [self[dims[1]].min(),self[dims[1]].max()]
        else:
            if isinstance(dim2range, ( type(()), type([]) ) ):
                if len(dim2range) != 2:
                    raise ValueError('dim2range must have two elements!')
            else: raise TypeError('dim2range must be a tuple or list!')

        if debug:
            print('%s range = %f, %f' % (dims[0], dim1range[0], dim1range[1]))
            print('%s range = %f, %f' % (dims[1], dim2range[0], dim2range[1]))

        # Now, Regrid.
        grid1 = np.arange(dim1range[0], dim1range[1]+cellsize, cellsize)
        grid2 = np.arange(dim2range[0], dim2range[1]+cellsize, cellsize)

        for key in self:
            # Skip grid-type entries.
            if key in (self['grid'].attrs['dims']+['grid']): continue
            self[key] = griddata(self[dims[0]], self[dims[1]],
                                 self[key], grid1, grid2)

        # Change grid, gridtype, gridsize, and npoints to match new layout.
        self['grid'].attrs['gtype'] = 'Regular'
        self['grid'].attrs['npoints'] = len(grid1) * len(grid2)
        self['grid'].attrs['resolution'] = cellsize
        self[dims[0]] = grid1; self[dims[1]] = grid2



    #############################
    # EXTRACTION/INTERPOLATION
    #############################
    def extract(self, x, y, vars='All'):
        '''
        For x, y of a 1D curve, extract values along that curve
        and return slice as a new data set.
        '''

        from spacepy.pybats import PbData, dmarray
        from spacepy.pybats.batsmath import interp_2d_reg

        # If our x, y locations are not numpy arrays, fix that.
        # Additionally, convert scalars to 1-element vectors.
        x, y = np.array(x), np.array(y)
        if not x.shape: x=x.reshape(1)
        if not y.shape: y=y.reshape(1)

        # Default: all variables are extracted except coordinates.
        if vars == 'All':
            vars = list(self.keys())
            for var in ('x','y','z','grid'):
                if var in vars: vars.remove(var)

        # Create data object for holding extracted values.
        # One vector for each value w/ same units as parent object.
        out = PbData(attrs=self.attrs)
        for var, val in zip(self['grid'].attrs['dims'], (x, y)):
            out[var] = dmarray(val, attrs=self[var].attrs)
        nPts = len(x)
        for var in vars:
            out[var]=dmarray(np.zeros(nPts), attrs=self[var].attrs)

        # Some helpers:
        xAll = self[self['grid'].attrs['dims'][0]]
        yAll = self[self['grid'].attrs['dims'][1]]

        # Navigate the quad tree, interpolating as we go.
        for k in self.qtree:
            # Only leafs, of course!
            if not self.qtree[k].isLeaf: continue
            # Find all points that lie in this block.
            # If there are none, just keep going.
            pts = \
                (x >= self.qtree[k].lim[0]) & (x <= self.qtree[k].lim[1]) & \
                (y >= self.qtree[k].lim[2]) & (y <= self.qtree[k].lim[3])
            if not pts.any(): continue
            locs = self.qtree[k].locs
            for var in vars:
                out[var][pts] = interp_2d_reg(x[pts], y[pts], xAll[locs], 
                                              yAll[locs], self[var][locs])

        return out

    ######################
    # TRACING TOOLS
    ######################
    def get_stream(self, x, y, xvar, yvar, method='rk4', style='mag',
                   maxPoints=40000):
        '''
        Trace a 2D streamline through the domain, returning a Stream
        object to the caller.

        x and y set the starting point for the tracing.

        xvar and yvar are string keys to self.data that define the
        vector field through which this function traces.  

        The method kwarg sets the numerical method to use for the
        tracing.  Default is Runge-Kutta 4 (rk4).
        '''

        stream = Stream(self, x, y, xvar, yvar, style=style, 
                        maxPoints=maxPoints, method=method)

        return stream

    ######################
    # VISUALIZATION TOOLS
    ######################
    def add_grid_plot(self, target=None, loc=111, DoLabel=True, DoShowNums=False,
                      title='BATS-R-US Grid Layout'):
        '''
        Create a plot of the grid resolution by coloring regions of constant
        resolution.  Kwarg "target" specifies where to place the plot and can
        be a figure, an axis, or None.  If target is a figure, a new subplot
        is created.  The subplot location is set by the kwarg "loc", which
        defaults to 111.  If target is an axis, the plot is placed into that
        axis object.  If target is None, a new figure and axis are created
        and used to display the plot. 

        Resolution labels can be disabled by setting kwarg DoLabel to False.

        Plot title is set using the 'title' kwarg, defaults to 'BATS-R-US
        Grid Layout'.

        Note that if target is not an axis object, the axis will automatically
        flip the direction of positive X GSM and turn on equal aspect ratio.
        In other words, if you want a customized axis, it's best to create
        one yourself.

        Figure and axis, even if none given, are returned.
        '''
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.ticker import MultipleLocator
        from numpy import linspace

        if self['grid'].attrs['gtype'] == 'Regular':
            raise ValueError('Function not compatable with regular grids')

        # Get dimensions over which we shall plot.
        xdim, ydim = self['grid'].attrs['dims'][0:2]

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=loc)
        ax.set_aspect('equal')

        # Set plot range based on quadtree.
        ax.set_xlim([self.qtree[1].lim[0],self.qtree[1].lim[1]])
        ax.set_ylim([self.qtree[1].lim[2],self.qtree[1].lim[3]])
        # Plot.

        for key in list(self.qtree.keys()):
            self.qtree[key].plotbox(ax)
        self.qtree.plot_res(ax, tagLeafs=DoShowNums)

        # Customize plot.
        ax.set_xlabel('GSM %s' % xdim.upper())
        ax.set_ylabel('GSM %s' % ydim.upper())
        ax.set_title(title)
        #if xdim=='x':
        #    ax.invert_xaxis()
        #if ydim=='y':
        #    ax.invert_yaxis()
        self.add_body(ax)

        return fig, ax

    def find_earth_lastclosed(self, tol=np.pi/360., method='rk4',
                              max_iter=100, debug=False):
        '''
        For Y=0 cuts, attempt to locate the last-closed magnetic field line
        for both day- and night-sides.  This is done using a bisection 
        approach to precisely locate the transition between open and closed
        geometries.  The method stops once this transition is found within
        a latitudinal tolerance of *tol*, which defaults to $\pi/360.$, or
        one-half degree.  The tracing *method* can be set via keyword and
        defaults to 'rk4' (4th order Runge Kutta, see 
        :class:`~spacepy.pybats.bats.Stream` for more information).
        The maximum number of iterations the algorithm will take is set
        by *max_iter*, which defaults to 100.

        This method returns 5 objects: 
        
        * The dipole tilt in radians
        * A tuple of the northern/southern hemisphere polar angle of 
          footpoints for the dayside last-closed field line.
        * A tuple of the northern/southern hemisphere polar angle of
          footpoints for the nightside last-closed field line.
        * The dayside last-closed field line as a 
          :class:`~spacepy.pybats.bats.Stream` object.
        * The nightside last-closed field line as a 
          :class:`~spacepy.pybats.bats.Stream` object.

        In each case, the angle is defined as elevation from the positive
        x-axis, in radians.

        '''

        import matplotlib.pyplot as plt
        from numpy import (cos, sin, pi, arctan, arctan2, sqrt)

        # Get the dipole tilt by tracing a field line near the inner
        # boundary.  Find the max radial distance; tilt angle == angle off
        # equator of point of min R (~=max |B|).
        x_small = self.attrs['rbody']*-1.2
        stream = self.get_stream(x_small, 0, 'bx', 'bz', method=method)
        r = stream.x**2 + stream.y**2
        loc = r==r.max()
        tilt = arctan(stream.y[loc]/stream.x[loc])[0]

        if debug:
            print('Dipole is tilted {} degress above the z=0 plane.'.format(
                tilt*180./pi))
        
        # Dayside- start by tracing from plane of min |B| and perp. to that: 
        R = self.attrs['rbody']*1.10
        s1 = self.get_stream(R*cos(tilt), R*sin(tilt), 'bx','bz', method=method)
        s2 = self.get_stream(R*cos(pi/2.+tilt), R*sin(pi/2.+tilt), 
                             'bx', 'bz', method=method)
        
        theta = tilt
        dTheta=np.pi/4. # Initially, search 90 degrees.
        nIter = 0
        while (dTheta>tol)or(s1.open):
            nIter += 1
            closed = not(s1.open)
            # Adjust the angle towards the open-closed boundary.
            theta += closed*dTheta 
            theta -= s1.open*dTheta
            # Trace at the new theta to further restrict angular range:
            s1 = self.get_stream(R*cos(theta), R*sin(theta), 'bx', 'bz', 
                                 method=method)
            dTheta /= 2.
            if nIter>max_iter:
                if debug: print('Did not converge before reaching max_iter')
                break

        # Use last line to get southern hemisphere theta:
        npts = s1.x.size/2
        r = sqrt(s1.x**2+s1.y**2)
        loc = np.abs(r-self.attrs['rbody'])==np.min(np.abs(
            r[:npts]-self.attrs['rbody']))
        xSouth, ySouth = s1.x[loc], s1.y[loc]
        # "+ 0" syntax is to quick-copy object.
        theta_day = [theta+0, 2*np.pi+arctan(ySouth/xSouth)[0]+0]
        day = s1

        # Nightside: 
        theta+=tol
        R = self.attrs['rbody']*1.0
        s2 = self.get_stream(R*cos(theta),R*sin(theta),'bx','bz',method=method)
        s1 = self.get_stream(R*cos(pi+tilt), R*sin(pi+tilt), 
                             'bx', 'bz', method=method)
        
        dTheta=(pi+tilt-theta)/2.
        theta = pi+tilt
        nIter = 0
        while (dTheta>tol)or(s1.open):
            nIter += 1
            closed = not(s1.open)
            theta -= closed*dTheta 
            theta += s1.open*dTheta
            #print("Theta = {:f} (dTheta={})".format(theta, dTheta))
            s1 = self.get_stream(R*cos(theta),R*sin(theta), 'bx','bz', 
                                 method=method)
            dTheta /= 2.
            if nIter>max_iter:
                if debug: print('Did not converge before reaching max_iter')
                break
            
        # Use last line to get southern hemisphere theta:
        npts = s1.x.size/2
        r = sqrt(s1.x**2+s1.y**2)
        loc = np.abs(r-self.attrs['rbody'])==np.min(np.abs(
            r[:npts]-self.attrs['rbody']))
        xSouth, ySouth = s1.x[loc], s1.y[loc]
        theta_night = [theta+0, 2*np.pi+arctan2(ySouth,xSouth)[0]+0]
        night = s1
        #plt.plot(s1.x, s1.y, 'r-')

        return tilt, theta_day, theta_night, day, night

    def add_b_magsphere_new(self, target=None, loc=111,  style='mag', 
                            DoLast=True, DoOpen=True, DoTail=True, 
                            method='rk4', tol=np.pi/360., DoClosed=True,
                            nOpen=5, nClosed=15, **kwargs):
        '''
        Create an array of field lines closed to the central body in the
        domain.  Add these lines to Matplotlib target object *target*.
        If no *target* is specified, a new figure and axis are created.

        A tuple containing the figure, axes, and LineCollection object
        is returned.  Any additional kwargs are handed to the LineCollection
        object.

        Note that this should currently only be used for GSM y=0 cuts
        of the magnetosphere.

        Algorithm:  This method, unlike its predecessor, starts by finding
        the last closed field lines via 
        :func:`~spacepy.pybats.bats.Bats2d.find_earth_lastclosed`.  It then
        fills the regions between the open and closed regions.  Currently, it
        does not treat purely IMF field lines.

        ========== ===========================================================
        Kwarg      Description
        ========== ===========================================================
        target     The figure or axes to place the resulting lines.
        style      The color coding system for field lines.  Defaults to 'mag'.
                   See :class:`spacepy.pybats.bats.Stream`.
        loc        The location of the subplot on which to place the lines.
        DoLast     Plot last-closed lines as red lines.  Defaults to **True**.
        DoOpen     Plot open field lines.  Defaults to **True**.
        DoClosed   Plot closed field lines.  Defaults to **True**.
        nOpen      Number of closed field lines to trace per hemisphere.  
                   Defaults to 5.
        nClosed    Number of open field lines to trace per hemisphere.
                   Defaults to 15.
        method     The tracing method; defaults to 'rk4'.   See 
                   :class:`spacepy.pybats.bats.Stream`.
        tol        Tolerance for finding open-closed boundary; see
                   :func:`~spacepy.pybats.bats.Bats2d.find_earth_lastclosed`.
        ========== ===========================================================
        
        Extra kwargs are passed to Matplotlib's LineCollection class.

        Three objects are returned: the figure and axes on which lines are
        placed and the LineCollection object containing the plotted lines.
        '''
        
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from numpy import (arctan, cos, sin, where, pi, log, 
                           arange, sqrt, linspace, array)
        
        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=111)
        self.add_body(ax)

        # Lines and colors:
        lines = []
        colors= []

        # Start by finding open/closed boundary.
        tilt, thetaD, thetaN, s1, s2 = self.find_earth_lastclosed(method=method)
        if DoLast:
            lines+=[array([s1.x,s1.y]).transpose(),
                    array([s2.x,s2.y]).transpose()]
            colors+=2*['r']

        # Useful parameters for the following traces:
        R = self.attrs['rbody']
        dTheta  = 1.5*np.pi/180.
        dThetaN = .05*np.abs(thetaN[0]-thetaD[0])
        dThetaS = .05*np.abs(thetaN[1]-thetaD[1])
        
        ## Do closed field lines ##
        if DoClosed:
            for tDay, tNit in zip(
                    linspace(0,     thetaD[0]-dTheta, nClosed),
                    linspace(np.pi, thetaN[1]-dTheta, nClosed)):
                x, y = R*cos(tDay), R*sin(tDay)
                sD   = self.get_stream(x,y,'bx','bz',method=method)
                x, y = R*cos(tNit), R*sin(tNit)
                sN   = self.get_stream(x,y,'bx','bz',method=method)
                # Append to lines, colors.
                lines.append(array([sD.x, sD.y]).transpose())
                lines.append(array([sN.x, sN.y]).transpose())
                colors.append(sD.style[0])
                colors.append(sN.style[0])

        ## Do open field lines ##
        if DoOpen:
            for tNorth, tSouth in zip(
                    linspace(thetaD[0]+dThetaN, thetaN[0]-dThetaN, nOpen),
                    linspace(thetaN[1]+dThetaS, thetaD[1]-dThetaS, nOpen)):
                x, y = R*cos(tNorth), R*sin(tNorth)
                sD   = self.get_stream(x,y,'bx','bz',method=method)
                x, y = R*cos(tSouth), R*sin(tSouth)
                sN   = self.get_stream(x,y,'bx','bz',method=method)
                # Append to lines, colors.
                lines.append(array([sD.x, sD.y]).transpose())
                lines.append(array([sN.x, sN.y]).transpose())
                colors.append(sD.style[0])
                colors.append(sN.style[0])  
                    
            
        # Create line collection & plot.
        collect = LineCollection(lines, colors=colors, **kwargs)
        ax.add_collection(collect)

        return fig, ax, collect

    def add_b_magsphere(self, target=None, loc=111, style='mag', 
                        DoImf=False, DoOpen=False, DoTail=False, 
                        DoDay=True, method='rk4', **kwargs):
        '''
        Create an array of field lines closed to the central body in the
        domain.  Add these lines to Matplotlib target object *target*.
        If no *target* is specified, a new figure and axis are created.

        A tuple containing the figure, axes, and LineCollection object
        is returned.  Any additional kwargs are handed to the LineCollection
        object.

        Note that this should currently only be used for GSM y=0 cuts
        of the magnetosphere.

        Method:
        First, the title angle is approximated by tracing a dipole-like
        field line and finding the point of maximum radial distance on
        that line.  This is used as the magnetic equator.  From this 
        equator, many lines are traced starting at the central body
        radius.  More lines are grouped together at higher magnetic
        latitude to give good coverage at larger L-shells.  Once an
        open field line is detected, no more lines are traced.  This
        is done on the day and night side of the central body.

        Because this method rarely captures the position of the night
        side x-line, more field lines are traced by marching radially
        outwards away from the furthest point from the last traced and
        closed field line.  This is repeated until open lines are found.
        '''
        
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from numpy import (arctan, cos, sin, where, pi, log, 
                           arange, sqrt, linspace, array)
        
        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=loc)

        lines = []
        colors= []

        # Approximate the dipole tilt of the central body.
        stream = self.get_stream(3.0, 0, 'bx', 'bz', method=method)
        r = stream.x**2 + stream.y**2
        loc, = where(r==r.max())
        tilt = arctan(stream.y[loc[0]]/stream.x[loc[0]])
        
        # Initial values:
        daymax   = tilt + pi/2.0
        nightmax = tilt + 3.0*pi/2.0

        # Day side:
        n = arange(25)+1.0
        angle = tilt + 5.0*pi*log(n)/(12.0*log(n[-1]))
        for theta in angle:
            x = self.attrs['rbody'] * cos(theta)
            y = self.attrs['rbody'] * sin(theta)
            stream = self.get_stream(x,y,'bx','bz', style=style, method=method)
            if (stream.y[0] > self.attrs['rbody']) or (stream.style[0]=='k'):
                daymax = theta
                break
            savestream = stream
            if DoDay:
                lines.append(array([stream.x, stream.y]).transpose())
                colors.append(stream.style[0])

        # Add IMF field lines.
        if DoImf:
            stream = savestream
            r = sqrt(stream.x**2 + stream.y**2)
            loc, = where(r==r.max())
            x_mp = stream.x[loc[0]]+0.15
            y_mp = stream.y[loc[0]]
            delx = 2.0
            for i, x in enumerate(arange(x_mp, 15.0, delx)):
                # From dayside x-line out and up:
                y =y_mp-x_mp+x
                stream = self.get_stream(x, y, 'bx', 'bz', style=style, 
                                         method=method)
                lines.append(array([stream.x, stream.y]).transpose())
                colors.append(stream.style[0])

                # From top of magnetosphere down:
                y =x_mp+15.0-x+delx/3.0
                stream = self.get_stream(x-delx/3.0, y, 'bx', 'bz', 
                                         method=method, style=style)
                lines.append(array([stream.x, stream.y]).transpose())
                colors.append(stream.style[0])

                # From bottom of mag'sphere down:
                y =x_mp-10.0-x+2.0*delx/3.0
                stream = self.get_stream(x-2.0*delx/3.0, y, 'bx', 
                                         'bz', style=style, method=method)
                lines.append(array([stream.x, stream.y]).transpose())
                colors.append(stream.style[0])

        # Night side:
        angle = pi + tilt + pi*log(n)/(2.5*log(n[-1]))
        for theta in angle:
            x = self.attrs['rbody'] * cos(theta)
            y = self.attrs['rbody'] * sin(theta)
            stream = self.get_stream(x,y,'bx','bz', style=style, method=method)
            if stream.open:
                nightmax = theta
                break
            savestream = stream
            lines.append(array([stream.x, stream.y]).transpose())
            colors.append(stream.style[0])

        
        # March down tail.
        stream = savestream
        r = sqrt(stream.x**2 + stream.y**2)
        loc, = where(r==r.max())
        x1 = stream.x[loc[0]]
        y1 = stream.y[loc[0]]
        x = x1
        y = y1
        while (x-1.5)>self['x'].min():
            #print "Closed extension at ", x-1.5, y
            #ax.plot(x-1.5, y, 'g^', ms=10)
            stream = self.get_stream(x-1.5, y, 'bx', 'bz', style=style, 
                                     method=method)
            r = sqrt(stream.x**2 + stream.y**2)
            if stream.open:
                break
            lines.append(array([stream.x, stream.y]).transpose())
            colors.append(stream.style[0])
            loc, = where(r==r.max())
            x = stream.x[loc[0]]
            y = stream.y[loc[0]]

        if x1 == x:
            stream = self.get_stream(x1+1.0, y1, 'bx', 'bz', method=method)
            r = sqrt(stream.x**2 + stream.y**2)
            loc, = where(r==r.max())
            x1 = stream.x[loc[0]]
            y1 = stream.y[loc[0]]

        ## Add more along neutral sheet.
        if DoTail:
            m = (y-y1)/(x-x1)
        #print "Slope = ", m
        #print "From y, y1 = ", y, y1
        #print "and x, x1 = ", x, x1
        #ax.plot([x,x1],[y,y1], 'wo', ms=10)
            xmore = arange(x, -100, -3.0)
            ymore = m*(xmore-x)+y
            for x, y  in zip(xmore[1:], ymore[1:]):
                stream = self.get_stream(x, y, 'bx', 'bz', style=style, 
                                         method=method)
                lines.append(array([stream.x, stream.y]).transpose())
                colors.append(stream.style[0])
        #ax.plot(xmore, ymore, 'r+')

        # Add open field lines.
        if DoOpen:
            for theta in linspace(daymax,0.99*(2.0*(pi+tilt))-nightmax, 15):
                x = self.attrs['rbody'] * cos(theta)
                y = self.attrs['rbody'] * sin(theta)
                stream = self.get_stream(x,y,'bx','bz', method=method)
                if stream.open: 
                    lines.append(array([stream.x, stream.y]).transpose())
                    colors.append(stream.style[0])
                x = self.attrs['rbody'] * cos(theta+pi)
                y = self.attrs['rbody'] * sin(theta+pi)
                stream = self.get_stream(x,y,'bx','bz', method=method)
                if stream.open:
                    lines.append(array([stream.x, stream.y]).transpose())
                    colors.append(stream.style[0])

        if 'colors' in kwargs:
            colors=kwargs['colors']
            kwargs.pop('colors')
        collect = LineCollection(lines, colors=colors, **kwargs)
        ax.add_collection(collect)

        return fig, ax, collect

    def add_planet(self, ax=None, rad=1.0, ang=0.0, **extra_kwargs):
        '''
        Creates a circle of radius=self.attrs['rbody'] and returns the
        MatPlotLib Ellipse patch object for plotting.  If an axis is specified
        using the "ax" keyword, the patch is added to the plot.

        Unlike the add_body method, the circle is colored half white (dayside)
        and half black (nightside) to coincide with the direction of the 
        sun. Additionally, because the size of the planet is not intrinsically
        known to the MHD file, the kwarg "rad", defaulting to 1.0, sets the
        size of the planet.

        Extra keywords are handed to the Ellipse generator function.
        '''

        from matplotlib.patches import Circle, Wedge

        if 'rbody' not in self.attrs:
            raise KeyError('rbody not found in self.attrs!')

        body = Circle((0,0), rad, fc='w', zorder=1000, **extra_kwargs)
        arch = Wedge((0,0), rad, 90.+ang, -90.+ang, fc='k', 
                     zorder=1001, **extra_kwargs)
        
        if ax != None:
            ax.add_artist(body)
            ax.add_artist(arch)

        return body, arch

    def add_body(self, ax=None, facecolor='lightgrey', DoPlanet=True, ang=0.0, 
                 **extra_kwargs):
        '''
        Creates a circle of radius=self.attrs['rbody'] and returns the
        MatPlotLib Ellipse patch object for plotting.  If an axis is specified
        using the "ax" keyword, the patch is added to the plot.
        Default color is light grey; extra keywords are handed to the Ellipse
        generator function.

        Because the body is rarely the size of the planet at the center of 
        the modeling domain, add_planet is automatically called.  This can
        be negated by using the DoPlanet kwarg.
        '''
        from matplotlib.patches import Ellipse

        if 'rbody' not in self.attrs:
            raise KeyError('rbody not found in self.attrs!')

        dbody = 2.0 * self.attrs['rbody']
        body = Ellipse((0,0),dbody,dbody,facecolor=facecolor, zorder=999,
                       **extra_kwargs)

        if DoPlanet:
            self.add_planet(ax, ang=ang)
        if ax != None:
            ax.add_artist(body)

    def add_pcolor(self, dim1, dim2, value, zlim=None, target=None, loc=111, 
                   title=None, xlabel=None, ylabel=None,
                   ylim=None, xlim=None, add_cbar=False, clabel=None,
                   add_body=True, dolog=False, *args, **kwargs):
        '''        
        Create a pcolor plot of variable **value** against **dim1** on the 
        x-axis and **dim2** on the y-axis.  Pcolor plots shade each
        computational cell with the value at the cell center.  Because no 
        interpolation or smoothin is used in the visualization, pcolor plots
        are excellent for examining the raw output.
        
        Simple example:

        >>> self.add_pcolor('x', 'y', 'rho')

        If kwarg **target** is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object, 
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        =========== ==========================================================
        Kwarg       Description
        =========== ==========================================================
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        title       Sets title of axes.  Default is no title.
        xlabel      Sets x label of axes.  Defaults to **dim1** and units.
        ylabel      Sets y label of axes.  Defaults to **dim2** and units.
        xlim        Sets limits of x-axes.  Defaults to whole domain.
        ylim        Sets limits of y-axes.  Defaults to whole domain.
        zlim        Sets color bar range.  Defaults to variable max/min.
        add_cbar    Adds colorbar to plot.  Defaults to *False*.
        clabel      Sets colorbar label.  Defaults to **var** and units.
        add_body    Places planetary body in plot.  Defaults to **True**.
        dolog       Sets use of logarithmic scale.  Defaults to **False**.
        =========== ==========================================================
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=loc)

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=self[value].min(); zlim[1]=self[value].max()
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Logarithmic scale?
        if dolog:
            z=np.where(self[value]>zlim[0], self[value], 1.01*zlim[0])
            norm=LogNorm()
        else:
            z=self[value]
            norm=None

        if self['grid'].attrs['gtype']=='Regular':
            pass
        else:
            # Indices corresponding to QTree dimensions:
            ix=self['grid'].attrs['dims'].index(dim1)
            iy=self['grid'].attrs['dims'].index(dim2)
            for k in self.qtree:
                # Plot only leafs of the tree.
                if not self.qtree[k].isLeaf: continue
                leaf=self.qtree[k]
                x=leaf.cells[ix]
                y=leaf.cells[iy]
                z_local=z[leaf.locs]
                pcol=ax.pcolormesh(x,y,z_local,vmin=zlim[0],vmax=zlim[1],
                                   norm=norm,**kwargs)

        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(pcol, pad=0.01)
            if clabel==None: 
                clabel="%s (%s)" % (value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        if ylabel==None: ylabel='%s ($R_{E}$)'%dim2.upper()
        if xlabel==None: xlabel='%s ($R_{E}$)'%dim1.upper()
        ax.set_ylabel(ylabel); ax.set_xlabel(xlabel)
        if type(xlim)==type([]) and len(xlim)==2:
            ax.set_xlim(xlim)
        if type(ylim)==type([]) and len(ylim)==2:
            ax.set_ylim(ylim)

        # Add body/planet.
        self.add_body(ax)

        return fig, ax, pcol, cbar                              

    def add_contour(self, dim1, dim2, value, nlev=30, target=None, loc=111, 
                    title=None, xlabel=None, ylabel=None,
                    ylim=None, xlim=None, add_cbar=False, clabel=None,
                    filled=True, add_body=True, dolog=False, zlim=None,
                    *args, **kwargs):
        '''
        Create a contour plot of variable **value** against **dim1** on the 
        x-axis and **dim2** on the y-axis.
        
        Simple example:

        >>> self.add_contour('x', 'y', 'rho')

        If kwarg **target** is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object, 
        the plot is placed into that axis.

        Four values are returned: the matplotlib Figure and Axes objects,
        the matplotlib contour object, and the matplotlib colorbar object
        (defaults to *False* if not used.)

        =========== ==========================================================
        Kwarg       Description
        =========== ==========================================================
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        nlev        Number of contour levels.  Defaults to 30.
        title       Sets title of axes.  Default is no title.
        xlabel      Sets x label of axes.  Defaults to **dim1** and units.
        ylabel      Sets y label of axes.  Defaults to **dim2** and units.
        xlim        Sets limits of x-axes.  Defaults to whole domain.
        ylim        Sets limits of y-axes.  Defaults to whole domain.
        zlim        Sets contour range.  Defaults to variable max/min.
        add_cbar    Adds colorbar to plot.  Defaults to *False*.
        clabel      Sets colorbar label.  Defaults to **var** and units.
        add_body    Places planetary body in plot.  Defaults to **True**.
        dolog       Sets use of logarithmic scale.  Defaults to **False**.
        =========== ==========================================================
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import (LogNorm, Normalize)
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)

        # Set ax and fig based on given target.
        fig, ax = set_target(target, figsize=(10,10), loc=loc)

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=self[value].min(); zlim[1]=self[value].max()
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        # Set contour command based on grid type.
        if self['grid'].attrs['gtype'] != 'Regular':  # Non-uniform grids.
            if filled:
                contour=ax.tricontourf   
            else:
                contour=ax.tricontour   
        else:   # Uniform grids.
            if filled:
                contour=ax.contourf
            else:
                contour=ax.contour

        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]), 
                                            np.log10(zlim[1]), nlev))
            z=np.where(self[value]>zlim[0], self[value], 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z=self[value]
            norm=None
            ticks=None
            fmt=None

        # Plot contour.
        cont=contour(self[dim1],self[dim2],np.array(z),
                     levs,*args, norm=norm, **kwargs)
        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel==None: 
                clabel="%s (%s)" % (value, self[value].attrs['units'])
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        if ylabel==None: ylabel='%s ($R_{E}$)'%dim2.upper()
        if xlabel==None: xlabel='%s ($R_{E}$)'%dim1.upper()
        ax.set_ylabel(ylabel); ax.set_xlabel(xlabel)
        if type(xlim)==type([]) and len(xlim)==2:
            ax.set_xlim(xlim)
        if type(ylim)==type([]) and len(ylim)==2:
            ax.set_ylim(ylim)

        # Add body/planet.  Determine where the sun is first.
        if dim1.lower()=='x':
            ang=0.0
        elif dim2.lower()=='x':
            ang=90.0
        else: ang=0.0
        self.add_body(ax, ang=ang)

        return fig, ax, cont, cbar

class Mag(PbData):
    '''
    A container for data from a single BATS-R-US virtual magnetometer.  These
    work just like a typical :class:`spacepy.pybats.PbData` object.  Beyond
    raw magnetometer data, additional values are calculated and stored,
    including total pertubations (the sum of all global and ionospheric 
    pertubations as measured by the magnetometer).  Users will be interested
    in methods :meth:`~spacepy.pybats.bats.Mag.add_plot` and 
    :meth:`~spacepy.pybats.bats.Mag.recalc`.

    Instantiation is best done through :class: `spacepy.pybats.MagFile`
    objects, which load and parse organize many virtual magnetometers from a 
    single output file into a single object.  However, they can be created
    manually, though painfully.  Users must instantiate by handing the 
    new object the number of lines that will be parsed (rather, the number
    of data points that will be needed), a time vector, and (optionally) 
    the list of variables coming from the GM and IE module.  While the 
    latter two are keyword arguments, at least one should be provided.
    Next, the arrays whose keys were given by the *gmvars* and *ievars*
    keyword arguments in the instantiation step can either be filled manually
    or by using the :meth:`~spacepy.pybats.bats.Mag.parse_gmline` and
    :meth:`~spacepy.pybats.bats.Mag.parse_ieline` methods to parse lines of 
    ascii data from a magnetometer output file.  Finally, the
    :meth:`~spacepy.pybats.bats.Mag.recalc` method should be called to 
    calculate total perturbation.
    '''


    def __init__(self, nlines, time, gmvars=(), ievars=(), *args, **kwargs):
        from numpy import zeros

        super(Mag, self).__init__(*args, **kwargs)  # Init as PbData.

        self['time']=time
        self.attrs['nlines']=nlines
        
        self['x']=np.zeros(nlines)
        self['y']=np.zeros(nlines)
        self['z']=np.zeros(nlines)

        # Create IE and GM specific containers.
        for key in gmvars:
            self['gm_'+key]=np.zeros(nlines)
        for key in ievars:
            self['ie_'+key]=np.zeros(nlines)
            
    def parse_gmline(self, i, line, namevar):
        '''
        Parse a single line from a GM_mag*.out file and put into
        the proper place in the magnetometer arrays.  The line should 
        have the same number of variables as was initially given to the
        :class:`~spacepy.pybats.bats.Mag` object.  This method is best
        used through the :class:`~spacepy.pybats.bats.MagFile` class interface.

        Usage: 

        >>> self.parse_gmline(i, line, namevar)

        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts=line.split()
        self['x'][i]=float(parts[9])
        self['y'][i]=float(parts[10])
        self['z'][i]=float(parts[11])
        for j, key in enumerate(namevar):
            self['gm_'+key][i]=float(parts[j+12])

    def parse_ieline(self, i, line, namevar):
        '''
        Parse a single line from a IE_mag*.out file and put into
        the proper place in the magnetometer arrays.  The line should 
        have the same number of variables as was initially given to the
        :class:`~spacepy.pybats.bats.Mag` object.  This method is best
        used through the :class:`~spacepy.pybats.bats.MagFile` class interface.


        Usage: 

        >>> self.parse_gmline(i, line, namevar)

        where i is the entry number, line is the raw ascii line, and
        namevar is the list of variable names.
        '''
        parts=line.split()
        for j, key in enumerate(namevar):
            self['ie_'+key][i]=float(parts[j+11])

    def recalc(self):
        '''
        Calculate total :math:`\Delta B` from GM and IE; store under object keys 
        *totaln*, *totale*, and *totald* (one for each component of the HEZ 
        coordinate system).
        '''
        from numpy import sqrt, zeros

        # New containers:
        self['totaln']=np.zeros(self.attrs['nlines'])
        self['totale']=np.zeros(self.attrs['nlines'])
        self['totald']=np.zeros(self.attrs['nlines'])

        for key in list(self.keys()):
            if key[-2:]=='Bn':
                self['totaln']=self['totaln']+self[key]
            if key[-2:]=='Be':
                self['totale']=self['totale']+self[key]
            if key[-2:]=='Bd':
                self['totald']=self['totald']+self[key]
        #self['totalh']=sqrt(self['totaln']**2.0 + self['totale']**2.0)

    def add_plot(self, value, style='-', target=None, loc=111, label=None, 
                 **kwargs):
        '''
        Plot **value**, which should be a key corresponding to a data vector
        stored in the :class:`~spacepy.pybats.bats.Mag` object,
        against the object's *time*.  The **target** kwarg specifies the
        destination of the plot.  If not set, **target** defaults to None and
        a new figure and axis will be created.  If target is a matplotlib
        figure, a new axis is created at subplot location 111 (which can be 
        changed
        using kwarg **loc**).  If target is a matplotlib Axes object, the line
        is added to the plot as if ``Axes.plot()`` was used.  The line label,
        which is used for setting labels on figure legends,
        defaults to the value key but can be customized with the "label"
        kwarg.  The **style** keyword accepts basic Matplotlib line style
        strings such as '--r' or '+g'.  This string will be passed on to the
        plot command to customize the line.

        All extra kwargs are handed to ``Axes.plot``, allowing the user to set
        any additional options (e.g., line color and style, etc.).

        Three values are returned: the Figure object, Axis object, and 
        newly created line object.  These can be used to further customize
        the figure, axis, and line as necessary.

        Example: Plot total :math:`\Delta B_N` onto an existing axis with line 
        color blue, line style dashed, and line label "Wow!":

        >>> self.plot('totaln', target='ax', label='Wow!', lc='b', ls='--')

        Example: Plot total :math:`\Delta B_N` on a new figure, save returned
        values and overplot additional values on the returned axis.  Default 
        labels and line styles are used in this example.
        
        >>> fig, ax, line = self.plot('totaln')
        >>> self.plot('gm_dBe', target = ax)

        '''

        import matplotlib.pyplot as plt
        from spacepy.pybats import apply_smart_timeticks

        if not label:
            label=value

        # Set figure and axes based on target:
        fig, ax = set_target(target, figsize=(10,4), loc=loc)

        line=ax.plot(self['time'], self[value], style, label=label, **kwargs)
        apply_smart_timeticks(ax, self['time'], dolabel=True)
        
        return fig, ax, line

    def add_comp_plot(self, direc, target=None, loc=111):
        '''
        Create a plot with,  or add to an existing plot, an illustration of 
        how the
        separate components sum together to make the total disturbance in a
        given orthongal direction (arg **direc**).  The three possible 
        components are 'n' (northwards, towards the magnetic pole), 'e'
        (eastwards), or 'd' (downwards towards the center of the Earth.)  The
        components of the total disturbance in any on direction are 
        magnetospheric currents ('gm_dB'), gap-region field-aligned currents
        ('gm_facdB'), and ionospheric Hall and Pederson currents ('ie_Jp' and
        'ie_Jh').  

        Example usage:
        
        >>> self.add_comp_plot('n')

        This will create a new plot with the total disturbance in the 'n' 
        direction along with line plots of each component that builds this
        total.  This method uses the familiar PyBats **target** kwarg system
        to allow users to add these plots to existing figures or axes.
        '''

        import matplotlib.pyplot as plt
        from spacepy.pybats import apply_smart_timeticks

        # Set plot targets.
        fig, ax = set_target(target, figsize=(10,4), loc=loc)

        # Use a dictionary to assign line styles, widths.
        styles={'gm':'--', 'ie':'-.', 'to':'-'}
        widths={'gm':1.0,  'ie':1.0,  'to':1.5}
        colors={'gm_dB':'#FF6600', 'gm_facdB':'r', 'ie_JhdB':'b', 
                'ie_JpB':'c','total':'k'}
        
        # Labels:
        labels={'gm_dB':r'$J_{Mag}$', 'gm_facdB':r'$J_{Gap}$', 
                'ie_JhdB':r'$J_{Hall}$', 'ie_JpB':r'$J_{Peder}$',
                'total':r'Total $\Delta B$'}

        # Plot.
        for k in sorted(self):
            if (k[-1]!=direc) or (k=='time'): continue
            ax.plot(self['time'], self[k], label=labels[k[:-1]],
                    ls=styles[k[:2]], lw=widths[k[:2]], c=colors[k[:-1]])

        # Ticks, zero-line, and legend:
        apply_smart_timeticks(ax, self['time'], True, True)
        ax.hlines(0.0, self['time'][0], self['time'][-1], 
                  linestyles=':', lw=2.0, colors='k')
        ax.legend(ncol=3, loc='best')
    
        # Axis labels:
        ax.set_ylabel(r'$\Delta B_{%s}$ ($nT$)'%(direc.upper()))


class MagFile(PbData):
    '''
    BATS-R-US magnetometer files are powerful tools for both research and
    operations.  :class:`~spacepy.pybats.bats.MagFile` objects open, parse, 
    and visualize such output.

    The $\delta B$ calculated by the SWMF requires two components: GM (BATSRUS)
    and IE (Ridley_serial).  The data is spread across two files: GM_mag*.dat
    and IE_mag*.dat.  The former contains $\delta B$ caused by gap-region 
    (i.e., inside the inner boundary) FACs and the changing global field.  
    The latter contains the $\delta B$ caused by Pederson and Hall 
    currents in the ionosphere.  :class:`~spacepy.pybats.bats.MagFile objects
    can open one or both of these files at a time; when both are opened, the
    total $\delta B$ is calculated and made available to the user.

    Usage: 

    >>> # Open up the GM magnetometer file only.
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag')
    >>>
    >>> # Open up both the GM and IE file.
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag', 'IE_file.mag')
    >>>
    >>> # Open up the GM magnetometer file; search for the IE file.
    >>> obj = spacepy.pybats.bats.MagFile('GM_file.mag', find_ie=True)

    Note that the **find_ie** kwarg uses a simple search assuming the data 
    remain in a typical SWMF-output organizational tree (i.e., if the results
    of a simulation are in folder *results*, the GM magnetometer file can be 
    found in *results/GM/* or *results/GM/IO2/* while the IE file can be found
    in *results/IE/* or *results/IE/ionosphere/*).  It will also search the
    present working directory.  This method is not robust; the user must take
    care to ensure that the two files correspond to each other.
    '''
    
    def __init__(self, filename, ie_name=None, find_ie=False, *args, **kwargs):
        from glob import glob
        super(MagFile, self).__init__(*args, **kwargs)  # Init as PbData.
        self.attrs['gmfile']=filename
        if(find_ie and not ie_name):
            # Try to find the IE file based on our current location.
            basedir = filename[0:filename.rfind('/')+1]
            if glob(basedir + '../IE/IE_mag_*.mag'):
                self.attrs['iefile']=glob(basedir + '../IE/IE_mag_*.mag')[-1]
            elif glob(basedir + '../../IE/ionosphere/IE_mag_*.mag'):
                self.attrs['iefile'] = \
                    glob(basedir + '../../IE/ionosphere/IE_mag_*.mag')[-1]
            elif glob(basedir + '/IE_mag_*.mag'):
                self.attrs['iefile']= glob(basedir + '/IE_mag_*.mag')[-1]
            else:
                self.attrs['iefile']=None
        else:
            self.attrs['iefile']=ie_name
        self.readfiles()

    def readfiles(self):
        '''
        Read and parse GM file and IE file (if name given.)
        '''
        import datetime as dt
        from numpy import zeros

        # Slurp lines.
        infile = open(self.attrs['gmfile'], 'r')
        lines=infile.readlines()
        infile.close()

        # Parse header.
        trash=lines.pop(0) # Get station names.
        nmags=int((trash.split(':')[0]).split()[0])
        names=trash.split(':')[1] # Remove number of magnetometers.
        namemag = names.split()
        self.attrs['namemag']=namemag
        # Check nmags vs number of mags in header.
        if nmags != len(namemag):
            raise BaseException(
                'ERROR: GM file claims %i magnetomers, lists %i'
                % (nmags, len(namemag)))
        trash=lines.pop(0)
        gm_namevar = trash.split()[12:] #Vars skipping time, iter, and loc.
        self.attrs['gm_namevar']=gm_namevar
        self.attrs['nmag']=len(namemag)
        nlines = len(lines)/nmags

        # If there is an IE file, Parse that header, too.
        if self.attrs['iefile']:
            infile=open(self.attrs['iefile'], 'r')
            ielns =infile.readlines()
            infile.close()
            trash=ielns.pop(0)
            nmags=int((trash.split(':')[0]).split()[0])
            iestats=(trash.split(':')[1]).split()
            # Check nmags vs number of mags in header.
            if nmags != len(iestats):
                raise BaseException(
                    'ERROR: IE file claims %i magnetomers, lists %i'
                    % (nmags, len(namemag)))
            if iestats != self.attrs['namemag']:
                raise RuntimeError("Files do not have matching stations.")
            ie_namevar=ielns.pop(0).split()[11:]
            self.attrs['ie_namevar']=ie_namevar
            if (len(ielns)/self.attrs['nmag']) != (nlines-1):
                print('Number of lines do not match: GM=%d, IE=%d!' %
                      (nlines-1, len(ielns)/self.attrs['nmag']))
                nlines=min(ielns, nlines-1)
        else:
            ie_namevar=()
            self.attrs['ie_namevar']=()

        # Build containers.
        self['time']=np.zeros(nlines, dtype=object)
        self['iter']=np.zeros(nlines, dtype=float)
        for name in namemag:
            self[name]=Mag(nlines, self['time'], gm_namevar, ie_namevar)

        # Read file data.
        for i in range(nlines):
            line = lines.pop(0)
            parts=line.split()
            self['iter'][i]=parts[0]
            self['time'][i]=dt.datetime(
                int(parts[1]), #year
                int(parts[2]), #month
                int(parts[3]), #day
                int(parts[4]), #hour
                int(parts[5]), #minute
                int(parts[6]), #second
                int(parts[7])*1000 #microsec
                )
            self[namemag[0]].parse_gmline(i, line, gm_namevar)
            for j in range(1, nmags):
                self[namemag[j]].parse_gmline(i, lines.pop(0), gm_namevar)
            if self.attrs['iefile'] and i>0:
                line=ielns.pop(0)
                self[namemag[0]].parse_ieline(i, line, ie_namevar)
                for j in range(1, nmags):
                    self[namemag[j]].parse_ieline(i, ielns.pop(0), ie_namevar)
        self.recalc()
        
        # Get time res.
        self.attrs['dt']=(self['time'][1]-self['time'][0]).seconds/60.0
        
    def recalc(self):
        '''
        Sum all contributions from all models/regions to get total 
        :math:`\Delta B`.
        '''
        for mag in self.attrs['namemag']:
            self[mag].recalc()

class GeoIndexFile(LogFile):
    '''
    Geomagnetic Index files are a specialized BATS-R-US output that contain
    geomagnetic indices calculated from simulated ground-based magnetometers.
    Currently, the only index instituted is Kp through the faKe_p setup.  
    Future work will expand the system to include Dst, AE, etc.

    GeoIndFiles are a specialized subclass of pybats.LogFile.  It includes
    additional methods to quickly visualize the output, perform data-model
    comparisons, and more.
    '''

    def __repr__(self):
        return 'GeoIndexFile object at %s' % (self.attrs['file'])

    def __init__(self, filename, *args, **kwargs):
        '''
        Load ascii file located at self.attrs['file'].  
        '''
        # Call super __init__:
        super(GeoIndexFile, self).__init__(filename, *args, **kwargs)

        # Re-parse the file header; look for key info.
        head  = (self.attrs['descrip']).replace('=',' ')
        parts = head.split()
        if 'DtOutput=' in head:
            self.attrs['dt'] = float(parts[parts.index('DtOutput=')+1])
        if 'SizeKpWindow' in head:
            self.attrs['window'] = float(parts[parts.index(
                        'SizeKpWindow(Mins)')+1])
        if 'Lat' in head:
            self.attrs['lat'] = float(parts[parts.index('Lat')+1])
        if 'K9' in head:
            self.attrs['k9']  = float(parts[parts.index('K9') +1])
                                         

    def add_kp_quicklook(self, target=None, loc=111, label=None, 
                         plot_obs=True, **kwargs):
        '''
        Similar to "dst_quicklook"-type functions, this method fetches observed
        Kp from the web and plots it alongside the Kp read from the GeoInd file.
        Usage:
        >>> obj.kp_quicklook(target=SomeMplTarget)
        The target kwarg works like in other PyBats plot functions: it can be
        a figure, an axes, or None, and it determines where the plot is placed.

        Other kwargs customize the line.  Label defaults to fa$K$e$_{P}$, extra
        kwargs are passed to pyplot.plot.
        '''
        import matplotlib.pyplot as plt
        from spacepy.pybats import apply_smart_timeticks

        # Set up plot target.
        fig, ax = set_target(target, figsize=(10,4), loc=loc)
        
        # Create label:
        if not(label):
            label = 'fa$K$e$_{P}$'
            if ('lat' in self.attrs) and ('k9' in self.attrs):
                label += ' (Lat=%04.1f$^{\circ}$, K9=%03i)' % \
                    (self.attrs['lat'], self.attrs['k9'])

        # Sometimes, the "Kp" varname is caps, sometimes not.
        kp = 'Kp'
        if kp not in self.keys(): kp='kp'

        ax.plot(self['time'], self['Kp'], label=label,**kwargs)
        ax.set_ylabel('$K_{P}$')
        ax.set_xlabel('Time from '+ self['time'][0].isoformat()+' UTC')
        ax.grid()
        apply_smart_timeticks(ax, self['time'])

        fig.tight_layout()

        if plot_obs:
            try:
                import spacepy.pybats.kyoto as kt
            except ImportError:
                print("kyotodst package unavailable.")
                return fig, ax
        
            try:
                stime = self['time'][0]; etime = self['time'][-1]
                if not hasattr(self, 'obs_kp'):
                    self.obs_kp = kt.fetch('kp', (stime.year, stime.month), 
                                           (etime.year, etime.month))
            except BaseException as args:
                print('WARNING! Failed to fetch Kyoto Kp: ' + args)
            else:
                self.obs_kp.add_histplot(target=ax, color='k', ls='--', lw=3.0)
                ax.legend(loc='best')
                apply_smart_timeticks(ax, self['time'])
                
        return fig, ax

class VirtSat(LogFile):
    '''
    A :class:`spacepy.pybats.LogFile` object tailored to virtual satellite
    output; includes special satellite-specific plotting methods.
    '''

    def __init__(self, *args, **kwargs):

        from re import findall
        from scipy.interpolate import interp1d
        from matplotlib.dates import date2num

        super(VirtSat, self).__init__(*args, **kwargs)
        # Attempt to extract the satellite's name and save it.
        try:
            s = self.attrs['file']
            a = findall('sat_([\-\w]+)_\d+_n\d+\.sat|sat_(\w+)_n\d+\.sat', s)[0]
            name = filter(None, a)[0]
        except IndexError:
            name = None
        self.attrs['name']=name

        # Create interpolation functions for position.
        self._interp={}
        tnum=date2num(self['time'])
        for x in 'xyz':
            self._interp[x] = interp1d(tnum, self[x])

    def __repr__(self):
        return 'Satellite object {} at {}'.format(self.attrs['name'],
                                                  self.attrs['file'])

    def calc_ndens(self):
        '''
        Calculate number densities for each fluid.  Species mass is 
        ascertained either by searching file parameters or by assuming
        mass density from the species name (e.g. OpRho is clearly oxygen).
        If neither can be found, an exception is raised.

        New values are saved using the keys *SpeciesN* (e.g. *OpN*).
        '''

        from spacepy.datamodel import dmarray

        mass={'hp':1.0, 'op':16.0, 'he':4.0, 'sw':1.0, 'o':16.0, 'h':1.0}
        species,names = [],[]

        # Find all species, the variable names end in "Rho".
        for k in self:
            if (k[-3:] == 'rho')   \
                    and (k!='rho') \
                    and (k[:-3]+'N' not in self.keys()):
                species.append(k)
                names.append(k[:-3])
            if (k[:3] == 'rho')   \
                    and (k!='rho') \
                    and (k[3:]+'N' not in self.keys()):
                species.append(k)
                names.append(k[3:])
        # Individual number density
        for s,n in zip(species, names):
            self[n+'N'] = dmarray(self[s] / mass[n], 
                                  attrs={'units':'$cm^{-3}$'})

        # Total N is sum of individual  number densities.
        self['N'] = dmarray(np.zeros(self['rho'].shape), 
                            attrs={'units':'$cm^{-3}$'}) 
        if names:
            for n in names: self['N']+=self[n+'N']
        else:
            self['N'] += self['rho']
                                

    def calc_temp(self, units='eV'):
        '''
        Calculate plasma temperature for each fluid.  Number density is
        calculated using *calc_ndens* if it hasn't been done so already.
        
        Temperature is obtained via density and pressure through the simple
        relationship P=nkT.

        Use the units kwarg to set output units.  Current choices are
        KeV, eV, and K.  Default is eV.
        '''

        from spacepy.datamodel import dmarray

        units = units.lower()

        # Create dictionary of unit conversions.
        conv  = {'ev' : 6241.50935,  # nPa/cm^3 --> eV.
                 'kev': 6.24150935,  # nPa/cm^3 --> KeV.
                 'k'  : 72429626.47} # nPa/cm^3 --> K.

        # Calculate number density if not done already.
        if not 'N' in self:
            self.calc_ndens()
        
        # Find all number density variables.
        for key in self.keys():
            # Next variable if not number density:
            if key[-1] != 'N':
                continue
            # Next variable if no matching pressure:
            if not key[:-1]+'p' in self:
                continue
            self[key[:-1]+'t'] = dmarray(
                conv[units] * self[key[:-1]+'p']/self[key],
                attrs = {'units':units})

    def calc_bmag(self):
        '''
        Calculates total magnetic field magnitude via:
        $|B|=\sqrt{B_X^2+B_Y^2+B_Z^2}$.  Results stored as variable name *b*.
        '''

        if 'b' not in self:
            self['b']=np.sqrt(self['bx']**2+
                              self['by']**2+
                              self['bz']**2)

        return True

    def calc_magincl(self, units='deg'):
        '''
        Magnetic inclination angle (a.k.a. inclination angle) is defined:
        $\Theta = sin^{-1}(B_Z/B)$.  It is a crucial value when examining
        magnetic dynamics about geosychronous orbit, the tail, and other
        locations.

        This function calculates the magnetic inclination for the Virtual 
        Satellite object and saves it as *b_incl*.  Units default to degrees;
        the keyword **units** can be changed to 'rad' to change this.
        '''

        if 'b' not in self: self.calc_bmag()

        incl = np.arcsin(self['bz']/self['b'])

        if units=='deg':
            self['b_incl'] = dmarray(incl*180./np.pi, {'units':'degrees'})
        elif units=='rad':
            self['b_incl'] = dmarray(incl, {'units':'radians'})
        else:
            raise ValueError('Unrecognized units.  Use "deg" or "rad"')

        return True
            
    def get_position(self, time):
        '''
        For an arbitrary time, *time*, return a tuple of coordinates for
        the satellite's position at that time in GSM coordinates.

        *time* should be type datetime.datetime, matplotlib.dates.date2num
        output (i.e., number of days (fraction part represents hours,
        minutes, seconds) since 0001-01-01 00:00:00 UTC, *plus* *one*), or
        a sequence of either.
        The satellite's position is interpolated (linearly) to *time*.
        '''
        
        from matplotlib.dates import date2num

        # Test if sequence.
        if isinstance(time, (list, np.ndarray)):
            testval=time[0]
        else:
            testval=time

        # Test if datetime or not.
        if type(testval) == type(self['time'][0]):
            time=date2num(time)

        # Interpolate, using "try" as to not pass time limits and extrapolate.
        try:
            loc = (self._interp['x'](time), 
                   self._interp['y'](time), 
                   self._interp['z'](time))
        except ValueError:
            loc = (None, None, None)

        return loc

    def add_sat_loc(self, time, target, plane='XY', dobox=False,
                    dolabel=False, size=12, c='k', **kwargs):
        '''
        For a given axes, *target*, add the location of the satellite at
        time *time* as a circle.  If kwarg *dolabel* is True, the satellite's 
        name will be used to label the dot.  Optional kwargs are any accepted by
        matplotlib.axes.Axes.plot.  The kwarg *plane* specifies the plane
        of the plot, e.g., 'XY', 'YZ', etc.
        '''
        
        # Get Axes' limits:
        xlim = target.get_xlim()
        ylim = target.get_ylim()
            
        plane=plane.lower()
        loc = self.get_position(time)
        if None in loc: return
        x=loc['xyz'.index(plane[0])]
        y=loc['xyz'.index(plane[1])]

        # Do not label satellite if outside axes bounds.
        if (x<min(xlim))or(x>max(xlim))or \
           (y<min(ylim))or(y>max(ylim)): dolabel=False

        target.plot(x, y, 'o', **kwargs)
        if dolabel:
            xoff = 0.03*(xlim[1]-xlim[0])
            if dobox:
                target.text(x+xoff,y,self.attrs['name'], 
                            bbox={'fc':'w','ec':'k'},size=size, va='center')
            else:
                target.text(x+xoff,y,self.attrs['name'], 
                            size=size, va='center', color=c)
                

        # Restore Axes' limits.
        target.set_ylim(ylim)
        target.set_xlim(xlim)

    def add_orbit_plot(self, plane='XY', target=None, loc=111, rbody=1.0,
                       title=None, trange=None):
        '''
        Create a 2D orbit plot in the given plane (e.g. 'XY' or 'ZY').
        '''
        
        from spacepy.pybats.ram import add_body, grid_zeros
        import matplotlib.pyplot as plt

        fig, ax = set_target(target, figsize=(5,5), loc=loc)

        plane=plane.upper()

        # Set time range of plot.
        if not trange: trange = [self['time'].min(), self['time'].max()]
        tloc = (self['time']>=trange[0])&(self['time']<=trange[-1])
        
        # Extract orbit X, Y, or Z.
        plane=plane.lower()
        if plane[0] in ['x','y','z']:
            x = self[plane[0]][tloc]
        else:
            raise ValueError('Bad dimension specifier: ' + plane[0])
        if (plane[1] in ['x','y','z']) and (plane[0]!=plane[1]):
            y = self[plane[1]][tloc]
        else:
            raise ValueError('Bad dimension specifier: ' + plane[1])
                

        ax.plot(x, y, 'g.')
        add_body(ax, rad=rbody, add_night=('X' in plane))

        # Finish customizing axis.
        ax.axis('equal')
        ax.set_xlabel('GSM %s'%(plane[0].upper()))
        ax.set_ylabel('GSM %s'%(plane[1].upper()))
        if title:
            ax.set_title(title)
        grid_zeros(ax)
        ax.grid()
        #set_orb_ticks(ax)

        return fig, ax
