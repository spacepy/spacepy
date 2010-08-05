'''
A module for reading, handling, and plotting RAM-SCB output.  
'''

import numpy as np
#from matplotlib.dates   import (date2num, DayLocator, MinuteLocator,
#                                DateFormatter, HourLocator, DayLocator,
#                                date2num, num2date)
#from matplotlib.ticker  import (MultipleLocator, LogLocator, LogFormatter,
#                                FuncFormatter, LogFormatterMathtext)
#from matplotlib.colors  import LogNorm
#from matplotlib.patches import Ellipse

############################################################################
#  A few useful functions:
def gen_egrid(nE=36, lb=0.1, ew=3E-2, power=1.27):
    '''
    Given number of points and a lower boundary (lb), build the RAM-SCB 
    energy grid.  Three arrays are returned: ebound(nE+1), the boundary
    of each bin, ecenter(nE), the center of each bin, and ewidth(nE),
    the width of each bin.  All output units are in KeV.

    Grid parameters can be set with kwargs but default to common RAM values.
    lb = lower boundary
    ew = first energy bin width
    power = power series exponent to determine grid spacing.
    '''
    
    ebound = np.zeros(nE+1)
    ecentr = np.zeros(nE)
    ewidth = np.zeros(nE)

    ebound[0] = lb
    ewidth[0] = ew
    ecentr[0] = ebound[0] + 0.5*ewidth[0]
    for i in range(nE-1):
        ewidth[i+1] = ewidth[i]*power
        ebound[i+1] = ebound[i] + ewidth[i]
        ecentr[i+1] = ecentr[i] + 0.5*(ewidth[i+1] + ewidth[i])
    ebound[-1] = ebound[-2] + ewidth[-1]

    return (ecentr, ebound, ewidth)

def viz_egrid(nE=36, lb=0.1, ew=3E-2, power=1.27):
    '''
    Vizualize the energy grid.  All kwargs correspond to those used by 
    gen_egrid.
    '''

    import matplotlib.pyplot as plt
    from numpy import zeros

    ecent, ebound, ewidth = gen_egrid(nE, lb, ew, power)
    y=zeros(len(ecent))

    fig = plt.figure(figsize=(12,4))
    fig.subplots_adjust(bottom=0.1,right=0.98, left=0.02)
    ax = fig.add_subplot(111)

    ax.errorbar(ecent, y, fmt='r', xerr=ewidth/2.0, label='Bin Widths',
                ecolor='r', elinewidth=3, capsize=30)
    ax.plot(ecent, y, 'g*', label='Bin Centers')
    ax.legend()
    ax.set_title('RAM-SCB Energy Grid - Power=%4.2f' % power)
    ax.set_xlabel('Energy ($KeV$)')
    ax.set_yticklabels([])
    ax.set_xlim([ebound[0]-20.0, ebound[-1]+20.0])
    ax.set_ylim([-1,1])

def gen_pgrid(nPa=72):
    '''
    Given number of points, build the RAM-SCB pitch angle grid.
    The return value is a nPax3 matrix containing the bin starts, centers,
    and ends in cosine of pitch angle.
    '''
    
    pgrid = np.zeros((nPa,3))

def read_t89file(filename):
    '''
    For some plots, we want T89 comparisons.  For now, we don't have that
    capability in Python (at least in a straight-forward manner.)  Hence,
    we've written some T89 files that need to be read on demand.

    Usage:
    time, b_dip, b_ext = read_t89file('somefile.t89')
    '''

    infile = open(filename, 'r')
    raw = infile.readlines()
    raw.pop(0) # skip header.
    
    nRec = len(raw)
    time = []
    bdip = np.zeros( (nRec, 3) )
    bext = np.zeros( (nRec, 3) )

    for i, line in enumerate(raw):
        parts = line.split()
        time.append(dt.datetime(
                int(parts[0]),
                int(parts[1]),
                int(parts[2]),
                int(parts[3]),
                int(parts[4]),
                int(parts[5]) ))
        bdip[i,:] = (float(parts[6]), float(parts[7]), float(parts[8]))
        bext[i,:] = (float(parts[9]), float(parts[10]),float(parts[11]))

    return time, bdip, bext

def grid_zeros(axis):
    '''
    Attempt to plot x=0 and y=0 gridlines w/o messing up 
    plot range.  This should be called last when creating
    a plot; after you have the range sorted out.
    '''
    xrng = axis.get_xlim()
    yrng = axis.get_ylim()
    axis.plot([-1*10^5,10^5], [0,0], 'k-')
    axis.plot([0,0], [-1*10^5,10^5], 'k-')
    axis.set_xlim(xrng)
    axis.set_ylim(yrng)
    
def set_orb_ticks(axis):
    '''
    Set major ticks to multiples of 1, minor ticks to 1/4.
    '''
    
    from matplotlib.ticker import MultipleLocator

    # Tick Locators:
    xMticks = MultipleLocator(1.00)
    xmticks = MultipleLocator(0.25)
    yMticks = MultipleLocator(1.00)
    ymticks = MultipleLocator(0.25)
    axis.xaxis.set_major_locator(xMticks)
    axis.xaxis.set_minor_locator(xmticks)
    axis.yaxis.set_major_locator(yMticks)
    axis.yaxis.set_minor_locator(ymticks)
    axis.grid(True, ls='k:')

def add_orbits(sat, title='', timelim=None):
    '''
    Add three orbit plots to the top of a plot.
    This is an internal function and should only be called by rampy.
    '''
    
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse

    if not timelim: 
        # Set default time limit if none given.
        timelim = [sat.time[0], sat.time[-1]]
        iMin=0
        iMax=-1
    else:
        # Use timelim to get indices that bound our plots.
        timediff = abs(sat.time - timelim[-1])
        iMax = np.nonzero(timediff == timediff.min())[0][0]
        timediff = abs(sat.time - timelim[0])
        iMin = np.nonzero(timediff == timediff.min())[0][0]

    # Here's a quicky to make arrows along our orbits.
    def make_arrow(index, coord1, coord2, plot):
        startx = sat.data['SM_xyz'][index,coord1]
        starty = sat.data['SM_xyz'][index,coord2]
        lengthx = sat.data['SM_xyz'][index+10,coord1] - startx
        lengthy = sat.data['SM_xyz'][index+10,coord2] - starty
        xlims = plot.get_xlim()
        dx = xlims[1] - xlims[0]
        arrow = plt.Arrow(startx, starty, lengthx, lengthy, 
                          width=dx/15.0, fc='g', ec='g')
        plot.add_patch(arrow)

    # x-y plane.
    orb1 = plt.subplot(4,3,1)
    # Show orbit:
    orb1.plot(sat.data['SM_xyz'][iMin:iMax,0], 
              sat.data['SM_xyz'][iMin:iMax,1], 'g')
    # Show Earth:
    ell1 = Ellipse((0,0),1,1,facecolor='b', alpha=0.2)
    orb1.add_artist(ell1)
    # Show axis:
    orb1.axis('equal')
    grid_zeros(orb1)
    # Labels, range, etc.
    orb1.set_xlabel('SM X')
    orb1.set_ylabel('SM Y')
    set_orb_ticks(orb1)
    # Annotate with arrows.
    #for j in range(0,len(sat.time)-10, 50):
    #make_arrow(0,0,1, orb1)
    #make_arrow(-11,0,1,orb1)
    
    # x-z plane.
    orb2 = plt.subplot(4,3,2)
    orb2.plot(sat.data['SM_xyz'][iMin:iMax,0],
              sat.data['SM_xyz'][iMin:iMax,2],'g')
    ell2 = Ellipse((0,0),1,1,facecolor='b', alpha=0.2)
    orb2.add_artist(ell2)
    # Show axis:
    orb2.axis('equal')
    grid_zeros(orb2)
    # Labels, range, etc.
    orb2.set_xlabel('SM X')
    orb2.set_ylabel('SM Z')
    set_orb_ticks(orb2)
    orb2.set_title(title)
    # Annotate with arrows.
    #for j in range(0,len(sat.time)-10, 50):
    #make_arrow(0,0,2,orb2)
    #make_arrow(-11,0,2,orb2)

    # y-z plane.
    orb3 = plt.subplot(4,3,3)
    orb3.plot(sat.data['SM_xyz'][iMin:iMax,1], 
              sat.data['SM_xyz'][iMin:iMax,2],'g')
    ell3 = Ellipse((0,0),1,1,facecolor='b', alpha=0.2)
    orb3.add_artist(ell3)
    # Show axis:
    orb3.axis('equal')
    grid_zeros(orb3)
    # Labels, range, etc.
    orb3.set_xlabel('SM y')
    orb3.set_ylabel('SM Z')
    set_orb_ticks(orb3)
    # Annotate with arrows.
    #for j in range(0,len(sat.time)-10, 50):
    #make_arrow(0,1,2, orb3)
    #make_arrow(-11,1,2, orb3)

def smart_timeticks(time):
    '''
    Given a list or array of time values, create intelligent timeticks based
    on the max and min of the input.
    Return three items: major tick locator, minor tick locator, and a 
    format string.

    Example:
    >>>Mtick, mtick, fmt = smart_timeticks([time1, time2])
    '''
    
    import matplotlib.dates as mdt

    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if nHours < 1:
        Mtick=mdt.MinuteLocator(byminute=[0,15,30,45])
        mtick=mdt.MinuteLocator(byminute=range(60), interval=5)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick=mdt.MinuteLocator(byminute=[0,30])
        mtick=mdt.MinuteLocator(byminute=range(60), interval=10)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick=mdt.HourLocator(byhour=range(24), interval=2)
        mtick=mdt.MinuteLocator(byminute=[0,15,30,45])
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 24:
        Mtick=mdt.HourLocator(byhour=[0,3,6,9,12,15,18,21])
        mtick=mdt.HourLocator(byhour=range(24))
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick=mdt.HourLocator(byhour=[0,6,12,18])
        mtick=mdt.HourLocator(byhour=range(24))
        fmt = mdt.DateFormatter('%H:%M UT')
    else:
        Mtick=mdt.DayLocator(bymonthday=range(1,32))
        mtick=mdt.HourLocator(byhour=[0,6,12,18])
        fmt = mdt.DateFormatter('%d %b')

    return(Mtick, mtick, fmt)

def apply_smart_timeticks(ax, time, dolimit=True):
    '''
    Given an axis 'ax' and a list/array of datetime objects, 'time', 
    use the smart_timeticks function to built smart time ticks and
    then immediately apply them to the give axis.

    The range of the 'time' input value will be used to set the limits
    of the x-axis as well.  Set kwarg 'dolimit' to False to override 
    this behavior.
    '''

    Mtick, mtick, fmt = smart_timeticks(time)
    ax.xaxis.set_major_locator(Mtick)
    ax.xaxis.set_minor_locator(mtick)
    ax.xaxis.set_major_formatter(fmt)
    if dolimit:
        ax.set_xlim([time[0], time[-1]])

    return True


############################################################################
class weimer(object):
    '''
    Read a Weimer electric field file.
    '''
    
    def __init__(self,filename):
        self.filename = filename
        self.read()

    def read(self):
        import datetime as dt

        infile = open(self.filename, 'r')
        
        # Get year, time:
        head = infile.readline()
        year = int(head.split()[-1])
        
        head = infile.readline()
        parts= head.split()
        doy = int(parts[0])
        hour= float(parts[1])
        
        self.time = dt.datetime(year, 1,1,0,0) + \
            dt.timedelta(days=doy-1) + \
            dt.timedelta(hours=hour)

        # skip rest of crap header.
        head = infile.readline()
        head = infile.readline()

        # Read data
        raw = infile.readlines()
        self.epot = np.zeros(len(raw)-1)
        self.l   = np.zeros(len(raw)-1)
        self.phi = np.zeros(len(raw)-1)

        for i, line in enumerate(raw[0:-1]):
            parts = line.split()
            self.l[i]    = float(parts[0])
            self.phi[i]  = float(parts[1])
            self.epot[i] = float(parts[2])

        # Make CPCP
        self.cpcp = self.epot.max() - self.epot.min()

        infile.close()

class ramsat(object):
    '''
    An object to handle NetCDF-formatted satellite output files.
    '''

    def __init__(self,filename):
        # Filename:
        self.filename=filename
        # Basic attributes:
        self.starttime=None
        self.data = {}
        self.global_attr = {}
        self.var_attr = {}
        self.dimensions = {}
        self.time = np.array([])

        # Read and parse NetCDF object.
        self.read()

    def read(self):
        import scipy.io
        import datetime as dt

        data = scipy.io.netcdf_file(self.filename, 'r')

        # Attempt to grab the start time.  Fall to default
        # if it is not found.
        try:
            stringtime = data.attributes['start_time']
        except KeyError:
            stringtime = '20000101000000'

        #parse start time and use it to create time object.
        self.starttime = dt.datetime( 
            int(stringtime[0:4]),   #year
            int(stringtime[4:6]),   #month
            int(stringtime[6:8]),   #dat
            int(stringtime[8:10]),  #hour
            int(stringtime[10:12]), #minute
            int(stringtime[12:14])) #second

        # Use start time to build time array.
        time = []
        for deltaT in data.variables['Time'][:]:
            time.append(self.starttime + 
                        dt.timedelta(seconds=float(deltaT)))
        # Convert time list to a numpy array.
        self.time = np.array(time)

        # Extract the variables; place in friendlier form.
        for key in data.variables:
            try:
                self.data[key] = data.variables[key][:]
            except IndexError:
                self.data[key] = data.variables[key].getValue()
            self.var_attr[key] = data.variables[key].attributes
            
        # Other properties to grab "as is":
        self.global_attr = data.attributes
        self.dimensions  = data.dimensions

        # Check for uniform time cadence.  Irregularity is very 
        # possible if the satellite leaves the domain.
        dT = np.zeros(len(time)-1)
        for i in range(1, len(self.time)):
            dT[i-1] = data.variables['Time'][i] - data.variables['Time'][i-1]
        # Do we have dT?  If not, generate it:
        if not self.data.has_key('DtWrite'):
            self.dt = min(dT)
        else:
            self.dt = self.data['DtWrite']

        self.uniform = True
        if any(dT - self.dt): self.uniform = False

        # Close and finish.
        data.close()
        
    def orbit_formatter(self, x, pos):
        '''
        A function that, when passed to the FuncFormatter class of 
        Matplotlib tick formatters, produces specialized ticks that give
        UT, LT, and inclination.
        '''
        
        import datetime as dt
        from matplotlib.dates import num2date

        nt = num2date(x, tz=None)
        nowtime = dt.datetime(  # Need "naive" dt object.
            nt.year, nt.month, nt.day, nt.hour, nt.minute, 
            nt.second, tzinfo=None)
        deltime = (self.time-nowtime)
        mintime = min(abs(deltime))
        try:
            index = (deltime.tolist()).index(mintime)
        except ValueError:
            # Build format string.
            fmtstring = \
                '%02d:%02d UT\n* MLT\n* MLat\nR=*$R_{E}$' % \
                (nowtime.hour, nowtime.minute)
            return(fmtstring)

        # Get local time from XY coords.
        x = self.data['SM_xyz'][index,0]
        y = self.data['SM_xyz'][index,1]
        z = self.data['SM_xyz'][index,2]
        R = np.sqrt(y**2 + x**2 + z**2)
        Mlat  = 180.0*np.arcsin(z/R)/np.pi
        theta = 180.0*y/abs(y)*np.arccos(x/np.sqrt(y**2+x**2))/np.pi + 180.0
        theta = np.mod(theta, 360.0)
        locHr = np.floor(theta/15.0)
        locMn = int(60.0 * (theta/15.0 - locHr))


        # Build format string.
        fmtstring = \
            '%02d:%02d UT\n%02d:%02d MLT\n%5.1f$^{\circ}$ MLat\nR=%4.2f $R_{E}$'\
            % (nowtime.hour, nowtime.minute, locHr, locMn, Mlat, R)

        return (fmtstring)

    def time_regrid(self):
        '''
        There is always the potential for missing data points in 
        a ramsat object: if the satellite is outside of the domain,
        the virtual sat routines will not dump data for that time step.
        This creates the difficulty of non-uniform time axis in the flux
        arrays, preventing quick and clean MatPlotLib plotting.

        This function checks to see if time regridding is necessary, and,
        if so, converts all numpy arrays to uniformly spaced (in time) but
        masked numpy arrays.  This is necessary to do before using many
        2D plot functions.
        
        Usage: simply call this method to regrid the data in a ramsat obect.
            >>> sat = rampy.ramsat('somefile.ncdf')
            >>> sat.time_regrid()
            
        The new time resolution is always set to the minimum dt found in the
        original file or the given dt in the original file.
        '''
        
        if(self.uniform): return

        npoints = 1 + (self.time[-1] - self.time[0]).seconds / self.dt
        npoints = int(npoints)

        # Uniform time array:
        time = [self.time[0]]
        for i in range(0,npoints-1):
            time.append(time[i] + dt.timedelta(seconds=self.dt))


        # Moar work to do.


        self.time = np.array(time)

    def create_omniflux(self):
        '''
        Integrate the flux(currently as a function of energy and 
        pitch angle) to get omnidirectional flux as a function of
        energy.  New units = (cm-2*s*keV)^-1

        Usage: just call this method to integrate all species.
        '''
        
        # Create new flux attributes:
        nTime = self.time.size
        nPa   = self.data['pa_grid'].size
        nEner = self.data['energy_grid'].size
        self.omniH = np.zeros((nTime, nEner))
        self.omniHe= np.zeros((nTime, nEner))
        self.omniO = np.zeros((nTime, nEner))

        # Create delta mu, where mu = cos(pitch angle)
        dMu = np.zeros(nPa)
        dMu[0] = self.data['pa_grid'][1]
        for i in range(1,nPa): # Factor of pi here so we don't need it later.
            dMu[i] = 4*np.pi*self.data['pa_grid'][i] -self.data['pa_grid'][i-1]

        # Integrate.
        for i in range(1, nPa):
            self.omniH =self.omniH + self.data['FluxH+'][:,i,:] * dMu[i]
            self.omniHe=self.omniHe+ self.data['FluxHe+'][:,i,:]* dMu[i]
            self.omniO =self.omniO + self.data['FluxO+'][:,i,:] * dMu[i]

        # Mask out bad values.
        self.omniH = np.ma.masked_where(self.omniH <=0, self.omniH)
        self.omniHe= np.ma.masked_where(self.omniHe<=0, self.omniHe)
        self.omniO = np.ma.masked_where(self.omniO <=0, self.omniO)

    def plot_b(self, timelim=None, style='b-'):
        '''
        Plot satellite location and magnetic field onto a plot 
        that is the size of a "Letter" piece of paper.
        
        Syntax:
        -C{ramsatobject.plot_b([style=(-I{Style String}])}
        where -I{Style String} is a MatPlotLib style string that
        would be accepted by a plot command.  This provides the user
        with custom line formats for the magnetic field plots.  
        Default is 'b-', which is a blue line.
        
        
        Example usage:
            >>> import rampy
            >>> file = rampy.ramsat('somesatfile.ncdf')
            >>> file.plot_b('k--')
        
        This produces the plot on screen with black dashed lines
        in the magnetic field plots.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter

        # Set default time limit if none given.
        if not timelim: timelim = [self.time[0], self.time[-1]] 

        # Check style keyword for correct type.
        # Switch to default if incorrect.
        if type(style) != type('b-'):
            style='b-'

        # Create plot window.
        fig = plt.figure(figsize=(8.38,10.88))
        plt.subplots_adjust(hspace=0.30, wspace=0.36, 
                            top=0.93, right=0.94, bottom=0.12)

        add_orbits(self, title = 'Orbit and Magnetic Field for ' + \
                       self.filename[0:self.filename.rfind('.')])

        # Time tick handlers.
        Mtick, mtick, fmt = smart_timeticks(self.time)

        #Bx
        bx = plt.subplot(4,1,2)
        bx.grid()
        bx.plot(self.time, self.data['B_xyz'][:,0], ls=style)
        bx.set_ylabel('B$_{x}$ ($nT$)')
        bx.set_xlim(timelim)
        bx.xaxis.set_major_locator(Mtick)
        bx.xaxis.set_minor_locator(mtick)
        bx.xaxis.set_major_formatter(fmt)
        
        #By
        by = plt.subplot(4,1,3)
        by.grid()
        by.plot(self.time, self.data['B_xyz'][:,1], style)
        by.set_ylabel('$B_{y} (nT)$')
        by.set_xlim(timelim)
        by.xaxis.set_major_locator(Mtick)
        by.xaxis.set_minor_locator(mtick)
        by.xaxis.set_major_formatter(fmt)
        
        # Bz.
        newfmt = FuncFormatter(self.orbit_formatter)
        bz = plt.subplot(4,1,4)
        bz.grid()
        bz.plot(self.time, self.data['B_xyz'][:,2], style)
        bz.set_ylabel('$B_{z} (nT)$')
        bz.set_xlim(timelim)
        bz.xaxis.set_major_locator(Mtick)
        bz.xaxis.set_minor_locator(mtick)
        bz.xaxis.set_major_formatter(newfmt)
        bz.set_xlabel('Universal Time from %s' % self.time[0].isoformat())
        
        #fig.show()

        return fig
        
    def plot_bext(self, timelim=None, style='b-', DoMask=False):
        '''
        Plot satellite location and external magnetic field onto a plot 
        that is the size of a "Letter" piece of paper.
        
        "External" magnetic field is B_total - B_dipole.

        Syntax:
        -C{ramsatobject.plot_b([style=(-I{Style String}])}
        where -I{Style String} is a MatPlotLib style string that
        would be accepted by a plot command.  This provides the user
        with custom line formats for the magnetic field plots.  
        Default is 'b-', which is a blue line.
        
        
        Example usage:
            >>> import rampy
            >>> file = rampy.ramsat('somesatfile.ncdf')
            >>> file.plot_b('k--')
        
        This produces the plot on screen with black dashed lines
        in the magnetic field plots.

        DoMask (default is False) can be set to mask values where the 
        satellite is within R=1.5RE.  These are regions where poor 
        resolution often creates bad results.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter

        # Set default time limit if none given.
        if not timelim: timelim = [self.time[0], self.time[-1]]
        
        # Check style keyword for correct type.
        # Switch to default if incorrect.
        if type(style) != type('b-'):
            style='b-'

        if DoMask:
            x = self.data['SM_xyz'][:,0]
            y = self.data['SM_xyz'][:,1]
            z = self.data['SM_xyz'][:,2]
            R = np.sqrt(y**2 + x**2 + z**2)
            Bx = np.ma.masked_where(R<2.5, self.data['Bext_xyz'][:,0])
            By = np.ma.masked_where(R<1.7, self.data['Bext_xyz'][:,1])
            Bz = np.ma.masked_where(R<1.7, self.data['Bext_xyz'][:,2])
        else:
            Bx = self.data['Bext_xyz'][:,0]
            By = self.data['Bext_xyz'][:,1]
            Bz = self.data['Bext_xyz'][:,2]

        # Create plot window.
        fig = plt.figure(figsize=(8.38,10.88))
        plt.subplots_adjust(hspace=0.30, wspace=0.36, 
                            top=0.93, right=0.94, bottom=0.12)

        add_orbits(self, title = 'Orbit and B$_{Total}$ - B$_{Dipole}$ for ' + \
                       self.filename[0:self.filename.rfind('.')])

        # Time tick handlers.
        Mtick, mtick, fmt = smart_timeticks(self.time)

        # Get T89 for comparison:
        try:
            t89file = '/Users/dwelling/ramstuff/RBSP_stuff/' + \
                self.filename.split('/')[-1][0:-4] + 't89'
            t89time, t89dip, t89ext = read_t89file(t89file)
        except:
            pass

        #Bx
        bx = plt.subplot(4,1,2)
        bx.grid()
        bx.plot(self.time, Bx, style)
        try:
            bx.plot(t89time, t89ext[:,0], 'k--')
        except:
            pass
        bx.set_ylabel('B$_{x}$ ($nT$)')
        bx.set_xlim(timelim)
        bx.xaxis.set_major_locator(Mtick)
        bx.xaxis.set_minor_locator(mtick)
        bx.xaxis.set_major_formatter(fmt)
        
        #By
        by = plt.subplot(4,1,3)
        by.grid()
        by.plot(self.time, By, style)
        try:
            by.plot(t89time, t89ext[:,1], 'k--')
        except:
            pass
        by.set_ylabel('$B_{y} (nT)$')
        by.set_xlim(timelim)
        by.xaxis.set_major_locator(Mtick)
        by.xaxis.set_minor_locator(mtick)
        by.xaxis.set_major_formatter(fmt)
        
        # Bz.
        newfmt = FuncFormatter(self.orbit_formatter)
        bz = plt.subplot(4,1,4)
        bz.grid()
        bz.plot(self.time, Bz, style)
        try:
            bz.plot(t89time, t89ext[:,2], 'k--')
        except:
            pass
        bz.set_ylabel('$B_{z} (nT)$')
        bz.set_xlim(timelim)
        bz.xaxis.set_major_locator(Mtick)
        bz.xaxis.set_minor_locator(mtick)
        bz.xaxis.set_major_formatter(newfmt)
        bz.set_xlabel('Universal Time from %s' % self.time[0].isoformat())
        
        #fig.show()
        
        return fig
        
    def plot_flux(self, timelim=None):
        '''
        blerg
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors  import LogNorm
        from matplotlib.ticker import (FuncFormatter, LogLocator, 
                                       LogFormatterMathtext)
        from matplotlib.dates import date2num

        # Set default time limit if none given.
        if not timelim: timelim = [self.time[0], self.time[-1]]

        fig = plt.figure(figsize=(8.38,10.88))
        plt.subplots_adjust(hspace=0.30, wspace=0.36, 
                            top=0.93, right=0.94)

        add_orbits(self, title='Flux by Energy for '+ \
                       self.filename[0:self.filename.rfind('.')],
                   timelim = timelim)

        # Flux (all pitch angles) vs energy
        self.create_omniflux()

        # Things shared by both plots.
        timelab = 'Universal Time from %s' % self.time[0].isoformat()
        energlab= 'Log$_{10}$ Energy ($KeV$)'
        fluxlab = '$cm^{-2}s^{-1}KeV^{-1}$'
        time = date2num(self.time)
        egrid= self.data['energy_grid']
        yrng = [0.5, egrid[26]]#max(egrid)]
        Mtick, mtick, fmt = smart_timeticks(timelim)
        maxFlux = max( [self.omniH.max(), self.omniO.max()] )


        fx1 = fig.add_subplot(3,1,2)
        flxH = fx1.pcolormesh(time, egrid, self.omniH[:-1,:-1].transpose(), 
                   norm=LogNorm(), vmin=1, vmax=maxFlux)
        fx1.set_yscale('log')
        # Nudge position
        pos = fx1.get_position()
        pos.x1=pos.x1+0.1; pos.y0=pos.y0+0.05; pos.y1=pos.y1+0.05
        fx1.set_position(pos)
        cbar1 = plt.colorbar(flxH, pad=0.01, shrink=.85, ticks=LogLocator(), 
                             format=LogFormatterMathtext())
        cbar1.set_label(fluxlab)
        fx1.set_xlim(timelim)
        fx1.set_ylim(yrng)
        fx1.xaxis.set_major_locator(Mtick)
        fx1.xaxis.set_minor_locator(mtick)
        fx1.xaxis.set_major_formatter(fmt)
        fx1.set_ylabel(energlab)
        fx1.set_title('Omnidirectional H$^{+}$ Flux')


        fx2 = plt.subplot(3,1,3)
        fx2.set_yscale('log')
        flxO = fx2.pcolormesh(time, egrid, self.omniO[:-1,:-1].transpose(), 
                   norm=LogNorm(), vmin=1, vmax=maxFlux)
        # Nudge position
        pos = fx2.get_position()
        pos.x1=pos.x1+0.1; pos.y0=pos.y0+0.05; pos.y1=pos.y1+0.05
        fx2.set_position(pos)
        fx2.set_ylim([0.5, 500])
        cbar2 = plt.colorbar(flxO, pad=0.01, shrink=.85, ticks=LogLocator(), 
                             format=LogFormatterMathtext())
        cbar2.set_label(fluxlab)
        fx2.set_xlim(timelim)
        fx2.set_ylim(yrng)
        fx2.xaxis.set_major_locator(Mtick)
        fx2.xaxis.set_minor_locator(mtick)
        newfmt = FuncFormatter(self.orbit_formatter)
        fx2.xaxis.set_major_formatter(newfmt)
        fx2.set_xlabel(timelab)
        fx2.set_ylabel(energlab)
        fx2.set_title('Omnidirectional O$^{+}$ Flux')

        #fig.show()
        return fig

############################################################################
class plasmaboundary(object):
    '''
    Opens an ascii-format boundary file written from 
    IM_wrapper.f90:IM_put_from_gm
    '''

    def __init__(self,file):
        import re
        import datetime as dt

        self.filename = file
        self.read()
        # Extract time from file name.
        m = re.search('_d(\d{4})(\d{2})(\d{2})_t(\d{2})(\d{2})(\d{2})',
                         self.filename)
        self.time = dt.datetime(int(m.group(1)), #year
                                int(m.group(2)), #month
                                int(m.group(3)), #day
                                int(m.group(4)), #hour
                                int(m.group(5)), #min
                                int(m.group(6))) #sec

    def read(self):
        '''
        Read and parse the ascii boundary file.
        '''
        # Read contents, slurp.
        infile = open(self.filename, 'r')
        raw = infile.readlines()
        infile.close()

        # Throw away header.
        raw.pop(0)
        raw.pop(0)
        raw.pop(0)

        # Parse remaining data into a dictionary:
        self.npoints = len(raw)
        self.namevar = ['lt', 'RhoH', 'RhoHe', 'RhoO', 'tH', 'tHe', 'tO']
        self.data = {}
        
        for name in self.namevar:
            self.data[name] = np.zeros(self.npoints)

        for i, line in enumerate(raw):
            vals = line.split()
            for j, name in enumerate(self.namevar):
                if vals[j] == '*************':
                    vals[j] = 100000.0
                self.data[name][i] = float(vals[j])

    def write(self):
        '''
        Write self to self.filename.
        '''
        # Open file:
        out = open(self.filename, 'w')
        # Write header:
        out.write(' Time tag forthcoming.\n')
        out.write(' Units are Hours, cm-3, and eV\n')
        out.write(' lt RhoH RhoHe RhoO pH pHe pO\n')
        # Write data and close.
        for i in range(self.npoints):
            out.write('%02i %13.7f %13.7f %13.7f %13.7f %13.7f %13.7f\n' %
                      (self.data['lt'][i],
                       self.data['RhoH'][i],
                       self.data['RhoHe'][i],
                       self.data['RhoO'][i],
                       self.data['tH'][i],
                       self.data['tHe'][i],
                       self.data['tO'][i]))
        out.close()

############################################################################
def read_ram_dst(infile):
    '''
    A function to quickly read a RAM-SCB dst file as prepared during
    the post-processing step of a simulation.
    It currently returns a list of datetime objects and a numpy array
    of dst values (species-specific dst values are discarded.)
    Example: (time, dst) = read_ram_dst('dst_d20050831_t090000.txt')

    This function is a candidate for removal as more efficient ways
    of saving Dst from RAM_SCB are devised.
    '''

    import datetime as dt

    f = open(infile, 'r')
    lines = f.readlines()
    f.close()

    lines.pop(0) # Remove header.
    npoints = len(lines)
    time = []
    dst = np.zeros(npoints)

    for i, line in enumerate(lines):
        parts = line.split()
        time.append(dt.datetime(int(parts[0][1:5]),   #year
                                int(parts[0][5:7]),   #month
                                int(parts[0][7:9]),   #day
                                int(parts[0][11:13]), #hour
                                int(parts[0][13:15]), #min
                                int(parts[0][15:17])) #sec
                    )
        dst[i] = float(parts[-1]) * 1.3

    return(time, dst)
