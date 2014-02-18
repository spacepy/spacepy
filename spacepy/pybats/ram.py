'''
A module for reading, handling, and plotting RAM-SCB output.  
'''

# Global imports
import datetime as dt

import numpy as np
from scipy.io import netcdf

from spacepy.pybats import PbData
from spacepy.datamodel import dmarray

############################################################################
#  A few useful functions:
def gen_rgrid(nR=20, rmin=1.75, rmax=6.5):
    '''
    Given the number of radial points, nR, the min and max of the radial
    extent, rmin and rmax, repsectively, create and return the radial grid
    as a vector of points in Earth radii.

    Note that the first and last cells are ghost cells.
    '''
    r=np.zeros(nR+1)
    dR=(rmax-rmin)/(nR-1)
    for i in range(nR+1):
        r[i]=rmin+i*dR
    return r

def gen_tgrid(nT=25):
    '''
    Given the number of local time grid points, return, as vectors, the grid
    points in local time in units of radians and then local time hours.
    '''
    lt=np.zeros(nT); dT=24.0/(nT-1)
    phi=np.zeros(nT);dp=np.pi/(nT-1) * 2.0
    for i in range(nT):
        lt[i]= dT*i
        phi[i]=dp*i
    return phi, lt

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

    Usage: 
    (ecentr, ebound, ewidth)=gen_egrid(nE=36, lb=0.1, ew=3E-2, power=1.27)
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

def young_comp(kp, f107):
    '''
    FracH, FracHe, FracO = young_comp(kp, f107)

    Determine plasma sheet composition using the Young et al. empirical 
    relationship based on Kp and F10.7.
    Young et al. (JGR, 1982, Vol. 87 No. A11)

    Returns fraction of total number density that is Hydrogen, 
    Helium, and Oxygen.
    '''
    from numpy import exp

    ratOH = 4.5E-2 * exp(0.17*kp + 0.01*f107) # Eq. 5, pg. 9088
    ratHeH= 0.618182*ratOH*exp(-0.24*kp - 0.011*f107) + 0.011*ratOH
    fracH = 1.0 / (1.0 + ratHeH + ratOH)
    fracHe= ratHeH * fracH
    fracO = ratOH  * fracH

    return fracH, fracHe, fracO

def viz_young_comp():
    '''
    Create a quick-look plot of % O+ as calculated by Young et al. for
    various values of Kp and F10.7.
    '''
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    import numpy as np
    
    kp = np.arange(0, 9.1, 0.1)
    f7 = np.arange(0, 301, 1.0)
    fracH =np.zeros( (kp.size, f7.size) )
    fracO =np.zeros( (kp.size, f7.size) )
    fracHe=np.zeros( (kp.size, f7.size) )

    for i,f in enumerate(f7):
        fracH[:,i], fracHe[:,i], fracO[:,i] = young_comp(kp, f)
        

    levs=np.linspace(0, 0.8, 50)
    fig=plt.figure()
    ax2=fig.add_subplot(111)
    cnt2=ax2.contourf(f7,kp, fracO, levs)
    ax2.set_title('Young et al. Fraction O$^{+}$')
    ax2.set_xlabel('F10.7 Flux')
    ax2.set_ylabel('Kp Index')
    cb2=plt.colorbar(cnt2, ticks=MultipleLocator(0.1))
    cb2.set_label('Fraction O$^{+}$')
    ax2.grid()

    if plt.isinteractive():
        plt.draw()
    else:
        plt.show()


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
    axis.axvline(0, ls='-', color='k')
    axis.axhline(0, ls='-', color='k')
    
def set_orb_ticks(axis):
    '''
    Set major ticks to multiples of 1, minor ticks to 1/4.
    '''
    
    from matplotlib.ticker import MultipleLocator

    # Tick Locators:
    xMticks = MultipleLocator(2.00)
    xmticks = MultipleLocator(0.5)
    yMticks = MultipleLocator(2.00)
    ymticks = MultipleLocator(0.5)
    axis.xaxis.set_major_locator(xMticks)
    axis.xaxis.set_minor_locator(xmticks)
    axis.yaxis.set_major_locator(yMticks)
    axis.yaxis.set_minor_locator(ymticks)
    axis.grid(True)#, ls='k:')

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
        mtick=mdt.MinuteLocator(byminute=list(range(60)), interval=5)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick=mdt.MinuteLocator(byminute=[0,30])
        mtick=mdt.MinuteLocator(byminute=list(range(60)), interval=10)
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick=mdt.HourLocator(byhour=list(range(24)), interval=2)
        mtick=mdt.MinuteLocator(byminute=[0,15,30,45])
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 24:
        Mtick=mdt.HourLocator(byhour=[0,3,6,9,12,15,18,21])
        mtick=mdt.HourLocator(byhour=list(range(24)))
        fmt = mdt.DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick=mdt.HourLocator(byhour=[0,6,12,18])
        mtick=mdt.HourLocator(byhour=list(range(24)))
        fmt = mdt.DateFormatter('%H:%M UT')
    elif deltaT.days < 15:
        Mtick=mdt.DayLocator(bymonthday=list(range(1,32)))
        mtick=mdt.HourLocator(byhour=[0,6,12,18])
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 32:
        Mtick=mdt.DayLocator(bymonthday=list(range(5,35,5)))
        mtick=mdt.HourLocator(byhour=[0,6,12,18])
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 60:
        Mtick=mdt.MonthLocator()
        mtick=mdt.DayLocator(bymonthday=list(range(5,35,5)))
        fmt = mdt.DateFormatter('%b %d')
    elif deltaT.days < 731:
        Mtick=mdt.MonthLocator()
        mtick=mdt.DayLocator(bymonthday=15)
        fmt = mdt.DateFormatter('%b %Y')
    else:
        Mtick=mdt.YearLocator()
        mtick=mdt.MonthLocator(bymonth=7)
        fmt = mdt.DateFormatter('%Y')
    return(Mtick, mtick, fmt)

def apply_smart_timeticks(ax, time, dolimit=True, dolabel=False):
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
    if dolabel:
        ax.set_xlabel('Time from %s' % time[0].isoformat())
    return True

def add_body(ax, rad=2.0, rotate=0.0, add_night=True,**extra_kwargs):
    '''
    Creates a circle of radius kwarg rad (default 2.0 RE) at the center of 
    axis ax.  Then, a planet of radius 1.0 RE is added to the center to 
    represent the Earth.  Extra kwargs are passed to the patch generators.

    Three patches are created by this function: the inner boundary (a light
    grey Circle patch), the planet (a white Circle patch) and the night side
    of the planet (a black wedge.)  All three are returned to the caller.

    The kwarge "rotate", which defaults to 0.0, can be set to any angle
    (in degrees) to rotate the location of local noon/midnight on the 
    planet patches.  The default is local noon points towards +X.  A
    rotate value of 90 will cause local noon to point towards +Y.

    Example
    >>>ax=pylab.subplot(111)
    >>>inner, planet, night = rampy.add_body(ax)

    If you are using a polar axis, try using add_body_polar.
    '''

    from matplotlib.patches import Circle, Wedge

    inner = Circle((0,0), rad, fc='lightgrey', ec='k', zorder=1000,
                   **extra_kwargs)
    planet= Circle((0,0), 1.0, fc='w', zorder=1001, **extra_kwargs)
    night = Wedge( (0,0), 1.0, 90+rotate, -90+rotate, fc='k', 
                   zorder=1002, **extra_kwargs)

    ax.add_artist(inner)
    ax.add_artist(planet)
    if add_night:
        ax.add_artist(night)
    
    return inner, planet, night

def add_body_polar(ax, noon=90.0):
    '''
    Add the "Earth" to the center of a polar dial plot located on 
    axis "ax".  Kwarg noon specifies the angle counter clockwise from the
    x-axis, in degrees, that is local noon (e.g. sun facing.)  Default is
    90 degrees, or in the +Y direction.
    '''
    from numpy import pi, zeros, linspace, array
    from matplotlib.patches import Polygon

    # Number of points:
    nPts = 101

    # Convert "noon" to radians.
    noon = pi*noon/180.0

    # First, a circle.
    x=linspace(0., 2.0*pi, nPts) # Polar angle.
    y=zeros(nPts) + 1.0          # Polar radius.
    circ = Polygon( array([x,y]).transpose(), fc='w', ec='k', zorder=1001)

    # Now, nightside.
    x=linspace(noon+pi/2.0, noon+3.0*pi/2.0, nPts)
    arc = Polygon( array([x,y]).transpose(), fc='k', ec='k', zorder=1002)

    # Add to plot.
    ax.add_artist(circ)
    ax.add_artist(arc)

def _adjust_dialplot(ax, rad, title='12',labelsize=15):
    '''
    Ram output is often visualized with equatorial dial plots.  This
    function quickly adjusts those plots to give a uniform, clean
    look and feel.
    '''
    from numpy import max, pi
    from matplotlib.ticker import MultipleLocator
    # Constrain range of plot:
    ax.set_ylim([0,max(rad)])
    # Set MLT labels:
    lt_labels = ['06', title, '18', '00']
    xticks    = [   0,   pi/2,   pi, 3*pi/2]
    ax.set_xticks(xticks)
    ax.set_xticklabels(lt_labels)
    ax.tick_params('x', labelsize=labelsize)
    # Set L labels and grid.  Turn off label at L=0.
    ax.yaxis.set_major_locator(MultipleLocator(2))
    labels=ax.get_yticklabels()
    labels[0].set_visible(False)
    ax.tick_params('y', labelcolor='w', labelsize=labelsize)
    ax.grid(True, c='w', lw=1.5, ls=':')
    # Change background color so labels stand out.
    ax.set_axis_bgcolor('gray')
    add_body_polar(ax)

def get_iono_cb(ct_name='bwr'):
    '''
    Several custom colorbars used by RIM and AMIE have become standard when
    visualizing data from these models.  These are 'blue_white_red' and 
    'white_red', used for data that have positive and negative values and
    for data that have only positive values, respectively.  This function
    builds and returns these colorbars when called with the initials of the
    color table name as the only argument.

    For RAM-SCB usage, the 'bwr' color map is useful for plotting ionospheric
    potential and will create plots that look similar to those produced by
    AMIE and RIM.  The 'wr' color map is useful when plotting anisotropy.

    Example - get each colormap:
    
    >>>bwr_map = get_iono_cb('bwr')
    >>>wr_map  = get_iono_cb('wr')
    '''
    from matplotlib.colors import LinearSegmentedColormap as lsc

    if ct_name=='bwr':
        table = {
            'red':  [(0.,0.,.0),(.34,0.,0.),(.5,1.,1.),(1.,1.,1.)],
            'green':[(0.,0.,0.),(.35,1.,1.),(.66,1.,1.),(1.,0.,0.)],
            'blue' :[(0.,1.,1.),(.5,1.,1.),(.66,0.,0.),(.85,0.,0.),(1.,.1,.1)]
            }
        cmap = lsc('blue_white_red',table)
    elif ct_name=='wr':
        table = {
            'red':  [(0.,1.,1.),(1.,1.,1.)],
            'green':[(0.,1.,1.),(1.,0.,0.)],
            'blue' :[(0.,1.,1.),(1.,.0,.0)]
            }
        cmap = lsc('white_red',table)

    return cmap

############################################################################
class weimer(object):
    '''
    Read a Weimer electric field file.
    '''
    
    def __init__(self,filename):
        self.filename = filename
        self.read()

    def _parse_head(self, headlines):
        parts= headlines.split()
        self.doy = int(parts[0])
        self.hour= float(parts[1])
        self.time = dt.datetime(year, 1,1,0,0) + \
            dt.timedelta(days=doy-1) + \
            dt.timedelta(hours=hour)

    def read(self):
        import datetime as dt

        infile = open(self.filename, 'r')
        
        # Get year, time:
        head = infile.readline()
        year = int(head.split()[-1])
        
        self._parse_head(infile.readline())

        # skip rest of crap header.
        head = infile.readline()
        head = infile.readline()

        # Read data
        raw = infile.readlines()
        self.epot = np.zeros(len(raw)-1)
        self.l   = np.zeros(len(raw)-1)
        self.lt  = np.zeros(len(raw)-1)

        for i, line in enumerate(raw[0:-1]):
            parts = line.split()
            self.l[i]    = float(parts[0])
            self.lt[i]   = float(parts[1])
            self.epot[i] = float(parts[2])

        # Make polar angle:
        self.phi=self.lt*np.pi/12.0
        # Make CPCP
        self.cpcp = self.epot.max() - self.epot.min()

        infile.close()


    def add_potplot(self, target=None, loc=111, zlim=50, n=31, add_cbar=True):
        '''
        Quickly add a potential plot to MPL object "target".
        '''
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.ticker import MultipleLocator
        from numpy import linspace

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc, polar=True)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(5,4))
            ax = fig.add_subplot(loc, polar=True)

        cmap=get_iono_cb('bwr')
        crange = Normalize(vmin=-1.*zlim, vmax=zlim)
        levs = linspace(-1*zlim, zlim, n)

        cont=ax.tricontourf(self.phi+np.pi/2.0, self.l, self.epot, levs, 
                            norm=crange, cmap=cmap)
        _adjust_dialplot(ax,self.l)
        cbar=False
        if add_cbar:
            cticks = MultipleLocator(25)
            cbar = plt.colorbar(cont, ticks=cticks, shrink=0.8, pad=0.08)
            cbar.set_label('kV')

        return fig, ax, cont, cbar

class WeqFile(weimer):
    '''
    Slight variation on weimer class to read weq_***.in files.
    '''
    def _parse_head(self, headlines):
        parts=headlines.split()
        self.UT=float(parts[0])

############################################################################
class RamSat(object):
    '''
    A class to handle and plot RAM-SCB virtual satellite files.  Instantiate
    an object by simply calling:
    >>>sat=rampy.RamSat('SatFile.nc')
    Simple satellite file info, such as the time vector, file name, and
    output frequency are found as direct object attributes:
    >>>sat.time
         array([2005-08-31 09:29:00, ...,2005-09-02 00:58:00],dtype=object)
    Data arrays are stored as object keys.  For example, obtain satellite 
    coordinates by typing:
    >>>sat['SM_xyz']
    RamSats work, for the most part, like dictionaries.  Add a new data entry
    (an array, constant, etc.) by simply calling:
    >>>sat['New Thing']=(0,1,2,3,4)
    Be sure to peruse the docstrings of the many object methods to discover
    the many plotting functions and canned plotting routines.
    '''

    # Make'em work like a dictionary:
    def __getitem__(self, key):
        if key in self.filedata:
            return self.filedata[key][...]
        elif key in self.data:
            return self.data[key]
        else:
            raise KeyError(key)
        
    def __setitem__(self, key, value):
        self.data[key] = value
        if not key in self.namevars:
            self.namevars.append(key)

    def keys(self):
        '''
        List all keys for data objects stored in the RamSat object.
        '''
        return list(self.filedata.keys()) + list(self.data.keys())

    def __contains__(self, key):
        return key in self.filedata or key in self.data

    def __init__(self, filename):
        # Filename:
        self.filename=filename
        self._read()

    def close(self):
        self.f.close()

    def __del__(self):
        self.close()

    def _read(self):
        '''
        Load the NCDF file.  Should only be called by __init__().
        '''
        self.f = netcdf.netcdf_file(self.filename, mode='r', mmap=False)
        self.namevars = self.f.variables.keys()
        self.attrs = {}
        #split off the netCDF attributes from the Python attributes
        for k in dir(self.f):
            if k[0] == '_' or k in ('dimensions', 'filename', 'fp', 'mode',
                                    'use_mmap', 'variables', 'version_byte'):
                continue
            self.attrs[k] = getattr(self.f, k)
        # self.filedata contains the raw netcdf_variable objects.
        self.filedata = self.f.variables
        # New values saved as self[key] are stored in self.data, not
        # self.filedata which cannot be changed.
        self.data={}

        # Get start time, set default if not found.
        try:
            stringtime = self.attrs['start_time']
        except KeyError:
            stringtime = '20000101000000'

        # parse start time and use it to create time object.
        self.starttime = dt.datetime.strptime(stringtime, '%Y%m%d%H%M%S')
        
        # Convert time in seconds to datetime.
        secs = self.f.variables['Time'][...]
        self.time = np.array(
            [self.starttime + dt.timedelta(seconds=float(s)) for s in secs])

        # Some other variables to keep handy:
        try:
            self.dt = self['DtWrite']
        except KeyError:
            # Old files do not have DtWrite.
            self.dt = np.diff(time).min()

    def create_omniflux(self):
        '''
        Integrate the flux(currently as a function of energy and 
        pitch angle) to get omnidirectional flux as a function of
        energy.  New units = (cm-2*s*keV)^-1

        Usage: just call this method to integrate all species.
        '''
        
        # Create new flux attributes:
        nTime = self.time.size
        nPa   = self['pa_grid'].size
        nEner = self['energy_grid'].size
        self['omniH'] = np.zeros((nTime, nEner))
        self['omniHe']= np.zeros((nTime, nEner))
        self['omniO'] = np.zeros((nTime, nEner))
        self['omnie'] = np.zeros((nTime, nEner))
        # Create delta mu, where mu = cos(pitch angle)
        dMu = np.zeros(nPa)
        if 'pa_width' in self:
            dMu=4*np.pi*self['pa_width']
        else:
            dMu[0] = self['pa_grid'][1]
            for i in range(1,nPa): 
                # Factor of pi here so we don't need it later.
                dMu[i] = 4*np.pi*self['pa_grid'][i] -self['pa_grid'][i-1]

        # Integrate.
        tempH=self['FluxH+']; tempO=self['FluxO+']
        tempe=self['Fluxe-']; tempHe=self['FluxHe+']
        for i in range(1, nPa):
            tempH[:,i,:] =tempH[:,i,:] * dMu[i]
            tempe[:,i,:] =tempe[:,i,:] * dMu[i]
            tempO[:,i,:] =tempO[:,i,:] * dMu[i]
            tempHe[:,i,:]=tempHe[:,i,:]* dMu[i]

        self.data['omniH']=tempH.sum(1)
        self.data['omnie']=tempe.sum(1)
        self.data['omniO']=tempO.sum(1)
        self.data['omniHe']=tempHe.sum(1)

        # Mask out bad values.
        self.data['omniH'] =np.ma.masked_where(self.data['omniH'] 
                                               <=0,self.data['omniH'])
        self.data['omniHe']=np.ma.masked_where(self.data['omniHe']
                                               <=0,self.data['omniHe'])
        self.data['omniO'] =np.ma.masked_where(self.data['omniO'] 
                                               <=0,self.data['omniO'])
        self.data['omnie'] =np.ma.masked_where(self.data['omnie'] 
                                               <=0,self.data['omnie'])

    #####RamSat Viz Routines#####
    def _orbit_formatter(self, x, pos):
        '''
        A function that, when passed to the FuncFormatter class of 
        Matplotlib tick formatters, produces specialized ticks that give
        UT, LT, and inclination.
        '''
        
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
        x = self['SM_xyz'][index,0]
        y = self['SM_xyz'][index,1]
        z = self['SM_xyz'][index,2]
        R = np.sqrt(y**2 + x**2 + z**2)
        Mlat  = 180.0*np.arcsin(z/R)/np.pi
        theta = 180.0*np.sign(y)*np.arccos(x/np.sqrt(y**2+x**2))/np.pi
        theta = np.mod(theta, 360.0)
        locHr = np.floor(theta/15.0)
        locMn = int(60.0 * (theta/15.0 - locHr))

        # Build format string.
        fmtstring = \
            '%02d:%02d UT\n%02d:%02d MLT\n%5.1f$^{\circ}$ MLat\nR=%4.2f $R_{E}$'\
            % (nowtime.hour, nowtime.minute, locHr, locMn, Mlat, R)

        return (fmtstring)

    def add_orbit_plot(self, plane='XY', target=None, timelim=False, loc=111, 
                       ls='g.', title=False, invertX=True):
        """
        Add a simple, 2D plot of the satellite orbit in a given (SM) plane.

        Other Parameters
        ================
        plane : string
             Plane to plot. Options are 'XY', 'XZ', and 'YZ'; defaults 'XY'.

        target : Figure or Axes
             If None (default), a new figure is generated from scratch.
             If a matplotlib Figure object, a new axis is created
             to fill that figure.
             If a matplotlib Axes object, the plot is placed
             into that axis.

        loc : int
            Use to specify the subplot placement of the axis
            (e.g. loc=212, etc.) Used if target is a Figure or None.
            Default 111 (single plot).

        invertX : boolean
            Reverse the SM X axis so "Sun is to the left." Default True.
            This reverses the SM axis, not the plot axis, so it has no effect
            if plane is YZ.
        """
        import matplotlib.pyplot as plt
        
        if not plane.upper() in ('XY','XZ','YZ'):
            raise ValueError("{0} is not a valid plot plane.".format(plane))

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(5,5))
            ax = fig.add_subplot(loc)

        # Variables to map plot plane to correct variables:
        plane = plane.upper()
        ijk = {'X':0, 'Y':1, 'Z':2}
        i = ijk[plane[0]]
        j = ijk[plane[1]]

        if not timelim: 
        # Set default time limit if none given.
            timelim = [self.time[0], self.time[-1]]
            iMin=0
            iMax=-1
        else:
            # Use timelim to get indices that bound our plots.
            timediff = abs(self.time - timelim[-1])
            iMax = np.nonzero(timediff == timediff.min())[0][0]
            timediff = abs(self.time - timelim[0])
            iMin = np.nonzero(timediff == timediff.min())[0][0]
            
        # Add orbit:
        ax.plot(self['SM_xyz'][iMin:iMax,i], self['SM_xyz'][iMin:iMax,j],ls)
        # Add body:
        add_body(ax,add_night=(plane!='YZ'))

        # Axis details:
        ax.axis('equal')
        if plane.upper() in ('XY','XZ') and invertX:
            xmin, xmax = ax.get_xlim()
            if xmin < xmax:
                ax.invert_xaxis()
        ax.set_xlabel('SM %s'%(plane[0]))
        ax.set_ylabel('SM %s'%(plane[1]))
        if title:
            ax.set_title(title)
        grid_zeros(ax)
        set_orb_ticks(ax)

        return fig, ax

    def add_omniflux_plot(self, nameflux, target=None, zlim=[1E4,1E9],
                          add_cbar=True, do_orbticks=False, title=False,
                          timelim=False, loc=111, no_xlabels=False):
        """
        Create and place a pcolor-type plot of omnidirectional flux.

        Omnidirectional fluxes are calculated by the object method
        "create_omniflux" and are saved into the object with keys such as
        'omniO' and 'omniHe', etc.  

        Other Parameters
        ================


        target : Figure or Axes
             If None (default), a new figure is generated from scratch.
             If a matplotlib Figure object, a new axis is created
             to fill that figure.
             If a matplotlib Axes object, the plot is placed
             into that axis.

        loc : int
            Use to specify the subplot placement of the axis
            (e.g. loc=212, etc.) Used if target is a Figure or None.
            Default 111 (single plot).

        do_orbticks : boolean
            Activate orbit ticks (default False): add complex but
            informational tick marks that list the satellite's coordinates
            along with the time on the x-axis.  However, using these will
            require some "massaging" of the plot to make them properly
            visible.  Be sure to add lots of padding to the bottom of the plot.

        no_xlabels : boolean
            (Default False) Completely omit xlabel and xticklabels.
            Convenient for stacking flux plots on each other.
        """

        import matplotlib.pyplot as plt
        from matplotlib.colors  import LogNorm
        from matplotlib.ticker import (FuncFormatter, LogLocator, 
                                       LogFormatterMathtext)
        from matplotlib.dates import date2num

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(loc)

        # Check for omni fluxes, calculate as necessary.
        if not nameflux in self:
            self.create_omniflux()
            if not nameflux in self:
                raise KeyError('%s is not a valid omnidirectional flux.' 
                               % nameflux)
        # Create a time vector that binds each pixel correctly.
        time=np.zeros(self.time.size+1)
        time[0]=date2num(self.time[0]-dt.timedelta(seconds=self.dt/2.0))
        time[1:]=date2num(self.time+dt.timedelta(seconds=self.dt/2.0))
        #egrid=self['energy_grid']
        ecenter, eboundary, ewidth=gen_egrid(nE=self['energy_grid'].size)
#        print("Need better energy grid setup for pcolormesh.")
        flx=ax.pcolormesh(time,eboundary,self[nameflux].transpose(),
                          norm=LogNorm(),vmin=zlim[0],vmax=zlim[1])
        ax.set_yscale('log')
        ax.set_ylim( [eboundary[0],eboundary[-1]] )
        if not timelim:
            timelim=[self.time[0],self.time[-1]]
        apply_smart_timeticks(ax,timelim,dolabel=True)
        if no_xlabels:
            ax.set_xlabel('')
            ax.set_xticklabels([''])
            do_orbticks=False
        ax.set_ylabel('E ($keV$)')
        if title:  #If title not set, use a default:
            ax.set_title(title)
        else:
            labels={'omniH':'H$^{+}$','omniHe':'He$^{+}$',
                    'omniO':'O$^{+}$','omnie':'e$^{-}$'}
            ax.set_title('Omnidirectional %s Flux' % (labels[nameflux]))
        if do_orbticks:
            ax.xaxis.set_major_formatter(FuncFormatter(self._orbit_formatter))
        if add_cbar:
            cbar = plt.colorbar(flx, pad=0.01, shrink=.85, ticks=LogLocator(), 
                                format=LogFormatterMathtext(), ax=ax)
            cbar.set_label('$cm^{-2}s^{-1}keV^{-1}$')
        else:
            cbar=False

        return fig, ax, flx, cbar

    def plot_omni_quicklook(self, flux_opts=None, eflux_opts=None,
                            hflux_opts=None, oflux_opts=None):
        """
        Create a quick-look plot of omnidirectional fluxes.

        Other Parameters
        ================
        flux_opts : dict
            dictionary of keyword arguments to pass to
            :meth:`add_omniflux_plot` for all flux plots

        eflux_opts : dict
            as flux_opts, but for the electron flux plot only

        hflux_opts : dict
            as flux_opts, but for the H+ flux plot only

        oflux_opts : dict
            as flux_opts, but for the O+ flux plot only

        Returns
        =======
        out : Figure
            Has 9 axes instances (axes attribute)
            0: XY orbit plot
            1: XZ orbit plot
            2: YZ orbit plot
            3: e- flux
            4: H+ flux
            5: O+ flux
            6: e- flux color bar
            7: H+ flux color bar
            8: O+ flux color bar
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        
        fig=plt.figure(figsize=(11,7))
        fig.subplots_adjust(left=0.07, right=0.99, bottom=0.19, 
                            top=0.94, wspace=0.4, hspace=0.25)
        gs=gridspec.GridSpec(3,3)

        # Do orbits first.
        a1=fig.add_subplot(gs[0,0])
        a2=fig.add_subplot(gs[1,0])
        a3=fig.add_subplot(gs[2,0])
        self.add_orbit_plot('XY', target=a1)
        self.add_orbit_plot('XZ', target=a2)
        self.add_orbit_plot('YZ', target=a3)

        # Add fluxes.
        a1=fig.add_subplot(gs[0,1:])
        a2=fig.add_subplot(gs[1,1:])
        a3=fig.add_subplot(gs[2,1:])
        if eflux_opts is None:
            eflux_opts = {}
        if hflux_opts is None:
            hflux_opts = {}
        if oflux_opts is None:
            oflux_opts = {}
        if flux_opts is None:
            flux_opts = {}
        for k in flux_opts:
            for d in (eflux_opts, hflux_opts, oflux_opts):
                if not k in d:
                    d[k] = flux_opts[k]
        self.add_omniflux_plot('omnie', target=a1, no_xlabels=True,
                               **eflux_opts)
        self.add_omniflux_plot('omniH', target=a2, no_xlabels=True,
                               **hflux_opts)
        self.add_omniflux_plot('omniO', target=a3, do_orbticks=True,
                               **oflux_opts)
        
        return fig

############################################################################
class PlasmaBoundary(PbData):
    '''
    Opens an ascii-format boundary file written from 
    IM_wrapper.f90:IM_put_from_gm
    '''

    def __init__(self,filename, time=None, *args, **kwargs):
        
        # Init base object.
        super(PlasmaBoundary, self).__init__(*args, **kwargs)

        # Set up this instance:
        self.attrs['file']=filename
        self.read(time)

    def read(self, time=None):
        '''
        Read and parse the ascii boundary file.
        '''
        import datetime as dt

        # Read contents, slurp.
        infile = open(self.attrs['file'], 'r')
        raw = infile.readlines()
        infile.close()

        # First line of header has date and time.
        if time:
            self.attrs['time']=time
            raw.pop(0)
        else:
            parts = raw.pop(0).split()
            self.attrs['elapsed']=float(parts[0])
            self.attrs['time'] = dt.datetime(int(parts[1]), #year
                                    int(parts[2]), #month
                                    int(parts[3]), #day
                                    int(parts[4]), #hour
                                    int(parts[5]), #min
                                    int(parts[6])) #sec
        # Save header info that is useful.
        raw.pop(0)
        self.attrs['namevar']=raw.pop(0).split()
        self.attrs['npoints']=len(raw)
        
        # Parse remaining data into a dictionary:
        units = {'l':'hours', 'R':'cm-3', 'p':'eV'}
        for name in self.attrs['namevar']:
            self[name] = dmarray(np.zeros(self.attrs['npoints']), 
                                 attrs={'units':units[name[0]]})

        for i, line in enumerate(raw):
            vals = line.split()
            for j, name in enumerate(self.attrs['namevar']):
                if vals[j] == '*************':
                    vals[j] = 100000.0
                self[name][i] = float(vals[j])

    def write(self):
        '''
        Write self to self.filename.
        '''
        # Open file:
        out = open(self.attrs['file'], 'w')
        # Write header:
        out.write('%10.5E %4i %2i %2i %2i %2i %2i\n'%
                  (self.attrs['elapsed'], self.attrs['time'].year,
                   self.attrs['time'].month, self.attrs['time'].day,
                   self.attrs['time'].hour, self.attrs['time'].minute,
                   self.attrs['time'].second)
                  )
        out.write(' Units are Hours, cm-3, and eV\n')
        out.write(' '+' '.join(self.attrs['namevar'])+'\n')
        # Write data and close.
        for i in range(self.attrs['npoints']):
            out.write('%02i %13.6E %13.6E %13.6E %13.6E %13.6E %13.6E\n' %
                      (self[self.attrs['namevar'][0]][i],
                       self[self.attrs['namevar'][1]][i],
                       self[self.attrs['namevar'][2]][i],
                       self[self.attrs['namevar'][3]][i],
                       self[self.attrs['namevar'][4]][i],
                       self[self.attrs['namevar'][5]][i],
                       self[self.attrs['namevar'][6]][i]))
        out.close()


############################################################################
class BoundaryGroup(PbData):
    '''
    A class that collects many PlasmaBoundary objects together to work as
    a coherent group.
    '''

    def __init__(self, path='.', rotate=True, *args, **kwargs):
        from glob import glob
        from matplotlib.dates import date2num

        # Init base object.
        super(BoundaryGroup, self).__init__(*args, **kwargs)
        
        # Load files and count them!
        files = glob(path+'/bound_plasma*.out')
        nfiles=len(files)
        if nfiles==0:
            raise ValueError('No files found in path='+indir)

        # Use info from first file to set up data structures.
        temp = PlasmaBoundary(files[0])
        namevar = temp.attrs['namevar']
        self['time'] = dmarray(np.zeros(nfiles, dtype=object))

        for name in namevar:
            self[name] = dmarray(np.zeros((25,nfiles)), 
                                 attrs={'units':temp[name].attrs['units']})

        for i, f in enumerate(files):
            temp = PlasmaBoundary(f)
            for name in namevar:
                self[name][:,i] = temp[name]
            self['time'][i]=(temp.attrs['time'])
        
        # Create some new values to plot:
        self['nAll'] = dmarray(self['RhoH'] + self[namevar[2]] + self['RhoO'],
                               attrs={'units':'cm-3'})
        self['ratO']= dmarray(100.0*self['RhoO']/self['nAll'], 
                              attrs={'units':'\%'})

        if rotate:
            for val in self.keys():
                if val == 'time': continue
                temparray = self[val].copy()
                temparray[0:12,  :] = self[val][13:25, :]
                temparray[12:25, :] = self[val][0:13,  :]
                self[val][:,:]   = temparray[:,:]

        # Put pressure/temp in keV.
        for val in self.keys():
            if (val[0] == 'p') and (self[val].attrs['units']=='eV'):
                self[val] /= 1000.0
                self[val].attrs['units'] = 'keV'

        # Internals for plotting:
        self._y         = np.arange(0,26) - 0.5
        self._lt_labels = ['Dusk', 'Midnight','Dawn']
        self._yticks    = [   6.0,       12.0,  18.0]
        self._dtime     = date2num(self['time']) # create decimal time.
        self._zlims     = {'p':[0,40], 'R':[0,8], 'r':[0,60], 'n':[0,8]}

    def add_ltut(self, var, target=None, loc=111, cmap='jet', zlim=None, 
                 add_cbar=True, clabel=None, xlabel='full', title=None,
                 grid=True, ntick=5):
        '''
        Plot variable *var* as a contour against local time (y-axis) and
        universal time (x-axis) using the PyBats *target* method of other
        standard plotting methods.  Four items are returned: the Matplotlib
        Figure, Axes, Mesh, and ColorBar objects used (if add_cbar is set to
        **False**, the returned ColorBar object is simply set to **False**.)

        ========== =======================================================
        Kwarg      Description
        ---------- -------------------------------------------------------
        target     Select plot destination.  Defaults to new figure/axis.
        loc        The location of any generated subplots.  Default is 111.
        add_cbar   Toggles the automatic colorbar.  Default is**True**.
        cmap       Selects Matplotlib color table.  Defaults to *jet*.
        zlim       Limits for z-axis.  Defaults to best-guess.
        clabel     Sets colorbar label.  Defaults to *var* units.
        xlabel     Sets x-axis limits, use 'full', 'ticks', or **None**.
        title      Sets axis title; defaults to **None**.
        grid       Show white dotted grid?  Defaults to **True**
        ntick      Number of attempted cbar ticks.  Defaults to 5.
        ========== =======================================================
        
        '''

        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator
        from spacepy.pybats import apply_smart_timeticks

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax  = target
            fig = ax.figure
        else:
            fig = plt.figure()
            ax  = fig.add_subplot(loc)

        # Set up z-limits; use variable-specific defaults.
        if zlim==None:
            try:
                zlim = self._zlims[var[0]]
            except ValueError:
                print("No default zlimits for ", var)
                zlim=None
            
        # Create plot:
        mesh = ax.pcolormesh(self._dtime, self._y, self[var], 
                             cmap=plt.get_cmap(cmap), 
                             vmin=zlim[0], vmax=zlim[-1])
        # Use LT ticks and markers on y-axis:
        ax.set_yticks(self._yticks)
        ax.set_yticklabels(self._lt_labels)
        ax.set_ylim([4,20])

        # White ticks, slightly thicker:
        ax.tick_params(axis='both', which='both', color='w', width=1.2)

        # Grid marks:
        if grid: ax.grid(c='w')

        if title: ax.set_title(title)
        if xlabel == 'full':
            # Both ticks and label.
            apply_smart_timeticks(ax, self['time'], dolabel=True)
        elif xlabel == 'ticks':
            # Ticks, but no date label.
            apply_smart_timeticks(ax, self['time'], dolabel=False)
        else:
            # A blank x-axis is often useful.
            apply_smart_timeticks(ax, self['time'], dolabel=False)
            ax.set_xticklabels('')
        # Add cbar as necessary:
        if add_cbar:
            lct = MultipleLocator(np.ceil(10*(zlim[1]-zlim[0])/(ntick-1))/10.0)
            cbar=plt.colorbar(mesh, ax=ax, pad=0.01, shrink=0.85, ticks=lct)
            cbar.set_label(clabel)
        else:
            cbar=None
            
        return fig, ax, mesh, cbar

############################################################################
class PressureFile(PbData):
    '''
    A class for reading and visualizing pressure_####.in files.
    '''

    def __init__(self, filename, *args, **kwargs):
        from numpy import array, pi, round
        from dateutil.parser import parse

        # Init base object.
        super(PressureFile, self).__init__(*args, **kwargs)

        # Set up object.
        self.attrs['file']=filename

        # Read and parse file.
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()

        # Create variables:
        nRec=len(lines[2:])
        self['L']    =dmarray(np.zeros(nRec),attrs={'units':'RE'})
        self['mlt']  =dmarray(np.zeros(nRec),attrs={'units':'Hours'})
        self['perH'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['parH'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['perO'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['parO'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['perHe']=dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['parHe']=dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['pere'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['pare'] =dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})
        self['total']=dmarray(np.zeros(nRec),attrs={'units':'keV/cm-3'})

        try:
            self.attrs['time']=parse(lines[0][5:28], fuzzy=True)
        except:
            self.attrs['time']='unknown'

        for i, line in enumerate(lines[2:]): #Skip header.
            parts = line.split()
            self['L'][i]     =float(parts[0])
            self['mlt'][i]   =float(parts[1])
            self['perH'][i]  =float(parts[2])
            self['parH'][i]  =float(parts[3]) 
            self['perO'][i]  =float(parts[4])
            self['parO'][i]  =float(parts[5]) 
            self['perHe'][i] =float(parts[6])
            self['parHe'][i] =float(parts[7])
            self['pere'][i]  =float(parts[8])
            self['pare'][i]  =float(parts[9])
            self['total'][i] =float(parts[10])

        # Theta is an angle used for polar plots.
        self['theta']=self['mlt']*pi/12.0 - pi/2.0
        self['theta'].attrs={'units':'rad'}
        
        # Grid spacing/size:
        self.attrs['dTheta'] = self['theta'][1]-self['theta'][0]
        self.attrs['nTheta'] = int(np.round(2.*pi/self.attrs['dTheta'])+1)
        self.attrs['nL'] = len(self['L'])/self.attrs['nTheta']
        self.attrs['dL'] = self['L'][self.attrs['nTheta']] - self['L'][0]

        # Calculate isotropic pressures and anisotropies.
        self['totH']  = (2./3.)*self['perH'] +(1./3.)*self['parH']
        self['totHe'] = (2./3.)*self['perHe']+(1./3.)*self['parHe']
        self['totO']  = (2./3.)*self['perO'] +(1./3.)*self['parO']
        self['tote']  = (2./3.)*self['pere'] +(1./3.)*self['pare']
        self['anie']  = self['pere']  / self['pare']  - 1.0
        self['aniH']  = self['perH']  / self['parH']  - 1.0
        self['aniHe'] = self['perHe'] / self['parHe'] - 1.0
        self['aniO']  = self['perO']  / self['parO']  - 1.0

        # Isotropy units.
        self['anie'].attrs={'units':''}
        self['aniH'].attrs={'units':''}
        self['aniHe'].attrs={'units':''}
        self['aniO'].attrs={'units':''}

        # Labels for plots.
        self.labels={}
        self['pere'].attrs['label'] = r'$\bot$ Pressure, e$^{-}$'
        self['pare'].attrs['label'] = r'$\parallel$ Pressure, e$^{-}$'
        self['tote'].attrs['label'] = r'Pressure, e$^{-}$'
        self['perH'].attrs['label'] = r'$\bot$ Pressure, H$^{+}$'
        self['parH'].attrs['label'] = r'$\parallel$ Pressure, H$^{+}$'
        self['totH'].attrs['label'] = r'Pressure, H$^{+}$'
        self['aniH'].attrs['label'] = r'Anisotropy, H$^{+}$'
        self['perHe'].attrs['label']= r'$\bot$ Pressure, He$^{+}$'
        self['parHe'].attrs['label']= r'$\parallel$ Pressure, He$^{+}$'
        self['totHe'].attrs['label']= r'Pressure, He$^{+}$'
        self['aniHe'].attrs['label']= r'Anisotropy, He$^{+}$'
        self['perO'].attrs['label'] = r'$\bot$ Pressure, O$^{+}$'
        self['parO'].attrs['label'] = r'$\parallel$ Pressure, O$^{+}$'
        self['totO'].attrs['label'] = r'Pressure, O$^{+}$'
        self['aniO'].attrs['label'] = r'Anisotropy, O$^{+}$'
        self['total'].attrs['label']= r'Total Pressure'


    def add_cont_press(self, var='total', n=31, target=None, maxz=1000.0, 
                       minz=1.0, loc=111, add_cbar=False, npa=False, 
                       labelsize=15,  title='auto', **kwargs):
        '''
        Create a polar log-axis contour plot of pressure and add it to axis
        'ax'.  For speedier plots, use plot_cont_press, which makes its
        own axis.
        '''
        from numpy import linspace, power, max, pi, log10
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.cm import get_cmap
        from matplotlib.pyplot import colorbar
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)

        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc, polar=True)
        elif type(target).__base__ in (plt.Axes, plt.PolarAxes):
            ax  = target
            fig = ax.figure
        else:
            fig = plt.figure()
            ax  = fig.add_subplot(loc, polar=True)

        p=self[var]
        mapname='spectral'
        label='$KeV/cm^-3$'
        if title=='auto':
            title=self[var].attrs['label']
        if npa:
            p=p*0.16 # Energy Density to Pressure in nPa.
            mapname='jet'
            label='nPa'
        # Set up color bar & levels.
        levs = power(10, linspace(log10(minz), log10(maxz), n))
        minz=0.01
        cont = ax.tricontour(self['theta'], self['L'], p, levs, 
                             norm=LogNorm(), cmap=get_cmap(mapname))
        cont = ax.tricontourf(self['theta'], self['L'], p, levs, 
                              norm=LogNorm(), cmap=get_cmap(mapname))
        _adjust_dialplot(ax, self['L'], title=title, labelsize=labelsize)
        if add_cbar:
            cbar = colorbar(cont, pad=0.1, ticks=LogLocator(), ax=ax,
                            format=LogFormatterMathtext(), shrink=0.8)
            #cbar.set_label(self.labels[var]+' ($KeV/cm^-3)$')
            cbar.set_label(label)
            
    def add_pcol_press(self, var='total', target=None, maxz=1000.0, minz=1.0,
                       add_cbar=False, loc=111, labelsize=15, title='auto',
                       **kwargs):
        '''
        Add a pcolor plot of the pressure object to target.
        '''
        from numpy import reshape, linspace, pi
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.cm import get_cmap
        from matplotlib.pyplot import colorbar
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)
        
        # Set ax and fig based on given target.
        if type(target) == plt.Figure:
            fig = target
            ax  = fig.add_subplot(loc, polar=True)
        elif type(target).__base__ in (plt.Axes,plt.PolarAxes):
            ax  = target
            fig = ax.figure
        else:
            fig = plt.figure()
            ax  = fig.add_subplot(loc, polar=True)

        # Set title.
        if title=='auto':
            title=self[var].attrs['label']

        # Set up grid centered on gridpoints.
        dT=self.attrs['dTheta']
        dL=self.attrs['dL']
        T=linspace(-1.0*dT/2.0, 2.*pi-dT/2.0, self.attrs['nTheta'])
        T=T-pi/2.0
        R=linspace(self['L'][0]-dL/2.0,self['L'][-1]+dL/2.0,self.attrs['nL']+1)
        p=reshape(self[var], [self.attrs['nL'], self.attrs['nTheta']])
        pcol = ax.pcolormesh(T, R, p[:,:-1], norm=LogNorm(),
                             vmin=minz, vmax=maxz, cmap=get_cmap('spectral'))
        _adjust_dialplot(ax, R, title=title, labelsize=15)
        if add_cbar:
            cbar = colorbar(pcol, pad=0.1, ticks=LogLocator(), ax=ax,
                            format=LogFormatterMathtext(), shrink=0.8)
            #cbar.set_label(self.labels[var]+' ($keV/cm^{-3}$)')
            cbar.set_label('$keV/cm^{-3}$')


############################################################################
class BoundaryFluxFile(object):
    '''
    Read, plot, and edit flux (*.swf) or (*.dat) files.
    *.dat files are output into the Dsbnd directory.
    '''

    def __init__(self, filename):
        self.filename=filename
        if filename[-3:]=='swf':
            self.read_swf()
        else:
            self.read_dsbnd()

    def read_dsbnd(self):
        from datetime import datetime
        from re import search
        from numpy import zeros, linspace
        f=open(self.filename, 'r')
        lines = f.readlines()
        f.close()

        parts=lines[0].split()
        self.rtime=float(parts[-1])

        # Try to get run time, species from filename.
        match = search('ds(\w{2})\_d(\d{4})(\d{2})(\d{2})\_t(\d{2})(\d{2})(\d{2})', 
                  self.filename)
        if match:
            # New filename format! Huzzah!
            t = match.groups()
            self.species=t[0].strip('_')
            self.time = datetime(int(t[1]), int(t[2]), int(t[3]), 
                                 int(t[4]), int(t[5]), int(t[6]))
        else:
            self.time = None
            self.species=None

        # Determine size of file and allocate array.
        self.nE = len(lines) - 1
        parts = lines[1].split()
        self.nLT = len(parts)-1
        self.flux=zeros([self.nLT, self.nE])
        self.E=zeros(self.nE)

        # Read and parse lines.
        for i, line in enumerate(lines[1:]):
            parts = line.split()
            self.E[i]=float(parts[0])
            for j in range(self.nLT):
                self.flux[j,i]=float(parts[j+1])

        # Set up local time grid.
        self.LT = linspace(0, 24, self.nLT)

    def read_swf(self):
        from numpy import zeros, linspace
        f=open(self.filename, 'r')
        lines = f.readlines()
        f.close()
        
        # Determine size of file and allocate array.
        self.nLT = len(lines) - 2
        parts = lines[1].split()
        self.nE = len(parts)-2
        self.flux=zeros([self.nLT, self.nE])
        
        # Read and parse lines.
        for i, line in enumerate(lines[1:-1]):
            parts = line.split()
            for j in range(self.nE):
                self.flux[i,j]=float(parts[j+2])

        # Set up local time grid.
        self.LT = linspace(0, 24, self.nLT)

    def quickplot(self, zlim=[10,1E7]):
        '''
        A quick plot of input flux on a fresh axis.
        '''
        from numpy import array
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        egrid = array(range(self.nE))
        flux = self.flux.transpose()
        flux[flux<0.01] = 0.01
        flx = ax.pcolor(self.LT, egrid, flux, norm=LogNorm(), vmin=zlim[0],
                        vmax=zlim[-1], cmap=plt.get_cmap('gnuplot2'))
        cbar = plt.colorbar(flx, pad=0.01, shrink=0.85, ticks=LogLocator(),
                            format=LogFormatterMathtext())
        cbar.set_label('$cm^{-2}s^{-1}ster^{-1}keV^{-1}$')
        ax.set_xlim([0,24])
        ax.set_ylim([0,self.nE-1])
        if self.time:
            ax.set_title(self.time.isoformat())
        ax.set_xlabel('Local Time Sector')
        if hasattr(self, 'E'):
            ax.set_ylabel('Energy (keV)')
            newlabs=[]
            for val in ax.get_yticks()[:-1]:
                newlabs.append('%6.2f' % self.E[int(val)])
            ax.set_yticklabels(newlabs)
        else:
            ax.set_ylabel('Energy Bin')

        return fig, ax

############################################################################
class LogFile(PbData):
    def __init__(self, file, *args, **kwargs):
        super(LogFile, self).__init__(*args, **kwargs)
        self.attrs['file'] = file
        self.read()

    def read(self):
        '''
        Load the ascii logfile located at self.filename.
        This method is automatically called upon instantiation.
        '''
        import numpy as np
        import datetime as dt

        # Slurp in entire file.
        infile = open(self.attrs['file'], 'r')
        raw = infile.readlines()
        infile.close()

        # Parse the header.
        self.attrs['descrip'] = raw.pop(0)
        namevar = (raw.pop(0)).split()
        nCols = namevar
        loc={}
        for i, name in enumerate(namevar):
            loc[name] = i

        # Check the last line for completeness.
        # Files read while simulation is running will have incomplete
        # lines.
        nPts=len(raw)
        checkline = raw[-1].split()
        if len(checkline) != nCols:
            nPts-=1
        self.attrs['npts'] = nPts

        # Pop time/date names off of Namevar.
        namevar.pop(namevar.index('year'))
        namevar.pop(namevar.index('mo'))
        namevar.pop(namevar.index('dy'))
        namevar.pop(namevar.index('hr'))
        namevar.pop(namevar.index('mn'))
        namevar.pop(namevar.index('sc'))
        namevar.pop(namevar.index('time'))
        
        # Create containers for data:
        self['runtime']=dmarray(np.zeros(nPts),attrs={'units':'s'})
        self['time'] = dmarray(np.zeros(nPts, dtype=object))
        self['iter'] = dmarray(np.zeros(nPts))
        for name in namevar:
            self[name] = dmarray(np.zeros(nPts))

        for i, line in enumerate(raw[:nPts]):
            vals = line.split()
            # Set time:
            self['time'][i]=(dt.datetime(
                int(vals[loc['year']]), # Year
                int(vals[loc['mo']  ]), # Month
                int(vals[loc['dy']]), # Day
                int(vals[loc['hr']]), # Hour
                int(vals[loc['mn']]), # Minute
                int(vals[loc['sc']]), # Second
                ))
            self['runtime'][i]=float(vals[loc['time']])
            # Collect data
            for j, name in enumerate(namevar):
                self[name][i] = float(vals[loc[name]])

    def add_dst_quicklook(self, target=None, loc=111, showObs=True,
                          showBiot=True):
        '''
        Create a quick-look plot of Dst.  If kyotodst module is
        installed, compare to observed Dst.

        If kwarg 'target' is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure.
        if target is a matplotlib Axes object, the plot is placed
        into that axis.

        Two values are returned: the matplotlib Figure and Axes objects
        used to generate the plot (in that order.)

        =========== ================================================
        Kwarg       Description
        ----------- ------------------------------------------------
        target      Set plot destination.  Defaults to new figure.
        loc         Set subplot location.  Defaults to 111.
        showObs     Show observed Dst?  Defaults to **True**
        showBiot    Show Biot-Savart Dst?  Defaults to **True**
        =========== ================================================

        '''
        
        import matplotlib.pyplot as plt
        
        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(loc)
        
        ax.plot(self['time'], self['dstRam'], label='RAM Dst (DPS)')
        if 'dstBiot' in self and showBiot:
            ax.plot(self['time'], self['dstBiot'], label='RAM Dst (Biot)')
        ax.hlines(0.0, self['time'][0], self['time'][-1], 
                  'k', ':', label='_nolegend_')
        apply_smart_timeticks(ax, self['time'])
        ax.set_ylabel('Dst ($nT$)')
        ax.set_xlabel('Time from '+ self['time'][0].isoformat()+' UTC')

        try:
            import spacepy.pybats.kyoto as kt
        except ImportError:
            return fig, ax
        
        if showObs:
            try:
                stime = self['time'][0]; etime = self['time'][-1]
                if not hasattr(self, 'obs_dst'):
                    self.obs_dst = kt.fetch('dst',stime,etime)
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

############################################################################
class IonoPotScb(object):
    '''
    The 3D equilibrium code produces NetCDF files that contain the 
    ionospheric potential on the polar cap as well as mapped to the 
    equatorial plane.  The IonoPotScb object can be used to parse and
    visualize the data quickly.

    For quick parsing of these files, PyNIO is required.
    '''

    def __getitem__(self, key):
        if key in self.filedata:
            return self.filedata[key].get_value()
        elif key in self.data:
            return self.data[key]
        else:
            raise KeyError('Key not found in object.')
    def __setitem__(self, key,value):
        self.data[key]=value

    def keys(self):
        '''
        List all keys for data objects stored in the IonoPotScb object.
        '''
        return list(self.filedata.keys()) + list(self.data.keys())

    def __init__(self, filename):
        try:
            from PyNIO import Nio
        except ImportError:
            try:
                import Nio
            except ImportError:
                raise ImportError('PyNIO required, not found.')

        self.filename=filename
        
        # Load file as Nio object.
        f = Nio.open_file(filename, 'r')
        self.NameVars = list(f.variables.keys())
        self.attrs = f.attributes
        # self.filedata contains the raw NioVariable objects.
        self.filedata = f.variables
        self.time = f.variables['time'].get_value()
        # New values saved as self[key] are stored in self.data, not
        # self.filedata which cannot be changed.
        self.data={}

        # Determine units of potential.
        if self['PhiIono'].max()/1000.0 > 10.0:
            self.units='V'
        else:
            self.units='kV'

        self._netcdf = f

    def calc_pot_drop(self):
        '''
        Calculate the electric potential drop across the inner magnetosphere
        for the whole time period.  This value is analogous to cross polar
        cap potential drop, but is only valid inside of the RAM-SCB domain.

        This value is simply the maximum potential minus the minimum potential.
        The new data entry is stored using key 'ceqp', which stands for the
        cross equatorial potential.
        '''
        from numpy import zeros
        
        self['ceqp'] = zeros(len(self.time))
        
        for i in range(len(self.time)):
            self['ceqp'][i] = self['PhiIono'][i,:,:].max() - \
                self['PhiIono'][i,:,:].min()
        if self.units=='V':
            self['ceqp']=self['ceqp']/1000.0


    def close(self):
        '''
        Close the NetCDF file.  This will make values stored inside of the
        file unavailable to the user.
        '''
        self._netcdf.close()

    def plot_eqPot(self, time, target=None, range=200, n=31, add_cbar=True):
        '''
        Plot the equatorial electric potential in kV to target, where
        target may be a matplotlib figure or axis or None.  If target 
        is a figure, a new subplot is created.  If target is None, a 
        new figure AND axis is created.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
        from matplotlib.ticker import MultipleLocator
        from numpy import linspace

        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(111)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        cmap=get_iono_cb('bwr')
        crange = Normalize(vmin=-1.*range, vmax=range)
        levs = linspace(-1*range, range, n)

        factor = 1.0
        if self.units=='V':
            factor = 1000.0
        
        ax.set_aspect('equal')
        ax.set_xlabel('Y (R$_{E}$)')
        ax.set_ylabel('X (R$_{E}$)')
        cont = ax.contourf(self['yEq'][time,:,:], self['xEq'][time,:,:],
                           self['PhiIono'][time,:,:]/factor, 
                           levs, norm=crange, cmap=cmap)
        add_body(ax, rotate=90.0)
        cbar=False
        if add_cbar:
            cticks = MultipleLocator(50)
            cbar = plt.colorbar(cont, ticks=cticks, shrink=0.8, pad=0.08)
            cbar.set_label('kV')
        return fig, ax, cont, cbar

############################################################################
class Currents(object):
    '''
    The 3D equilibrium code produces NetCDF files that contain the 
    electric currents throughout the domain.  Currents objects parse
    and visualize this data.

    For quick parsing of these files, PyNIO is required.
    '''

    def __getitem__(self, key):
        if key in self.filedata:
            return self.filedata[key].get_value()
        elif key in self.data:
            return self.data[key]
        else:
            raise KeyError('Key not found in object.')
    def __setitem__(self, key,value):
        self.data[key]=value

    def keys(self):
        '''
        List all keys for data objects stored in the IonoPotScb object.
        '''
        return list(self.filedata.keys()) + list(self.data.keys())

    def __init__(self, filename):
        try:
            from PyNIO import Nio
        except ImportError:
            try:
                import Nio
            except ImportError:
                raise ImportError('PyNIO required, not found.')

        self.filename=filename
        
        # Load file as Nio object.
        f = Nio.open_file(filename, 'r')
        self.NameVars = list(f.variables.keys())
        self.attrs = f.attributes
        # self.filedata contains the raw NioVariable objects.
        self.filedata = f.variables
        self.time = f.variables['time'].get_value()
        # New values saved as self[key] are stored in self.data, not
        # self.filedata which cannot be changed.
        self.data={}

        self._netcdf = f

    def close(self):
        '''
        Close the NetCDF file.  This will make values stored inside of the
        file unavailable to the user.
        '''
        self._netcdf.close()

    def plot_jpar(self, ax, time,maxz=0.5):
        from spacepy.pybats.rim import get_iono_cb
        from matplotlib.colors import Normalize
        cmap=get_iono_cb('bwr')
        crange = Normalize(vmin=-1.*maxz, vmax=maxz)
        cnt1=ax.contourf(self['xPole'][time,:,:], self['yPole'][time,:,:],
                         self['jParVas'][time,:,:], norm=crange, cmap=cmap)

############################################################################
class ParamFile(object):
    '''
    It is sometimes necessary to read and parse a PARAM.in file.  ParamFile
    objects quickly parse the info and store it in a convenient data structure.
    '''

    def __repr__(self):
        return 'PARAM file at %s with %i commands.' % \
            (self.path+self.file, len(self.cmd))

    # Make 'er work like a dict.
    def __getitem__(self,key):
        return self.cmd[key]
    def keys(self):
        return list(self.cmd.keys())

    def __contains__(self, key):
        return key in self.cmd

    def __init__(self, infile):
        '''
        Read the file, parse each command.
        '''
        # Grab the full path and save it.
        if infile.rfind('/')>-1:
            self.path=infile[0:infile.rfind('/')+1]+'/'
        else:
            self.path='./'
        self.file=infile
        # Slurp entire file.
        f=open(infile,'r')
        lines=f.readlines()
        f.close()

        # Some commands have specialized parse functions.
        # Those are connected here.
        self.known_cmds={
            'STARTTIME':self._parse_starttime,
            'STOP':self._parse_stop}
        self.cmd={}
        self._parse(lines)

    # Command specific parsers:
    def _add_cmd(self,cmd_name,lines):
        self.cmd[cmd_name]=[]; npop=0
        for l in lines:
            if l.strip()=='' or l[0]=='#':
                break
            self.cmd[cmd_name].append(l.strip())
            npop+=1
        return npop    
    def _parse_stop(self,lines):
        self.iters=int(lines[0].split()[0])
        self.dur  =float(lines[1].split()[0])
        self.stoptime=self.starttime+dt.timedelta(seconds=self.dur)
        return 2
    def _parse_starttime(self, lines):
        self.starttime=dt.datetime(
            int(lines[0].split()[0]), #year
            int(lines[1].split()[0]), #month
            int(lines[2].split()[0]), #day
            int(lines[3].split()[0]), #hour
            int(lines[4].split()[0]), #minute
            int(lines[5].split()[0]), #second
            int(float(lines[6].split()[0])*1E6)) #millisec
        return 7

    # Quick function to pop multiple lines.
    def _parse(self, lines):
        # A quick function for multi-popping.
        def pop(arr, n):
            for i in range(n):
                arr.pop(0)
            return arr
        # Loop throug all lines!
        while len(lines)>0:
            l=lines.pop(0)
            if l[0]!='#': continue # Find first/next #COMMAND statement.
            l=l.strip()

            # Skip begin/end comps.
            if l.find('BEGIN_COMP')>-1 or l.find('END_COMP')>-1:
                continue

            if l[1:] in self.known_cmds:
                # If recognized command, process and pop arguments.
                npop=self.known_cmds[l[1:]](lines)
                pop(lines, npop)
            else:
                npop=self._add_cmd(l[1:],lines)
                pop(lines,npop)


############################################################################
class GeoMltFile(object):
    '''
    GeoMltFile is a class to open, read, manipulate and write files that 
    contain LANL geosynchronous multi-satellite averaged, MLT-interpolated
    fluxes.
    '''
    
    def __init__(self, filename=None, scrub=True):
        # Create empty arrays.
        self.flux=np.zeros([288,24,36])
        self.time=np.zeros(288, dtype=object)
        self.egrid=np.zeros(36)
        self.nsats=np.zeros( (288, 24, 36), dtype=int )
        self.header=[]
        self.particle='proton'

        if filename:
            self.filename=filename
            self.read()

    def read(self):
        '''
        Load contents from *self.filename*.
        '''

        ### Open file. ###
        f=open(self.filename, 'r')

        ### Parse header information. ###
        # Read header
        for i in range(3):
            self.header.append(f.readline())
        # Set measurement type:
        if (self.header[0]).find('electron')>0: self.particle='electron'
        # Parse energy grid.
        parts=f.readline().split()
        parts.pop(0); parts.pop(0); parts.pop(0)
        self.egrid[:]=parts[-36:]

        ### Parse fluxes. ###
        # Slurp rest of file and close.
        lines = f.readlines()
        f.close()
        # Grab L-grid.
        lt=[]
        lt.append(float(lines[0].split()[1]))
        lt.append(float(lines[1].split()[1]))
        i=2
        while lt[0]!=lt[-1]:
            lt.append(float(lines[i].split()[1]))
            i+=1
        self.lgrid=np.array(lt[:-1])
        nL=self.lgrid.size
        # Parse file one epoch at a time.
        for i in range(288):
            # Get one epoch worth of data.
            sublines=lines[nL*i:nL*i+nL]
            # Get time for this epoch.
            t=sublines[0].split()[0]
            self.time[i]=dt.datetime(
                int(t[0:4]),   int(t[5:7]),   int(t[8:10]),
                int(t[11:13]), int(t[14:16]), int(t[17:19]), 1000*int(t[20:23]))
            # Parse rest of lines.
            for j,l in enumerate(sublines):
                # Use some string comprehension magic.
                self.nsats[i,j,:]=[l[32 + 2*x: 32+ 2*x+ 2] for x in range(36)]
                self.flux[i,j,:] =[l[104+18*x:104+18*x+18] for x in range(36)]
                
    def scrub(self, lastflux=None):
        '''
        GeoMlt files often have bad data througout in the form of negative
        or zero fluxes.  In RAM-SCB, these are replaced with the last good
        data value.  This function "scrubs" the data in the same manner as
        RAM-SCB.  The kwarg *lastflux* may be given as a nLT x nE numpy array
        of the last good fluxes, e.g. the correct values at the epoch directly
        preceeding the first epoch of this file.  These will be used if the
        first time entry of the file contains bad values.  If not given, bad
        values at the first time entry will be set to zero.
        '''
        
        # Handle first line:
        if lastflux:
            # Check lastflux.
            if lastflux.shape != (self.lgrid.size, self.egrid.size):
                raise ValueError("Shape of lastflux is incorrect.")
            self.flux[0,self.flux[0,:,:]<=0.0]=lastflux[self.flux[0,:,:]<=0.0]
        else:
            self.flux[0,self.flux[0,:,:]<=0.0]=0.0

        for i in range(1,288):
            self.flux[i,self.flux[i,:,:]<=0.0] = \
                self.flux[i-1,self.flux[i,:,:]<=0.0]
                               

    def __iadd__(self, other):
        '''
        Append another **GeoMltFile** object to this one.
        '''

        if type(self) != type(other):
            raise TypeError(
                'object can only be combined with one of same type.')

        # Grid checks:
        if any(self.egrid != other.egrid):
            raise ValueError('both items must have identical energy grids!')
        if any(self.lgrid != other.lgrid):
            raise ValueError('both items must have identical MLT grids!')
        if self.particle != other.particle:
            raise ValueError('both items must have identical partilces!')

        # Append data together:
        self.time = np.append(self.time, other.time)
        self.flux = np.append(self.flux, other.flux, 0)
        self.nsats= np.append(self.nsats, other.nsats, 0)

        # Append file names:
        if type(self.filename)==str:
            self.filename = [self.filename]
        self.filename.append(other.filename)

        return self


    def plot_epoch_flux(self, epoch=0, target=None, loc=111):
        '''
        Plot fluxes for a single file epoch.
        '''

        from numpy import array
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)


        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
        else:
            fig = plt.figure()
            ax = fig.add_subplot(loc)

        egrid = array(range(len(self.egrid)+1))
        lgrid = array(range(len(self.lgrid)+1))
        flux = self.flux[epoch,:,:].transpose()
        flux[flux<0.01] = 0.01
        flx = ax.pcolormesh(lgrid, egrid, flux, norm=LogNorm(),
                            vmin=0.01, vmax=1e10, 
                            cmap=plt.get_cmap('gnuplot2'))
        cbar = plt.colorbar(flx, pad=0.01, shrink=0.85, ticks=LogLocator(),
                            format=LogFormatterMathtext())
        cbar.set_label('$cm^{-2}s^{-1}ster^{-1}keV^{-1}$')
        ax.set_xlim([0,24])
        ax.set_ylim([0,len(egrid)])
        ax.set_title('Flux at Boundary - %s'%self.filename)
        ax.set_xlabel('Local Time Sector')
        # Use actual energy channels.
        ax.set_ylabel('Energy (keV)')
        newlabs=[]
        for val in ax.get_yticks()[:-1]:
            newlabs.append('%6.2f' % self.egrid[int(val)])
        ax.set_yticklabels(newlabs)
        

        return fig, ax
