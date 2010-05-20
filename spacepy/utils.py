# -*- coding: utf-8 -*-
"""
Set of generic utilities.

--++-- By Steve Morley --++--

smorley@lanl.gov/morley_steve@hotmail.com,
Los Alamos National Laboratory, ISR-1,
PO Box 1663, Los Alamos, NM 87545

"""
from spacepy import help
print "Warning (utils): This module is deprecated. Use toolbox instead."


def doy2md(year, doy):
    """Convert day-of-year to month and day
    
    Inputs:
    =======
    Year, Day of year (Jan 1 = 001)
    
    Returns:
    ========
    Month, Day (e.g. Oct 11 = 10, 11)
    
    Note: Implements full and correct leap year rules.
    
    Modification history:
    Created by Steve Morley (in July '05), rewritten for
    Python in October 2009
    """
    
    import numpy as np
    year, doy = np.array(year, dtype=float), np.array(doy, dtype=float)
    
    mn_arr = np.array([1,32,60,91,121,152,182,213,244,274,305,335], dtype=float)
    leapyr = (year % 4) == 0
    if leapyr == True:
        if (year % 100 == 0) and (year % 400) != 0:
            pass
        else:
            mn_arr[2:] = mn_arr[2:]+1
    [a] = np.where(doy >= mn_arr)
    month = a[-1:]+1
    day = (doy-mn_arr[a[-1:]])+1
    
    return [month,day]


def t_overlap(ts1,ts2):
    """Finds the overlapping elements in two lists of datetime objects
    
    Returns:
    ========
     - indices of 1 within interval of 2, & vice versa
    
    Example:
    ========
     - Given two series of datetime objects, event_dates and omni['Time']:
    
    >>> import spacepy.utils as utils
    >>> [einds,oinds] = utils.t_overlap(event_dates, omni['Time'])
    >>> omni_time = omni['Time'][oinds[0]:oinds[-1]+1]
    >>> print omni_time
    [datetime.datetime(2007, 5, 5, 17, 57, 30), datetime.datetime(2007, 5, 5, 18, 2, 30),
    ... , datetime.datetime(2007, 5, 10, 4, 57, 30)]
    """
    
    import datetime as dt
    import numpy as np
    from matplotlib.dates import date2num
    
    tn1, tn2 = date2num(ts1),date2num(ts2)
    
    dum = abs(tn1-tn2[0])
    in1st = np.where(min(dum) == dum)
    dum = abs(tn1-tn2[-1])
    in1en = np.where(min(dum) == dum)
    inds1 = range(int(in1st[0][0]),int(in1en[0][0])+1)
    
    dum = abs(tn2-tn1[0])
    in2st = np.where(dum == min(dum))
    dum = abs(tn2-tn1[-1])
    in2en = np.where(dum == min(dum))
    inds2 = range(int(in2st[0][0]),int(in2en[0][0])+1)
    
    return inds1, inds2

def t_common(ts1, ts2, mask_only=True):
    """Finds the elements in a list of datetime objects present in another
    
    Returns:
    ========
     - Two element tuple of truth tables (of 1 present in 2, & vice versa)
    """
    
    import datetime as dt
    import numpy as np
    from matplotlib.dates import date2num, num2date
    
    tn1, tn2 = date2num(ts1),date2num(ts2)
    
    el1in2 = np.setmember1d(tn1,tn2) #makes mask of present/absent
    el2in1 = np.setmember1d(tn2,tn1)
    
    if mask_only:
        return el1in2, el2in1
    else:
        truemask1 = np.abs(np.array(el1in2)-1)
        truemask2 = np.abs(np.array(el2in1)-1)
        time1 = np.ma.masked_array(tn1, mask=truemask1)
        time2 = np.ma.masked_array(tn2, mask=truemask2)
        dum1 = num2date(time1.compressed())
        dum2 = num2date(time2.compressed())
        if type(ts1)==np.ndarray:
            dum1 = np.array(dum1)
            dum2 = np.array(dum2)
            
        return dum1, dum2

def ShueMP(P,Bz):
    """Calculates the Shue et al. (1997) subsolar magnetopause radius
    
    Ported from Drew Turner's (LASP) MatLab script
    
    Inputs:
    =======
    SW ram pressure [nPa], IMF Bz (GSM) [nT]
    
    Output:
    =======
    Magnetopause (sub-solar point) standoff distance [Re]
    """
    
    import numpy as np
    
    # Initialize r0 and make P and Bz numpy arrays
    r0 = np.zeros((len(P)),dtype=float)
    Bz = np.array(Bz)
    P = np.array(P)
    
    # Find where Bz >= 0 and where it is < 0
    iBzPos = np.where(Bz >= 0)
    iBzNeg = np.where(Bz < 0)
    
    # Calculate r0
    r0[iBzPos] = (11.4 + 0.013*Bz[iBzPos])*P[iBzPos]**(-1/6.6)
    r0[iBzNeg] = (11.4 + 0.140*Bz[iBzNeg])*P[iBzNeg]**(-1/6.6)
    
    return r0


def MoldwinPP(Kp):
    """Calculates the Moldwin et al. (2002) subsolar magnetopause radius
    
    Inputs: Kp index (uses max Kp from prev 12 hrs)
    Output: Plasmapause radius [Re]
    """
    
    import numpy as np
    
    Kp = np.array(Kp)
    Lpp = np.zeros(len(Kp))
    for i in range(3,len(Kp)):
        Kpmax = np.max(Kp[i-3:i+1])
        Lpp[i] = 5.99 - 0.382*Kpmax
    Lpp[:3] = np.nan
    
    return Lpp , Kpmax


def readOMNI(fname):
    """Read amalgamated OMNI2 file into dictionary of data arrays
    
    Input: filename (str)
    
    Output: dictionary with datetime objects and data arrays
    """
    
    import numpy as np
    import datetime as dt
    from utils import doy2md
    
    dum = np.fromfile(fname, dtype=float, count=-1, sep=' ')
    dum = dum.reshape(-1,54) #reshape to 54 column (reads as flat array)
    OMNItime = []
    for i, year in enumerate(dum[0:,0]):
        [mon, day] = doy2md(year, dum[i,1])
        d = dt.datetime(int(year), int(mon), int(day), int(dum[i,2]), 30)
        OMNItime.extend([d]) 
    
    dum[dum[0:,16]==999.9, 16] = np.nan
    dum[dum[0:,25]==9999., 25] = np.nan
    dum[dum[0:,22]==9999999., 22] = np.nan
    dum[dum[0:,24]==9999., 24] = np.nan
    dum[dum[0:,23]==999.9, 23] = np.nan
    dum[dum[0:,28]==99.99, 28] = np.nan
    dum[dum[0:,40]==99999., 40] = np.nan
    dum[dum[0:,38]==99., 38] = np.nan
            
    out_dict = {'Time': OMNItime, 'Bz': dum[:,16], 'Vx': dum[:,24], 'Vy': dum[:,25], 'n': dum[:,23], \
        'T': dum[:,22], 'P': dum[:,28], 'Dst': dum[:,40], 'Kp': dum[:,38]}
    return out_dict


def readOMNIhi(fname):
    """Read high-resOMNI file into dictionary of data arrays
    
    Input: filename (str)
    
    Output: dictionary with datetime objects and data arrays
    """
    import numpy as np
    import datetime as dt
    from utils import doy2md
    
    dum = np.loadtxt(fname, dtype=float)
    OMNItime = []
    for i, year in enumerate(dum[0:,0]):
        [mon, day] = doy2md(year, dum[i,1])
        d = dt.datetime(int(year), int(mon), int(day), int(dum[i,2]), \
            int(dum[i,3])+2, 30)
        OMNItime.extend([d]) 
    
    dum[dum[0:,16]>=999.9, 16] = np.nan
    dum[dum[0:,21]>=9999., 21] = np.nan
    dum[dum[0:,23]>=9999., 23] = np.nan
    dum[dum[0:,25]>=999.9, 25] = np.nan
    dum[dum[0:,26]>=9999999.,26] = np.nan
    dum[dum[0:,27]>=99.99, 27] = np.nan
    #dum[dum[0:,41]>=??, 41] = np.nan  #what is bad val for SymH?

    out_dict = {'Time': OMNItime, 'Bz': dum[:,16], 'Vx': dum[:,21], \
        'Vy': dum[:,23], 'n': dum[:,25], 'P': dum[:,27], 'SymH': dum[:,41], \
        'T': dum[:,26]}

    return out_dict


def mu_ai(energy, b=100.e-9, rme=.511):
    """Calculate 1st adiabatic invariant given energy in [MeV/G]
    
    Input: energy (req'd) in MeV
    
    Uses E = Ek + E0, E = sqrt(p^2c^2 + E0^2);
    Then uses mu = p^2/2mB to arrive at
    mu = (p^2c^2)/(2*E0*B)
    """

    p2c2 = energy*(energy+(2*rme))
    bg = b*1.e4
    mu = (p2c2)/(2.*rme*bg)
    return mu


def dm_ll(kp=2.7,l=6.6):
    """Magnetic field diffusion coefficient from Brautigam and Albert"""

    fac = 0.506*kp-9.325
    dmll = (10**fac)*(l**10.)
    
    return dmll


def windowmean(data, time=[], winsize=0, overlap=0, pts=True):
    """Windowing mean function, window overlap is user defined
    
    Inputs:
    data - 1D series of points;
    time - series of timestamps, optional (format as numeric or datetime);
    For non-overlapping windows set overlap to zero.
    e.g.,
    
    >>> wsize, olap = datetime.timedelta(1), datetime.timedelta(0,3600)
    
    >>> outtime, outdata = windowmean(data, time, winsize=wsize, overlap=olap)
    
    where the time, winsize and overlap are either numberic or datetime objects,
    in this example the window size is 1 day and the overlap is 1 hour.
    
    Caveats: This is a quick and dirty function - it is NOT optimised, at all.
    """
    import numpy as np
    import datetime as dt
    from matplotlib.dates import date2num,num2date
    
    #check inputs and initialize
    #Set resolution to 1 if no times supplied
    if len(time)==0:
        startpt, res = 0, 1
        time = range(len(data))
        dateflag = False
    else:
        try:
            assert len(data) == len(time)
        except:
            return 'windowmean error: data and time must have same length'
        #First check if datetime objects: yes? - convert to serial
        if (type(time[0]) == dt.datetime):
            if not winsize:
                return 'windowmean error: winsize must be set for datetime input'
            else:
                try:
                    assert type(winsize) == dt.timedelta
                    assert type(overlap) == dt.timedelta
                except:
                    return 'windowmean error: winsize/overlap must be timedeltas'
            #get window and overlap in serial format
            winsize = date2num(time[0]+winsize) - date2num(time[0])
            overlap = date2num(time[0]+overlap) - date2num(time[0])
            time = date2num(time) #put time in serial format
            pts = False #force time-based averaging
            dateflag = True #set flag so we return outtime as datetime objects
        else:
            try:
                assert type(winsize) == dt.timedelta
                assert type(overlap) == dt.timedelta
            except:
                return 'windowmean error: winsize/overlap must be timedeltas'
            dateflag = False
            pts = False
            startpt = time[0]
            #res = time[1]-time[0]
    
    #now actually do windowing mean
    outdata, outtime = [], []
    if pts:
        #loop for fixed number of points in window
        if winsize % 1:
            winsize = round(winsize)
            print 'windowmean error: non-integer windowsize, rounding to %d' \
            % winsize
        if winsize < 1:
            winsize = 1
            print 'windowmean error: window length < 1, defaulting to 1'
        if overlap >= winsize:
            overlap = winsize - 1
            print '''windowmean error: overlap longer than window, truncated to
            %d''' % overlap
        lastpt = winsize-1 #set last point to end of window size
        while lastpt < len(data):
            datwin = np.ma.masked_where(np.isnan(data[startpt:startpt+winsize]), \
                data[startpt:startpt+winsize])
            getmean = np.mean(datwin.compressed()) #mean of window, excl. NaNs
            print len(time),startpt,winsize
            gettime = (time[startpt+winsize] - time[startpt])/2. \
                + time[startpt]#new timestamp
            startpt = startpt+winsize-overlap
            lastpt = startpt+winsize
            outdata.append(getmean) #construct output arrays
            outtime.append(gettime)
    else:
        #loop with time-based window
        print '''time-based averaging'''
        lastpt = time[0] + winsize
        startpt = time[0]
        res = time[1]-time[0]
        if overlap >= winsize:
            overlap = winsize - res
            print '''windowmean error: overlap longer than window, truncated to
            %d''' % overlap
        while lastpt < time[-1]:
            keep0 = np.where(time < startpt+winsize)
            keep1 = np.where(time >= startpt)
            getinds = np.intersect1d(keep0[0],keep1[0])
            if getinds.any():
                getdata = np.ma.masked_where(np.isnan(data[getinds]),data[getinds])
                getmean = np.mean(getdata.compressed()) #find mean excluding NaNs
            else:
                getmean = np.nan
            gettime = startpt + winsize/2. #new timestamp
            startpt = startpt + winsize - overlap #advance window start
            lastpt = startpt + winsize
            outdata.append(getmean) #construct output arrays
            outtime.append(gettime)
            
    if dateflag: #return times as datetime objects if flag set
        outtime = num2date(outtime)
        
    return outdata, outtime
    
    
def medabsdev(series):
    """Calculate median absolute deviation of a given input series
    
    Median absolute deviation (MAD) is a robust and resistant measure of
    the spread of a sample (same purpose as standard deviation). The
    MAD is preferred to the interquartile range as the interquartile
    range only shows 50% of the data whereas the MAD uses all data but
    remains robust and resistant. See e.g. Wilks, Statistical methods
    for the Atmospheric Sciences, 1995, Ch. 3.
    
    This implementation is robust to presence of NaNs
    
    Example:
    Find the median absolute deviation of a data set. Here we use the log-
    normal distribution fitted to the population of sawtooth intervals, see
    Morley and Henderson, Comment, Geophysical Research Letters, 2009.
    
    >>> data = numpy.random.lognormal(mean=5.1458, sigma=0.302313, size=30)
    >>> print data
    array([ 181.28078923,  131.18152745, ... , 141.15455416, 160.88972791])
    >>> utils.medabsdev(data)
    28.346646721370192
    """
    
    import numpy as np
    #ensure input is numpy array (and make 1-D)
    series = (np.array(series, dtype=float)).ravel()
    #mask for NaNs
    series = np.ma.masked_where(np.isnan(series),series)
    #get median absolute deviation of unmasked elements
    perc50 = np.median(series.compressed())
    mad = np.median(abs(series.compressed()-perc50))
    
    return mad
    
    
def makePoly(x, y1, y2, face = 'blue', alpha=0.5):
    """Make filled polygon for plotting
    
    Can be replaced by built-in matplotlib function fill_between
    
    >>> poly0c = makePoly(x, ci_low, ci_high, face='red', alpha=0.8)
    >>> ax0.add_patch(poly0qc)
    """
    
    import numpy as np
    import matplotlib as mpl
    x2, y1 = x[-1::-1], y1[-1::-1]
    polyx = np.concatenate((x,x2))
    polyy = np.concatenate((y2,y1)) 
    xy = np.empty((len(polyy),2))
    xy[:,0], xy[:,1] = polyx, polyy 
    madePoly = mpl.patches.Polygon(xy, facecolor = face, alpha = alpha)
    
    return madePoly
    
def binHisto(data):
    """Calculates bin width and number of bins for histogram using Freedman-Diaconis rule
    """
    import numpy as np
    from matplotlib.mlab import prctile
    pul = prctile(data, p=(25,75)) #get confidence interval
    ql, qu = pul[0], pul[1]
    iqr = qu-ql
    binw = 2.*iqr/(len(data)**(1./3.))
    nbins = round((np.max(data)-np.min(data))/binw)
    
    return (binw, nbins)