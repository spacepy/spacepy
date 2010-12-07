#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Toolbox of various functions and generic utilities.

"""

from __future__ import division
try:
    from spacepy import help
except ImportError:
    pass
except:
    pass
__version__ = "$Revision: 1.63 $, $Date: 2010/12/07 17:02:38 $"
__author__ = 'S. Morley and J. Koller'


def tOverlap(ts1, ts2):
    """Finds the overlapping elements in two lists of datetime objects

    @param ts1: first set of datetime object
    @type ts1: datetime
    @param ts2: datatime object
    @type ts2: datetime
    @return: indices of ts1 within interval of ts2, & vice versa
    @rtype: list
    
    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com
    
    @version: 09-Jun-2010: Total rewrite for clarity, efficiency and brevity.
        
    Given two series of datetime objects, event_dates and omni['Time']:
    
    >>> import spacepy.toolbox as tb
    >>> [einds,oinds] = tb.tOverlap(event_dates, omni['Time'])
    >>> omni_time = omni['Time'][oinds[0]:oinds[-1]+1]
    >>> print omni_time
    [datetime.datetime(2007, 5, 5, 17, 57, 30), datetime.datetime(2007, 5, 5, 18, 2, 30),
    ... , datetime.datetime(2007, 5, 10, 4, 57, 30)]
    
    """
    
    #elements of ts2 within bounds of ts1
    t_lower, t_upper = min(ts1), max(ts1)
    bool1 = [v > t_lower for v in ts2]
    bool2 = [v < t_upper for v in ts2]
    mask2in1 = [b1 and b2 for b1, b2 in zip(bool1, bool2)]
    inds2in1 = [i for i, val in enumerate(mask2in1) if val==True]
    
    #elements of ts1 within bounds of ts2
    t_lower, t_upper = min(ts2), max(ts2)
    bool1 = [v > t_lower for v in ts1]
    bool2 = [v < t_upper for v in ts1]
    mask1in2 = [b1 and b2 for b1, b2 in zip(bool1, bool2)]
    inds1in2 = [i for i, val in enumerate(mask1in2) if val==True]
    
    if len(inds2in1) == 0:
        inds2in1 = None
    if len(inds1in2) == 0:
        inds1in2 = None
    
    return inds1in2, inds2in1
    

def tCommon(ts1, ts2, mask_only=True):
    """Finds the elements in a list of datetime objects present in another

    @param ts1: first set of datetime object
    @type ts1: datetime
    @param ts2: datatime object
    @type ts2: datetime
    @return: Two element tuple of truth tables (of 1 present in 2, & vice versa)
    @rtype: tuple

    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com
    
    """
    import numpy as np
    from matplotlib.dates import date2num, num2date
    
    tn1, tn2 = date2num(ts1), date2num(ts2)
    
    v_test = np.__version__.split('.')
    if v_test[0] == 1 and v_test[1] <= 3:
        el1in2 = np.setmember1d(tn1, tn2) #makes mask of present/absent
        el2in1 = np.setmember1d(tn2, tn1)
    else:
        el1in2 = np.in1d(tn1, tn2, assume_unique=True) #makes mask of present/absent
        el2in1 = np.in1d(tn2, tn1, assume_unique=True)
    
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

def loadpickle(fln):
    """
    load a pickle and return content as dictionary

    @param fln: filename
    @type fln: string
    @return: dictionary with content from file
    @rtype: dictionary
    
    @see: savepickle

    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010

    >>> d = loadpickle('test.pbin')
    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle

    fh = open(fln, 'rb')
    content = pickle.load(fh)
    fh.close()

    return content


# -----------------------------------------------
def savepickle(fln, dict):
    """
    save dictionary variable dict to a pickle with filename fln 

    @param fln: filename
    @type fln: string
    @param dict:  container with stuff
    @type dict: dictionary
    
    @see: loadpickle

    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010
    
    >>> d = {'grade':[1,2,3], 'name':['Mary', 'John', 'Chris']}
    >>> savepickle('test.pbin', d)
    
    """
    try:
        import cPickle as pickle
    except:
        import pickle

    fh = open(fln, 'wb')
    pickle.dump(dict, fh, 2) # 2 ... fast binary
    fh.close()

    return

# -----------------------------------------------
def assemble(fln_pattern, outfln, sortkey='ticks'):
    """
    assembles all pickled files matching fln_pattern into single file and 
    save as outfln. Pattern may contain simple shell-style wildcards *? a la fnmatch
    file will be assembled along time axis given by Ticktock (key: 'ticks') in dictionary
    If sortkey = None, then nothing will be sorted

    @param fln_pattern: pattern to match filenames
    @type fln_pattern: string
    @param outfln: filename to save combined files to
    @type outfln: string
    @return: dcomb - dictionary with combined values
    @rtype: dict
    
    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010

    >>> assemble('input_files_*.pkl', 'combined_input.pkl')
    adding input_files_2001.pkl
    adding input_files_2002.pkl
    adding input_files_2004.pkl
    writing: combined_input.pkl
    """

    import glob
    import numpy as n
    from . import time as t
    from . import coordinates as c

    filelist = glob.glob(fln_pattern)
    filelist = human_sort(filelist)
    # read all files
    d = {}
    for fln in filelist:
        print("adding ", fln)
        d[fln] = loadpickle(fln)
    
    # combine them
    dcomb = d[filelist[0]]  # copy the first file over
    for fln in filelist[1:]:
       # check if sortkey is actually available
        assert (sortkey in d[fln] or sortkey==None), 'provided sortkey ='+sortkey+' is not available'
        
        if sortkey:
            TAIcount = len(d[fln][sortkey])
        else:
            TAIcount = len(d[fln][ list(d[fln].keys())[0] ])
        for key in d[fln]:
            #print fln, key
            dim = n.array(n.shape(d[fln][key]))
            ax = n.where(dim==TAIcount)[0]
            if len(ax) == 1: # then match with TAI length is given (jump over otherwise like for 'parameters')
                if isinstance(dcomb[key], t.Ticktock):
                    dcomb[key] = dcomb[key].append(d[fln][key])
                elif isinstance(dcomb[key], c.Coords):
                    dcomb[key] = dcomb[key].append(d[fln][key])
                else:
                    dcomb[key] = n.append(dcomb[key], d[fln][key], axis=ax)

    if sortkey:    #  then sort     
        if isinstance(dcomb[sortkey], t.Ticktock):
            idx = n.argsort(dcomb[sortkey].RDT)
        else:
            idx = n.argsort(dcomb[sortkey])
        TAIcount = len(dcomb[sortkey])
        for key in dcomb.keys():
            dim = n.array(n.shape(dcomb[key]))
            ax = n.where(dim==TAIcount)[0]
            if len(ax) == 1: # then match with length of TAI
                dcomb[key] = dcomb[key][idx] # resort
    else:
        # do nothing
        pass

    print('\n writing: ', outfln)
    savepickle(outfln, dcomb)

    return dcomb

# -----------------------------------------------    
def human_sort( l ): 
    """ Sort the given list in the way that humans expect. 
    http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
    """ 
    import re 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key ) 

    return l

# -----------------------------------------------
def feq(x, y, precision=0.0000005):
    """
    compare two floating point values if they are equal
    after: http://www.lahey.com/float.htm
    
    further info at::
        - http://docs.python.org/tut/node16.html
        - http://www.velocityreviews.com/forums/t351983-precision-for-equality-of-two-floats.html
        - http://www.boost.org/libs/test/doc/components/test_tools/floating_point_comparison.html
        - http://howto.wikia.com/wiki/Howto_compare_float_numbers_in_the_C_programming_language

    @param x: a number
    @type x: float
    @param y: float  or array of floats
    @type y: float
    @keyword precision: precision for equal (default 0.0000005)
    @type precision: float
    @return: True (equal) or False (not equal)
    @rtype: boolean

    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010
    @version: V2: 18-May-2010: User-specified precision added
       
    >>> index = where( feq(Lpos,Lgrid) ) # use float point comparison
 
    """

    boolean = abs(x-y) <= (abs(x+y)*precision)

    return boolean


# -----------------------------------------------
def dictree(in_dict, verbose=False, spaces=None, levels=True):
    """ pretty print a dictionary tree

    @param in_dict: a complex dictionary (with substructures)
    @type in_dict: dict
    @keyword verbose: print more info
    @type verbose: boolean
    @keyword spaces: string will added for every line
    @type spaces: string
    @keyword levels: number of levels to recurse through (True means all)
    @type levels: integer
    
    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010
    @version: V1.1: 24-Feb-2010 S. Morley, added verbose option
    @version: v1.2: 17-May-2010 S. Morley, added levels option

    >>> d = {'grade':{'level1':[4,5,6], 'level2':[2,3,4]}, 'name':['Mary', 'John', 'Chris']}
    >>> dictree(d)
    +
    |____grade
        |____level1
        |____level2
    |____name

    """

    try:
        assert isinstance(in_dict, dict)
    except AssertionError:
        raise TypeError('dictree: Input must be dictionary')

    if not spaces:
        spaces = ''
        print('+')
    
    if levels:
        try:
            assert levels is True
        except AssertionError:
            levels -= 1
            if levels == 0:
                levels = None

    for key in sorted(in_dict.keys()):
        if verbose:
            typestr = str(type(in_dict[key])).split("'")[1]
            if type(in_dict[key]) != dict:
                try:
                    dimstr = in_dict[key].shape
                    dimstr = ' ' + str(dimstr)
                except AttributeError:
                    try:
                        dimstr = len(in_dict[key])
                        dimstr = ' [' + str(dimstr) + ']'
                    except:
                        dimstr = ''
                print(spaces + '|____' + key + ' ('+ typestr + dimstr + ')')
            else:
                print(spaces + '|____' + key)
        else:
            print(spaces + '|____' + key)
        if isinstance(in_dict[key], dict) and levels:
            dictree(in_dict[key], spaces = spaces + '     ', verbose = verbose, levels = levels)

    return None

# -----------------------------------------------
def printfig(fignum, saveonly=False, pngonly=False, clean=False, filename=None):
    """save current figure to file and call lpr (print).

    This routine will create a total of 3 files (png, ps and c.png) in the 
    current working directory with a sequence number attached. Also, a time 
    stamp and the location of the file will be imprinted on the figure. The 
    file ending with c.png is clean and no directory or time stamp are 
    attached (good for powerpoint presentations).

    @param fignum: matplotlib figure number
    @type fignum: integer
    @keyword saveonly:  True (don't print and save only to file)  False (print and save)
    @type saveonly: boolean
    @keyword pngolny: True (only save png files and print png directly) False (print ps file, and generate png, ps; can be slow)
    @type pngonly: boolean
    @keyword clean: True (print and save only clean files without directory info) False (print and save directory location as well)
    @type clean: boolean
    @keyword filename: None (If specified then the filename is set and code does not use the sequence number)

    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010
    @version: V2: 19-Feb-2010: added pngonly and clean options, array/list support (JK)
    @version: V3: 21-Jul-2010: added filename keyword (BAL)

    >>> pyplot.plot([1,2,3],[2,3,2])
    >>> spacepy.printfig(1)
    """

    import os, glob, datetime
    import pylab as p
    
    try:
        nfigs = len(fignum)
    except:
        nfigs = 1
        fignum = [fignum]

    for ifig in fignum:
        # active this figure
        p.figure(ifig)

        if filename == None:
            # create a filename for the figure
            cwd = os.getcwd()
            num = len(glob.glob('*.png'))
            fln = cwd+'/figure_'+str(num)
        else:
            fln = filename
        # truncate fln if too long
        if len(fln) > 60: 
            flnstamp = '[...]'+fln[-60:]
        else:
            flnstamp = fln


        # save a clean figure without timestamps
        if clean == True:
            p.savefig(fln+'_clean.png')

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")
        # add the filename to the figure for reference
        p.figtext(0.01, 0.01, timestamp+flnstamp+'.png', rotation='vertical', va='bottom', size=8)

        # now save the figure to this filename
        if pngonly == False:
            p.savefig(fln+'.ps')
        
        p.savefig(fln+'.png')

        # send it to the printer
        if saveonly != True:
            if pngonly == False:
                os.popen('lpr '+fln+'.ps')
            else:
                os.popen('lpr '+fln+'.png')
            
    return

# -----------------------------------------------
def update(all=True, omni=False, leapsecs=False):
    """
    Download and update local database for omni, leapsecs etc

    @keyword all: if True, update all of them
    @type all: boolean
    @keyword omni: if True. update only onmi
    @type omni: boolean 
    @keyword leapsecs:  if True, update only leapseconds
    @type leapsecs: boolean
    @return: data directory where things are saved
    @rtype: string


    @author: Josef Koller
    @organization: Los Alamos National Lab
    @contact: jkoller@lanl.gov
    
    @version: V1: 20-Jan-2010
    @version: V1.1: 24-May-2010 Minor modification to return data directory (BAL)
    @version: V1.2: 11-Jun-2010 moved pickle_omni in here and added Qbits (JK)

    >>> update(omni=True)
    """
    
    try:
        from spacepy.time import Ticktock, doy2date
    except:
        exec(open('spacepy/time.py').read())
        #import sys
        #sys.path.append('./spacepy')
        #from spacepy import time as st
        
    
    import os, sys
    import zipfile
    import datetime
    #import spacepy.time as st
    from spacepy import savepickle, DOT_FLN, OMNI_URL, LEAPSEC_URL
    import numpy as n
    if sys.version_info[0]<3:
        import urllib as u
    else:
        import urllib.request as u
    #import time
    
    datadir = DOT_FLN+'/data'
    
    #leapsec_url ='ftp://maia.usno.navy.mil/ser7/tai-utc.dat'
    leapsec_fname = DOT_FLN+'/data/tai-utc.dat'

    # define location for getting omni
    #omni_url = 'ftp://virbo.org/QinDenton/hour/merged/latest/WGhour-latest.d.zip'
    omni_fname_zip = DOT_FLN+'/data/WGhour-latest.d.zip'
    omni_fname_pkl = DOT_FLN+'/data/omnidata.pkl'

    if all == True:
        omni = True
        leapsecs = True

    if omni == True:
        # retrieve omni, unzip and save as table
        print("Retrieving omni file ...")
        u.urlretrieve(OMNI_URL, omni_fname_zip)
        fh_zip = zipfile.ZipFile(omni_fname_zip)
        data = fh_zip.read(fh_zip.namelist()[0])
        A = n.array(data.split('\n'))
        #dd = data.split('\n')
        # save data as ascii file
        #fh = open(omni_fname_pkl, 'w')
        #fh.writelines(data)
        #fh.flush()
        #fh.close
        print("Now pickling (this will take a few minutes) ...")
        
        # create a keylist
        keys = A[0].split()
        keys.remove('8')
        keys.remove('6')
        keys[keys.index('status')] = '8_status'
        keys[keys.index('stat')] = '6_status'
        keys[keys.index('dst')] = 'Dst'
        keys[keys.index('kp')] = 'Kp'
        #keys[keys.index('Hr')] = 'Hr'
        keys[keys.index('V_SW')] = 'velo'
        keys[keys.index('Den_P')] = 'dens'
        keys[keys.index('Day')] = 'DOY'
        keys[keys.index('Year')] = 'Year'
        
        # remove keyword lines and empty lines as well
        idx = n.where(A != '')[0]
        A = A[idx[1:]]
        
        # put it into a 2D table
        #tab = n.zeros((len(A),len(keys)))
        tab = [['']*len(keys)]*len(A)
        stat8 = ['']*(len(A))
        stat6 = ['']*(len(A))
        for i in n.arange(len(A)):
            tab[i] = A[i].split()
            stat8[i] = A[i].split()[11]
            stat6[i] = A[i].split()[27]
        
        tab = n.array(tab, dtype='float') 
        # take out where Dst not available ( = 99999) or year == 0
        idx = n.where((tab[:,12] !=99.0) & (tab[:,0] != 0))[0]
        tab = tab[idx,:]
        stat8 = n.array(stat8)[idx]
        stat6 = n.array(stat6)[idx]
        
        omnidata = {} 
        # sort through and make an omni dictionary
        # extract keys from line above
        for ikey, i  in zip(keys,range(len(keys))):
            omnidata[ikey] = tab[:,i]
        
        # add TAI to omnidata
        nTAI = len(omnidata['DOY'])
        omnidata['UTC'] = ['']*nTAI
        omnidata['RDT'] = n.zeros(nTAI)
        
        
        #t1 = time.time()
        # add interpolation quality flags
        omnidata['Qbits'] = {}
        for ik, key in enumerate(['ByIMF', 'BzIMF', 'velo', 'dens', 'Pdyn', 'G1', 'G2', 'G3']):
            arr = n.array(list(n.array(stat8).tostring()), dtype=int).reshape((8,nTAI))
            omnidata['Qbits'][key] = arr[ik,:]
        for ik, key in enumerate(['W1', 'W2', 'W3', 'W4', 'W5', 'W6']):
            arr = n.array(list(n.array(stat6).tostring()), dtype=int).reshape((6,nTAI))
            omnidata['Qbits'][key] = arr[ik,:]
            
        #remove string status keys
        foo = omnidata.pop('6_status')
        foo = omnidata.pop('8_status')
        
        # add time information to omni pickle (long loop)
        for i in range(nTAI):
            year = int(omnidata['Year'][i])
            doy = int(omnidata['DOY'][i])
            month, day = doy2date(year,doy)
            UT_hr = omnidata['Hr'][i]
            hour, minute = divmod(UT_hr*60., 60)
            minute, second = divmod(minute*60., 60)  
            omnidata['UTC'][i] = datetime.datetime(year, month, day, int(hour), int(minute), int(second))
        
        omnidata['ticks'] = Ticktock(omnidata['UTC'], 'UTC')
        omnidata['RDT'] = omnidata['ticks'].RDT
        
        #t2 = time.time()
        #print t2-t1
        # save as pickle
        savepickle(omni_fname_pkl, omnidata)
            
        # delete left-overs
        os.remove(omni_fname_zip)

    if leapsecs == True:
        print("Retrieving leapseconds file ... ")
        u.urlretrieve(LEAPSEC_URL, leapsec_fname)
        
    return datadir

def windowMean(data, time=[], winsize=0, overlap=0, st_time=None):
    """Windowing mean function, window overlap is user defined

    @param data: 1D series of points;
    @type data: 
    @keyword time:  series of timestamps, optional (format as numeric or datetime)
       For non-overlapping windows set overlap to zero.
    @type time: 
    @keyword winsize: window size
    @type winsize: 
    @keyword overlap: amount of window overlap
    @type overlap: 
    @keyword st_time: for time-based averaging, a start-time other than the first 
        point can be specified
    @type st_time:
    @return: the windowed mean of the data, and an associated reference time vector
    
    @todo: Finish documentation

    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com

    @version: V1 pre 8-Jul-2010


    For non-overlapping windows set overlap to zero.
    e.g. (time-based averaging),
    
    >>> wsize, olap = datetime.timedelta(1), datetime.timedelta(0,3600)
    >>> outdata, outtime = windowmean(data, time, winsize=wsize, overlap=olap)
    
    in this example the window size is 1 day and the overlap is 1 hour.
    
    e.g. (pointwise averaging),
    
    >>> outdata, outtime = windowmean(data, winsize=10, overlap=9)
    
    where winsize and overlap are numeric,
    in this example the window size is 10 points and the overlap is 9 points.
    The output vectors start at winsize/2 and end at N-(winsize/2), the output time vector
    is basically a reference to the nth point in the original series.
    
    @note: This is a quick and dirty function - it is NOT optimized, at all.

    """
    import numpy as np
    import datetime as dt
    
    #check inputs and initialize
    #Set resolution to 1 if no times supplied
    if len(time) == 0:
        startpt, res = 0, 1
        time = list(range(len(data)))
        pts = True
    else:
        try:
            assert len(data) == len(time)
        except:
            return 'windowmean error: data and time must have same length'
        #First check if datetime objects
        if (type(time[0]) == dt.datetime):
            if not winsize:
                return 'windowmean error: winsize must be set for datetime input'
            else:
                try:
                    assert type(winsize) == dt.timedelta
                    assert type(overlap) == dt.timedelta
                except:
                    return 'windowmean error: winsize/overlap must be timedeltas'
            pts = False #force time-based averaging
        else:
            try:
                assert type(winsize) == dt.timedelta
                assert type(overlap) == dt.timedelta
            except:
                return 'windowmean error: winsize/overlap must be timedeltas'
            pts = False
            startpt = time[0]
    
    #now actually do windowing mean
    outdata, outtime = [], []
    data = np.array(data)
    if pts:
        #loop for fixed number of points in window
        if winsize % 1:
            winsize = round(winsize)
            print('windowmean error: non-integer windowsize, rounding to %d' \
            % winsize)
        if winsize < 1:
            winsize = 1
            print ('windowmean error: window length < 1, defaulting to 1')
        if overlap >= winsize:
            overlap = winsize - 1
            print ('''windowmean error: overlap longer than window, truncated to
            %d''' % overlap)
        lastpt = winsize-1 #set last point to end of window size
        while lastpt < len(data):
            datwin = np.ma.masked_where(np.isnan(data[startpt:startpt+winsize]), \
                data[startpt:startpt+winsize])
            getmean = np.mean(datwin.compressed()) #mean of window, excl. NaNs
            gettime = (time[startpt+winsize] - time[startpt])/2. \
                + time[startpt]#new timestamp
            startpt = startpt+winsize-overlap
            lastpt = startpt+winsize
            outdata.append(getmean) #construct output arrays
            outtime.append(gettime)
    else:
        #loop with time-based window
        lastpt = time[0] + winsize
        if st_time:
            startpt = st_time
        else:
            startpt = time[0]
        if overlap >= winsize:
            raise ValueError
        while lastpt < time[-1]:
            [getinds,dum] = tOverlap(time, [startpt,startpt+winsize])
            if getinds: #if not None
                getdata = np.ma.masked_where(np.isnan(data[getinds]), data[getinds])
                getmean = np.mean(getdata.compressed()) #find mean excluding NaNs
            else:
                getmean = np.nan
            gettime = startpt + winsize//2 #new timestamp -floordiv req'd with future division
            startpt = startpt + winsize - overlap #advance window start
            lastpt = startpt + winsize
            outdata.append(getmean) #construct output arrays
            outtime.append(gettime)
        
    return outdata, outtime
    
    
def medAbsDev(series):
    """Calculate median absolute deviation of a given input series
    
    Median absolute deviation (MAD) is a robust and resistant measure of
    the spread of a sample (same purpose as standard deviation). The
    MAD is preferred to the inter-quartile range as the inter-quartile
    range only shows 50% of the data whereas the MAD uses all data but
    remains robust and resistant. See e.g. Wilks, Statistical methods
    for the Atmospheric Sciences, 1995, Ch. 3.
    
    @param series: the input data series
    @type series: TODO
    @return: the median absolute deviation
    @rtype: float

    @todo: finish documentation
        
    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com
    
    @version: V1 pre 8-Jul-2010

    Find the median absolute deviation of a data set. Here we use the log-
    normal distribution fitted to the population of sawtooth intervals, see
    Morley and Henderson, Comment, Geophysical Research Letters, 2009.
    
    >>> data = numpy.random.lognormal(mean=5.1458, sigma=0.302313, size=30)
    >>> print data
    array([ 181.28078923,  131.18152745, ... , 141.15455416, 160.88972791])
    >>> toolbox.medabsdev(data)
    28.346646721370192

    @note: This implementation is robust to presence of NaNs
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

    @param x: List of x start and stop values
    @param y1: List of y lower bounds for fill
    @param y2: List of y upper bounds for fill
    @keyword face: color of the fill (default blue)
    @keyword alpha: alpha of the fill (default 0.5)

    @deprecated: Equivalent functionality to built-in matplotlib function fill_between

    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com

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
    
def binHisto(data, verbose=False):
    """
    Calculates bin width and number of bins for histogram using Freedman-Diaconis rule
    if rule fails, defaults to square-root method

    @param data: list/array of data values
    @type data: arraylike
    @keyword verbose: print out some more information
    @type verbose: boolean
    @return: calculated width of bins using F-D rule, number of bins (nearest integer) to use for histogram
    @rtype: tuple
    
    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov/morley_steve@hotmail.com
    @version: V1 SM
    @version: V2: 6Jul2010 Steve Morley added fallback rule

    >>> import numpy, spacepy
    >>> import matplotlib.pyplot as plt
    >>> data = numpy.random.randn(100)
    >>> binw, nbins = spacepy.toolbox.binHisto(data)
    >>> plt.hist(data, bins=nbins, histtype='step', normed=True)

    """
    import numpy as np
    from matplotlib.mlab import prctile
    pul = prctile(data, p=(25,75)) #get confidence interval
    ql, qu = pul[0], pul[1]
    iqr = qu-ql
    binw = 2.*iqr/(len(data)**(1./3.))
    nbins = round((np.max(data)-np.min(data))/binw)
    # if nbins is 0 or inf dont use the F-D rule just use sqrt(num) rule
    if nbins == 0 or nbins == np.inf:
        nbins = round(np.sqrt(len(data)))
        binw = len(data)/nbins
        if verbose:
            print("Used sqrt rule")
    else:
        if verbose:
            print("Used F-D rule")
    return (binw, nbins)
    
def smartTimeTicks(time):
    """Returns major ticks, minor ticks and format for time-based plots
    
    smartTimeTicks takes a list of datetime objects and uses the range
    to calculate the best tick spacing and format.  Returned to the user
    is a tuple containing the major tick locator, minor tick locator, and
    a format string -- all necessary to apply the ticks to an axis.

    It is suggested that, unless the user explicitly needs this info,
    to use the convenience function applySmartTimeTicks to place the
    ticks directly on a given axis.

    @param time: list of datetime objects
    @type time: list
    @return: tuple of Mtick - major ticks, mtick - minor ticks, fmt - format
    @rtype: tuple

    @author: Dan Welling
    @organization: Los Alamos National Lab
    @contact: dwelling@lanl.gov/dantwelling@gmail.com
    """
    from matplotlib.dates import (MinuteLocator, HourLocator, 
                                  DayLocator, DateFormatter)
    
    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if nHours < 1:
        Mtick = MinuteLocator(byminute = [0,15,30,45])
        mtick = MinuteLocator(byminute = list(range(60)), interval = 5)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick = MinuteLocator(byminute = [0,30])
        mtick = MinuteLocator(byminute = list(range(60)), interval = 10)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick = HourLocator(byhour = list(range(24)), interval = 2)
        mtick = MinuteLocator(byminute = [0,15,30,45])
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 24:
        Mtick = HourLocator(byhour = [0,3,6,9,12,15,18,21])
        mtick = HourLocator(byhour = list(range(24)))
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick = HourLocator(byhour = [0,6,12,18])
        mtick = HourLocator(byhour = list(range(24)))
        fmt = DateFormatter('%H:%M UT')
    else:
        Mtick = DayLocator(bymonthday = list(range(1,32)))
        mtick = HourLocator(byhour = [0,6,12,18])
        fmt = DateFormatter('%d %b')

    return (Mtick, mtick, fmt) 

def applySmartTimeTicks(ax, time, dolimit = True):
    '''
    Given an axis 'ax' and a list/array of datetime objects, 'time', 
    use the smartTimeTicks function to build smart time ticks and
    then immediately apply them to the given axis.  The first and
    last elements of the time list will be used as bounds for the
    x-axis range.

    The range of the 'time' input value will be used to set the limits
    of the x-axis as well.  Set kwarg 'dolimit' to False to override 
    this behavior.

    @param ax: A matplotlib Axis object.
    @type ax: matplotlib.pyplot.Axes
    @param time: list of datetime objects
    @type time: list
    @keyword dolimit: The range of the 'time' input value will be used to set the limits
        of the x-axis as well. Setting this overrides this behavior.
    @type dolimit: bool
    @return: None
    @rtype: None

    @author: Dan Welling
    @organization: Los Alamos National Lab
    @contact: dwelling@lanl.gov/dantwelling@gmail.com
    '''

    Mtick, mtick, fmt = smartTimeTicks(time)
    ax.xaxis.set_major_locator(Mtick)
    ax.xaxis.set_minor_locator(mtick)
    ax.xaxis.set_major_formatter(fmt)
    if dolimit:
        ax.set_xlim([time[0], time[-1]])


def logspace(min, max, num, **kwargs):
    """Returns log spaced bins.  Same as numpy logspace except the min and max are the ,min and max 
    not log10(min) and log10(max)

    @param min: minimum value 
    @param max: maximum value
    @param num: number of log spaced bins
    @return: log spaced bins from min to max in a numpy array
    @rtype: numpy.array

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)

    >>> logspace(1, 100, 5)
    Out[2]: array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])

    """
    
    from numpy import logspace, log10
    return logspace(log10(min), log10(max), num, **kwargs)

def arraybin(array, bins):
    """Split a sequence into subsequences based on value.
    
    Given a sequence of values and a sequence of values representing the
    division between bins, return the indices grouped by bin.

    @param array: the input sequence to slice, must be sorted in ascending
                  order
    @type array: sequence
    @param bins: dividing lines between bins. Number of bins is len(bins)+1,
                 value that exactly equal a dividing value are assigned
                 to the higher bin
    @type bins: sequence
    @return: indices for each bin
    @rtype: list of lists

    @author: Brian Larsen, Jonathan Niehof
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov, jniehof@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)
    @version: V2: 07-Dec-2010 (JTN)
   
    >>> arraybin(range(10), [4.2])
    Out[4]: [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
    """
    value = 0
    bin_it = (i for i in range(len(array)) if array[i] >= value)
    splits = [next(bin_it, len(array)) for value in bins]
    return [list(range(start_idx, stop_idx)) for (start_idx, stop_idx)
            in zip([0] + splits, splits + [len(array)])]

def mlt2rad(mlt, midnight = False):
    """
    Convert mlt values to radians for polar plotting
    transform mlt angles to radians from -pi to pi
    referenced from noon by default

    @param mlt:  array of mlt values
    @type mlt: numpy.array
    @keyword midnight: reference to midnight instead of noon
    @type midnight: boolean
    @return: array of radians
    @rtype: numpy.array

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)
    
    >>> mlt2rad(array([3,6,9,14,22]))
    Out[9]: array([-2.35619449, -1.57079633, -0.78539816,  0.52359878,  2.61799388])
    """
    import numpy as np
    if midnight:
        mlt_arr =  mlt + 12
    else:
        mlt_arr = mlt
    rad_arr=(mlt_arr-12)*np.pi/12
    return rad_arr

def rad2mlt(rad, midnight=False):
    """
    Convert radian values to mlt 
    transform radians from -pi to pi to mlt
    referenced from noon by default

    @param rad:  array of radian values
    @type rad: numpy.array
    @keyword midnight: reference to midnight instead of noon
    @type midnight: boolean
    @return: array of mlt values
    @rtype: numpy.array

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)

    >>> rad2mlt(array([0,pi, pi/2.]))
    Out[8]: array([ 12.,  24.,  18.])
    """
    import numpy as np
    if midnight:
        rad_arr = rad + np.pi
    else:
        rad_arr = rad
    mlt_arr=rad_arr*(12/np.pi) + 12
    return mlt_arr


def leap_year(year, numdays=False, nobool=False):
    """
    return an array of boolean leap year, 
    a lot faster than the mod method that is normally seen

    @param year: <iterable> of years
    @type year: iterable
    @keyword numdays: optionally return the number of days in the year
    @type numdays: boolean
    @return: an array of boolean leap year, or array of number of days
    @rtype: numpy.array

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)
    @version: V2: 22-Nov-2010 (BAL) accepts lists not just arrays
    @version: V3: 07-Dec-2010 (BAL) cleanup
    >>> leap_year(numpy.arange(15)+1998)
    Out[10]: 
    array([False, False,  True, False, False, False,  True, False, False,
    ... False,  True, False, False, False,  True], dtype=bool)
    """
    mask400 = [(val % 400) == 0 for val in year]   # this is a leap year
    mask100 = [(val % 100) == 0 for val in year ]   # these are not leap years
    mask4   = [(val % 4) == 0 for val in year ]   # this is a leap year
    if numdays or nobool:
        if numdays:
            numdays=365
        else:
            numdays = 0
        return [numdays + ((val[0] | val[2]) & (~val[1] | val[0])) for val in zip(mask400, mask100, mask4)]
    else:
        return [bool(((val[0] | val[2]) & (~val[1] | val[0]))) for val in zip(mask400, mask100, mask4)]
        
leapyear = leap_year

    
def pmm(a, *b):
    """
    print min and max of input arrays

    @param a: input array
    @type a: numpy.array
    @keyword *b: some additional number of arrays
    @type *b: numpy.array
    @return: list of min, max for each array
    @rtype: list

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)
    @version: V2: 04-Sep-2010 (BAL) - changed to min() form a.min() works on lists
    
    >>> pmm(arange(10), arange(10)+3)
    Out[12]: [(0, 9), (3, 12)]
    """
    import numpy as np
    ans= [[ np.min(a), np.max(a) ]]
    for val in b:
        ans.append( [np.min(val), np.max(val)] )
    return ans

def timestamp(position=[1.003, 0.01], size='xx-small', draw=True, **kwargs):
    """
    print a timestamp on the current plot, vertical lower right 

    @keyword position: position for the timestamp
    @type position: list
    @keyword size: text size
    @type size: string
    @keyword draw: call draw to make sure it appears
    @type draw: boolean
    @keyword kwargs: other keywords to axis.annotate
    @type kwargs: keywords

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)

    >>> plot(arange(11))
    Out[13]: [<matplotlib.lines.Line2D object at 0x49072b0>]
    timestamp()
    """
    from datetime import datetime
    from matplotlib.pyplot import gca, draw
    now = datetime.now()
    strnow = now.strftime("%d%b%Y %H:%M")
    ax=gca()
    ax.annotate(strnow, position, xycoords='axes fraction', rotation='vertical', size=size, va='bottom',  **kwargs)
    if draw:
        draw()
        
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".

    @param question: string that is the question to ask
    @type question: string
    @keyword default: the default answer (yes)
    @type default: string
    @return: answer ('yes' or 'no')
    @rtype: string

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 14-Jun-2010 (BAL)

    >>> query_yes_no('Ready to go?')
    Ready to go? [Y/n] y
    Out[17]: 'yes'
    """
    import sys
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        if sys.version_info[0]==2:
            choice = raw_input().lower()
        elif sys.version_info[0]>2:
            choice = input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")

def interpol(newx, x, y, wrap=None, **kwargs):
    """1-D linear interpolation with interpolation of hours/longitude
    
    @return: interpolated data values for new abscissa values
    @rtype: numpy.masked_array
    
    @keyword hour: wraps interpolation at 24 (e.g. for local times)
    @type hour: string
    @keyword long: wraps interpolation at 360 (e.g. longitude)
    @type long: string
    @keyword sect: wraps interpolation based on user specified max    
    i.e. days=True is equivalent to sect=24
    @type sect: integer
    
    @author: Steve Morley
    @organization: Los Alamos National Lab
    @contact: smorley@lanl.gov
    """
    import scipy as sci
    import numpy as np
    
    if 'baddata' in kwargs:
        x = np.ma.masked_where(y==kwargs['baddata'], x)
        y = np.ma.masked_where(y==kwargs['baddata'], y)
        kwargs.__delitem__('baddata')
    else:
        tst = np.ma.core.MaskedArray
        if type(x)!=tst or type(y)!=tst or type(newx)!=tst:
            x = np.ma.masked_array(x)
            y = np.ma.masked_array(y)
            newx = np.ma.masked_array(newx)
    
    def wrap_interp(xout, xin, yin, sect):
        dpsect=360/sect
        yc = np.cos(np.deg2rad(y*dpsect))
        ys = np.sin(np.deg2rad(y*dpsect))
        new_yc = sci.interp(newx, x.compressed(), yc.compressed(), **kwargs)
        new_ys = sci.interp(newx, x.compressed(), ys.compressed(), **kwargs)
        try:
            new_bad = sci.interp(newx, x, y.mask)
        except ValueError:
            new_bad = np.zeros((len(newx)))
        newy = np.rad2deg(np.arctan(new_ys/new_yc))/dpsect
        #1st quadrant is O.K
        #2nd quadrant
        idx = [n for n in range(len(new_yc)) if new_yc[n]<0 and new_ys[n]>0]
        newy[idx] = sect/2 + newy[idx]
        #print('Sector 2 inds: %s' % idx)
        #3rd quadrant
        idx = [n for n in range(len(new_yc)) if new_yc[n]<0 and new_ys[n]<0]
        newy[idx] = sect/2 + newy[idx]
        #print('Sector 3 inds: %s' % idx)
        #4th quadrant
        idx = [n for n in range(len(new_yc)) if new_yc[n]>0 and new_ys[n]<0]
        newy[idx] = sect + newy[idx]
        #print('Sector 4 inds: %s' % idx)
        new_bad = np.ma.make_mask(new_bad)
        newy = np.ma.masked_array(newy, mask=new_bad)
        return newy
    
    if wrap=='hour':
        newy = wrap_interp(newx, x.compressed(), y.compressed(), 24)
    elif wrap=='lon':
        newy = wrap_interp(newx, x.compressed(), y.compressed(), 360)
    elif type(wrap)==int:
        newy = wrap_interp(newx, x.compressed(), y.compressed(), wrap)
    else:
        newy = sci.interp(newx, x.compressed(), y.compressed(), **kwargs)
        try:
            new_bad = sci.interp(newx, x, y.mask)
            new_bad = np.ma.make_mask(new_bad)
            newy = np.ma.masked_array(newy, mask=new_bad)
        except:
            pass
    return newy

# -----------------------------------------------


def normalize(vec):
    """
    Given an input vector normalize the vector

    @param vec: input vector to normalize
    @type vec: listlike
    @return: normalized vector
    @rtype: listlike
    """
    import numpy as N
    # check to see if vec is numpy array, this is fastest
    if isinstance(vec, N.ndarray):
        out = (vec - vec.min())/N.ptp(vec)
    else:
        vecmin = N.min(vec)
        ptp = N.ptp(vec)
        out = [(val -  vecmin)/ptp for val in vec]
    return out


def listUniq(inVal):
    """
    Given an input iterable (list, deque) return a list of the unique elements.
    Maintains order (keeps the first of non-unique elements

    @param inVal: Input iterable
    @type inVal: iterable
    @return: list of unique elements from iterable
    @rtype: list

    @author: Brian Larsen
    @organization: Los Alamos National Lab
    @contact: balarsen@lanl.gov
    
    @version: V1: 16-Nov-2010 (BAL)
    """
    seen = set()
    return [ x for x in inVal if x not in seen and not seen.add(x)]
