#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Toolbox of various functions and generic utilities.

"""
from __future__ import division
#from spacepy import help # not needed b/c help() is already defined in __init__.py 
__version__ = "$Revision: 1.20 $, $Date: 2010/06/15 16:30:37 $"
__author__ = 'S. Morley and J. Koller'


def tOverlap(ts1,ts2):
    """Finds the overlapping elements in two lists of datetime objects
    
    Returns:
    ========
     - indices of 1 within interval of 2, & vice versa
    
    Example:
    ========
     - Given two series of datetime objects, event_dates and omni['Time']:
    
    >>> import spacepy.toolbox as tb
    >>> [einds,oinds] = tb.tOverlap(event_dates, omni['Time'])
    >>> omni_time = omni['Time'][oinds[0]:oinds[-1]+1]
    >>> print omni_time
    [datetime.datetime(2007, 5, 5, 17, 57, 30), datetime.datetime(2007, 5, 5, 18, 2, 30),
    ... , datetime.datetime(2007, 5, 10, 4, 57, 30)]
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
    
    Modifications:
    ==============
    09-Jun-2010: Total rewrite for clarity, efficiency and brevity.
    """
    
    #elements of ts2 within bounds of ts1
    t_lower, t_upper = min(ts1), max(ts1)
    bool1 = [v > t_lower for v in ts2]
    bool2 = [v < t_upper for v in ts2]
    mask2in1 = [b1 and b2 for b1,b2 in zip(bool1,bool2)]
    inds2in1 = [i for i, val in enumerate(mask2in1) if val==True]
    
    #elements of ts1 within bounds of ts2
    t_lower, t_upper = min(ts2), max(ts2)
    bool1 = [v > t_lower for v in ts1]
    bool2 = [v < t_upper for v in ts1]
    mask1in2 = [b1 and b2 for b1,b2 in zip(bool1,bool2)]
    inds1in2 = [i for i, val in enumerate(mask1in2) if val==True]
    
    if len(inds2in1) == 0:
        inds2in1 = None
    if len(inds1in2) == 0:
        inds1in2 = None
    
    return inds1in2, inds2in1
    

def tCommon(ts1, ts2, mask_only=True):
    """Finds the elements in a list of datetime objects present in another
    
    Returns:
    ========
     - Two element tuple of truth tables (of 1 present in 2, & vice versa)
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
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

def loadpickle(fln):
    """
    load a pickle and return content as dictionary

    Input:
    ======
        - fln (string) : filename

    Returns:
    ========
        - d (dictionary) : dictionary with content from file

    Example:
    ========

    >>> d = loadpickle('test.pbin')

    See also:
    =========
    savepickle

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 20-Jan-2010
    """
    import cPickle

    fh = open(fln, 'rb')
    content = cPickle.load(fh)
    fh.close()

    return content


# -----------------------------------------------
def savepickle(fln, dict):
    """
    save dictionary variable dict to a pickle with filename fln 
    Author: Josef Koller, jkoller@lanl.gov

    Inputs:
    =======
        - fln (string) : filename
        - dict (dictionary) : container with stuff

    Example:
    ========
    >>> d = {'grade':[1,2,3], 'name':['Mary', 'John', 'Chris']}
    >>> savepickle('test.pbin', d)

    See also:
    =========
    loadpickle

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)

    Version:
    ========
    V1: 20-Jan-2010
    """

    import cPickle

    fh = open(fln, 'wb')
    cPickle.dump(dict, fh, 2) # 2 ... fast binary
    fh.close()

    return
# -----------------------------------------------
def assemble(fln_pattern, outfln):
    """
    assembles all pickled files matching fln_pattern into single file and 
    save as outfln. Pattern may contain simple shell-style wildcards *? a la fnmatch
    file will be assembled along time axis TAI, etc in dictionary

    Inputs:
    =======
        - fln_pattern (string) : pattern to match filenames
        - outfln (string) : filename to save combined files to

    Outputs:
    ========
        - dcomb (dict) : dictionary with combined values

    Example:
    ========
    
    >>> assemble('input_files_*.pbin', 'combined_input.pbin')
    adding input_files_2001.pbin
    adding input_files_2002.pbin
    adding input_files_2004.pbin
    writing: combined_input.pbin

    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 20-Jan-2010
    """

    import glob
    import numpy as n
    import time as t

    filelist = glob.glob(fln_pattern)

    # read all files
    d = {}
    for fln in filelist:
       print "adding ", fln
       d[fln] = loadpickle(fln)
    
    # combine them
    dcomb = d[filelist[0]]  # copy the first file over
    for fln in filelist[1:]:
       TAIcount = len(d[fln]['TAI'])
       for key in d[fln].keys():
          #print fln, key
          dim = n.array(n.shape(d[fln][key]))
          ax = n.where(dim==TAIcount)[0]
          if len(ax) == 1: # then match with TAI length is given (jump over otherwise)
             dcomb[key] = n.append(dcomb[key], d[fln][key], axis=ax)
             
    # sort in time
    idx = n.argsort(dcomb['TAI'])
    TAIcount = len(dcomb['TAI'])
    for key in dcomb.keys():
       dim = n.array(n.shape(dcomb[key]))
       ax = n.where(dim==TAIcount)[0]
       if len(ax) == 1: # then match with length of TAI
          dcomb[key] = dcomb[key][idx] # resort

    print '\n writing: ', outfln
    savepickle(outfln, dcomb)
    
    return dcomb

# -----------------------------------------------
def feq(x,y, precision=0.0000005):
    """
    compare two floating point values if they are equal
    after: http://www.lahey.com/float.htm
    
    further info at::
        http://docs.python.org/tut/node16.html
        http://www.velocityreviews.com/forums/t351983-precision-for-equality-of-two-floats.html
        http://www.boost.org/libs/test/doc/components/test_tools/floating_point_comparison.html
        http://howto.wikia.com/wiki/Howto_compare_float_numbers_in_the_C_programming_language

    Input:
    ======
        - x (float) : a number
        - y (array of floats) : float value or array of floats

    Returns:
    ========
        - boolean value: True or False

    Example:
    ========
    >>> index = where( feq(Lpos,Lgrid) ) # use float point comparison
    
    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 20-Jan-2010
    V2: 18-May-2010: User-specified precision added
    """

    boolean = abs(x-y) <= (abs(x+y)*precision)

    return boolean


# -----------------------------------------------
def dictree(in_dict, verbose=False, spaces=None, levels=True):
    """ pretty print a dictionary tree

    Input:
    ======
        - in_dict (dictionary) : a complex dictionary (with substructures)
        - spaces (string) : string will added for every line
        - levels (default False) : number of levels to recurse through (True means all)
        - boolean value: True or False

    Example:
    ========
    >>> d = {'grade':{'level1':[4,5,6], 'level2':[2,3,4]}, 'name':['Mary', 'John', 'Chris']}
    >>> dictree(d)
    +
    |____grade
        |____level1
        |____level2
    |____name

    
    Author:
    =======
    Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

    Version:
    ========
    V1: 20-Jan-2010
    V1.1: 24-Feb-2010 S. Morley, added verbose option
    v1.2: 17-May-2010 S. Morley, added levels option
    """

    import numpy as n

    if not spaces:
        spaces = ''
        print '+'
    
    if levels:
        try:
            assert levels is True
        except AssertionError:
            levels -= 1
            if levels == 0:
                levels = None

    for key in n.sort(in_dict.keys()):
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
                print spaces + '|____' + key + ' ('+ typestr + dimstr + ')'
            else:
                print spaces + '|____' + key
        else:
            print spaces + '|____' + key
        if type(in_dict[key]) == dict and levels:
            dictree(in_dict[key], spaces = spaces+ '     ', verbose = verbose, levels = levels)

    return None

# -----------------------------------------------
def printfig(fignum, saveonly=False, pngonly=False, clean=False):
   """save current figure to file and call lpr (print).
   
   This routine will create a total of 3 files (png, ps and c.png) in the 
   current working directory with a sequence number attached. Also, a time 
   stamp and the location of the file will be imprinted on the figure. The 
   file ending with c.png is clean and no directory or time stamp are 
   attached (good for powerpoint presentations).

   Input:
   ======
       - fignum (integer or array/list of integer) : matplotlib figure number
       - optional 
           - saveonly (boolean) : True (don't print and save only to file)
                                 False (print and save)
           - pngonly (boolean) : True (only save png files and print png directly)
           	                     False (print ps file, and generate png, ps; can be slow)
           - clean (boolean) : True (print and save only clean files without directory info)
                               False (print and save directory location as well)

   Example:
   ========
   >>> pylab.plot([1,2,3],[2,3,2])
   >>> spacepy.printfig(1)

   Author:
   =======
   Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

   Version:
   ========
   V1: 20-Jan-2010
   V2: 19-Feb-2010: added pngonly and clean options, array/list support (JK)
   """

   import pylab, os, sys, re, glob, datetime

   try:
   		nfigs = len(fignum)
   except:
   		nfigs = 1
   		fignum = [fignum]
   
   for ifig in fignum:
   		# active this figure
   		pylab.figure(ifig)

   		# create a filename for the figure
   		cwd = os.getcwd()
   		num = len(glob.glob('*.png'))
   		fln = cwd+'/figure_'+str(num)
   		# truncate fln if too long
   		if len(fln) > 60: 
   			flnstamp = '[...]'+fln[-60:]
   		else:
   			flnstamp = fln

	    # save a clean figure without timestamps
   		if clean == True:
	   		pylab.savefig(fln+'_clean.png')

   		timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")
   		# add the filename to the figure for reference
   		pylab.figtext(0.01, 0.01, timestamp+flnstamp+'.png', rotation='vertical', size=8)

   		# now save the figure to this filename
   		if pngonly == False:
	   		pylab.savefig(fln+'.ps')
	   
   		pylab.savefig(fln+'.png')

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

    Input:
    ======
       - all (bool) : if True, update all of them
       - omni (bool) : if True. update only onmi
       - leapsecs (bool) : if True, update only leapseconds


     Example:
     ========
     >>> update(omni=True)

     Author:
     =======
     Josef Koller, Los Alamos National Lab, jkoller@lanl.gov

     Version:
     ========
     V1: 20-Jan-2010
     V1.1: 24-May-2010 Minor modification to return data directory (BAL)
     V1.2: 11-Jun-2010 moved pickle_omni in here and added Qbits (JK)
     """

    import urllib as u
    import os
    import zipfile
    import re
    import datetime
    import spacepy.time as st
    from spacepy import savepickle, DOT_FLN, OMNI_URL, LEAPSEC_URL
    import numpy as n
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
        print "Retrieving omni file ..."
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
        print "Now pickling (this will take a few minutes) ..."
        
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
        tab = n.zeros((len(A),len(keys)))
        stat8 = ['']*(len(A))
        stat6 = ['']*(len(A))
        for i in n.arange(len(A)):
            tab[i,:] = A[i].split()
            stat8[i] = A[i].split()[11]
            stat6[i] = A[i].split()[27]
    
        tab = n.reshape(tab, n.shape(tab))
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
            month, day = st.doy2date(year,doy)
            UT_hr = omnidata['Hr'][i]
            hour, minute = divmod(UT_hr*60., 60)
            minute, second = divmod(minute*60., 60)  
            omnidata['UTC'][i] = datetime.datetime(year, month, day, int(hour), int(minute), int(second))
        
        omnidata['ticktock'] = st.Ticktock(omnidata['UTC'], 'UTC')
        omnidata['RDT'] = omnidata['ticktock'].RDT
        
        #t2 = time.time()
        #print t2-t1
        # save as pickle
        savepickle(omni_fname_pkl, omnidata)
            
        # delete left-overs
        os.remove(omni_fname_zip)

    if leapsecs == True:
        print "Retrieving leapseconds file ... "
        u.urlretrieve(LEAPSEC_URL, leapsec_fname)
        
    return datadir

def windowMean(data, time=[], winsize=0, overlap=0, st_time=None):
    """Windowing mean function, window overlap is user defined
    
    Inputs:
    data - 1D series of points;
    time - series of timestamps, optional (format as numeric or datetime);
    For non-overlapping windows set overlap to zero.
    e.g.,
    
    >>> wsize, olap = datetime.timedelta(1), datetime.timedelta(0,3600)
    
    >>> outdata, outtime = windowmean(data, time, winsize=wsize, overlap=olap)
    
    where the time, winsize and overlap are either numberic or datetime objects,
    in this example the window size is 1 day and the overlap is 1 hour.
    
    Caveats: This is a quick and dirty function - it is NOT optimised, at all.
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
    """
    import numpy as np
    import datetime as dt
    from matplotlib.dates import date2num,num2date
    
    #check inputs and initialize
    #Set resolution to 1 if no times supplied
    if len(time)==0:
        startpt, res = 0, 1
        time = range(len(data))
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
        if st_time:
            startpt = st_time
        else:
            startpt = time[0]
        if overlap >= winsize:
            raise ValueError
        while lastpt < time[-1]:
            [getinds,dum] = tOverlap(time, [startpt,startpt+winsize])
            if getinds: #if not None
                getdata = np.ma.masked_where(np.isnan(data[getinds]),data[getinds])
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
    >>> toolbox.medabsdev(data)
    28.346646721370192
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
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
    
    Equivalent functionality to built-in matplotlib function fill_between
    
    >>> poly0c = makePoly(x, ci_low, ci_high, face='red', alpha=0.8)
    >>> ax0.add_patch(poly0qc)
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
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
    
    Inputs:
    =======
    data - list/array of data values
    
    Outputs:
    ========
    binw - calculated width of bins using F-D rule
    nbins - number of bins (nearest integer) to use for histogram
    
    Example:
    ========
    >>> import numpy, spacepy
    >>> import matplotlib.pyplot as plt
    >>> data = numpy.random.randn(100)
    >>> binw, nbins = spacepy.toolbox.binHisto(data)
    >>> plt.hist(data, bins=nbins, histtype='step', normed=True)
    
    Author:
    =======
    Steve Morley, Los Alamos National Lab, smorley@lanl.gov/morley_steve@hotmail.com
    """
    import numpy as np
    from matplotlib.mlab import prctile
    pul = prctile(data, p=(25,75)) #get confidence interval
    ql, qu = pul[0], pul[1]
    iqr = qu-ql
    binw = 2.*iqr/(len(data)**(1./3.))
    nbins = round((np.max(data)-np.min(data))/binw)
    
    return (binw, nbins)
    
def smartTimeTicks(time):
    """Returns major ticks, minor ticks and format for time-based plots
    
    smartTimeTicks takes a list of datetime objects and uses the range
    to calculate the best tick spacing and format.
    
    Inputs:
    =======
    time - list of datetime objects
    
    Outputs:
    ========
    Mtick - major ticks
    mtick - minor ticks
    fmt - format
    
    Example:
    ========
    ? Meh.
    
    Author:
    =======
    Dan Welling, Los Alamos National Lab, dwelling@lanl.gov/dantwelling@gmail.com
    """
    from matplotlib.dates import MinuteLocator, HourLocator, DayLocator, DateFormatter
    
    deltaT = time[-1] - time[0]
    nHours = deltaT.days * 24.0 + deltaT.seconds/3600.0
    if nHours < 1:
        Mtick=MinuteLocator(byminute=[0,15,30,45])
        mtick=MinuteLocator(byminute=range(60), interval=5)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 4:
        Mtick=MinuteLocator(byminute=[0,30])
        mtick=MinuteLocator(byminute=range(60), interval=10)
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 12:
        Mtick=HourLocator(byhour=range(24), interval=2)
        mtick=MinuteLocator(byminute=[0,15,30,45])
        fmt = DateFormatter('%H:%M UT')
    elif nHours < 48:
        Mtick=HourLocator(byhour=[0,6,12,18])
        mtick=HourLocator(byhour=range(24))
        fmt = DateFormatter('%H:%M UT')
    else:
        Mtick=DayLocator(bymonthday=range(1,32))
        mtick=HourLocator(byhour=[0,6,12,18])
        fmt = DateFormatter('%d %b')

    return (Mtick, mtick, fmt) 


#function aliases to maintain compatibility with existing scripts using Toolbox
t_common = tCommon
t_overlap = tOverlap
smart_timeticks = smartTimeTicks


def logspace(min, max, num, **kwargs):
    """Returns log spaced bins.  Same as numpy logspace except the min and max are the ,min and max 
    not log10(min) and log10(max)

    logspace(min, max, num)

    Inputs:
    =======
    min - minimum value 
    max - maximum value
    num - number of log spaced bins

    Outputs:
    ========
    num log spaced bins from min to max in a numpy array
    
    Example:
    ========
    logspace(1, 100, 5)
    Out[2]: array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    
    from numpy import logspace, log10
    return logspace(log10(min), log10(max), num, **kwargs)

def arraybin(array, bins):
    """
    arraybin(array, bins)

    Given an array and a set of bins return the indices that are less than the 
    smallest bin between each set and larger than the largest bin
    bins should be sorted

    Inputs:
    =======
    array - the input array to slice
    bins - the bins to slice along (may be array or list)

    Outputs:
    ========
    list of indices 
    first element is less than fisrt bin and last bin is larger than last bin

    Example:
    ========
    arraybin(arange(10), [4.2])
    Out[4]: [(array([0, 1, 2, 3, 4]),), (array([5, 6, 7, 8, 9]),)]
 
    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    from pylab import where
    bin_ind=[]
    for i,b in enumerate(bins):
        if i == 0:
            ind=where( (array<bins[i]) )
        else:
            ind=where( (array>=bins[i-1]) & (array<bins[i]) )
        bin_ind.append(ind)
    ind=where( (array>=bins[i]) )        
    bin_ind.append(ind)
    return bin_ind

def mlt2rad(mlt, midnight=False):
    """
    mlt2rad(mlt, midnight=False)

    Convert mlt values to radians for polar plotting
    transform mlt angles to radians from -pi to pi
    referenced from noon by default

    Inputs:
    =======
    mlt - array of mlt values
    midnight=False - reference to midnioght instead of noon

    Outputs:
    ========
    array of radians 

    Example:
    ========
    mlt2rad(array([3,6,9,14,22]))
    Out[9]: array([-2.35619449, -1.57079633, -0.78539816,  0.52359878,  2.61799388])

    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    import numpy as np
    if midnight:
        mlt_arr = mlt + 12
    else:
        mlt_arr = mlt
    rad_arr=(mlt_arr-12)*np.pi/12
    return rad_arr

def rad2mlt(rad, midnight=False):
    """
    rad2mlt(rad, midnight=False)

    Convert radian values to mlt 
    transform radians from -pi to pi to mlt
    referenced from noon by default

    Inputs:
    =======
    rad - array of rad values
    midnight=False - reference to midnioght instead of noon

    Outputs:
    ========
    array of mlt 
 
    Example:
    ========
    rad2mlt(array([0,pi, pi/2.]))
    Out[8]: array([ 12.,  24.,  18.])

    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
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
    leap_year(year, numdays=False)

    return an array of boolean leap year, 
    a lot faster than the mod method that is normally seen
    
    Inputs:
    =======
    year - array of years
    numdays=False - optionally return the number of days in the year
    
    Outputs:
    ========
    an array of boolean leap year, or array of number of days

    Example:
    ========
    leap_year(arange(15)+1998)
    Out[10]: 
    array([False, False,  True, False, False, False,  True, False, False,
    ... False,  True, False, False, False,  True], dtype=bool)

    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    mask400 = (year % 400) == 0   # this is a leap year
    mask100 = (year % 100) == 0   # these are not leap years
    mask4   = (year %   4) == 0   # this is a leap year
    if numdays:
        numdays=365
        return numdays + ((mask400 | mask4) & (~mask100 | mask400))
    else:
        if nobool:
            return 0 + ((mask400 | mask4) & (~mask100 | mask400))
        return ((mask400 | mask4) & (~mask100 | mask400))
    
def pmm(a, *b):
    """
    pmm(a, *b)

    print min and max of input arrays

    Inputs:
    =======
    a - input array
    *b - some additional number of arrays

    Outputs:
    ========
    list of min, max for each array

    Example:
    ======== 
    pmm(arange(10), arange(10)+3)
    Out[12]: [(0, 9), (3, 12)]

    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    import numpy as np
    ans=[]
    ans.append( (a.min(), a.max()) )
    for val in b:
        ans.append( (val.min(), val.max()) )
    return ans

def timestamp(position=[1.003, 0.01], size='xx-small', draw=True, **kwargs):
    """
    timestamp(position=[1., 0.01], size='xx-small', **kwargs)

    print a timestamp on the current plot, vertical lower right 

    Inputs:
    =======
    (all optional)
    position - position for the timestamp
    size - text size
    draw - call draw to make sure it appears
    kwargs - other keywords to axis.annotate

    Outputs:
    ========
    timestamp written to the current plot

    Example:
    ======== 
    plot(arange(11))
    Out[13]: [<matplotlib.lines.Line2D object at 0x49072b0>]
    timestamp()

    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
    """
    from datetime import datetime
    from matplotlib.pyplot import annotate, gca, draw
    now = datetime.now()
    strnow = now.strftime("%d%b%Y %H:%M")
    ax=gca()
    ax.annotate(strnow, position, xycoords='axes fraction', rotation='vertical', size=size, **kwargs)
    if draw:
        draw()
        
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.
    
    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".

    Inputs:
    =======
    question - string that is the question to ask
    default - the default answer (yes)

    Outputs:
    ========
    answer ('yes' or 'no')
    
    Example:
    ======== 
    query_yes_no('Ready to go?')
    Ready to go? [Y/n] y
    Out[17]: 'yes'


    Author:
    =======
    Brian Larsen, Los Alamos National Lab, balarsen@lanl.gov
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
        choice = raw_input().lower()
        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


# -----------------------------------------------
def test():
    """
    test all spacepy routines

    Returns:
    ========
        - nFAIL (int) : number of failures

    Example:
    ========

    >>> test()
    test_ticktock: PASSED TEST 1
    test_ticktock: PASSED TEST 2
    0

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)

    Version:
    ========
    V1: 24-Jan-2010
    """
    
    import os
    import spacepy.toolbox as tb
    import pylab as p
    nFAIL = 0

    # savepickle
    D = {}
    D['names'] = ['John', 'Joe', 'Joyce']
    D['TAI'] = [1,2,3]
    tb.savepickle('test_pickle_1.pbin', D)

    # loadpickle
    DD = tb.loadpickle('test_pickle_1.pbin')
    if DD != D:
        print "test_toolbox: FAILED TEST load/save pickle"
        nFAIL =+ 1

    # assemble
    D2 = {}
    D2['names'] = ['Harry', 'Hans']
    D2['TAI'] = [4,5] 
    tb.savepickle('test_pickle_2.pbin', D2)
    try:
        Dcomb = tb.assemble('test_pickle_*.pbin', 'test_pickle_assembled.pbin')
        os.remove('test_pickle_1.pbin')
        os.remove('test_pickle_2.pbin')
        os.remove('test_pickle_assembled.pbin')
    except:
        print "test_toolbox: FAILED TEST assemble"
        nFAIL =+ 1
    
    # feq
    if not tb.feq(1.2, 1.2):
        print "test_toolbox: FAILED TEST assemble"
        nFAIL =+ 1

    # dictree
    try:
        tb.dictree(Dcomb)
    except:
        print "test_toolbox: FAILED TEST dictionary tree"
        nFAIL =+ 1

    # printfig
    p.plot(D['TAI'], D['TAI'])
    try:
        tb.printfig(1, saveonly=True)
        p.close(1)
        os.remove('figure_0.png')
        os.remove('figure_0.ps')
    except:
        print "testing toolbox: FAILED TEST printfig"
        nFAIL =+ 1

    # update
    #try:
    #    tb.update()
    #except:
    #    print "testing toolbox: FAILED TEST update"
    #    nFAIL =+ 1

    return nFAIL






