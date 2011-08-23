#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Toolbox of various functions and generic utilities.

Authors: Steve Morley, Jon Niehof, Brian Larsen, Josef Koller, Dan Welling
Institution: Los Alamos National Laboratory
Contact: smorley@lanl.gov, jniehof@lanl.gov, balarsen@lanl.gov, jkoller@lanl.gov, dwelling@lanl.gov
Los Alamos National Laboratory

Copyright ©2010 Los Alamos National Security, LLC.
"""
from __future__ import division

import math
import os
import sys
import zipfile
import datetime
import glob

import numpy as np

try:
    from spacepy import help
except ImportError:
    pass
except:
    pass

def tOverlap(ts1, ts2, *args, **kwargs):
    """
    Finds the overlapping elements in two lists of datetime objects

    Parameters
    ==========
    ts1 : datetime
        first set of datetime object
    ts2: datetime
        datatime object
    args:
        additional arguments passed to tOverlapHalf

    Returns
    =======
    out : list
        indices of ts1 within interval of ts2, & vice versa

    Examples
    ========
    Given two series of datetime objects, event_dates and omni['Time']:

    >>> import spacepy.toolbox as tb
    >>> from spacepy import omni
    >>> import datetime
    >>> event_dates = st.tickrange(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 10, 1), deltadays=3)
    >>> onni_dates = st.tickrange(datetime.datetime(2000, 1, 1), datetime.datetime(2000, 10, 1), deltadays=0.5)
    >>> omni = omni.get_omni(onni_dates)
    >>> [einds,oinds] = tb.tOverlap(event_dates, omni['ticks'])
    >>> omni_time = omni['ticks'][oinds[0]:oinds[-1]+1]
    >>> print omni_time
    [datetime.datetime(2000, 1, 1, 0, 0), datetime.datetime(2000, 1, 1, 12, 0),
    ... , datetime.datetime(2000, 9, 30, 0, 0)]
    """
    idx_1in2 = tOverlapHalf(ts2, ts1, *args, **kwargs)
    idx_2in1 = tOverlapHalf(ts1, ts2, *args, **kwargs)
    if len(idx_2in1) == 0:
        idx_2in1 = None
    if len(idx_1in2) == 0:
        idx_1in2 = None
    return idx_1in2, idx_2in1

def tOverlapHalf(ts1, ts2, presort=False):
    """
    Find overlapping elements in two lists of datetime objects

    This is one-half of tOverlap, i.e. it finds only occurances where
    ts2 exists within the bounds of ts1, or the second element
    returnd by tOverlap.

    Parameters
    ==========
    ts1 : list
        first set of datetime object
    ts2 : list
        datatime object
    presort : bool
        Set to use a faster algorithm which assumes ts1 and
                   ts2 are both sorted in ascending order. This speeds up
                   the overlap comparison by about 50x, so it is worth sorting
                   the list if one sort can be done for many calls to tOverlap

    Returns
    =======
    out : list
        indices of ts2 within interval of ts1

        **note:** Returns empty list if no overlap found
    """
    if presort:
        import bisect
        t_lower, t_upper = ts1[0], ts1[-1]
        return range(bisect.bisect_left(ts2, t_lower),
                     bisect.bisect_right(ts2, t_upper))
    else:
        t_lower, t_upper = min(ts1), max(ts1)
        return [i for i in range(len(ts2))
                if ts2[i] >= t_lower and ts2[i] <= t_upper]

def tCommon(ts1, ts2, mask_only=True):
    """Finds the elements in a list of datetime objects present in another

    Parameters
    ----------
    ts1 : list
        first set of datetime objects
    ts2 : list
        second set of datetime objects

    Returns
    -------
    out : tuple
        Two element tuple of truth tables (of 1 present in 2, & vice versa)
    """
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
        dum1 = [val.replace(tzinfo=None) for val in dum1]
        dum2 = [val.replace(tzinfo=None) for val in dum2]
        if type(ts1)==np.ndarray or type(ts2)==np.ndarray:
            dum1 = np.array(dum1)
            dum2 = np.array(dum2)
        return dum1, dum2

def loadpickle(fln):
    """
    load a pickle and return content as dictionary

    Parameters
    ==========
    fln : string
        filename

    Returns
    =======
    out : dict
        dictionary with content from file

    See Also
    ========
    savepickle

    Examples
    ========
    **note**: If fln is not found, but the same filename with '.gz'
           is found, will attempt to open the .gz as a gzipped file.

    >>> d = loadpickle('test.pbin')
    """
    try:
        import cPickle as pickle
    except ImportError:
        import pickle
    import os.path

    if not os.path.exists(fln) and os.path.exists(fln + '.gz'):
        try:
            import zlib
            with open(fln + '.gz', 'rb') as fh:
                return pickle.loads(
                    zlib.decompress(fh.read(), 16 + zlib.MAX_WBITS))
        except MemoryError:
            import gzip
            with open(fln + '.gz') as fh:
                gzh = gzip.GzipFile(fileobj=fh)
                contents = pickle.load(gzh)
                gzh.close()
            return contents
    else:
        with open(fln, 'rb') as fh:
            return pickle.load(fh)


# -----------------------------------------------
def savepickle(fln, dict, compress=None):
    """
    save dictionary variable dict to a pickle with filename fln

    Parameters
    ----------
    fln : string
        filename
    dict : dict
        container with stuff
    compress : bool
        write as a gzip-compressed file
                     (.gz will be added to L{fln}).
                     If not specified, defaults to uncompressed, unless the
                     compressed file exists and the uncomprssed does not.

    See Also
    ========
    loadpickle

    Examples
    ========
    >>> d = {'grade':[1,2,3], 'name':['Mary', 'John', 'Chris']}
    >>> savepickle('test.pbin', d)
    """
    try:
        import cPickle as pickle
    except:
        import pickle
    if compress == None:
        import os
        if not os.path.exists(fln) and os.path.exists(fln + '.gz'):
            compress = True
        else:
            compress = False
    if compress:
        import gzip
        with open(fln + '.gz', 'wb') as fh:
            gzh = gzip.GzipFile(fln, 'wb', compresslevel=3, fileobj=fh)
            pickle.dump(dict, gzh, 2)
            gzh.close()
    else:
        with open(fln, 'wb') as fh:
            pickle.dump(dict, fh, 2) # 2 ... fast binary
    return


# -----------------------------------------------
def assemble(fln_pattern, outfln, sortkey='ticks', verbose=True):
    """
    assembles all pickled files matching fln_pattern into single file and
    save as outfln. Pattern may contain simple shell-style wildcards *? a la fnmatch
    file will be assembled along time axis given by Ticktock (key: 'ticks') in dictionary
    If sortkey = None, then nothing will be sorted

    Parameters
    ----------
    fln_pattern : string
        pattern to match filenames
    outfln : string
        filename to save combined files to

    Returns
    =======
    out : dict
        dictionary with combined values

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> a, b, c = {'ticks':[1,2,3]}, {'ticks':[4,5,6]}, {'ticks':[7,8,9]}
    >>> tb.savepickle('input_files_2001.pkl', a)
    >>> tb.savepickle('input_files_2002.pkl', b)
    >>> tb.savepickle('input_files_2004.pkl', c)
    >>> a = tb.assemble('input_files_*.pkl', 'combined_input.pkl')
    ('adding ', 'input_files_2001.pkl')
    ('adding ', 'input_files_2002.pkl')
    ('adding ', 'input_files_2004.pkl')
    ('\\n writing: ', 'combined_input.pkl')
    >>> print(a)
    {'ticks': array([1, 2, 3, 4, 5, 6, 7, 8, 9])}
    """
    # done this way so it works before install
    from . import time as t
    from . import coordinates as c

    filelist = glob.glob(fln_pattern)
    filelist = human_sort(filelist)
    # read all files
    d = {}
    for fln in filelist:
        if verbose: print("adding ", fln)
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
            dim = np.array(np.shape(d[fln][key]))
            ax = np.where(dim==TAIcount)[0]
            if len(ax) == 1: # then match with TAI length is given (jump over otherwise like for 'parameters')
                if isinstance(dcomb[key], t.Ticktock):
                    dcomb[key] = dcomb[key].append(d[fln][key])
                elif isinstance(dcomb[key], c.Coords):
                    dcomb[key] = dcomb[key].append(d[fln][key])
                else:
                    dcomb[key] = np.append(dcomb[key], d[fln][key], axis=ax)

    if sortkey:    #  then sort
        if isinstance(dcomb[sortkey], t.Ticktock):
            idx = np.argsort(dcomb[sortkey].RDT)
        else:
            idx = np.argsort(dcomb[sortkey])
        TAIcount = len(dcomb[sortkey])
        for key in dcomb.keys():
            dim = np.array(np.shape(dcomb[key]))
            ax = np.where(dim==TAIcount)[0]
            if len(ax) == 1: # then match with length of TAI
                dcomb[key] = dcomb[key][idx] # resort
    else:
        # do nothing
        pass

    if verbose: print('\n writing: ', outfln)
    savepickle(outfln, dcomb)
    return dcomb

def human_sort( l ):
    """ Sort the given list in the way that humans expect.
    http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html

    Parameters
    ==========
    l : list
        list of objects to human sort

    Returns
    =======
    out : list
        sorted list

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> dat = ['r1.txt', 'r10.txt', 'r2.txt']
    >>> dat.sort()
    >>> print dat
    ['r1.txt', 'r10.txt', 'r2.txt']
    >>> tb.human_sort(dat)
    ['r1.txt', 'r2.txt', 'r10.txt']
    """
    import re
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )
    return l

def feq(x, y, precision=0.0000005):
    """
    compare two floating point values if they are equal
    after: http://www.lahey.com/float.htm

    See Also
    ========
        - http://docs.python.org/tut/node16.html
        - http://www.velocityreviews.com/forums/t351983-precision-for-equality-of-two-floats.html
        - http://www.boost.org/libs/test/doc/components/test_tools/floating_point_comparison.html
        - http://howto.wikia.com/wiki/Howto_compare_float_numbers_in_the_C_programming_language

    Parameters
    ==========
    x : float
        a number
    y : float or array of flaots
        otehr numbers to compare
    precision : float (optional)
        precision for equal (default 0.0000005)

    Returns
    =======
    out : bool
        True (equal) or False (not equal)

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> x = 1 + 1e-4
    >>> y = 1 + 2e-4
    >>> tb.feq(x, y)
    False
    >>> tb.feq(x, y, 1e-3)
    True
    """
    boolean = abs(x-y) <= (abs(x+y)*precision)
    return boolean

def dictree(in_dict, verbose=False, spaces=None, levels=True, attrs=False, **kwargs):
    """
    pretty print a dictionary tree

    Parameters
    ==========
    in_dict : dict
        a complex dictionary (with substructures)
    verbose : boolean (optional)
        print more info
    spaces : string (optional)
        string will added for every line
    levels : integer (optional)
        number of levels to recurse through (True means all)
    attrs : boolean (optional)
        display information for attributes

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> d = {'grade':{'level1':[4,5,6], 'level2':[2,3,4]}, 'name':['Mary', 'John', 'Chris']}
    >>> tb.dictree(d)
    +
    |____grade
         |____level1
         |____level2
    |____name

    More complicated example using a datamodel:

    >>> from spacepy import datamodel
    >>> counts = datamodel.dmarray([2,4,6], attrs={'units': 'cts/s'})
    >>> data = {'counts': counts, 'PI': 'Dr Zog'}
    >>> tb.dictree(data)
    +
    |____PI
    |____counts
    >>> tb.dictree(data, attrs=True, verbose=True)
    +
    |____PI (str [6])
    |____counts (spacepy.datamodel.dmarray (3,))
        :|____units (str [5])

    Attributes of, e.g., a CDF or a datamodel type object (obj.attrs)
    are denoted by a colon.
    """
    try:
        assert hasattr(in_dict, 'keys')
    except AssertionError:
        try:
            assert hasattr(in_dict, 'attrs')
        except:
            raise TypeError('dictree: Input must be dictionary-like')

    if not spaces:
        spaces = ''
        print('+')

    if 'toplev' in kwargs:
        toplev = kwargs['toplev']
    else:
        toplev = True
    try:
        if toplev and attrs:
            dictree(in_dict.attrs, spaces = ':', verbose = verbose, levels = levels, attrs=attrs, toplev=True)
            toplev = False
    except:
        pass

    # TODO, if levels is True why check again?
    if levels:
        try:
            assert levels is True
        except AssertionError:
            levels -= 1
            if levels == 0:
                levels = None

    try:
        for key in sorted(in_dict.keys()):
            bar = '|____' + str(key)
            if verbose:
                typestr = str(type(in_dict[key])).split("'")[1]
                #check entry for dict-like OR .attrs dict
                try:
                    dimstr = in_dict[key].shape
                    dimstr = ' ' + str(dimstr)
                except AttributeError:
                    try:
                        dimstr = len(in_dict[key])
                        dimstr = ' [' + str(dimstr) + ']'
                    except:
                        dimstr = ''
                print(spaces + bar + ' ('+ typestr + dimstr + ')')
            else:
                print(spaces + bar)
            if hasattr(in_dict[key], 'attrs') and attrs:
                dictree(in_dict[key].attrs, spaces = spaces + '    :', verbose = verbose, levels = levels, attrs=attrs, toplev=False)
            if hasattr(in_dict[key], 'keys') and levels:
                dictree(in_dict[key], spaces = spaces + '     ', verbose = verbose, levels = levels, attrs=attrs, toplev=False)
    except:
        pass
    return None

def printfig(fignum, saveonly=False, pngonly=False, clean=False, filename=None):
    """save current figure to file and call lpr (print).

    This routine will create a total of 3 files (png, ps and c.png) in the
    current working directory with a sequence number attached. Also, a time
    stamp and the location of the file will be imprinted on the figure. The
    file ending with c.png is clean and no directory or time stamp are
    attached (good for powerpoint presentations).

    Parameters
    ==========
    fignum : integer
        matplotlib figure number
    saveonly : boolean (optional)
        True (don't print and save only to file)  False (print and save)
    pngolny : boolean (optional)
        True (only save png files and print png directly) False (print ps file, and generate png, ps; can be slow)
    clean : boolean (optional)
        True (print and save only clean files without directory info) False (print and save directory location as well)
    filename : string (optional)
        None (If specified then the filename is set and code does not use the sequence number)

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> import matplotlib.pyplot as plt
    >>> p = plt.plot([1,2,3],[2,3,2])
    >>> tb.printfig(1, pngonly=True, saveonly=True)
    """
    import matplotlib.pyplot as plt

    try:
        nfigs = len(fignum)
    except:
        nfigs = 1
        fignum = [fignum]

    for ifig in fignum:
        # active this figure
        plt.figure(ifig)

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
            plt.savefig(fln+'_clean.png')
            plt.savefig(fln+'_clena.ps')

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S  ")
        # add the filename to the figure for reference
        plt.figtext(0.01, 0.01, timestamp+flnstamp+'.png', rotation='vertical', va='bottom', size=8)

        # now save the figure to this filename
        if pngonly == False:
            plt.savefig(fln+'.ps')

        plt.savefig(fln+'.png')

        # send it to the printer
        if saveonly != True:
            if pngonly == False:
                os.popen('lpr '+fln+'.ps')
            else:
                os.popen('lpr '+fln+'.png')
    return

def update(all=True, omni=False, leapsecs=False, PSDdata=False):
    """
    Download and update local database for omni, leapsecs etc

    Parameters
    ==========
    all : boolean (optional)
        if True, update all of them
    omni : boolean (optional)
        if True. update only onmi
    leapsecs : boolean (optional)
        if True, update only leapseconds

    Returns
    =======
    out : string
        data directory where things are saved

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.update(omni=True)
    """
    from spacepy.time import Ticktock, doy2date
    from spacepy import savepickle, DOT_FLN, OMNI_URL, LEAPSEC_URL, PSDDATA_URL

    if sys.version_info[0]<3:
        import urllib as u
    else:
        import urllib.request as u

    datadir = DOT_FLN+'/data'

    #leapsec_url ='ftp://maia.usno.navy.mil/ser7/tai-utc.dat'
    leapsec_fname = DOT_FLN+'/data/tai-utc.dat'

    # define location for getting omni
    #omni_url = 'ftp://virbo.org/QinDenton/hour/merged/latest/WGhour-latest.d.zip'
    omni_fname_zip = DOT_FLN+'/data/WGhour-latest.d.zip'
    omni_fname_pkl = DOT_FLN+'/data/omnidata.pkl'

    PSDdata_fname = DOT_FLN+'/data/psd_dat.sqlite'

    if all == True:
        omni = True
        leapsecs = True

    if omni == True:
        # retrieve omni, unzip and save as table
        print("Retrieving omni file ...")
        u.urlretrieve(OMNI_URL, omni_fname_zip, reporthook=progressbar)
        fh_zip = zipfile.ZipFile(omni_fname_zip)
        data = fh_zip.read(fh_zip.namelist()[0])
        A = np.array(data.split('\n'))
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
        idx = np.where(A != '')[0]
        # put it into a 2D table
        tab = [val.split() for val in A[idx[1:]]]
        stat8 = [val[11] for val in tab]
        stat6 = [val[27] for val in tab]

        tab = np.array(tab, dtype='float32')
        # take out where Dst not available ( = 99999) or year == 0
        idx = np.where((tab[:,12] !=99.0) & (tab[:,0] != 0))[0]
        tab = tab[idx,:]
        stat8 = np.array(stat8)[idx]
        stat6 = np.array(stat6)[idx]

        omnidata = {}
        # sort through and make an omni dictionary
        # extract keys from line above
        for ikey, i  in zip(keys,range(len(keys))):
            if ikey in ('Year', 'DOY', 'Hr', 'Dst'):
                omnidata[ikey] = np.array(tab[:, i], dtype='int16')
            else:
                omnidata[ikey] = tab[:,i]

        # add TAI to omnidata
        nTAI = len(omnidata['DOY'])

        # add interpolation quality flags
        omnidata['Qbits'] = {}
        arr = np.array(list(np.array(stat8).tostring()), dtype=np.byte).reshape((8,nTAI))
        for ik, key in enumerate(['ByIMF', 'BzIMF', 'velo', 'dens', 'Pdyn', 'G1', 'G2', 'G3']):
            omnidata['Qbits'][key] = arr[ik,:]
        arr = np.array(list(np.array(stat6).tostring()), dtype=np.byte).reshape((6,nTAI))
        for ik, key in enumerate(['W1', 'W2', 'W3', 'W4', 'W5', 'W6']):
            omnidata['Qbits'][key] = arr[ik,:]

        #remove string status keys
        foo = omnidata.pop('6_status')
        foo = omnidata.pop('8_status')

        # add time information to omni pickle (long loop)
        omnidata['UTC'] = [datetime.datetime(int(omnidata['Year'][i]), 1, 1) +
                 datetime.timedelta(days=int(omnidata['DOY'][i]) - 1,
                                    hours=int(omnidata['Hr'][i]))
                 for i in range(nTAI)]

        omnidata['ticks'] = Ticktock(omnidata['UTC'], 'UTC')
        omnidata['RDT'] = omnidata['ticks'].RDT
        del omnidata['ticks'] #Can be quickly regenerated on import
        del omnidata['Year']
        del omnidata['Hr']

        # save as pickle
        savepickle(omni_fname_pkl, omnidata)

        # delete left-overs
        os.remove(omni_fname_zip)

    if leapsecs == True:
        print("Retrieving leapseconds file ... ")
        u.urlretrieve(LEAPSEC_URL, leapsec_fname)

    if PSDdata == True:
        print("Retrieving PSD sql database")
        u.urlretrieve(PSDDATA_URL, PSDdata_fname, reporthook=progressbar)
    return datadir

def progressbar(count, blocksize, totalsize):
    """
    print a progress bar with urllib.urlretrieve reporthook functionality

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> import urllib
    >>> urllib.urlretrieve(PSDDATA_URL, PSDdata_fname, reporthook=tb.progressbar)
    """
    percent = int(count*blocksize*100/totalsize)
    sys.stdout.write("\rDownload Progress " + "...%d%%" % percent)
    if percent == 100: print('\n')
    sys.stdout.flush()


def windowMean(data, time=[], winsize=0, overlap=0, st_time=None):
    """
    Windowing mean function, window overlap is user defined

    Parameters
    ==========
    data : array_like
        1D series of points;
    time : list (optional)
        series of timestamps, optional (format as numeric or datetime)
        For non-overlapping windows set overlap to zero.
    winsize : integer or datetime.timedelta (optional)
        window size
    overlap : integer or datetime.timedelta (optional)
        amount of window overlap
    st_time : datetime.datetime (optional)
        for time-based averaging, a start-time other than the first
        point can be specified

    Returns
    =======
    out : tuple
        the windowed mean of the data, and an associated reference time vector

    Examples
    ========
    For non-overlapping windows set overlap to zero.
    e.g. (time-based averaging)
    Given a data set of 100 points at hourly resolution (with the time tick in
    the middle of the sample), the daily average of this, with half-overlapping
    windows is calculated:

    >>> import spacepy.toolbox as tb
    >>> from datetime import datetime, timedelta
    >>> wsize = datetime.timedelta(days=1)
    >>> olap = datetime.timedelta(hours=12)
    >>> data = [10, 20]*50
    >>> time = [datetime.datetime(2001,1,1) + datetime.timedelta(hours=n, minutes = 30) for n in range(100)]
    >>> outdata, outtime = tb.windowMean(data, time, winsize=wsize, overlap=olap, st_time=datetime.datetime(2001,1,1))
    >>> outdata, outtime
    ([15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0],
     [datetime.datetime(2001, 1, 1, 12, 0),
      datetime.datetime(2001, 1, 2, 0, 0),
      datetime.datetime(2001, 1, 2, 12, 0),
      datetime.datetime(2001, 1, 3, 0, 0),
      datetime.datetime(2001, 1, 3, 12, 0),
      datetime.datetime(2001, 1, 4, 0, 0),
      datetime.datetime(2001, 1, 4, 12, 0)])

    When using time-based averaging, ensure that the time tick corresponds to
    the middle of the time-bin to which the data apply. That is, if the data
    are hourly, say for 00:00-01:00, then the time applied should be 00:30.
    If this is not done, unexpected behaviour can result.

    e.g. (pointwise averaging),

    >>> outdata, outtime = tb.windowMean(data, winsize=24, overlap=12)
    >>> outdata, outtime
    ([15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0], [12.0, 24.0, 36.0, 48.0, 60.0, 72.0, 84.0])

    where winsize and overlap are numeric,
    in this example the window size is 24 points (as the data are hourly) and
    the overlap is 12 points (a half day). The output vectors start at
    winsize/2 and end at N-(winsize/2), the output time vector
    is basically a reference to the nth point in the original series.

    **note** This is a quick and dirty function - it is NOT optimized, at all.
    """
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
            raise ValueError('windowmean error: data and time must have same length')
        #First check if datetime objects
        try:
            assert type(winsize) == datetime.timedelta
            assert type(overlap) == datetime.timedelta
        except AssertionError:
            raise TypeError('windowmean error: winsize/overlap must be timedeltas')
        pts = False #force time-based averaging
        if (type(time[0]) != datetime.datetime):
            startpt = time[0]

    #now actually do windowing mean
    outdata, outtime = [], []
    data = np.array(data)
    if pts:
        #loop for fixed number of points in window
        if not isinstance(winsize, (int, long)):
            winsize = int(round(winsize))
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
            raise ValueError('Overlap requested greater than size of window')
        while lastpt < time[-1]:
            getinds = tOverlapHalf([startpt,startpt+winsize], time, presort=True)
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
    """
    Calculate median absolute deviation of a given input series

    Median absolute deviation (MAD) is a robust and resistant measure of
    the spread of a sample (same purpose as standard deviation). The
    MAD is preferred to the inter-quartile range as the inter-quartile
    range only shows 50% of the data whereas the MAD uses all data but
    remains robust and resistant. See e.g. Wilks, Statistical methods
    for the Atmospheric Sciences, 1995, Ch. 3.

    Parameters
    ==========
    series : array_like
        the input data series

    Returns
    =======
    out : float
        the median absolute deviation

    Examples
    ========
    Find the median absolute deviation of a data set. Here we use the log-
    normal distribution fitted to the population of sawtooth intervals, see
    Morley and Henderson, Comment, Geophysical Research Letters, 2009.

    >>> import numpy
    >>> import spacepy.toolbox as tb
    >>> numpy.random.seed(8675301)
    >>> data = numpy.random.lognormal(mean=5.1458, sigma=0.302313, size=30)
    >>> print data
    array([ 181.28078923,  131.18152745, ... , 141.15455416, 160.88972791])
    >>> tb.medabsdev(data)
    28.346646721370192

    **note** This implementation is robust to presence of NaNs
    """
    #ensure input is numpy array (and make 1-D)
    series = (np.array(series, dtype=float)).ravel()
    #mask for NaNs
    series = np.ma.masked_where(np.isnan(series),series)
    #get median absolute deviation of unmasked elements
    perc50 = np.median(series.compressed())
    mad = np.median(abs(series.compressed()-perc50))
    return mad

def makePoly(x, y1, y2, face = 'blue', alpha=0.5):
    """
    Make filled polygon for plotting

    Parameters
    ==========
    x : list
        List of x start and stop values
    y1 : list
        List of y lower bounds for fill
    y2 : list
        List of y upper bounds for fill
    face : string (optional)
        color of the fill (default blue)
    alpha : float (optional)
        alpha of the fill (default 0.5)

    .. deprecated:: vesion 0.1
    Equivalent functionality to built-in matplotlib function fill_between

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> poly0c = tb.makePoly(x, ci_low, ci_high, face='red', alpha=0.8)
    >>> ax0.add_patch(poly0qc)
    """
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
    Calculates bin width and number of bins for histogram using Freedman-Diaconis rule, if rule fails, defaults to square-root method

    The Freedman-Diaconis method is detailed in:
        Freedman, D., and P. Diaconis (1981), On the histogram as a density estimator: L2 theory, Z. Wahrscheinlichkeitstheor. Verw. Geb., 57, 453–476

    and is also described by:
        Wilks, D. S. (2006), Statistical Methods in the Atmospheric Sciences, 2nd ed.

    Parameters
    ==========
    data : array_like
        list/array of data values
    verbose : boolean (optional)
        print out some more information

    Returns
    =======
    out : tuple
        calculated width of bins using F-D rule, number of bins (nearest integer) to use for histogram

    Examples
    ========
    >>> import numpy, spacepy
    >>> import matplotlib.pyplot as plt
    >>> numpy.random.seed(8675301)
    >>> data = numpy.random.randn(1000)
    >>> binw, nbins = spacepy.toolbox.binHisto(data)
    >>> print(nbins)
    19.0
    >>> p = plt.hist(data, bins=nbins, histtype='step', normed=True)
    """
    from matplotlib.mlab import prctile
    pul = prctile(data, p=(25,75)) #get confidence interval
    ql, qu = pul[0], pul[1]
    iqr = qu-ql
    binw = 2.*iqr/(len(data)**(1./3.))
    nbins = round((max(data)-min(data))/binw)
    # if nbins is 0, NaN or inf dont use the F-D rule just use sqrt(num) rule
    if nbins==0 or not np.isfinite(nbins) or binw==0:
        nbins = round(np.sqrt(len(data)))
        binw = len(data)/nbins
        if verbose:
            print("Used sqrt rule")
    else:
        if verbose:
            print("Used F-D rule")
    return (binw, nbins)

def smartTimeTicks(time):
    """
    Returns major ticks, minor ticks and format for time-based plots

    smartTimeTicks takes a list of datetime objects and uses the range
    to calculate the best tick spacing and format.  Returned to the user
    is a tuple containing the major tick locator, minor tick locator, and
    a format string -- all necessary to apply the ticks to an axis.

    It is suggested that, unless the user explicitly needs this info,
    to use the convenience function applySmartTimeTicks to place the
    ticks directly on a given axis.

    Parameters
    ==========
    time : list
        list of datetime objects

    Returns
    =======
    out : tuple
        tuple of Mtick - major ticks, mtick - minor ticks, fmt - format
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
    """
    Given an axis 'ax' and a list/array of datetime objects, 'time',
    use the smartTimeTicks function to build smart time ticks and
    then immediately apply them to the given axis.  The first and
    last elements of the time list will be used as bounds for the
    x-axis range.

    The range of the 'time' input value will be used to set the limits
    of the x-axis as well.  Set kwarg 'dolimit' to False to override
    this behavior.

    Parameters
    ==========
    ax : matplotlib.pyplot.Axes
        A matplotlib Axis object.
    time : list
        list of datetime objects
    dolimit : boolean (optional)
        The range of the 'time' input value will be used to set the limits
        of the x-axis as well. Setting this overrides this behavior.
    """
    Mtick, mtick, fmt = smartTimeTicks(time)
    ax.xaxis.set_major_locator(Mtick)
    ax.xaxis.set_minor_locator(mtick)
    ax.xaxis.set_major_formatter(fmt)
    if dolimit:
        ax.set_xlim([time[0], time[-1]])

def logspace(min, max, num, **kwargs):
    """
    Returns log-spaced bins. Same as numpy.logspace except the min and max are the min and max
    not log10(min) and log10(max)

    Parameters
    ==========
    min : float
        minimum value
    max : float
        maximum value
    num : integer
        number of log spaced bins

    Other Parameters
    ================
    kwargs : dict
        additional keywords passed into matplotlib.dates.num2date

    Returns
    =======
    out : array
        log-spaced bins from min to max in a numpy array

    Notes
    =====
    This function works on both numbers and datetime objects

    Examples
    ========
    >>> import sapcepy.toolbox as tb
    >>> tb.logspace(1, 100, 5)
    array([   1.        ,    3.16227766,   10.        ,   31.6227766 ,  100.        ])

    See Also
    ========
    toolbox.geomspace
    toolbox.linspace
    """
    from numpy import logspace, log10
    if isinstance(min, datetime.datetime):
        from matplotlib.dates import date2num, num2date
        ans = num2date(logspace(log10(date2num(min)), log10(date2num(max)), num, **kwargs))
        ans = [val.replace(tzinfo=None) for val in ans]
        return np.array(ans)
    else:
        return logspace(log10(min), log10(max), num, **kwargs)

def linspace(min, max, num=50, endpoint=True, retstep=False):
    """
    Returns linearly spaced numbers.  Same as numpy.linspace except
    allows for support of datetime objects

    Parameters
    ==========
    start : float
        The starting value of the sequence.
    stop : float
        The end value of the sequence, unless `endpoint` is set to False.
        In that case, the sequence consists of all but the last of ``num + 1``
        evenly spaced samples, so that `stop` is excluded.  Note that the step
        size changes when `endpoint` is False.
    num : int (optional)
        Number of samples to generate. Default is 50.
    endpoint : bool, optional
        If True, `stop` is the last sample. Otherwise, it is not included.
        Default is True.
    retstep : bool (optional)
        If True, return (`samples`, `step`), where `step` is the spacing
        between samples.

    Returns
    =======
    samples : array
        There are `num` equally spaced samples in the closed interval
        ``[start, stop]`` or the half-open interval ``[start, stop)``
        (depending on whether `endpoint` is True or False).
    step : float (only if `retstep` is True)
        Size of spacing between samples.

    See Also
    ========
    toolbox.geomspace
    toolbox.logspace
    """
    from numpy import linspace, log10
    if isinstance(min, datetime.datetime):
        from matplotlib.dates import date2num, num2date
        ans = num2date(linspace(date2num(min), date2num(max),
                                 num=num, endpoint=endpoint, retstep=retstep))
        ans = [val.replace(tzinfo=None) for val in ans]
        return np.array(ans)
    else:
        return linspace(min, max,
                        num=num, endpoint=endpoint, retstep=retstep)

def geomspace(start, ratio=None, stop=False, num=50):
    """
    Returns geometrically spaced numbers.

    Parameters
    ==========
    start : float
        The starting value of the sequence.
    ratio : float (optional)
        The ratio between subsequent points
    stop: float (optional)
        End value, if this is selected `num' is overridden
    num : int (optional)
        Number of samples to generate. Default is 50.

    Returns
    =======
    seq : array
        geometrically spaced sequence

    Examples
    ========
    To get a geometric progression between 0.01 and 3 in 10 steps

    >>> import spacepy.toolbox as tb
    >>> tb.geomspace(0.01, stop=3, num=10)
    [0.01,
     0.018846716378431192,
     0.035519871824902655,
     0.066943295008216955,
     0.12616612944575134,
     0.23778172582285118,
     0.44814047465571644,
     0.84459764235318191,
     1.5917892219322083,
     2.9999999999999996]

     To get a geometric progression with a specified ratio, say 10

    >>> import spacepy.toolbox as tb
    >>> tb.geomspace(0.01, ratio=10, num=5)
     [0.01, 0.10000000000000001, 1.0, 10.0, 100.0]

    See Also
    ========
    toolbox.linspace
    toolbox.logspace
    """
    if not ratio and stop != False:
        ratio = (stop/start)**(1/(num-1))
    seq = []
    seq.append(start)
    if stop == False:
        for j in range(1, num):
            seq.append(seq[j-1]*ratio)
        return seq
    else:
        val, j = start, 1
        while val <= stop or feq(val, stop, ):
            val = seq[j-1]*ratio
            seq.append(val)
            j+=1
        return seq[:-1]

def arraybin(array, bins):
    """
    Split a sequence into subsequences based on value.

    Given a sequence of values and a sequence of values representing the
    division between bins, return the indices grouped by bin.

    Parameters
    ==========
    array : array_like
        the input sequence to slice, must be sorted in ascending order
    bins : array_like
        dividing lines between bins. Number of bins is len(bins)+1,
            value that exactly equal a dividing value are assigned
            to the higher bin

    Returns
    =======
    out : list
        indices for each bin (list of lists)

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.arraybin(range(10), [4.2])
    [[0, 1, 2, 3, 4], [5, 6, 7, 8, 9]]
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

    Parameters
    ==========
    mlt : numpy array
        array of mlt values
    midnight : boolean (optional)
        reference to midnight instead of noon

    Returns
    =======
    out : numpy array
        array of radians

    Examples
    ========
    >>> from numpy import array
    >>> mlt2rad(array([3,6,9,14,22]))
    array([-2.35619449, -1.57079633, -0.78539816,  0.52359878,  2.61799388])
    """
    if midnight:
        try:
            mlt_arr =  [val + 12 for val in mlt]
        except TypeError:
            mlt_arr = mlt + 12
    else:
        mlt_arr = mlt
    try:
        rad_arr = [(val-12)*np.pi/12. for val in mlt_arr]
    except TypeError:
        rad_arr = (mlt_arr-12)*np.pi/12.
    return rad_arr

def rad2mlt(rad, midnight=False):
    """
    Convert radian values to mlt
    transform radians from -pi to pi to mlt
    referenced from noon by default

    Parameters
    ==========
    rad : numpy array
        array of radian values
    midnight : boolean (optional)
        reference to midnight instead of noon

    Returns
    =======
    out : numpy array
        array of mlt values

    Examples
    ========
    >>> rad2mlt(array([0,pi, pi/2.]))
    array([ 12.,  24.,  18.])
    """
    if midnight:
        rad_arr = rad + np.pi
    else:
        rad_arr = rad
    mlt_arr=rad_arr*(12/np.pi) + 12
    return mlt_arr

def leap_year(year, numdays=False):
    """
    return an array of boolean leap year,
    a lot faster than the mod method that is normally seen

    Parameters
    ==========
    year : array_like
        array of years
    numdays : boolean (optional)
        optionally return the number of days in the year

    Returns
    =======
    out : numpy array
        an array of boolean leap year, or array of number of days

    Examples
    ========
    >>> import numpy
    >>> import spacepy.toolbox as tb
    >>> leap_year(numpy.arange(15)+1998)
    array([False, False,  True, False, False, False,  True, False, False,
    ... False,  True, False, False, False,  True], dtype=bool)
    """
    if not isinstance(year, (tuple, np.ndarray, list)):
        year = [year]
    mask400 = [(val % 400) == 0 for val in year]   # this is a leap year
    mask100 = [(val % 100) == 0 for val in year ]   # these are not leap years
    mask4   = [(val % 4) == 0 for val in year ]   # this is a leap year
    if numdays:
        numdays=365
        ans = [numdays + ((val[0] | val[2]) & (~val[1] | val[0])) for val in zip(mask400, mask100, mask4)]
        if len(ans) == 1:
            return ans[0]
        else:
            return ans
    else:
        ans = [bool(((val[0] | val[2]) & (~val[1] | val[0]))) for val in zip(mask400, mask100, mask4)]
        if len(ans) == 1:
            return ans[0]
        else:
            return ans

leapyear = leap_year

def pmm(a, *b):
    """
    print min and max of input arrays

    Parameters
    ==========
    a : numpy array
        input array
    b : list argumants
        some additional number of arrays

    Returns
    =======
    out : list
        list of min, max for each array

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> from numpy import arange
    >>> tb.pmm(arange(10), arange(10)+3)
    [[0, 9], [3, 12]]
    """
    ans= [[ np.min(a), np.max(a) ]]
    for val in b:
        ans.append( [np.min(val), np.max(val)] )
    return ans

def timestamp(position=[1.003, 0.01], size='xx-small', draw=True, **kwargs):
    """
    print a timestamp on the current plot, vertical lower right

    Parameters
    ==========
    position : list
        position for the timestamp
    size : string (optional)
        text size
    draw : bollean (optional)
        call draw to make sure it appears
    kwargs : keywords
        other keywords to axis.annotate

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> from pylab import plot, arange
    >>> plot(arange(11))
    [<matplotlib.lines.Line2D object at 0x49072b0>]
    >>> tb.timestamp()
    """
    from matplotlib.pyplot import gca, draw
    now = datetime.datetime.now()
    strnow = now.strftime("%d%b%Y %H:%M")
    ax=gca()
    ax.annotate(strnow, position, xycoords='axes fraction', rotation='vertical', size=size, va='bottom',  **kwargs)
    if draw:
        draw()

def query_yes_no(question, default="yes"):
    """
    Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
    It must be "yes" (the default), "no" or None (meaning
    an answer is required of the user).

    The "answer" return value is one of "yes" or "no".

    Parameters
    ==========
    question : string
        the question to ask
    default : string (optional)

    Returns
    =======
    out : string
        answer ('yes' or 'no')

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.query_yes_no('Ready to go?')
    Ready to go? [Y/n] y
    'yes'
    """
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
    """
    1-D linear interpolation with interpolation of hours/longitude

    Parameters
    ==========
    newx : array_like
        x values where we want the interpolated values
    x : array_like
        x values of the origional data
    y : array_like
        7 values of the origional data
    wrap : string, optional
        for continous x data that wraps in y at 'hours' (24), 'longitude' (360),
        or arbitary value (int, float)
    kwargs : dict
        additional keywords, currently accepts baddata that sets baddata for
        masked arrays

    Returns
    =======
    out : numpy.masked_array
        interpolated data values for new abscissa values

    Examples
    ========
    For a simple interpolation

    >>> import spacepy.toolbox as tb
    >>> import numpy
    >>> x = numpy.arange(10)
    >>> y = numpy.arange(10)
    >>> tb.interpol(numpy.arange(5)+0.5, x, y)
    array([ 0.5,  1.5,  2.5,  3.5,  4.5])

    To use the wrap functionality, without the wrap keyword you get the wrong answer

    >>> y = range(24)*2
    >>> x = range(len(y))
    >>> tb.interpol([1.5, 10.5, 23.5], x, y, wrap='hour').compressed() # compress removed the masked array
    array([  1.5,  10.5,  23.5])
    >>> tb.interpol([1.5, 10.5, 23.5], x, y)
    array([  1.5,  10.5,  11.5])
    """
    import scipy as sci

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

    Parameters
    ==========
    vec : array_like
        input vector to normalize

    Returns
    =======
    out : array_like
        normalized vector

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.normalize([1,2,3])
    [0.0, 0.5, 1.0]
    """
    # check to see if vec is numpy array, this is fastest
    if isinstance(vec, np.ndarray):
        out = (vec - vec.min())/np.ptp(vec)
    else:
        vecmin = np.min(vec)
        ptp = np.ptp(vec)
        out = [(val -  vecmin)/ptp for val in vec]
    return out

def listUniq(inVal):
    """
    Given an input iterable (list, deque) return a list of the unique elements.
    Maintains order (keeps the first of non-unique elements

    Parameters
    ==========
    inVal : iterable
        Input iterable

    Returns
    =======
    out : list
        list of unique elements from iterable

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> a = [1,1,2,3,3,4,5,5]
    >>> tb.listUniq(a)
    [1, 2, 3, 4, 5]
    """
    seen = set()
    return [ x for x in inVal if x not in seen and not seen.add(x)]

def intsolve(func, value, start=None, stop=None, maxit=1000):
    """
    Find the function input such that definite integral is desired value.

    Given a function, integrate from an (optional) start point until the
    integral reached a desired value, and return the end point of the
    integration.

    Parameters
    ==========
    func : callable
        function to integrate, must take single parameter
    value : float
        desired final value of the integral
    start : float (optional)
        value at which to start integration, default -Infinity
    stop : float (optional)
        value at which to stop integration, default +Infinity
    maxit : integer
        maximum number of iterations

    Returns
    =======
    out : float
        x such that the integral of L{func} from L{start} to x is L{value}

    **Note:** Assumes func is everywhere positive, otherwise solution may
           be multi-valued.
    """
    from scipy import inf
    from scipy.integrate import quad
    from warnings import warn
    if start is None:
        start = -inf
    if stop is None:
        stop = inf
    lower_bound = start
    upper_bound = stop
    it = 0
    while it < maxit:
        it += 1
        if upper_bound == inf:
            if lower_bound == -inf:
                test_bound = 0
            elif lower_bound < 1.0:
                test_bound = 1.0
            else:
                test_bound = lower_bound * 2
        else:
            if lower_bound == -inf:
                if upper_bound > -1.0:
                    test_bound = -1.0
                else:
                    test_bound = upper_bound * 2
            else:
                test_bound = (lower_bound + upper_bound) / 2.0
        (test_value, err) = quad(func, start, test_bound)
        if abs(value - test_value) <= err:
            break
        elif value < test_value:
            upper_bound = test_bound
        else:
            lower_bound = test_bound

    if abs(value - test_value) > err:
        warn('Difference between desired value and actual is ' +
             str(abs(value - test_value)) +
             ', greater than integral error ' +
             str(err), UserWarning, stacklevel=2)
    return test_bound

def dist_to_list(func, length, min=None, max=None):
    """
    Convert a probability distribution function to a list of values

    This is a deterministic way to produce a known-length list of values
    matching a certain probability distribution. It is likely to be a closer
    match to the distribution function than a random sampling from the
    distribution.

    Parameters
    ==========
    func : callable
        function to call for each possible value, returning
            probability density at that value (does not need to be
            normalized.)
    length : int
        number of elements to return
    min : float
        minimum value to possibly include
    max : float
        maximum value to possibly include

    Examples
    ========
    >>> import matplotlib
    >>> import numpy
    >>> import spacepy.toolbox as tb
    >>> gauss = lambda x: math.exp(-(x ** 2) / (2 * 5 ** 2)) / \
                          (5 * math.sqrt(2 * math.pi))
    >>> vals = tb.dist_to_list(gauss, 1000, -numpy.inf, numpy.inf)
    >>> print vals[0]
    -16.45263...
    >>> p1 = matplotlib.pyplot.hist(vals, bins=[i - 10 for i in range(21)], \
                               facecolor='green')
    >>> matplotlib.pyplot.hold(True)
    >>> x = [i / 100.0 - 10.0 for i in range(2001)]
    >>> p2 = matplotlib.pyplot.plot(x, [gauss(i) * 1000 for i in x], 'red')
    >>> matplotlib.pyplot.draw()
    """
    from scipy import inf
    from scipy.integrate import quad
    from warnings import warn
    if min is None:
        min = -inf
    if max is None:
        max = inf
    total = quad(func, min, max)[0]
    step = float(total) / length
    return [intsolve(func, (0.5 + i) * step, min, max) for i in range(length)]

def bin_center_to_edges(centers):
    """
    Convert a list of bin centers to their edges

    Given a list of center values for a set of bins, finds the start and
    end value for each bin. (start of bin n+1 is assumed to be end of bin n).
    Useful for e.g. matplotlib.pyplot.pcolor.

    Edge between bins n and n+1 is arithmetic mean of the center of n
    and n+1; edge below bin 0 and above last bin are established to make
    these bins symmetric about their center value.

    Parameters
    ==========
    centers : list
        list of center values for bins

    Returns
    =======
    out : list
        list of edges for bins

    **note:** returned list will be one element longer than centers

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.bin_center_to_edges([1,2,3])
    [0.5, 1.5, 2.5, 3.5]
    """
    return [1.5 * centers[0] - 0.5 * centers[1] if i == 0 else
            1.5 * centers[-1] - 0.5 * centers[-2] if i == len(centers) else
            (centers[i - 1] + centers[i]) / 2.0
            for i in range(len(centers) + 1)]

def hypot(*vals):
    """
    Compute sqrt(vals[0] **2 + vals[1] **2 ...), ie. n-dimensional hypoteneuse

    Parameters
    ==========
    vals : float (arbitary number)
        arbitary number of float values

    Returns
    =======
    out : float
        the Euclidian distance of the points ot the origin

    Examples
    ========
    >>> import spacepy.toolbox as tb
    >>> tb.hypot(3,4)
    5.0
    >>> a = [3, 4]
    >>> tb.hypot(*a)
    5.0
    >>> tb.hypot(*range(10))
    16.88194...
    """
    return math.sqrt(sum((v ** 2 for v in vals)))

def thread_job(job_size, thread_count, target, *args, **kwargs):
    """
    Split a job into subjobs and run a thread for each

    Each thread spawned will call L{target} to handle a slice of the job.

    This is only useful if a job:
      1. Can be split into completely independent subjobs
      2. Relies heavily on code that does not use the Python GIL, e.g.
         numpy or ctypes code
      3. Does not return a value. Either pass in a list/array to hold the
         result, or see L{thread_map}

    Examples
    ========
    squaring 100 million numbers:

    >>> import numpy
    >>> import spacepy.toolbox as tb
    >>> numpy.random.seed(8675301)
    >>> a = numpy.random.randint(0, 100, [100000000])
    >>> b = numpy.empty([100000000], dtype='int64')
    >>> def targ(in_array, out_array, start, count): \
             out_array[start:start + count] = in_array[start:start + count] ** 2
    >>> tb.thread_job(len(a), 0, targ, a, b)
    >>> print(b[0:5])
    [2704 7225  196 1521   36]

    This example:
      - Defines a target function, which will be called for each thread.
        It is usually necessary to define a simple "wrapper" function
        like this to provide the correct call signature.
      - The target function receives inputs C{in_array} and C{out_array},
        which are not touched directly by C{thread_job} but are passed
        through in the call. In this case, C{a} gets passed as
        C{in_array} and C{b} as C{out_array}
      - The target function also receives the start and number of elements
        it needs to process. For each thread where the target is called,
        these numbers are different.

    Parameters
    ==========
    job_size : int
        Total size of the job. Often this is an array size.
    thread_count : int
        Number of threads to spawn. If =0 or None, will
            spawn as many threads as there are cores available on
            the system. (Each hyperthreading core counts as 2.)
            Generally this is the Right Thing to do.
            If NEGATIVE, will spawn abs(thread_count) threads,
            but will run them sequentially rather than in
            parallel; useful for debugging.
    target : callable
        Python callable (generally a function, may also be an
            imported ctypes function) to run in each thread. The
            *last* two positional arguments passed in will be a
            "start" and a "subjob size," respectively;
            frequently this will be the start index and the number
            of elements to process in an array.
    args : sequence
        Arguments to pass to L{target}. If L{target} is an instance
            method, self must be explicitly pssed in. start and
            subjob_size will be appended.
    kwargs : dict
        keyword arguments to pass to L{target}.
    """
    try:
        import threading
    except:
        return (target(*(args +  (0, job_size)), **kwargs), )
    if thread_count == None or thread_count == 0:
        try:
            import multiprocessing
            thread_count = multiprocessing.cpu_count()
        except: #Do something not too stupid
            thread_count = 8
    if thread_count < 0:
        thread_count *= -1
        seq = True
    else:
        seq = False
    count = float(job_size) / thread_count
    starts = [int(count * i + 0.5) for i in range(thread_count)]
    subsize = [(starts[i + 1] if i < thread_count - 1 else job_size) -
               starts[i]
               for i in range(thread_count)]
    threads = []
    for i in range(thread_count):
        t = threading.Thread(
            None, target, None,
            args + (starts[i], subsize[i]), kwargs)
        t.start()
        if seq:
            t.join()
        else:
            threads.append(t)
    if not seq:
        for t in threads:
            t.join()

def thread_map(target, iterable, thread_count=None, *args, **kwargs):
    """
    Apply a function to every element of a list, in separate threads

    Interface is similar to multiprocessing.map, except it runs in threads

    Examples
    ========
    find totals of several arrays

    >>> import numpy
    >>> from spacepy import toolbox
    >>> inputs = range(100)
    >>> totals = toolbox.thread_map(numpy.sum, inputs)
    >>> print(totals[0], totals[50], totals[99])
    (0, 50, 99)

    Parameters
    ==========
    target : callable
        Python callable to run on each element of iterable.
            For each call, an element of iterable is appended to
            args and both args and kwargs are passed through.
            Note that this means the iterable element is always the
            *last* positional argument; this allows the specification
            of self as the first argument for method calls.
    iterable : iterable
        elements to pass to each call of L{target}
    args : sequence
        arguments to pass to target before each element of
            iterable
    thread_count : integer
        Number of threads to spawn; see L{thread_job}.
    kwargs : dict
        keyword arguments to pass to L{target}.

    Returns
    =======
    out : list
        return values of L{target} for each item from L{iterable}
    """
    try:
        jobsize = len(iterable)
    except TypeError:
        iterable = list(iterable)
        jobsize = len(iterable)
    def array_targ(function, it, retvals, arglist, kwarglist, start, size):
        for i in range(start, start + size):
            retvals[i] = function(*(arglist + (it[i],)), **kwarglist)
    retvals = [None] * jobsize
    thread_job(jobsize, thread_count, array_targ,
               target, iterable, retvals, args, kwargs)
    return retvals


if __name__ == "__main__":
    import doctest
    doctest.testmod()
