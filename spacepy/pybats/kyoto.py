#!/usr/bin/env python
'''
kyoto is a tool set for obtaining and handling geomagnetic indices stored at the
Kyoto World Data Center (WDC) website.  Indices can be loaded from file or
fetched from the web.

Instantiation of objects from this module should be done through the constructor
functions fetch and load.  Use help on these objects for more
information.
'''

# Module wide imports.
import numpy as np
import datetime as dt
from pybats import PbData
from spacepy.datamodel import dmarray, SpaceData

#############################################################################
# DATA CLASSES
#############################################################################

class KyotoDst(PbData):
    '''
    Handles the hourly dst index from Kyoto WDC.  

    Use the specialized constructor functions fetch and load to
    instantiate from the web or from file.
    '''

    def __init__(self, lines=None, *args, **kwargs):
        # Init as PbData:
        super(KyotoDst, self).__init__(*args, **kwargs)
        if lines:
            self._parse(lines)
        else:
            pass

    def _parse(self, lines):
        '''
        Given raw ascii input as a list of lines, parse into object.
        '''

        if lines[0][0:7]==' Format': #IAGA-2002 formatted file.
            data=parse_iaga(lines, iagacode='DST')
            for key in data:
                self[key]=data[key]
            self['dst'].attrs['units']='nT'
            self.attrs['npts']=data['time'].size
            return
    
        self.attrs['npts']=len(lines)
        time = []
        dst  = np.zeros(24*self.attrs['npts'])
        for i,line in enumerate(lines):
            # Get year, month, day.
            try:
                yy = int(line[14:16]) * 100
            except:
                yy = 1900
            yy = yy + int(line[3:5])
            dd = int(line[8:10])
            mm = int(line[5:7 ])
        
        # Parse the rest of the data.
            for j in range(0,24):
                time.append(dt.datetime(yy, mm, dd, j))
                loc = 20 + 4*j
                dst[24*i + j] = float(line[loc:loc+4])

        self['time']= dmarray(time)
        self['dst'] = dmarray(dst, attrs={'units':'nT'})

#===========================================================================
class KyotoKp(PbData):
    '''
    Handles global Kp index from Kyoto WDC.
    
    Use the specialized constructor functions fetch and load to
    instantiate from the web or from file.
    '''
    def __init__(self, lines=None, *args, **kwargs):
        # Init as PbData:
        super(KyotoKp, self).__init__(*args, **kwargs)
        if lines:
            self._parse(lines)
        else:
            pass

    def _parse(self, lines):
        '''
        Given raw ascii input as a list of lines, parse into object.
        '''
        
        # Kill header as necessary.
        if lines[0][1:5]=='HTML':
            lines.pop(0);lines.pop(0);lines.pop(-1)
        if lines[0][0:4]=='YYYY':
            lines.pop(0)

        npoints = len(lines)
        time = np.zeros(8*npoints, dtype=object)
        t1 = np.zeros(8*npoints, dtype=object)
        t2 = np.zeros(8*npoints, dtype=object)
        kp   = np.zeros(8*npoints)
        hr1  = [0,3,6,9,12,15,18,21,0]
        hrs  = [1,4,7,10,13,16,19,22]
        frac = {' ':0.0, '-':-1./3., '+':1./3.}

        # Parse lines.
        for i,line in enumerate(lines):
            yy = int(line[0:4])
            mm = int(line[4:6])
            dd = int(line[6:8])
            for j in range(8):
                kp[8*i+j]  = float(line[9+2*j])+frac[line[10+2*j]]
                time[8*i+j]= dt.datetime(yy,mm,dd,hrs[j],30,0)
                t1[8*i+j]= dt.datetime(yy,mm,dd,hr1[j],0,0)
                t2[8*i+j]= dt.datetime(yy,mm,dd,hr1[j+1],0,0) \
                    +dt.timedelta(days=j/7)
        # Store data into object.
        self['time']=dmarray(time)
        self['binstart']=dmarray(t1)
        self['binstop']=dmarray(t2)
        self['kp']=dmarray(kp, attrs={'units':'None'})
        self.attrs['npts']=npoints

    def add_histplot(self, target=False, loc=111, label='Kyoto $K_{p}$', 
                     **kwargs):
        '''
        Make a quick histogram-style plot of the Kp data.
        '''

        import matplotlib.pyplot as plt
        from pybats import apply_smart_timeticks

        # Shortcuts for the lazy.
        bstart=self['binstart']
        bstop =self['binstop']
        npts  =self.attrs['npts']

        # Reformulate time to get histogram-type look.
        newtime=np.zeros(npts*24, dtype=object)
        newkp  =np.zeros(npts*24)
        for i in range(npts*8):
            newtime[3*i  ] = bstart[i]
            newtime[3*i+1] = self['time'][i]
            newtime[3*i+2] = bstop[i]
            newkp[3*i:3*i+3] = self['kp'][i], self['kp'][i], self['kp'][i]
    
        if type(target) == plt.Figure:
            fig = target
            ax = fig.add_subplot(loc)
            line=ax.plot(newtime, newkp, label=label, **kwargs)
            apply_smart_timeticks(ax, newtime, dolabel=True)
        elif type(target).__base__ == plt.Axes:
            ax = target
            fig = ax.figure
            line=ax.plot(newtime, newkp, label=label, **kwargs)
        else:
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(loc)
            line=ax.plot(newtime, newkp, label=label, **kwargs)
            apply_smart_timeticks(ax, newtime, dolabel=True)
        return fig, ax

#############################################################################
# CONSTRUCTOR FUNCTIONS
#############################################################################

def fetch(dtype, tstart, tstop, debug=False):
    '''
    Fetch Kyoto data of type dtype from the web given a start time, tstart,
    and a stop time, tstop.

    dtype can be any of the following strings: dst, kp
    Start and stop times can be year-month tuples or datetime objects.

    The correct, filled data class will be returned.

    Example 1: Get Kp for the whole year of 2005:
    kp=fetch('kp', (2005,1), (2005,12))
    
    Example 2: Get Dst for the a few months from 1981:
    import datetime
    t1=datetime.datetime(1981,9, 1)
    t2=datetime.datetime(1981,11, 31)
    dst=fetch('dst', t1, t2)
    '''

    objts={'dst':KyotoDst,'kp':KyotoKp}
    funcs={'dst':dstfetch,'kp':kpfetch}

    time=[]
    # Check and parse input times.
    for t in ([tstart, tstop]):
        if type(t)==dt.datetime:
            time.append( (t.year, t.month) )
        elif type(t)==type( () ) or type(t)==type( [] ):
            time.append(t)
        else:
            raise ValueError, 'Unrecognized type for input time.'

    if debug:
        print 'Fetching from datetuple ', t[0], ' to ', t[1]
        print 'Using functions for ', dtype
        print 'Calling with ', time[0][0], time[0][1], time[1][0], time[1][1]
    lines=funcs[dtype](time[0][0], time[0][1], time[1][0], time[1][1])
    obj=objts[dtype]()
    obj._parse(lines)
    
    return obj

#############################################################################
# WEBFETCH FUNCTIONS
#############################################################################

def inttomonth(inmonth):
    '''
    Turn integer month into string month.
    Necessary for webfetch.
    '''
    if inmonth > 12: inmonth = 12
    if inmonth < 1:  inmonth =  1
    strmonth = [
        'Jan', 'Feb', 'Mar', 'Apr',
        'May', 'Jun', 'Jul', 'Aug',
        'Sep', 'Oct', 'Nov', 'Dec'
        ]
    return(strmonth[inmonth-1])

#===========================================================================

def dstfetch(yrstart, mostart, yrstop, mostop):
    '''
    A function to fetch Kyoto Dst directly from the Kyoto WDC website.
    Returns raw ascii lines.
    '''
    import urllib as url
    # Set up all form info.
    forminfo = {} # Date-time stuff - requires special formatting.
    forminfo2= {} # Other form stuff.

    # Start time:
    if (yrstart-2000) >= 0:
        forminfo['SCent'] = 20
    else:
        forminfo['SCent'] = 19
    forminfo['STens'] = int(('%d'%yrstart)[2])
    forminfo['SYear'] = int(('%d'%yrstart)[3])
    forminfo['SMonth'] = '%02i' % mostart
    # End Time:
    if (yrstop-2000) >= 0:
        forminfo['ECent'] = 20
    else:
        forminfo['ECent'] = 19
    forminfo['ETens'] = int(('%d'%yrstop)[2])
    forminfo['EYear'] = int(('%d'%yrstop)[3])
    forminfo['EMonth'] = '%02i' % mostop

    # Build and open URL.  Note that the forminfo must be
    # in a certain order for the request to work...
    email = url.quote('spacepy@lanl.gov')
    target = 'http://wdc.kugi.kyoto-u.ac.jp/cgi-bin/dstae-cgi?' + \
        '%s=%2i&%s=%i&%s=%i&%s=%s&%s=%2i&%s=%i&%s=%i&%s=%s&' % \
        ('SCent', forminfo['SCent'], \
             'STens', forminfo['STens'], \
             'SYear', forminfo['SYear'], \
             'SMonth',forminfo['SMonth'],\
             'ECent', forminfo['ECent'], \
             'ETens', forminfo['ETens'], \
             'EYear', forminfo['EYear'], \
             'EMonth',forminfo['EMonth'])
    # Form format updated 03/14/2011
    target = target + \
        'Image+Type=GIF&COLOR=COLOR&AE+Sensitivity=0' + \
        '&Dst+Sensitivity=0&Output=DST&Out+format=WDC&Email=' +\
        email

    # Fetch data from web.
    f = url.urlopen(target)
    lines = f.readlines()

    return(lines)

#===========================================================================

def kpfetch(yrstart, mostart, yrstop, mostop):
    '''
    A function to fetch Kyoto Kp directly from the Kyoto WDC website.
    Returns raw ascii lines.
    '''
    import urllib as url
    # Set up all form info.
    forminfo = {} # Date-time stuff - requires special formatting.
    forminfo2= {} # Other form stuff.

    # Start time:
    if (yrstart-2000) >= 0:
        forminfo['SCent'] = 20
    else:
        forminfo['SCent'] = 19
    forminfo['STens'] = int(('%d'%yrstart)[2])
    forminfo['SYear'] = int(('%d'%yrstart)[3])
    forminfo['SMonth'] = '%02i' % mostart
    # End Time:
    if (yrstop-2000) >= 0:
        forminfo['ECent'] = 20
    else:
        forminfo['ECent'] = 19
    forminfo['ETens'] = int(('%d'%yrstop)[2])
    forminfo['EYear'] = int(('%d'%yrstop)[3])
    forminfo['EMonth'] = '%02i' % mostop

    # Build and open URL.  Note that the forminfo must be
    # in a certain order for the request to work...
    email = url.quote('dantwelling@lanl.gov')
    data=url.urlencode(forminfo)
    target='http://wdc.kugi.kyoto-u.ac.jp/cgi-bin/kp-cgi'#
    data='%s=%2i&%s=%i&%s=%i&%s=%s&%s=%2i&%s=%i&%s=%i&%s=%s&%s=%s' % \
        ('SCent', forminfo['SCent'], \
             'STens', forminfo['STens'], \
             'SYear', forminfo['SYear'], \
             'SMonth',forminfo['SMonth'],\
             'ECent', forminfo['ECent'], \
             'ETens', forminfo['ETens'], \
             'EYear', forminfo['EYear'], \
             'EMonth',forminfo['EMonth'],\
             'Email' ,email)

    # Fetch data from web.
    f = url.urlopen(target, data) 
    lines = f.readlines()
    return(lines)

#===========================================================================

def parse_iaga(lines, iagacode=None):
    '''
    KyotoWDC uses two format types: WDC, which is data specific, and 
    IAGA-2002, which is general for all data types.  This function is
    a general reader for this format.  It returns a dictionary of vectors,
    each corresponding to a column from the file.

    'lines' is simply a list of lines from the IAGA-formatted file.
    'iagacode', if given, should be a string containing the IAGA code for the
    file contents.  If given, this function will raise an
    exception if iagacode does not match the file's code.  This is 
    useful for ensure the correct data values are located in this file.
    '''

    from dateutil.parser import parse

    # Begin by parsing header; ensuring the correct file format.
    fmt=(lines.pop(0)).split()
    if (fmt[0]!='Format') or (fmt[1]!='IAGA-2002'):
        raise Exception, 'Data is not in IAGA-2002 format.'

    # Parse mandatory IAGA header lines.
    source=(lines.pop(0)).split()[1]
    lines.pop(0)
    code=(lines.pop(0)).split()[2]
    for i in range(8):
        lines.pop(0)
    
    # Check Iaga Code as necessary.
    if iagacode:
        if iagacode != code:
            raise Exception, "IAGA Code does not match required code."

    # Loop through and count optional header lines.
    nHead=12
    while True:
        line=lines.pop(0)
        if line[:2]!=' #': break
        nHead+=1

    # Parse column labels.  We don't need time or DOY.
    parts=line.lower().split()[3:-1]
    data={'time':[], 'doy':[]}
    for name in parts:
        data[name]=[]

    # Read all data.
    for l in lines:
        if l[-2]=='|':continue # skip repeat headers.

        p=l.split()
        data['time'].append(parse(' '.join(p[0:2])))
        data['doy'].append(int(p[2]))
        
        for i,name in enumerate(parts):
            data[name].append(float(p[i+3]))

    # Convert to dmarrays.
    for name in data:
        data[name]=dmarray(data[name])

    return data
