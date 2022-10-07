#!/usr/bin/env python
'''
kyoto is a tool set for obtaining and handling geomagnetic indices stored at the
`Kyoto World Data Center (WDC) website <http://wdc.kugi.kyoto-u.ac.jp/wdc/Sec3.html>`_.  
Indices can be loaded from file or fetched from the web.

Instantiation of objects from this module should be done through the constructor
functions :func:`fetch` and :func:`load`.  Use help on these objects
for more information.
'''

# Module wide imports.
import numpy as np
import datetime as dt
from spacepy.plot import set_target, applySmartTimeTicks, levelPlot
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray, SpaceData

#############################################################################
# DATA CLASSES
#############################################################################

class KyotoDst(PbData):
    '''
    Handle hourly dst index from Kyoto WDC.  

    Use the specialized constructor functions :function:`fetch` and 
    :function:`load` to instantiate from the web or from file.
    '''

    def __repr__(self):
        return 'KyotoDst object for handling observed Dst data.'

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
class KyotoAe(PbData):
    '''
    Handle Ae index from Kyoto WDC.  

    Use the specialized constructor functions :function:`fetch` and 
    :function:`load` to instantiate from the web or from file.
    '''

    def __repr__(self):
        return 'KyotoDst object for handling auroral index data.'

    def __init__(self, lines=None, *args, **kwargs):
        # Init as PbData:
        super(KyotoAe, self).__init__(*args, **kwargs)
        if lines:
            self._parse(lines)
        else:
            pass

    def _parse(self, lines):
        # Decode binary data:
        if type(lines[0])==bytes:
            for i in range(len(lines)):
                lines[i] = lines[i].decode()
                
        # Call parse_iaga to parse lines:
        data = parse_iaga(lines)
        for k in data:
            self[k] = data[k]
            
#===========================================================================
class KyotoSym(PbData):
    '''
    Handle Sym-H and related indices from Kyoto WDC.  

    Use the specialized constructor functions :function:`fetch` and 
    :function:`load` to instantiate from the web or from file.
    '''

    def __repr__(self):
        return 'KyotoDst object for handling observations-based '+\
            'SYM- and ASY- data.'

    def __init__(self, lines=None, *args, **kwargs):
        # Init as PbData:
        super(KyotoSym, self).__init__(*args, **kwargs)
        if lines:
            self._parse(lines)
        else:
            pass

    def _parse(self, lines):
        # Decode binary data:
        if type(lines[0])==bytes:
            for i in range(len(lines)):
                lines[i] = lines[i].decode()
        # Call parse_iaga to parse lines:
        data = parse_iaga(lines, 'ASY/SYM')
        for k in data:
            self[k] = data[k]

#===========================================================================
class KyotoKp(PbData):
    '''
    Handles global Kp index from Kyoto WDC.
    
    Use the specialized constructor functions :function:`fetch` and 
    :function:`load` to instantiate from the web or from file.
    '''
    def __repr__(self):
        return 'KyotoDst object for handling observations-based global '+\
            'Kp index data'

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

        # Decode binary data:
        if type(lines[0])==bytes:
            for i in range(len(lines)):
                lines[i] = lines[i].decode()
        
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
                t2[8*i+j]= t1[8*i+j]+dt.timedelta(hours=3)
        # Store data into object.
        self['time']=dmarray(time)
        self['binstart']=dmarray(t1)
        self['binstop']=dmarray(t2)
        self['kp']=dmarray(kp, attrs={'units':'None'})
        self.attrs['npts']=npoints

    def add_histplot(self, target=False, loc=111, label='Kyoto $K_{p}$', 
                     time_range=None, filled=False, level_kwargs={},
                     **kwargs):
        '''
        Make a quick histogram-style plot of the Kp data.

        Returns
        =======
        fig : matplotlib figure object
        ax  : matplotlib axes object

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

        label : string
            The label applied to the line when a legend is added to the axes.
            Defaults to 'Kyoto $K_{p}$'.

        time_range : tuple of datetimes
            The time range to plot.  Only the first and last values in the
            tuple (or list) are used if the number of elements is greater than
            two.  Defaults to **None**, meaning that the full time range
            available is used.

        filled : Boolean
            If True, make a filled 'traffic light' plot of Kp using
            spacepy.plot.levelPlot. Extra keyword arguments for levelPlot can
            be supplied via level_kwargs.

        Extra keyword arguments are passed to :function:`matplotlib.pyplot.plot`
        to customize the line style (and are ignored if filled is True).

        Examples
        ========
        >>> import matplotlib.pyplot as plt
        >>> import spacepy.pybats.kyoto as kt
         >>> kp = kt.fetch('kp', (1981, 11), (1981, 11))
        >>> kp.add_histplot(lw=2.0, color='r', label='Example Kp')
        >>> ax = plt.gca()
        >>> kp.add_histplot(filled=True)

        '''
        import matplotlib.pyplot as plt

        if not time_range:
            time_range = self['time']
    
        fig, ax = set_target(target, figsize=(10,4),  loc=loc)
        if not filled:
            tvar = self['binstart'].tolist()
            tvar.append(tvar[-1]+dt.timedelta(hours=3))
            kpvar = self['kp'].tolist()
            kpvar.append(kpvar[-1])
            line = ax.plot(tvar, kpvar, label=label, drawstyle='steps-post',
            **kwargs)
        else:
            ax = levelPlot(self, time='binstart', var='kp', target=ax, **level_kwargs)
        applySmartTimeTicks(ax, time_range, dolabel=True)

        return fig, ax

#############################################################################
# CONSTRUCTOR FUNCTIONS
#############################################################################

def fetch(dtype, tstart, tstop, debug=False):
    '''
    Fetch Kyoto data of type dtype from the web given a start time, tstart,
    and a stop time, tstop.

    dtype can be any of the following strings: 

    1. dst -- The hourly Dst index.
    2. kp  -- The 3-hour global Kp index.
    3. ae  -- The AE, AL, and AU minute-resolution auroral indices.
    4. sym -- The minute resolution SYM and ASY indices for H and D directions.

    Start and stop times can be year-month tuples or datetime objects.

    The filled corresponding data class will be returned.

    Example 1: Get Kp for the whole year of 2005:
    kp=fetch('kp', (2005,1), (2005,12))
    
    Example 2: Get Dst for the a few months from 1981:
    import datetime
    t1=datetime.datetime(1981,9, 1)
    t2=datetime.datetime(1981,11, 31)
    dst=fetch('dst', t1, t2)
    '''

    from calendar import monthrange as mrng
    
    objts={'dst':KyotoDst,'kp':KyotoKp,'ae':KyotoAe,'sym':KyotoSym}
    funcs={'dst':dstfetch,'kp':kpfetch,'ae':aefetch,'sym':symfetch}

    # Convert times into datetime objects:
    # Start time:
    if type(tstart)==list or type(tstart)==tuple:
        # Convert tuples into datetimes:
        tstart = dt.datetime(tstart[0], tstart[1], 1, 0, 0)
    elif type(tstart)!=dt.datetime and type(tstop)!=dt.date:
        # Not a date time or list/tuple?  Fail.
        raise ValueError('Unrecognized type for input time.')
    # Stop time:
    if type(tstop)==list or type(tstop)==tuple:
        # Convert tuples into datetime:
        tstop  = dt.datetime(tstop[0], tstop[1], mrng(tstop[0],tstop[1])[1],0,0)
    elif type(tstop)!=dt.datetime and type(tstop)!=dt.date:
        # Not a date time or list/tuple?  Fail.
        raise ValueError('Unrecognized type for input time.')
    
    if debug:
        print('Using functions for ', dtype)
        print('Fetching from {0:%Y-%m-%d %H:%M} to {1:%Y-%m-%d %H:%M}'.format(
            tstart,tstop))

    # Call correct function:
    if dtype=='dst' or dtype=='kp':
        lines=funcs[dtype](tstart.year, tstart.month, tstop.year, tstop.month)
    else:
        lines=funcs[dtype](tstart, tstop)

    # Instantiate proper object and parse lines:    
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
    try:
        import urllib.parse, urllib.request
    except ImportError: #glue the python 3 names onto python 2
        import urllib
        urllib.parse = urllib
        urllib.request = urllib
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
    email = urllib.parse.quote('spacepy@lanl.gov')
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
    f = urllib.request.urlopen(target)
    lines = f.readlines()

    return(lines)

#===========================================================================

def kpfetch(yrstart, mostart, yrstop, mostop):
    '''
    A function to fetch Kyoto Kp directly from the Kyoto WDC website.
    Returns raw ascii lines.
    '''
    try:
        import urllib.parse, urllib.request
    except ImportError: #glue the python 3 names onto python 2
        import urllib
        urllib.parse = urllib
        urllib.request = urllib
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
    email = urllib.parse.quote('spacepy@lanl.gov')
    data=urllib.parse.urlencode(forminfo)
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
    f = urllib.request.urlopen(target, data.encode()) 
    lines = f.readlines()
    return(lines)
#===========================================================================

def symfetch(t_start, t_stop):
    '''
    A function to fetch Kyoto SYM-H directly from the Kyoto WDC website.
    Returns raw ascii lines obtained from the website representing the data
    in IAGA 2002 format.  

    Note that, because of the higher density of the data, full datetimes
    are required for the start and stop time.  This is in contrast to
    other fetch functions, where hourly or three-hourly resolution allows
    for only the start and stop years and months are required as integers.

    Parameters
    ----------
    t_start : datetime.datetime
       Start time for range of requested data as a datetime object.
    t_stop : datetime.datetime
       Stop time for range of requested data as a datetime object.

    '''
    try:
        import urllib.parse, urllib.request
    except ImportError: #glue the python 3 names onto python 2
        import urllib
        urllib.parse = urllib
        urllib.request = urllib

    # Get difference between start and stop time:
    deltaT = t_stop - t_start
        
    # Set up all form info.
    forminfo = {} # Date-time stuff - requires special formatting.

    # GENERATE FORM INFO:
    # First three digits of year, last digit of year:
    forminfo['Tens'] = int(t_start.year/10)
    forminfo['Year'] = int(t_start.year - 10*forminfo['Tens'])
    # Month:
    forminfo['Month']= t_start.month
    # Day (right now, locked in to one month of data):
    forminfo['Day_Tens'] = int(t_start.day/10)
    forminfo['Days']     = int(t_start.day - 10*forminfo['Day_Tens'])
    # Hour/minute:
    forminfo['Hour'] = t_start.hour
    forminfo['min']  = t_start.minute
    # Duration:
    forminfo['Dur_Day_Tens'] = int(deltaT.days/10)
    forminfo['Dur_Day']      = int(deltaT.days - 10*forminfo['Dur_Day_Tens'])
    forminfo['Dur_Hour']     = int(deltaT.seconds/3600.)
    forminfo['Dur_Min']      = int(deltaT.seconds/60.-forminfo['Dur_Hour']*60.)

    # Email and other information:
    forminfo['Output']     = 'ASY'
    forminfo['Out+format'] = 'IAGA2002'
    forminfo['Email']      =  urllib.parse.quote('spacepy@lanl.gov')

    # Build and open URL.  Note that the forminfo must be
    # in a certain order for the request to work...
    target='http://wdc.kugi.kyoto-u.ac.jp/cgi-bin/aeasy-cgi'
    data='Tens={0[Tens]:03d}&Year={0[Year]:d}&'.format(forminfo)
    data+='Month={0[Month]:02d}&Day_Tens={0[Day_Tens]:d}&'.format(forminfo)
    data+='Days={0[Days]:d}&Hour={0[Hour]:02d}&min={0[min]:02d}'.format(forminfo)
    data+='&Dur_Day_Tens={0[Dur_Day_Tens]:02d}&'.format(forminfo)
    data+='Dur_Day={0[Dur_Day]:d}&Dur_Hour={0[Dur_Hour]:02d}&'.format(forminfo)
    data+='Dur_Min={0[Dur_Min]:02d}&Image+Type=GIF&COLOR='.format(forminfo)
    data+='COLOR&AE+Sensitivity=0&ASY%2FSYM++Sensitivity=0&'
    data+='Output={0[Output]}&Out+format='.format(forminfo)
    data+='IAGA2002&Email={0[Email]}'.format(forminfo)
    #print good
    #print target+'?'+data
    
    # Fetch data from web.
    f = urllib.request.urlopen(target+'?'+data) 
    lines = f.readlines()
    return(lines)

#===========================================================================

def aefetch(t_start, t_stop):
    '''
    A function to fetch Kyoto AE directly from the Kyoto WDC website.
    Returns raw ascii lines obtained from the website representing the data
    in IAGA 2002 format.  

    Note that, because of the higher density of the data, full date times
    are required for the start and stop time.  This is in contrast to
    other fetch functions, where hourly or three-hourly resolution allows
    for only the start and stop years and months are required as integers.

    Parameters
    ----------
    t_start : datetime.datetime
       Start time for range of requested data as a datetime object.
    t_stop : datetime.datetime
       Stop time for range of requested data as a datetime object.

    '''
    try:
        import urllib.parse, urllib.request
    except ImportError: #glue the python 3 names onto python 2
        import urllib
        urllib.parse = urllib
        urllib.request = urllib

    # Get difference between start and stop time:
    deltaT = t_stop - t_start
        
    # Set up all form info.
    forminfo = {} # Date-time stuff - requires special formatting.

    # GENERATE FORM INFO:
    # First three digits of year, last digit of year:
    forminfo['Tens'] = int(t_start.year/10)
    forminfo['Year'] = int(t_start.year - 10*forminfo['Tens'])
    # Month:
    forminfo['Month']= t_start.month
    # Day (right now, locked in to one month of data):
    forminfo['Day_Tens'] = int(t_start.day/10)
    forminfo['Days']     = int(t_start.day - 10*forminfo['Day_Tens'])
    # Hour/minute:
    forminfo['Hour'] = t_start.hour
    forminfo['min']  = t_start.minute
    # Duration:
    forminfo['Dur_Day_Tens'] = int(deltaT.days/10)
    forminfo['Dur_Day']      = int(deltaT.days - 10*forminfo['Dur_Day_Tens'])
    forminfo['Dur_Hour']     = int(deltaT.seconds/3600.)
    forminfo['Dur_Min']      = int(deltaT.seconds/60.-forminfo['Dur_Hour']*60.)

    # Email and other information:
    forminfo['Output']     = 'AE'
    forminfo['Out+format'] = 'IAGA2002'
    forminfo['Email']      =  urllib.parse.quote('spacepy@lanl.gov')

    # Keep this here for reference.  After testing, remove.
    #good='http://wdc.kugi.kyoto-u.ac.jp/cgi-bin/aeasy-cgi?Tens=201&Year=0&Month=04&Day_Tens=0&Days=0&Hour=00&min=00&Dur_Day_Tens=03&Dur_Day=0&Dur_Hour=00&Dur_Min=00&Image+Type=GIF&COLOR=COLOR&AE+Sensitivity=0&ASY%2FSYM++Sensitivity=0&Output=AE&Out+format=IAGA2002&Email=dog%40fart.crap'

    # Build and open URL.  Note that the forminfo must be
    # in a certain order for the request to work...
    target='http://wdc.kugi.kyoto-u.ac.jp/cgi-bin/aeasy-cgi'
    data='Tens={0[Tens]:03d}&Year={0[Year]:d}&'.format(forminfo)
    data+='Month={0[Month]:02d}&Day_Tens={0[Day_Tens]:d}&'.format(forminfo)
    data+='Days={0[Days]:d}&Hour={0[Hour]:02d}&min={0[min]:02d}'.format(forminfo)
    data+='&Dur_Day_Tens={0[Dur_Day_Tens]:02d}&'.format(forminfo)
    data+='Dur_Day={0[Dur_Day]:d}&Dur_Hour={0[Dur_Hour]:02d}&'.format(forminfo)
    data+='Dur_Min={0[Dur_Min]:02d}&Image+Type=GIF&COLOR='.format(forminfo)
    data+='COLOR&AE+Sensitivity=0&ASY%2FSYM++Sensitivity=0&'
    data+='Output=AE&Out+format='.format(forminfo)
    data+='IAGA2002&Email={0[Email]}'.format(forminfo)
    #print good
    #print target+'?'+data
    
    # Fetch data from web.
    f = urllib.request.urlopen(target+'?'+data) 
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
        raise Exception('Data is not in IAGA-2002 format.')

    # Parse mandatory IAGA header lines.
    source=(lines.pop(0)).split()[1]
    lines.pop(0)
    code=(lines.pop(0)).split()[2]
    for i in range(8):
        lines.pop(0)
    
    # Check Iaga Code as necessary.
    if iagacode:
        if iagacode != code:
            raise Exception("IAGA Code does not match required code.")

    # Loop through and count optional header lines.
    nHead=12
    while True:
        line=lines.pop(0)
        if line[:2]!=' #': break
        nHead+=1

    # Parse column labels.  We don't need time or DOY.
    parts=line.lower().split()[3:-1]

    # 
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
