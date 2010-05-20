#!/usr/bin/env python
"""
Spacecraft data functions
"""

__version__ = "$Revision: 1.1 $, $Date: 2010/05/20 17:19:44 $"
__author__ = 'J. Koller, Los Alamos National Lab (jkoller@lanl.gov)'

# load omni file during import
from spacepy import __path__
from spacepy import loadpickle
PSDfln = __path__[0]+'/data/PSDdb.pbin'
PSDdb = loadpickle(PSDfln)

# -----------------------------------------------
# PSDdata class
# -----------------------------------------------    
class PSDdataClass(object):
    """

    Author:
    =======
    Josef Koller, Los Alamos National Lab (jkoller@lanl.gov)
    
    Version:
    ========
    V1: 17-Mar-2010 (JK)
    
    """
    
    def __init__(self, MUval, Kval, satname):
        self.MU = MUval
        self.K = Kval
        self.satname = satname
        
        return
        
    # -----------------------------------------------    
    def __str__(self):
        """
        """
 
        return '<PSDdataClass instance satname='+self.satname+'; (MU,K)=('+str(self.MU)+';'+str(self.K)+')>'
    __repr__ = __str__
    

        
# -----------------------------------------------    
def get_obs(tickstart, tickend, MU, K, satlist=None):
    """
    return all observations between tickstart and tickend
    possible MU values: 167, 462, 1051, 2083
    possible K values: 0.01, 0.03, 0.1, 0.3, 1.0
    
    """
    
    import numpy as n
    import spacepy.time as st
    import spacepy.satdata as sat
    
    
    if isinstance(satlist, str):
        satlist = [satlist]
    
    satkeys = sat.PSDdb.keys()
    if satlist is None:
        satlist = satkeys
       
    PSD = []
    Lstar = []
    
    satstring = ''
    for i, key in enumerate(satlist):
        idx1 = n.where((sat.PSDdb[key]['MUKgrid'] == [MU,K]).all(axis=1))[0]
        print idx1, len(idx1)
        if len(idx1) > 0: # then no matching MU,K combination found for this s/c
            ticks = sat.PSDdb[key]['PSDdataClass'][idx1].ticktock
            idx2 = n.where( (tickstart < ticks) & (ticks < tickend) )[0]
            PSD = n.append(PSD, sat.PSDdb[key]['PSDdataClass'][idx1].PSD[idx2])
            Lstar = n.append(Lstar, sat.PSDdb[key]['PSDdataClass'][idx1].Lstar[idx2])
            if len(idx2) > 0:
                satstring = ' '.join([satstring, key])
        
    if len(Lstar) is 0:
        Lstar = None
        PSD = None
        
        
    # add observational error    
    obserr = n.array([0.31]*len(PSD))
    
    # setup PSDdataClass class instance
    satstring = satstring.strip()
    dataobj = PSDdataClass(MU, K, satstring)
    dataobj.Lstar = Lstar
    dataobj.PSD = PSD
    dataobj.PSDerr = PSD*0.31
    dataobj.ticktock = st.Ticktock([tickstart.UTC, tickend.UTC], 'UTC')        
    
    return dataobj
# -----------------------------------------------    
def convert2PSDclass(fln, MUval, Kval, satname):
    """
    """
    
    import numpy as n
    import spacepy.time as st
    import spacepy.satdata as sat
    import datetime, time
    
    # read file
    fh = open(fln)
    alltext = fh.readlines()
    fh.close()
    
    nlines = len(alltext)
    for line in alltext:
        if line.startswith('\n') or line.startswith(';'):
            nlines -= 1
            
    ISOdate = ['']*nlines
    UTC = ['']*nlines
    Lstar = n.zeros(nlines)
    PSD = n.zeros(nlines)
    i = 0
    for line in alltext:
        if (not line.startswith('\n')) and (not line.startswith(';')):            
            iTAI, ISO, iLstar, iPSD = line.split()   
            ISOdate[i] =  ISO
            Lstar[i] = float(iLstar)
            PSD[i] = float(iPSD)
            i += 1
            
    # convert time information
    for i, ISO in enumerate(ISOdate):
        foo = time.strptime(ISO, '%Y%m%d%H%M%S')
        UTC[i] = datetime.datetime(foo.tm_year, foo.tm_mon, foo.tm_mday, \
            foo.tm_hour, foo.tm_min, foo.tm_sec)
    # create PSD class instance and add stuff
    obj = sat.PSDdataClass(MUval, Kval, satname)
    obj.ticktock = st.Ticktock(UTC, 'UTC')
    obj.Lstar = Lstar
    obj.PSD = PSD
    
    return obj
    
    
# -----------------------------------------------    
def updatePSDdb(path, satlist=None):
    """
    calls makePSDdb and updates all
	example: updatePSDdb('/Users/jkoller/Research/input_data/')
    """
    
    import spacepy.toolbox as tb
    import numpy as n
    import spacepy.time as st
    import glob

    if satlist is None:
        satlist = ['Polar', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A', '1990-095', '1991-080']

    if isinstance(satlist, str):
        satlist = [satnlist]
  
    db = {}
    for satname in satlist:
        db[satname] = {}
        flist = glob.glob(path+satname+'*.dat')
        obj = ['']*len(flist)
        MUKgrid = n.zeros((len(flist),2))   
        # find unique  MU, K
        for i, fln in enumerate(flist):
            MU = int(fln.split('_')[4].lstrip('MU'))
            Kstring = fln.split('_')[5].lstrip('K').rstrip('.dat').replace('P','.')
            K = float(Kstring)
            MUKgrid[i] = [MU, K]    
        uniqMU = n.unique(MUKgrid[:,0])
        uniqK = n.unique(MUKgrid[:,1])        
        # work through all unique combinations
        i = -1
        for iMU in uniqMU:
            for iK in uniqK:
                Kstring = str(iK).rstrip('0').replace('.','P')
                MUstring = str(int(iMU))
                fln = glob.glob(path+satname+'*_MU'+MUstring+'_K'+Kstring+'.dat')
                print "reading "+str(fln)
                i += 1            
                MUKgrid[i] = [iMU, iK]
                if len(fln) is 1:
                    obj[i] = convert2PSDclass(fln[0], iMU, iK, satname)                
                elif len(fln) > 1:
                    obj[i] = convert2PSDclass(fln[0], iMU, iK, satname) 
                    for ifln in fln[1:]:
                        foo = convert2PSDclass(ifln, iMU, iK, satname)
                        obj[i].PSD = n.append(obj[i].PSD, foo.PSD)
                        obj[i].Lstar = n.append(obj[i].Lstar, foo.Lstar)
                        obj[i].ticktock = obj[i].ticktock.append(foo.ticktock)    
        # sort in time
        MUKgrid = MUKgrid[:i+1]
        obj = obj[:i+1]        
        # sort in time
        for i in range(len(obj)):
            idx = obj[i].ticktock.UNX.argsort()
            obj[i].Lstar = obj[i].Lstar[idx]
            obj[i].PSD = obj[i].PSD[idx]
            obj[i].ticktock = obj[i].ticktock[idx]
        db[satname]['MUKgrid'] = MUKgrid
        db[satname]['PSDdataClass'] = obj
        

    tb.savepickle(PSDfln, db)
    
    return
