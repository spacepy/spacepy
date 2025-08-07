#!/usr/bin/env python
'''
PyBats submodule for handling input/output for the CIMI Model
, one of the choices for the IM module in the SWMF.
'''

# Global imports:
import numpy as np
import datetime as dt
import spacepy.plot.apionly
from spacepy.plot import set_target
from spacepy.pybats import PbData
from spacepy.datamodel import dmarray




class CimiEq(PbData):
    '''
    Class for loading a CIMIeq.out file.
    At instantiation time, user may wish to set the start date and time of 
    the simulation using the starttime kwarg.  If not given, start time
    will default to Jan. 1st, 2000, 00:00UT.
    '''
    def __init__(self, filename, starttime=None, *args, **kwargs):
        super(CimiEq, self).__init__(*args, **kwargs) # Init as PbData.
        self.attrs['file']=filename
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(starttime)

    def __repr__(self):
        return 'CimiEq single output file %s' % (self.attrs['file'])

    def _read(self, starttime):
        '''
        Read CimiEq file; should only be called upon instantiation.
        '''

        # Reads file assuming it is ascii
        # Slurp whole file.
        f=open(self.attrs['file'], 'r')
        lines=f.readlines()
        f.close()
        
        # Determine size of file.
        nTimes=lines.count(lines[0])
        nLat =int(lines[2].split()[0])
        nLon =int(lines[2].split()[1]) 

        self.attrs['nLat'] =nLat
        self.attrs['nLon'] =nLon
        self.attrs['nTime']=nTimes

        nDataLength = nLat*nLon
        
        # Start building time array.
        self['time']=np.zeros(nTimes, dtype=object)
        
        # Get variable names; pop g and rbody
        var=(lines[4].split())[0:-2]
        self._rawvar=var
        
        # initialize arrays to hold data
        for v in var:
            self[v]=dmarray(np.zeros((nTimes, nDataLength)))
            
            # Set units if possible (for now just set none
            self[v].attrs['units']=None
                    
        # Loop through rest of data to fill arrays.
        for i in range(nTimes):
            t=float((lines[i*(nDataLength+5)+1].split())[1])
            self['time'][i]=starttime+dt.timedelta(seconds=t)
            for j, l in enumerate(lines[i*(nDataLength+5)+5:(i+1)*(nDataLength+5)]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][i,j]=float(parts[k])
                            

    def add_slice(self, var, time, nlev=31, zlim=None, xlim=None, ylim=None,
                  target=None, 
                  loc=111, title=None,
                  add_cbar=True, clabel=None,
                  show_pts=False, dolog=False, add_body=True,  
                  figsize=(8.34,7),
                  *args, **kwargs):
        '''
        Create a plot of variable *var* at a given time.

        If kwarg **target** is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object, 
        the plot is placed into that axis.
        '''

        import matplotlib.pyplot as plt
        from matplotlib.colors import (LogNorm, Normalize)
        from matplotlib.patches import Circle
        from matplotlib.ticker import (LogLocator, LogFormatter, 
                                       LogFormatterMathtext, MultipleLocator)
        from matplotlib.patches import Wedge
        
        # Grab the slice of data that we want:
        if type(time) == type(self['time'][0]):
            if time not in self['time']:
                raise ValueError('Time not in object')
            time = np.arange(self.attrs['nTime'])[self['time']==time][0]

        fig, ax = set_target(target, loc=loc, figsize=figsize)
        ax.set_aspect('equal')

        # Grab values from correct time/location.
        x = self['x'][time, :]
        y = self['y'][time, :]
        value = self[var][time, :]

        # Get max/min if none given.
        if zlim==None:
            zlim=[0,0]
            zlim[0]=value.min(); zlim[1]=value.max()
            if dolog and zlim[0]<=0:
                zlim[0] = np.min( [0.0001, zlim[1]/1000.0] )

        #add points at "infinity" with values set zmin so that color
        #background fills in smoothly
        x=x.append(x,-15.0)
        x=x.append(x,-15.0)
        x=x.append(x,15.0)
        x=x.append(x,15.0)

        y=y.append(y,-15.0)
        y=y.append(y,15.0)
        y=y.append(y,-15.0)
        y=y.append(y,15.0)

        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])
        value=value.append(value,zlim[0])
        
        # Create levels and set norm based on dolog.
        if dolog:
            levs = np.power(10, np.linspace(np.log10(zlim[0]), 
                                            np.log10(zlim[1]), nlev))
            z=np.where(value>zlim[0], value, 1.01*zlim[0])
            norm=LogNorm()
            ticks=LogLocator()
            fmt=LogFormatterMathtext()
        else:
            levs = np.linspace(zlim[0], zlim[1], nlev)
            z=value
            norm=None
            ticks=None
            fmt=None
        
        # Create contour plot.
        cont=ax.tricontourf(np.asarray(x), np.asarray(y), np.asarray(z), \
                            np.asarray(levs), *args, norm=norm, **kwargs)
        if show_pts:
            ax.plot(x, y, '+w')

        if xlim != None:
            ax.set_xlim(xlim[0],xlim[1])

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])
            
        # Add cbar if necessary.
        if add_cbar:
            cbar=plt.colorbar(cont, ticks=ticks, format=fmt, pad=0.01)
            if clabel==None: 
                clabel="%s" % (var)
            cbar.set_label(clabel)
        else:
            cbar=None # Need to return something, even if none.
 
        # Set title, labels, axis ranges (use defaults where applicable.)
        if title: ax.set_title(title)
        #ax.set_yticks([]), ax.set_xticks([])
        ax.set_xlabel('X [Re]')
        ax.set_ylabel('Y [Re]')

        if add_body:
            angle=-90.0
            radius=1.0
            colors=('w','k')
            theta1, theta2 = angle, angle + 180
            center = (0,0)
            w1 = Wedge(center, radius, theta1, theta2, fc=colors[0], ec='k', **kwargs)
            w2 = Wedge(center, radius, theta2, theta1, fc=colors[1], ec='k', **kwargs)
            for wedge in [w1,w2]:
                ax.add_artist(wedge)

            
        return fig, ax, cont, cbar
    
    def interpolate(self, var, xpt, ypt, time):
        '''
        extract solution at given point and time
        '''
        from scipy.interpolate import griddata
        
        #find closest time to do the interpolation
        dt=self['time']-time 
        iTime=np.where(abs(dt) == np.min(abs(dt))) [0][0]
        
        #get number of points
        nPoints=len(self['x'][iTime,:]) 
        
        #create point pairs
        Points=np.zeros([nPoints,2]) 
        for i in range(nPoints): 
            Points[i][0]=self['x'][iTime,i] 
            Points[i][1]=self['y'][iTime,i]

        Values=self[var][iTime,:]
            
        InterpValue = griddata(Points, Values, (xpt, ypt),
                               method='linear',fill_value=0.0)
        
        return InterpValue

class CimiLog(PbData):
    '''
    A class for reading and plotting a complete set of CIMI log files.
    Open using a glob string that encompasses all of the log files that are
    intended to be read, e.g. 'CIMI*.log'.  Note that unix wildcards are
    accepted.
    '''

    def __init__(self, filepattern, starttime=None, *args, **kwargs):
        super(CimiLog, self).__init__(*args, **kwargs) # Init as PbData.
        if not starttime:
            starttime=dt.datetime(2000,1,1,0,0,0)
        self._read(filepattern, starttime)
        
    def __repr__(self):
        return 'CIMI log files:', self['files']

    def _read(self, filepattern, starttime):
        '''
        Read all ascii line files; should only be called upon instantiation.
        '''

        from glob import glob
        
        self['files']=sorted(glob(filepattern))



        #get total length of combined data
        nDataLength = 0
        for file in self['files']:
            # Reads file assuming it is ascii
            # Slurp whole file.
            f=open(file, 'r')
            lines=f.readlines()
            f.close()
            nDataLength = nDataLength+len(lines)-2
        
        
        iCounter = 0
        offset = 0
        for file in self['files']:
            # Reads file assuming it is ascii
            # Slurp whole file.
            f=open(file, 'r')
            lines=f.readlines()
            f.close()
            
            if iCounter == 0:
                # Get variable names; pop g and rbody
                var=(lines[1].split())
                self._rawvar=var
                # initialize arrays to hold data
                for v in var:
                    self[v]=dmarray(np.zeros((nDataLength)))

            #read in data
            for j,l in enumerate(lines[2:]):
                parts=l.split()
                for k,v in enumerate(var):
                    self[v][offset+j]=float(parts[k])

            iCounter=iCounter+1
            offset = offset +len(lines)-2
            
        # Start building time array.
        self['time']=np.zeros(nDataLength, dtype=object)

        for i in range(nDataLength):
            t= self['t'][i]
            self['time'][i]=starttime+dt.timedelta(seconds=t)
            
    def plotlog(self, vars, xlim=None, ylim=None,
                  target=None, 
                  loc=111, title=None,
                  dolog=False, 
                  figsize=(8.34,7),
                  *args, **kwargs):
        '''
        Create a plot of variable *var* at a given time.

        If kwarg **target** is None (default), a new figure is 
        generated from scratch.  If target is a matplotlib Figure
        object, a new axis is created to fill that figure at subplot
        location **loc**.  If **target** is a matplotlib Axes object, 
        the plot is placed into that axis.
        '''

        import matplotlib.pyplot as plt
        

        fig, ax = set_target(target, loc=loc, figsize=figsize)
        #ax.set_aspect('equal')

        # Grab values from correct time/location.
        for var in vars:
            x = self['time'][:]
            y = self[var][:]
            ax.plot(x,y,label=var)

        
        # Create levels and set norm based on dolog.
        if dolog:
            plt.yscale('log')

        #add legend
        ax.legend()



        if xlim != None:
            ax.set_xlim(xlim[0],xlim[1])

        if ylim != None:
            ax.set_ylim(ylim[0],ylim[1])
            

            
        return fig, ax
