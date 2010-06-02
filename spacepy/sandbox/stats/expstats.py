#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A class for exploratory statistics

"""

from spacepy import help
__version__ = "$Revision: 1.1 $, $Date: 2010/06/02 17:14:29 $"
__log__ = """
$Log: expstats.py,v $
Revision 1.1  2010/06/02 17:14:29  balarsen
initial addition of exploratory stats class, this is the NIST exploratory suggestions


"""
__author__ = 'Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)'


# -----------------------------------------------
# expstats class
# -----------------------------------------------    
class Expstats():
    """
    a = Expstats(x, y)
    
    A class to perform exploratory statistics 
        
    Input:
    ======
        - x : independent coordinate points 
        - y : dependent coordinate points 
        
    Returns:
    ========
        - instance with a.lagplot, etc
    
    Example:
    ========  
    
       
    See Also:
    =========
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
    Version:
    ========
    V1: 19-May-2010 (BAL)
    
    """
    
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def lagplot(self,fignum=1, subplot=111, **pltkeywords):
        import matplotlib.pyplot as plt
        fig=plt.figure(fignum)
        ax = fig.add_subplot(111)
        pp = ax.plot(self.y[:-1], self.y[1:], 'x', **pltkeywords)
        ax.set_xlabel('Y(i-1)')
        ax.set_ylabel('Y(i)')
        ax.set_title('Lag Plot')
        return pp
        
    def hist(self,fignum=1, subplot=111,  **pltkeywords):
        import matplotlib.pyplot as plt
        import spacepy.toolbox as tb
        bins = tb.binHisto(self.y)
        fig=plt.figure(fignum)
        ax = fig.add_subplot(111)
        pp = ax.hist(self.y, bins[1], **pltkeywords)
        ax.set_xlabel('Y')
        ax.set_ylabel('Counts')
        ax.set_title('Histogram')
        return pp

    def runSequence(self, fignum=1, subplot=111, **pltkeywords):
        import matplotlib.pyplot as plt
        fig=plt.figure(fignum)
        ax = fig.add_subplot(111)
        pp = ax.plot(self.y, 'x', **pltkeywords)
        ax.set_xlabel('Index')
        ax.set_ylabel('Y')
        ax.set_title('Run Sequence Plot')
        return pp

    def probplot(self, sparams=(), dist='norm', fit=1, fignum=1, subplot=111, **pltkeywords):
        import scipy.stats as ss
        import numpy as np
        import matplotlib.pyplot as plt
        ## so far only functions on normal distributions
        prob = ss.probplot(self.y, sparams=sparams, dist=dist, fit=fit)
        ## make a line from the regression best fit
        xfit = np.arange(np.floor(np.min(prob[0][0])), np.ceil(np.max(prob[0][1])), 0.25)
        yfit = xfit*[prob[1][0]] + prob[1][1]
        fig = plt.figure(fignum)
        ax = fig.add_subplot(111)
        pp = ax.plot(prob[0][0], prob[0][1], 'x', **pltkeywords)
        ax.set_xlabel('Statistic medians (' + dist + ')')
        ax.set_ylabel('Ordered response')
        ax.set_title('Probability plot (' + dist + ')')
        pp1 = ax.plot(xfit, yfit)
        return pp, pp1

    def fourplot(self, fignum=1, clear=True):
        import matplotlib.pyplot as plt
        fig=plt.figure(fignum)
        if clear:
            plt.clf()
        ax=fig.add_subplot(self.runSequence(subplot=221))
        ax=fig.add_subplot(self.lagplot(subplot=222))
        ax=fig.add_subplot(self.hist(subplot=223))
        ax=fig.add_subplot(self.probplot(subplot=224))

