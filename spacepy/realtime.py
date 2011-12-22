#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The realtime module fetches realtime data and forecasts from other sites

Authors: Josef Koller

Institution: Los Alamos National Laboratory

Contact: jkoller@lanl.gov

Copyright ©2012 Los Alamos National Security, LLC.


About realtime
--------------

bla bla

Examples
--------

badf asdfasdf asdf

    >>> import spacepy.realtime as rt

"""

__contact__ = 'Josef Koller, jkoller@lanl.gov'

# -------------------------------------
def now():

    return


# -------------------------------------	
def fcst():

    from spacepy.time import Ticktock, Tickdelta
    from spacepy import savepickle, DOT_FLN
    import os, sys
 
    RICE_URL_1 = 'http://mms.rice.edu/realtime/Predictions_1.last'
    RICE_URL_3 = 'http://mms.rice.edu/realtime/Predictions_3.last'

 
    if sys.version_info[0]<3:
        import urllib as u
    else:
        import urllib.request as u

    # retrieve file
    rc1_fln = os.path.join(DOT_FLN, 'data', 'rice1.tmp')
    rc3_fln = os.path.join(DOT_FLN, 'data', 'rice3.tmp')
    u.urlretrieve(RICE_URL_1, rc1_fln)
    u.urlretrieve(RICE_URL_3, rc3_fln)

    # parse file
    dd={}
    dd['fcst_hrs'] = [1,3]
    yr, mo, day, hrmin, Kp, Dst, AE = open(rc1_fln).readline().split()
    t1 = Ticktock(yr+'-'+mo+'-'+day+'T'+hrmin[:2]+':'+hrmin[2:]+':00', 'ISO')
    dd['Kp'] = [float(Kp)]
    dd['Dst'] = [float(Dst)]
    dd['AE'] = [float(AE)]
    
    yr, mo, day, hrmin, Kp, Dst, AE = open(rc3_fln).readline().split()
    t3 = Ticktock(yr+'-'+mo+'-'+day+'T'+hrmin[:2]+':'+hrmin[2:]+':00', 'ISO')
    dd['tick_calc'] = [t1, t3]
    dd['tick_fcst'] = [t1+Tickdelta(hours=1), t3+Tickdelta(hours=3)]
    dd['Kp'].append(float(Kp))
    dd['Dst'].append(float(Dst))
    dd['AE'].append(float(AE))
    
        
    return dd
    


