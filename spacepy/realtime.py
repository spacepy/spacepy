#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The realtime module fetches realtime data and forecasts from other sites

Authors: Josef Koller

Institution: Los Alamos National Laboratory

Contact: jkoller@lanl.gov

Copyright 2012 Los Alamos National Security, LLC.


About realtime
--------------

bla bla

Examples
--------

badf asdfasdf asdf

    >>> import spacepy.realtime as rt

"""

__contact__ = 'Josef Koller, jkoller@lanl.gov'

import datetime
import urllib.request as u

import spacepy.datamodel as dm
import dateutil.parser as dup



# -------------------------------------
def now():
    """
    Return the current values for the parameters

    Returns
    -------
    out : datamodel.SpaceData
        A SpaceData containing the current values of the indices
    """
    raise(NotImplementedError('Not yet implemented (nudge nudge)'))


# -------------------------------------	
def nextForecast():
    """
    Return only the next value
    
    Returns
    -------
    out : datamodel.SpaceData
        A SpaceData containing the next forcast for the indices
    """
    dd = forecast()
    data = dm.SpaceData()
    data.attrs = dd.attrs
    data['forecast_hours'] = [1,3]
    data['AE'] = [dd['AE_1hr'][-1], dd['AE_3hr'][-1]]
    data['Kp'] = [dd['Kp_1hr'][-1], dd['Kp_3hr'][-1]]
    data['Dst'] = [dd['Dst_1hr'][-1], dd['Dst_3hr'][-1]]
    data['Calc'] = [dd['Calc_1hr'][-1], dd['Calc_3hr'][-1]]
    data['Epoch'] = [dd['Epoch_1hr'][-1], dd['Epoch_3hr'][-1]]
    return data


def _parseRICE(data, hours):
    """
    parse the data read form the Rice website
    
    .. warning: This is an internal function do not call directly

    Parameters
    ----------
    data : list
        data from the Rice website
    hours : str
        the forecast hours, used int h dict keys
    
    Returns
    -------
    out : datamodel.SpaceData
        datemnodel object of the parsed data
    
    """
    dd = dm.SpaceData()
    # strip the newlines
    data = [val.rstrip() for val in data]
    header = data.pop(0) # grab the header
    header = header.split()[4:]
    units = []
    for i, val in enumerate(header):
        if '(' in val:
            val = val.split('(')
            units.append(val[-1].split(')')[0])
            header[i] = val[0]
        else:
            units.append(None)
    data = dm.dmarray([val.split() for val in data])
    times = dm.dmarray([dup.parse('{0}{1:02d}{2:02d}T{3:04d}'.format(val[0], int(val[1]), int(val[2]), int(val[3]))) for val in data[:, 0:4]])
    data = dm.dmarray(data[:,4:], dtype=float)
    for i, (key, unit) in enumerate(zip(header, units)):
        dd[key+'_' + hours + 'hr'] = dm.dmarray(data[:,i], attrs={'units':unit})
    dd['Calc_' + hours + 'hr'] = times
    dd['Epoch_' + hours + 'hr'] = dd['Calc_' + hours + 'hr'] + datetime.timedelta(hours=int(hours))
    return dd


def forecast():
    RICE_URL_1_last = 'http://mms.rice.edu/realtime/Predictions_1.last'
    RICE_URL_3_last = 'http://mms.rice.edu/realtime/Predictions_3.last'
    RICE_Boyle_all = 'http://mms.rice.edu/realtime/File1.txt'
    RICE_1hr_Kp_Dst = 'http://mms.rice.edu/realtime/File2.txt'
    RICE_3hr_Kp_Dst = 'http://mms.rice.edu/realtime/File3.txt'

    # grab all the 1 hour data
    hr1 = u.urlopen(RICE_1hr_Kp_Dst)
    data = hr1.readlines()
    hr1.close()
    dd1 = _parseRICE(data, '1')

    # grab all the 3 hour data
    hr3 = u.urlopen(RICE_3hr_Kp_Dst)
    data = hr3.readlines()
    hr3.close()
    dd3 = _parseRICE(data, '3')

    dd = dm.SpaceData()
    dd.attrs['URL_1hr'] = RICE_1hr_Kp_Dst
    dd.attrs['URL_3hr'] = RICE_3hr_Kp_Dst
    dd.attrs['retrive_time'] = datetime.datetime.now()
    
    for key1, key3 in zip(dd1, dd3):
        dd[key1] = dd1[key1]
        dd[key3] = dd3[key3]
            
    return dd
  


