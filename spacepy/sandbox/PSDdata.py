#!/usr/bin/env python
# -*- coding: utf-8 -*-


def get_PSD(ticks, MU=1051, K=0.005, sats=None, _query=None):
    """
    get_PSD(ticks, mu, k, sats)
    
    Get the PSD from an input mu, k, date range, and list of sats from the psd database

    ** Warning ** With the sqlite backend this is not screaming fast.  Be a little patient
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - ticks - a list of start and stop ticks exclusive (Ticktock objects)
    - MU (optional) - a single mu value (default 1051)
    - K (optional) - a single k value (default 0.005)
    - sats (optional) - a list of sats to return data from 
        possible values: =['1990-095', '1991-080', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A']
    - _query (optional) - specify the query to take place (for the helper functions)
    
    Returns:
    ========
    - a dictionary of the form {'Ticks':time, 'Lstar':lstar, 'PSD':psd, 'Sat':sats}
    
    Example:
    ========
    
    See Also:
    =========
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
    Version:
    ========
    V1: 02-Jun-2010 (BAL)
    """
    try:
        from sqlalchemy.sql import select
    except ImportError:
        raise(ImportError("SQLAlchemy not installed, install it and then try again"))
    from sqlalchemy import create_engine, func
    from sqlalchemy import MetaData, Column, Table, ForeignKey
    from sqlalchemy import Integer, String, Float
    import sqlalchemy
    from sqlalchemy.orm import mapper, sessionmaker
    from sqlalchemy import Column, Integer, String, DateTime, BigInteger, Boolean, Date, Float, Table
    from sqlalchemy import desc, and_
    import os.path
    import datetime
    import subprocess      
    import os
    import spacepy.toolbox as tb
    import spacepy.time
    import numpy as n

    #####################################
    ###  Input checking #################
    #####################################
    # ticks a ticktock object with 2 elements?
    
    try:
        if len(ticks) != 2:
            raise(ValueError('ticks must be 2 element list of start stop ticks'))
    except:
        raise(ValueError('ticks must be 2 element list of start stop dates')) 
  
    ## check to see if the db exists
    file_local = os.path.isfile("psd_dat.sqlite")
    file_dot = os.path.isfile(os.environ['HOME']+'/.spacepy/data/psd_dat.sqlite')
    
    if not (file_local or file_dot):
        raise(Exception("Must have the database either in the current directory or your ~/.spacepy/data directory. \n \
        If you are at LANL, you can download a version of the database on the LANL only site:\n \
        https://sf4.lanl.gov/sf/docman/do/downloadDocument/projects.spacepy/docman.root/doc142771/1"))

    if file_local:
        db_loc = 'psd_dat.sqlite'
    else:
        db_loc = os.environ['HOME']+'/.spacepy/data/psd_dat.sqlite'

    ## setup the engine to the database
    engine = sqlalchemy. create_engine('sqlite:///' + db_loc, echo=False)
    ## this holds the metadata, tables names, attributes, etc
    metadata = sqlalchemy.MetaData(bind=engine)
    ## a session is what you use to actually talk to the DB, set one up with the current engine
    Session = sessionmaker(bind=engine)
    session = Session()

    ## ask for the table names form the database (does not grab views)
    table_names = engine.table_names()

    ## create a dictionary of all the table names that will be used as calss names.
    ## this uses the db table name as the tabel name and a cap 1st letter as the class
    ## when interacting using python use the class
    table_dict = {}
    for val in table_names:
        table_dict[val[0].upper() + val[1:]] = val

    ##  dynamincally create all the classes (c1)
    ##  dynamicallly create all the tables in the db (c2)
    ##  dynaminically create all the mapping between class and table (c3)
    ## this just saves a lot of typing and is equilivant to:
    ##     class Missions(object):
    ##         pass
    ##     missions = Table('missions', metadata, autoload=True)
    ##     mapper(Missions, missions)
    for val in table_dict:
        c1 = compile("""class %s(object):\n\tpass""" % (val), '', 'exec')
        c2 = compile("%s = Table('%s', metadata, autoload=True)" % (str(table_dict[val]), table_dict[val]) , '', 'exec')
        c3 = compile("mapper(%s, %s)" % (val, str(table_dict[val])), '', 'exec')
        exec(c1)
        exec(c2)
        exec(c3)



    obstime = []
    psd = []
    lstar = []
    sats_l = []
    mu_l = []
    k_l = []


    if _query == 'availablesats':
        s = select([func.distinct(Psd.sat)], ticks.UTC[1] > Psd.time > ticks.UTC[0])
        ans = s.execute()
        for val in ans:
            sats_l.append(val)
        ret = {'sat': sats_l}
        
    if _query == 'availablemu':
        s = select([func.distinct(Psd.mu)], ticks.UTC[1] > Psd.time > ticks.UTC[0])
        ans = s.execute()
        for val in ans:
            mu_l.append(val)
        ret = {'MU': mu_l}
        
    if _query == 'availablek':
        s = select([func.distinct(Psd.k)], ticks.UTC[1] > Psd.time > ticks.UTC[0])
        ans = s.execute()
        for val in ans:
            k_l.append(val)
        ret = {'K': k_l}
            
    if _query == None:
        ## session.query is what does the query
        ans = session.query(Psd.time, Psd.lstar, Psd.psd, Psd.sat).filter_by(mu = MU) \
            .filter_by(k=K).filter(Psd.time > ticks.UTC[0]).filter(Psd.time < ticks.UTC[1]) \
            .order_by(Psd.time).all()
        for val in ans:
            obstime.append(val[0])
            lstar.append(val[1])
            psd.append(val[2])
            sats_l.append(val[3])
            mu_l.append(MU)
            k_l.append(K)
        
        if len(obstime) == 0:
           # no obs, return empy dictionary
           return {}
        
        Ticks = spacepy.time.Ticktock(obstime, 'UTC')
        ret = {'Ticks':Ticks, 'Lstar':n.array(lstar), 'PSD':n.array(psd), \
            'sat':n.array(sats_l), 'MU':n.array(mu_l), 'K':n.array(k_l)}

        # filter out unwanted s/c
        if sats: # then specific s/c provided, remove unwanted ones
            idx = n.array([], dtype='int')
            for satid in sats:
                idx = n.append(idx, n.where(ret['sat'] == satid)[0])

            for key in list(ret.keys()):
                ret[key] = ret[key][idx]

        # sort in time
        idx = ret['Ticks'].argsort()
        for key in list(ret.keys()):
            ret[key] = ret[key][idx]
    
    ## annoyingly in the interest of getting this done not sure how to select on the sats specified so do that on this dict
    return ret

def get_PSD_availablesats(ticks):
    """
    get_PSD_availablesats(ticks)
    
    Get a list of sats with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - ticks - start and stop dates as a Ticktock object of len=2
    
    Returns:
    ========
    - a list of the sats that have data for thje gioven time period
    
    Example:
    ========
    
    See Also:
    =========
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
    Version:
    ========
    V1: 02-Jun-2010 (BAL)
    """
    import numpy as np
    ans = get_PSD(ticks, _query='availablesats')
    return ans['sat']

def get_PSD_availablemu(ticks):
    """
    get_PSD_availablemu(ticks)
    
    Get a list of mu values with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - ticks - start and stop dates as a Ticktock object of len=2
    
    Returns:
    ========
    - a list of the mu values that have data for the given time period
    
    Example:
    ========
    
    See Also:
    =========
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
    Version:
    ========
    V1: 02-Jun-2010 (BAL)
    """
    import numpy as np
    ans = get_PSD(ticks, _query='availablemu')
    return ans['MU']

def get_PSD_availablek(ticks):
    """
    get_PSD_availablek(ticks)
    
    Get a list of k values with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - ticks - start and stop dates as a Ticktock object of len=2
    
    Returns:
    ========
    - a list of the k values that have data for the given time period
    
    Example:
    ========
    
    See Also:
    =========
    
    Author:
    =======
    Brian Larsen, Los Alamos National Lab (balarsen@lanl.gov)
    
    Version:
    ========
    V1: 02-Jun-2010 (BAL)
    """
    import numpy as np
    ans = get_PSD(ticks, _query='availablek')
    return ans['K']


if __name__ == "__main__":
    from datetime import datetime
    import spacepy.time
    from pylab import *
    import numpy as np
    ticks = spacepy.time.Ticktock([datetime(2002, 8, 1), datetime(2002, 8, 2)], 'UTC')
    mu = 462
    k = 0.03
    # sats=['1990-095', '1991-080', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A']
    ans = get_PSD_availablek(ticks)
    print(ans)
    ans = get_PSD_availablemu(ticks)
    print(ans)
    ans = get_PSD_availablesats(ticks)
    print(ans)
    
    ans = get_PSD(ticks, mu, k, sats=['GPS-ns41', 'LANL-02A'])
    
    semilogy(ans['Ticks'].UTC, ans['PSD'])
    ax = gca()
    ax.set_ylabel('PSD')
    print('Data from these sats:')
    print(np.unique(ans['sat']))



