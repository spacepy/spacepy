

def get_PSD(dates, mu=1051, k=0.005,
            sats=['1990-095', '1991-080', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A'],
            _query=None):
    """
    get_PSD(dates, mu, k, sats=['1990-095', '1991-080', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A'])
    
    Get the PSD from an input mu, k, date range, and list of sats from the psd database

    ** Warning ** With the sqlite backend this is not screaming fast.  Be a little patient
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - dates - a list of start and stop dates exclusive (Ticktock objects)
    - mu (optional) - a single mu value (default 1051)
    - k (optional) - a single k value (default 0.005)
    - sats (optional) - a list of sats to returen data from **Currently not implemented**
    - _query (optional) - specify the query to take place (for the helper functions)
    
    Returns:
    ========
    - a dictionary of the form {'time':time, 'lstar':lstar, 'psd':psd, 'sat':sats}
    
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

    #####################################
    ###  Input checking #################
    #####################################
    # dates a ticktock object with 2 elements?
    try:
        if len(dates) != 2:
            raise(ValueError('Dates must be 2 element list of start stop dates'))
    except:
        raise(ValueError('Dates must be 2 element list of start stop dates')) 
  
    ## check to see if the db exists
    file_local = os.path.isfile("psd_dat.sqlite")
    file_dot = os.path.isfile(os.environ['HOME']+'/.spacepy/data/psd_dat.sqlite')
    
    if not (file_local or file_dot):
        raise(Exception("Must have the database either in the current directory or your ~/.spacepy/data directory"))

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



    time = []
    psd = []
    lstar = []
    sats_l = []
    mu_l = []
    k_l = []


    if _query == 'availablesats':
        s = select([func.distinct(Psd.sat)], dates[1] > Psd.time > dates[0])
        ans = s.execute()
        for val in ans:
            sats_l.append(val)
    if _query == 'availablemu':
        s = select([func.distinct(Psd.mu)], dates[1] > Psd.time > dates[0])
        ans = s.execute()
        for val in ans:
            mu_l.append(val)
    if _query == 'availablek':
        s = select([func.distinct(Psd.k)], dates[1] > Psd.time > dates[0])
        ans = s.execute()
        for val in ans:
            k_l.append(val)
    if _query == None:
        ## session.query is what does the query
        ans = session.query(Psd.time, Psd.lstar, Psd.psd, Psd.sat).filter_by(mu = mu).filter_by(k=k).filter(Psd.time > dates[0]).filter(Psd.time < dates[1]).order_by(Psd.time).all()    
        for val in ans:
            time.append(val[0])
            lstar.append(val[1])
            psd.append(val[2])
            sats_l.append(val[3])
            mu_l.append(mu)
            k_l.append(k)

    ## annoyingly in the interest of getting this done not sure how to select on the sats specified so do that on this dict
    
    ret = {'time':time, 'lstar':lstar, 'psd':psd, 'sat':sats_l, 'mu':mu_l, 'k':k_l}
    return ret

def get_PSD_availablesats(dates):
    """
    get_PSD_availablesats(dates)
    
    Get a list of sats with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - dates - a list of start and stop dates (datetime objects)
    
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
    ans = get_PSD(dates, _query='availablesats')
    return ans['sat']

def get_PSD_availablemu(dates):
    """
    get_PSD_availablemu(dates)
    
    Get a list of mu values with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - dates - a list of start and stop dates (datetime objects)
    
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
    ans = get_PSD(dates, _query='availablemu')
    return ans['mu']

def get_PSD_availablek(dates):
    """
    get_PSD_availablek(dates)
    
    Get a list of k values with data for the given time range
    (Databease created 2-Jun-2010 by Brian Larsen)
    
    Input:
    ======
    - dates - a list of start and stop dates (datetime objects)
    
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
    ans = get_PSD(dates, _query='availablek')
    return ans['k']


if __name__ == "__main__":
    from datetime import datetime
    from pylab import *
    import numpy as np
    dates = [datetime(2005, 1, 1), datetime(2005, 7, 1)]
    mu = 462
    k = 0.03
    ## sats=['1990-095', '1991-080', 'GPS-ns41', 'LANL-01A', 'LANL-02A', 'LANL-97A']
    ans = get_PSD_availablek(dates)
    print ans
    ans = get_PSD_availablemu(dates)
    print ans
    ans = get_PSD_availablesats(dates)
    print ans
    
    ans = get_PSD(dates, mu, k)

    semilogy(ans['time'], ans['psd'])
    ax = gca()
    ax.set_ylabel('PSD')
    print('Data from these sats:')
    print(np.unique(ans['sat']))



