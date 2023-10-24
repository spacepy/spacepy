# -*- coding: utf-8 -*-

"""
module wrapper for irbem_lib
Reference for this library
https://sourceforge.net/projects/irbem/
D. Boscher, S. Bourdarie, P. O'Brien, T. Guild, IRBEM library V4.3, 2004-2008


Most functions in this module use an options list to define the models used and the settings
that define the quality level of the result. The options list is a 5-element list and is defined
as follows.

Options
-------
   - 1st element: 0 - don't compute L* or phi ;  1 - compute L*; 2- compute phi
   - 2nd element: 0 - initialize IGRF field once per year (year.5);  
       n - n is the  frequency (in days) starting on January 1st of each year 
       (i.e. if options(2nd element)=15 then IGRF will be updated on the following days of the 
       year: 1, 15, 30, 45 ...)
   - 3rd element: resolution to compute L* (0 to 9) where 0 is the recomended value to ensure a 
       good ratio precision/computation time
       (i.e. an error of ~2% at L=6). The higher the value the better will be the precision, the 
       longer will be the computing time. Generally there is not much improvement for values 
       larger than 4. Note that this parameter defines the integration step (theta) 
       along the field line such as dtheta=(2pi)/(720*[options(3rd element)+1])
   - 4th element: resolution to compute L* (0 to 9). The higher the value the better will be 
       the precision, the longer will be 
       the computing time. It is recommended to use 0 (usually sufficient) unless L* is not 
       computed on a LEO orbit. For LEO orbit higher values are recommended. Note that this 
       parameter defines the integration step (phi) along the drift shell such as 
       dphi=(2pi)/(25*[options(4th element)+1])
   - 5th element: allows to select an internal magnetic field model (default is set to IGRF)
       - 0 = IGRF
       - 1 = Eccentric tilted dipole
       - 2 = Jensen&Cain 1960
       - 3 = GSFC 12/66 updated to 1970
       - 4 = User-defined model (Default: Centred dipole + uniform [Dungey open model] )
       - 5 = Centred dipole

The routines also require specification of the external magnetic field model. The default is the
Tsyganenko 2001 storm-time model. The external model is always specified using the extMag keyword
and the following options exist.

extMag
------
    - '0' = No external field model
    - 'MEAD' = Mead and Fairfield
    - 'T87SHORT' = Tsyganenko 1987 short (inner magnetosphere)
    - 'T87LONG' = Tsyganenko 1987 long (valid in extended tail region)
    - 'T89' = Tsyganenko 1989
    - 'OPQUIET' = Olsen-Pfitzer static model for quiet conditions
    - 'OPDYN' = Olsen-Pfitzer static model for active conditions
    - 'T96' = Tsyganenko 1996
    - 'OSTA' = Ostapenko and Maltsev
    - 'T01QUIET' = Tsyganenko 2001 model for quiet conditions
    - 'T01STORM' = Tsyganenko 2001 model for active conditions
    - 'T05' = Tsyganenko and Sitnov 2005 model
    - 'ALEX' = Alexeev model
    - 'TS07' = Tsyganenko and Sitnov 2007 model

Many of these models have limits placed on the valid range of input parameters,
and outside these limits invalid (NaN) values will be returned.

    - MEAD	: Mead & Fairfield [1975] (uses 0<=Kp<=9 - Valid for rGEO<=17. Re)
    - T87SHORT: Tsyganenko short [1987] (uses 0<=Kp<=9 - Valid for rGEO<=30. Re)
    - T87LONG : Tsyganenko long [1987] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
    - T89	 : Tsyganenko [1989] (uses 0<=Kp<=9 - Valid for rGEO<=70. Re)
    - OPQUIET : Olson & Pfitzer quiet [1977] (default - Valid for rGEO<=15. Re)
    - OPDYN   : Olson & Pfitzer dynamic [1988] (uses 5.<=dens<=50., 300.<=velo<=500.,
        -100.<=Dst<=20. - Valid for rGEO<=60. Re)
    - T96	 : Tsyganenko [1996] (uses -100.<=Dst (nT)<=20., 0.5<=Pdyn (nPa)<10.,
        \\|ByIMF\\| (nT)<=10., \\|BzIMF\\| (nT)<=10. - Valid for rGEO<=40. Re)
    - OSTA	: Ostapenko & Maltsev [1997] (uses dst,Pdyn,BzIMF, Kp)
        T01QUIET: Tsyganenko [2002a,b] (uses -50.<Dst (nT)<20., 0.5<Pdyn (nPa)<=5.,
        \\|ByIMF\\| (nT)<=5., \\|BzIMF\\| (nT)<=5., 0.<=G1<=10., 0.<=G2<=10. - Valid for xGSM>=-15. Re)
    - T01STORM: Tsyganenko, Singer & Kasper [2003] storm  (uses Dst, Pdyn, ByIMF, BzIMF, G2, G3 -
        there is no upper or lower limit for those inputs - Valid for xGSM>=-15. Re)
    - T05	 : Tsyganenko & Sitnov [2005] storm  (uses Dst, Pdyn, ByIMF, BzIMF,
        W1, W2, W3, W4, W5, W6 - no upper or lower limit for inputs - Valid for xGSM>=-15. Re)
    - TS07   : Tsyganenko and Sitnov [2007] model. Uses specially calculated coefficient files.


Authors
-------
Josef Koller, Steve Morley 

Copyright 2010 Los Alamos National Security, LLC.
"""


from .irbempy import *

