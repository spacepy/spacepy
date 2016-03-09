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


Authors
-------
Josef Koller, Steve Morley 

Copyright 2010 Los Alamos National Security, LLC.
"""


from .irbempy import *

