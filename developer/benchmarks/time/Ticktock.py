
import datetime
import time
import timeit

import numpy as np
import spacepy.time as spt

number = 200

gps = np.linspace(6.96556813e+08, 6.97334413e+08, 1000)

tt = spt.Ticktock(gps, 'GPS')

print timeit.timeit(stmt= "tt = spt.Ticktock(gps, 'GPS')",
              setup='import numpy as np; import spacepy.time as spt; gps = np.linspace(6.96556813e+08, 6.97334413e+08, 1000)', number=number)/float(number)

# vanilla   0.0531055498123
# ldatetime=datetime.datetime 0.0437979602814  1.2x
# bisect 0.00696686983109 7.6x









