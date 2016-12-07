
import datetime
import time
import timeit

import numpy as np
import spacepy.time as spt

number = 200

gps = np.linspace(6.96556813e+08, 6.97334413e+08, 10)

tt = spt.Ticktock(gps, 'GPS')
print(tt)

print tt







