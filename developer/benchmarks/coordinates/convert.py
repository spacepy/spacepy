import datetime as dt
import warnings
import numpy as np
import spacepy.coordinates as spc
import spacepy.time as spt

warnings.filterwarnings("ignore")

npos = 20000
pos = np.ones((npos, 3))
st_tai = 2000000000
tt_same  = spt.Ticktock(np.ones(npos) + st_tai, 'TAI')
tt_close = spt.Ticktock(np.arange(npos) + st_tai, 'TAI')
tt_space = spt.Ticktock(np.arange(npos) * 240 + st_tai, 'TAI')

print('Each test uses {} positions to transform'.format(npos))

t0 = dt.datetime.now()
# convert same-time coords, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_same)
cc_close.convert('GSE', 'car')
print('SAME-TIME GEO->GSE (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert same-time coords, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_same, use_irbem=False)
cc_close.convert('GSE', 'car')
print('SAME-TIME GEO->GSE (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert same-time coords, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_same)
cc_close.convert('MAG', 'sph')
print('SAME-TIME GEO->MAG (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert same-time coords, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_same, use_irbem=False)
cc_close.convert('MAG', 'sph')
print('SAME-TIME GEO->MAG (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert closely spaced times, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close)
cc_close.convert('GSE', 'car')
print('CLOSE GEO->GSE (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert closely spaced times, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close, use_irbem=False)
cc_close.convert('GSE', 'car')
print('CLOSE GEO->GSE (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert closely spaced times, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close)
cc_close.convert('GSM', 'car')
print('CLOSE GEO->GSM (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert closely spaced times, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close, use_irbem=False)
cc_close.convert('GSM', 'car')
print('CLOSE GEO->GSM (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close)
cc_close.convert('GSE', 'car')
print('SPACED GEO->GSE (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, no magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close, use_irbem=False)
cc_close.convert('GSE', 'car')
print('SPACED GEO->GSE (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close)
cc_close.convert('GSM', 'car')
print('SPACED GEO->GSM (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, magnetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close, use_irbem=False)
cc_close.convert('GSM', 'car')
print('SPACED GEO->GSM (SPA): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, geodetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close)
cc_close.convert('GDZ', 'sph')
print('SPACED GEO->GDZ (IRB): {}'.format(dt.datetime.now()-t0))

t0 = dt.datetime.now()
# convert widely spaced times, geodetic systems
cc_close = spc.Coords(pos, 'GEO', 'car', ticks=tt_close, use_irbem=False)
cc_close.convert('GDZ', 'sph')
print('SPACED GEO->GDZ (SPA): {}'.format(dt.datetime.now()-t0))