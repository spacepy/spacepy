#!/usr/bin/env python

"""Calculate various approaches to UTC - TAI during rubber second era"""

import datetime

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy
import spacepy.pycdf
import spacepy.time

# These are values from USNO tai-utc during rubber-second era
# Piece-wise linear, calculated as a constant plus a rate times
# days since a particular day (so seconds/day). Days are UTC days
# The rate can change on 1 Jan and the constant can change any month;
# the days-since day is relative to the last rate change, not the
# last record
# Records are: UTC of start of validity, TAI-UTC at start, reference UTC-MJD,
# rate in seconds/MJD
rubbers = numpy.array([
    (datetime.datetime(1961, 1, 1), 1.422818, 37300.0, 0.001296),
    (datetime.datetime(1961, 8, 1), 1.372818, 37300.0, 0.001296),
    (datetime.datetime(1962, 1, 1), 1.845858, 37665.0, 0.0011232),
    (datetime.datetime(1963, 11, 1), 1.945858, 37665.0, 0.0011232),
    (datetime.datetime(1964, 1, 1), 3.24013, 38761.0, 0.001296),
    (datetime.datetime(1964, 4, 1), 3.34013, 38761.0, 0.001296),
    (datetime.datetime(1964, 9, 1), 3.44013, 38761.0, 0.001296),
    (datetime.datetime(1965, 1, 1), 3.54013, 38761.0, 0.001296),
    (datetime.datetime(1965, 3, 1), 3.64013, 38761.0, 0.001296),
    (datetime.datetime(1965, 7, 1), 3.74013, 38761.0, 0.001296),
    (datetime.datetime(1965, 9, 1), 3.84013, 38761.0, 0.001296),
    (datetime.datetime(1966, 1, 1), 4.31317, 39126.0, 0.002592),
    (datetime.datetime(1968, 2, 1), 4.21317, 39126.0, 0.002592),
    (datetime.datetime(1972, 1, 1), 10, 41317.0, 0.0),
    (datetime.datetime(1972, 7, 1), 11, 41317.0, 0.0),
    (datetime.datetime(1973, 1, 1), 12, 41317.0, 0.0),
    ],
    dtype=[('dt', 'O'), ('delta', 'f'), ('mjd', 'f'), ('rate', 'f')])


# Observed values of ut1 - tai
# Annual from 1955 (start of atomic time) through 1961 from
# https://stjarnhimlen.se/comp/time.html
# Every 3 mo start 1 Jan 1962, from
# https://www.iers.org/IERS/EN/Science/EarthRotation/UT1-TAI.html
ut1tai = numpy.array([
    1.114,
    0.834,
    0.504,
    0.004,
    -0.496,
    -0.966,
    -1.406,
    -1.813, -1.936, -2.058, -2.142,
    -2.289, -2.407, -2.551, -2.659,
    -2.847, -3.043, -3.217, -3.350,
    -3.558, -3.761, -3.964, -4.139,
    -4.360, -4.586, -4.817, -5.006,
    -5.248, -5.476, -5.700, -5.871,
    -6.111, -6.344, -6.573, -6.775,
    -7.021, -7.274, -7.521, -7.734,
    -7.997, -8.271, -8.526, -8.722,
    -8.985, -9.237, -9.503, -9.741,
    -10.045, -10.347, -10.638, -10.888,
    -11.189
])
ut1tai_dt = [datetime.datetime(1955 + i, 1, 1) for i in range(7)] \
            + [datetime.datetime(1962 + i // 4, (i % 4) * 3 + 1, 1)
               for i in range(len(ut1tai) - 7)]


def fake_taiutc(dt, leaplist):
    """TAI-UTC based on a list of leap second times"""
    return numpy.searchsorted(leaplist, dt)

# Calculate a list of leap seconds based on UT1
# Introduce leap second when leapsecondUTC-UT1 > 0.4
# NIST adds ls when UTC-UT1 > 0.4 to keep UTC-UT1 < 0.9
ut1_ls = []
for i, dt in enumerate(ut1tai_dt):
    if dt.month in (4, 10):
        continue # Can't do leapsecond
    if dt.year >= 1972: # 1972 is end of rubber second
        continue
    if -ut1tai[i] - len(ut1_ls) > 0.4:
        ut1_ls.append(dt)
# The non-rubber portions
if len(ut1_ls) < 10:
    ut1_ls.append(datetime.datetime(1972, 1, 1))
if len(ut1_ls) < 11:
    ut1_ls.append(datetime.datetime(1972, 7, 1))
if len(ut1_ls) < 12:
    ut1_ls.append(datetime.datetime(1973, 1, 1))
print('Leapseconds based on UTC - UT1 > 0.4s:')
print('\n'.join([str(t) for t in ut1_ls]))


def rubber_taiutc(dt):
    """Calculate TAI-UTC via rubber second for a UTC date dt"""
    idx = numpy.searchsorted(rubbers['dt'], dt, side='right') - 1
    if idx < 0:
        return 0.
    # assumes input time is start of day...
    mjd = (dt - datetime.datetime(1858, 11, 17)).days
    return rubbers[idx]['delta'] \
        + (mjd - rubbers[idx]['mjd'])  * rubbers[idx]['rate']


# Fake leapseconds, introduce leap second when leapsecondUTC - rubberUTC > 0.5
utc_ls = []
for dt in [datetime.datetime(1958 + i // 2, 1 + (i % 2) * 6, 1)
           for i in range(28)]:
    if rubber_taiutc(dt) - len(utc_ls) > 0.5:
        utc_ls.append(dt)
# The non-rubber portions
if len(utc_ls) < 10:
    utc_ls.append(datetime.datetime(1972, 1, 1))
if len(utc_ls) < 11:
    utc_ls.append(datetime.datetime(1972, 7, 1))
if len(utc_ls) < 12:
    utc_ls.append(datetime.datetime(1973, 1, 1))


# Fake leapseconds, current SpacePy
def spacepy_taiutc(dt):
    spacepy_ls = [
        datetime.datetime(int(y), int(m), 1)
        for y, m in zip(spacepy.time.year, spacepy.time.mon)]
    idx = numpy.searchsorted(spacepy_ls, dt, side='right') - 1
    if idx < 0:
        return 0
    return spacepy.time.secs[idx]

# And CDF
def cdf_taiutc(dt):
    # CDF tt2000 when UTC and TAI set equal
    tt2000_1958 = -1325419167816000000
    naive_seconds = (dt - datetime.datetime(1958, 1, 1)).total_seconds()
    actual_seconds = (spacepy.pycdf.lib.datetime_to_tt2000(dt) - tt2000_1958) \
                     / 1.e9
    return actual_seconds - naive_seconds


times = spacepy.time.tickrange('1955-1-1', '1973-1-1', 1).UTC


fig = plt.figure(figsize=(11, 8.5))
ax = fig.add_subplot(111)

ax.plot(times, [rubber_taiutc(t) for t in times],
        ls='-', marker=None, label='TAI-UTC')
ax.plot(ut1tai_dt, -ut1tai, ls='', marker='X', ms=5, label='TAI-UT1')
ax.plot(times, [fake_taiutc(t, ut1_ls) for t in times],
        ls='--', marker=None, label='Proposed (UT1)')
ax.plot(times, [fake_taiutc(t, utc_ls) for t in times],
        ls='-.', marker=None, label='Proposed (UTC)')
ax.plot(times, [spacepy_taiutc(t) for t in times],
        ls=':', marker=None, label='Current')
# Revert to <=0.2.2 leapsecond treatment
spacepy.time._read_leaps(oldstyle=True)
ax.plot(times, [spacepy_taiutc(t) for t in times],
        ls=':', marker=None, label='0.2.2')
ax.plot(times, [cdf_taiutc(t) for t in times],
        ls='--', marker=None, label='CDF TT2000')
ax.set_ylabel('seconds')
ax.set_xlabel('UTC date')
ax.legend(loc='best')
plt.savefig('rubber_seconds.pdf', dpi=200)
plt.savefig('rubber_seconds.png', dpi=200)
