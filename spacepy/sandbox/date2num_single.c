/*************************************************************************
# date2num
# take a single datetime object and compute the matplotlib.dates.date2num value
# but do it using an extension module so hat it is a heck of a lot faster
# 
# Brian Larsen & Jon Niehof
# balarsen@lanl.gov
# Copyright ©2010 - 2011 Los Alamos National Security, LLC.
*************************************************************************/
#include <Python.h>
#include <datetime.h>
#include <stdio.h>

#define HOURS_PER_DAY 24
#define MINUTES_PER_DAY (60*HOURS_PER_DAY)
#define SECONDS_PER_DAY (60*MINUTES_PER_DAY)
#define MUSECONDS_PER_DAY (1e6*SECONDS_PER_DAY)
#define SEC_PER_MIN 60
#define SEC_PER_HOUR 3600
#define SEC_PER_DAY (SEC_PER_HOUR*24)
#define SEC_PER_WEEK (SEC_PER_DAY*7)


// pulled from python source
static int is_leap(int year)
{
    /* Cast year to unsigned.  The result is the same either way, but
     * C can generate faster code for unsigned mod than for signed
     * mod (especially for % 4 -- a good compiler should just grab
     * the last 2 bits when the LHS is unsigned).
     */
    const unsigned int ayear = (unsigned int)year;
    return ayear % 4 == 0 && (ayear % 100 != 0 || ayear % 400 == 0);
}

// pulled from python source
static int _days_before_month[] = {
    0, /* unused; this vector uses 1-based indexing */
    0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
};

// pulled from python source
/* year, month -> number of days in year preceeding first day of month */
static int
days_before_month(int year, int month)
{
    int days;

    assert(month >= 1);
    assert(month <= 12);
    days = _days_before_month[month];
    if (month > 2 && is_leap(year))
        ++days;
    return days;
}

// pulled from python source
static int days_before_year(int year)
{
    int y = year - 1;
    /* This is incorrect if year <= 0; we really want the floor
     * here.  But so long as MINYEAR is 1, the smallest year this
     * can see is 0 (this can happen in some normalization endcases),
     * so we'll just special-case that.
     */
    assert (year >= 0);
    if (y >= 0)
        return y*365 + y/4 - y/100 + y/400;
    else {
        assert(y == -1);
        return -366;
    }
}

// pulled from python source
/* year, month, day -> ordinal, considering 01-Jan-0001 as day 1. */
static int ymd_to_ord(int year, int month, int day)
{
    return days_before_year(year) + days_before_month(year, month) + day;
}


/*==============================================================================
# date2num
# actually do the calculation, totoally stolen from matplotlib in how it works
#
# raise a ValueError for a non datetime input
==============================================================================*/
static double date2num(PyDateTime_DateTime *inval) {
    int microsecond, second, minute, hour, day, month, year;
    double ord; // the ordinal
    microsecond = PyDateTime_DATE_GET_MICROSECOND(inval);
    second = PyDateTime_DATE_GET_SECOND(inval);
    minute = PyDateTime_DATE_GET_MINUTE(inval);
    hour = PyDateTime_DATE_GET_HOUR(inval);
    day = PyDateTime_GET_DAY(inval);
    month = PyDateTime_GET_MONTH(inval);
    year = PyDateTime_GET_YEAR(inval);
    ord = (double)ymd_to_ord(year, month, day); // this is from datetimemodule.c
    // make sure we dont have any int division    
    ord += ((double)hour/HOURS_PER_DAY + 
            (double)minute/MINUTES_PER_DAY + 
            (double)second/SECONDS_PER_DAY + 
            (double)microsecond/MUSECONDS_PER_DAY);
    return (ord);
}


/*==============================================================================
 This handles the type checking and interface to then pass the data on to date2num()
==============================================================================*/
static PyObject *date2num_common(PyObject *self, PyObject *args) {
  PyDateTime_DateTime *inval;
  double outval;
  
  if (!PyArg_ParseTuple(args, "O", &inval)) 
    return NULL;
  // If it is not a datetime object raise ValueError
  if (!PyDate_Check(inval)) {
    PyErr_SetString(PyExc_ValueError, "Must be a datetime object");
    return NULL; 
  }
  outval = date2num(inval);

  /*Giving away our reference to the caller*/
  return Py_BuildValue("d", outval);
}

static PyMethodDef date2num_methods[] = {
   { "date2num_single", (PyCFunction)date2num_common, METH_VARARGS,
     "Return value is a floating point number (or sequence of floats) which \n"
     "gives the number of days (fraction part represents hours, minutes, seconds)\n"
     "since 0001-01-01 00:00:00 UTC, plus one. The addition of one here is a \n"
     "historical artifact. Also, note that the Gregorian calendar is assumed; \n"
     "this is not universal practice. For details, see the module docstring.\n"},
   { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initdate2num_single(void) {
   Py_InitModule3("date2num_single", date2num_methods, 
                    "Module for computing matplotlib.dates.date2num a lot faster");
   PyDateTime_IMPORT;
}

