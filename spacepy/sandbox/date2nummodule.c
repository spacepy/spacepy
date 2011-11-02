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
#include <numpy/arrayobject.h>

#define HOURS_PER_DAY 24
#define MINUTES_PER_DAY (60*HOURS_PER_DAY)
#define SECONDS_PER_DAY (60*MINUTES_PER_DAY)
#define MUSECONDS_PER_DAY (1e6*SECONDS_PER_DAY)
#define SEC_PER_MIN 60
#define SEC_PER_HOUR 3600
#define SEC_PER_DAY (SEC_PER_HOUR*24)
#define SEC_PER_WEEK (SEC_PER_DAY*7)


// pulled from python source
static npy_int is_leap(npy_int year)
{
    /* Cast year to unsigned.  The result is the same either way, but
     * C can generate faster code for unsigned mod than for signed
     * mod (especially for % 4 -- a good compiler should just grab
     * the last 2 bits when the LHS is unsigned).
     */
    const npy_uint ayear = (npy_uint)year;
    return ayear % 4 == 0 && (ayear % 100 != 0 || ayear % 400 == 0);
}

// pulled from python source
static npy_int _days_before_month[] = {
    0, /* unused; this vector uses 1-based indexing */
    0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334
};

// pulled from python source
/* year, month -> number of days in year preceeding first day of month */
static npy_int days_before_month(npy_int year, npy_int month)
{
    npy_int days;

    assert(month >= 1);
    assert(month <= 12);
    days = _days_before_month[month];
    if (month > 2 && is_leap(year))
        ++days;
    return days;
}

// pulled from python source
static npy_int days_before_year(npy_int year)
{
    npy_int y = year - 1;
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
static npy_int ymd_to_ord(npy_int year, npy_int month, npy_int day)
{
    return days_before_year(year) + days_before_month(year, month) + day;
}


/*==============================================================================
# date2num
# actually do the calculation, totoally stolen from matplotlib in how it works
#
# raise a ValueError for a non datetime input
==============================================================================*/
static npy_double date2num(PyDateTime_DateTime *inval) {
    npy_int microsecond, second, minute, hour, day, month, year;
    npy_double ord; // the ordinal
    // If inval is not a datetime object raise ValueError
    microsecond = PyDateTime_DATE_GET_MICROSECOND(inval);
    second = PyDateTime_DATE_GET_SECOND(inval);
    minute = PyDateTime_DATE_GET_MINUTE(inval);
    hour = PyDateTime_DATE_GET_HOUR(inval);
    day = PyDateTime_GET_DAY(inval);
    month = PyDateTime_GET_MONTH(inval);
    year = PyDateTime_GET_YEAR(inval);
    ord = (npy_double)ymd_to_ord(year, month, day); // this is from datetimemodule.c
    // make sure we dont have any int division    
    ord += ((npy_double)hour/HOURS_PER_DAY + 
            (npy_double)minute/MINUTES_PER_DAY + 
            (npy_double)second/SECONDS_PER_DAY + 
            (npy_double)microsecond/MUSECONDS_PER_DAY);
    return (ord);
}


/*==============================================================================
 This handles the type checking and interface to then pass the data on to date2num()
==============================================================================*/
static PyArrayObject *date2num_common(PyObject *self, PyObject *args) {
    PyDateTime_DateTime *inval;
    Py_ssize_t inval_len;
    Py_ssize_t ind;
    PyDateTime_DateTime *item;
    PyObject *outval;

    if (!PyArg_ParseTuple(args, "O", &inval)) 
        return NULL;

    // if inval is a datetime then call the conversion and return
    if (PyDate_Check((PyObject *)inval)) {
        // these were uninitialized, do I still need to DECREF?        
        //Py_DECREF(item);  // clean up the objects
        //Py_DECREF(outval);
        return (PyArrayObject *)PyFloat_FromDouble(date2num(inval));  
    } else {
    // if it is not a datetime then iterate over it doing the conversion for each
        if (!PySequence_Check((PyObject *)inval)) { // is the input a sequace of sorts
            PyErr_SetString(PyExc_ValueError, "Must be a datetime object or iterable of datetime objects");
            // these were uninitialized, do I still need to DECREF?        
            //Py_DECREF(item); // clean up the objects
            //Py_DECREF(outval);
            return NULL;         
        }
        // this is an iterable, do something with it
        inval_len = PySequence_Length((PyObject *)inval);
        // check the zeroth element just to be sure it is a datetime before making outarray, could be huge
        item = (PyDateTime_DateTime*)PySequence_GetItem((PyObject *)inval, 0);
        if (!PyDate_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "Iterable must contain datetime objects");
            Py_DECREF(item); // clean up the objects
            // this was uninitialized here does it still need decref?
            //Py_DECREF(outval);
            return NULL; 
        }
        Py_DECREF(item);

/****************
Have tried nearly every combination of PyArray_SimpleNew and PyArray_SETITEM so giving up
Moving to the less good solution of list then change to array
preallocating the list the populating it then change to array
***************/
        outval = PyList_New(inval_len);  // make a full sized list
        // step thru all the datetimes and  convert them, putting ans in the list
        for (ind=0;ind<inval_len;ind++) {
            item = (PyDateTime_DateTime*)PySequence_GetItem((PyObject *)inval, ind);
            // If this isn't a datetime error
            if (!PyDate_Check(item)) {
                PyErr_SetString(PyExc_ValueError, "Iterable must contain datetime objects");
                Py_DECREF(item);
                return NULL; 
            }

            PyList_SET_ITEM(outval, ind, (PyObject *)PyFloat_FromDouble(date2num(item)));
            Py_DECREF(item);
        }
    }
    /*Giving away our reference to the caller*/
    return (PyArrayObject *)PyArray_Return((PyArrayObject*)PyArray_FROM_O(outval));
}

static PyMethodDef date2num_methods[] = {
   { "date2num", (PyCFunction)date2num_common, METH_VARARGS,
     "Return value is a floating point number (or sequence of floats) which \n"
     "gives the number of days (fraction part represents hours, minutes, seconds)\n"
     "since 0001-01-01 00:00:00 UTC, plus one. The addition of one here is a \n"
     "historical artifact. Also, note that the Gregorian calendar is assumed; \n"
     "this is not universal practice. For details, see the module docstring.\n"},
   { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initdate2num(void) {
   Py_InitModule3("date2num", date2num_methods, 
                    "Module for computing matplotlib.dates.date2num a lot faster");
   PyDateTime_IMPORT;
   import_array();
}

