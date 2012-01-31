/*************************************************************************
# date2num
# take a single datetime object and compute the matplotlib.dates.date2num value
# but do it using an extension module so hat it is a heck of a lot faster
#
# Brian Larsen & Jon Niehof
# balarsen@lanl.gov
# Copyright 2010 - 2011 Los Alamos National Security, LLC.
*************************************************************************/
// the normal include for python extension modules
#include <Python.h>
// to use datetimes, it is separate from Python.h
#include <datetime.h>
// to use numpy arrays, there are many things in numpy
#include <numpy/arrayobject.h>

#define PyNumber_AsDouble(y) PyFloat_AsDouble(PyNumber_Float(y))


#define HOURS_PER_DAY 24
#define MINUTES_PER_DAY (60*HOURS_PER_DAY)
#define SECONDS_PER_DAY (60*MINUTES_PER_DAY)
#define MUSECONDS_PER_DAY (1e6*SECONDS_PER_DAY)
#define SEC_PER_MIN 60
#define SEC_PER_HOUR 3600
#define SEC_PER_DAY (SEC_PER_HOUR*24)
#define SEC_PER_WEEK (SEC_PER_DAY*7)
/* Number of days in 4, 100, and 400 year cycles.  That these have
 * the correct values is asserted in the module init function.
 */
#define DI4Y	1461	/* days_before_year(5); days in 4 years */
#define DI100Y	36524	/* days_before_year(101); days in 100 years */
#define DI400Y	146097	/* days_before_year(401); days in 400 years  */




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

/* For each month ordinal in 1..12, the number of days in that month,
 * and the number of days before that month in the same year.  These
 * are correct for non-leap years only.
 */
static int _days_in_month[] = {
	0, /* unused; this vector uses 1-based indexing */
	31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31
};

// pulled from python source
/* year, month -> number of days in that month in that year */
static int
days_in_month(int year, int month)
{
	assert(month >= 1);
	assert(month <= 12);
	if (month == 2 && is_leap(year))
		return 29;
	else
		return _days_in_month[month];
}

// pulled from python source
/* ordinal -> year, month, day, considering 01-Jan-0001 as day 1. */
static void
ord_to_ymd(npy_int ordinal, npy_int *year, npy_int *month, npy_int *day)
{
	npy_int n, n1, n4, n100, n400, leapyear, preceding;

	/* ordinal is a 1-based index, starting at 1-Jan-1.  The pattern of
	 * leap years repeats exactly every 400 years.  The basic strategy is
	 * to find the closest 400-year boundary at or before ordinal, then
	 * work with the offset from that boundary to ordinal.  Life is much
	 * clearer if we subtract 1 from ordinal first -- then the values
	 * of ordinal at 400-year boundaries are exactly those divisible
	 * by DI400Y:
	 *
	 *    D  M   Y            n              n-1
	 *    -- --- ----        ----------     ----------------
	 *    31 Dec -400        -DI400Y       -DI400Y -1
	 *     1 Jan -399         -DI400Y +1   -DI400Y      400-year boundary
	 *    ...
	 *    30 Dec  000        -1             -2
	 *    31 Dec  000         0             -1
	 *     1 Jan  001         1              0          400-year boundary
	 *     2 Jan  001         2              1
	 *     3 Jan  001         3              2
	 *    ...
	 *    31 Dec  400         DI400Y        DI400Y -1
	 *     1 Jan  401         DI400Y +1     DI400Y      400-year boundary
	 */
	assert(ordinal >= 1);
	--ordinal;
	n400 = ordinal / DI400Y;
	n = ordinal % DI400Y;
	*year = n400 * 400 + 1;

	/* Now n is the (non-negative) offset, in days, from January 1 of
	 * year, to the desired date.  Now compute how many 100-year cycles
	 * precede n.
	 * Note that it's possible for n100 to equal 4!  In that case 4 full
	 * 100-year cycles precede the desired day, which implies the
	 * desired day is December 31 at the end of a 400-year cycle.
	 */
	n100 = n / DI100Y;
	n = n % DI100Y;

	/* Now compute how many 4-year cycles precede it. */
	n4 = n / DI4Y;
	n = n % DI4Y;

	/* And now how many single years.  Again n1 can be 4, and again
	 * meaning that the desired day is December 31 at the end of the
	 * 4-year cycle.
	 */
	n1 = n / 365;
	n = n % 365;

	*year += n100 * 100 + n4 * 4 + n1;
	if (n1 == 4 || n100 == 4) {
		assert(n == 0);
		*year -= 1;
		*month = 12;
		*day = 31;
		return;
	}

	/* Now the year is correct, and n is the offset from January 1.  We
	 * find the month via an estimate that's either exact or one too
	 * large.
	 */
	leapyear = n1 == 3 && (n4 != 24 || n100 == 3);
	assert(leapyear == is_leap(*year));
	*month = (n + 50) >> 5;
	preceding = (_days_before_month[*month] + (*month > 2 && leapyear));
	if (preceding > n) {
		/* estimate is too large */
		*month -= 1;
		preceding -= days_in_month(*year, *month);
	}
	n -= preceding;
	assert(0 <= n);
	assert(n < days_in_month(*year, *month));

	*day = n + 1;
}



// pulled from python source
/* year, month -> number of days in year preceding first day of month */
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
# actually do the calculation, totally stolen from matplotlib in how it works
#
# raise a ValueError for a non datetime input
==============================================================================*/
// this is the c code that does all the work, does NOT have to be separate, just made sense to me
// PyDateTime_DateTime the C name for a datetime type
static npy_double date2num(PyDateTime_DateTime *inval) {
    // using the numpy types for variables in a platform independence try
    npy_int microsecond, second, minute, hour, day, month, year;
    npy_double ord; // the ordinal
    // If inval is not a datetime object raise ValueError
    // these function are all in datetime.h
    microsecond = PyDateTime_DATE_GET_MICROSECOND(inval);
    second = PyDateTime_DATE_GET_SECOND(inval);
    minute = PyDateTime_DATE_GET_MINUTE(inval);
    hour = PyDateTime_DATE_GET_HOUR(inval);
    day = PyDateTime_GET_DAY(inval);
    month = PyDateTime_GET_MONTH(inval);
    year = PyDateTime_GET_YEAR(inval);
    // this is in datetimemodule.c in python, reproduced here, working to change this
    ord = (npy_double)ymd_to_ord(year, month, day); // this is from datetimemodule.c

//    ord = (npy_double)ymd_to_ord(year, month, day); // this is from datetimemodule.c
    // make sure we don't have any int division
    ord += ((npy_double)hour/HOURS_PER_DAY +
            (npy_double)minute/MINUTES_PER_DAY +
            (npy_double)second/SECONDS_PER_DAY +
            (npy_double)microsecond/MUSECONDS_PER_DAY);
    return (ord);
}


/*==============================================================================
 This handles the type checking and interface to then pass the data on to date2num()

 excruciating documentation so that others can use this as an example
==============================================================================*/
// static: All of these are static so that they don't pollute the python name space
// name: doesn't matter what you call the function as the name is handled below
// self, args: this is the standard way to use a function that takes no kwargs
// return type: just whatever type you are going to return
static PyArrayObject *date2num_common(PyObject *self, PyObject *args) {
    // PyDateTime_DateTime the C name for a datetime type
    PyObject *inval;
    // using the python variable types in a platform independence try
    Py_ssize_t inval_len;
    Py_ssize_t ind;
    PyDateTime_DateTime *item;
    // the PyArrayObject that we will return is declared as a PyObject for ease later
    //    then casted to a PyObject* before return
    PyObject *outval;
    // this is the easy way to access the data of the output array, grab a double*
    //    to the data of the numpy array that we are creating
    npy_double *outval_dat;

    // This is required in each extension module, the "O" is an abject and tells
    //    python that there can only be one argument of type object
    //    the ! that is common can't be used on datetimes (no auto type checking)
    if (!PyArg_ParseTuple(args, "O", &inval))
        // this is what you do to let the exception actually be raised, PyArg_ParseTuple
        //     sets an exception flag then the return NULL lets it work and be raised
        return NULL;

    // if inval is a datetime then call the conversion and return
    if (PyDate_Check(inval))
        // if we got a datetime just do he calculation an return the answer
        // PyArrayObject*: this is casting it to the expected return type
        // PyFloat_FromDouble: changes a C double to Python float
        return (PyArrayObject *)PyFloat_FromDouble(date2num((PyDateTime_DateTime*)inval));
    // if it is not a datetime then iterate over it doing the conversion for each
    // PySequence_Check is magic, it tests if the input is a sequence (iterable)
    // PyObject *: cast inval to the type expected by PySequence_Check
    if (!PySequence_Check(inval)) { // is the input a sequence of sorts
        // PyErr_SetString: easy way to raise an exception
        PyErr_SetString(PyExc_ValueError, "Must be a datetime object or iterable of datetime objects");
        // return NULL: again this is how to get the exception to actually work
        //    be sure to do any required DECREF before return, none are required
        //    here as we have not initialized any python objects
        return NULL;
    }
    // this is an iterable, do something with it
    // PySequence_Length: how long s the iterable?
    // PyObject *: cast inval to the type expected by PySequence_Length
    inval_len = PySequence_Length(inval);
    // check the zeroth element just to be sure it is a datetime before making outarray, could be huge
    // PySequence_GetItem: gets an object from the sequence at the requested element
    // PyDateTime_DateTime*: cast the object to the type we want
    // PyObject *: cast inval to the type expected by PySequence_GetItem
    item = (PyDateTime_DateTime*)PySequence_GetItem(inval, 0);
    // same check as above, is it a datetime object?
    if (!PyDate_Check(item)) {
        // again raise an exception, really this sets the exception counter
        //     and the return is what lets python actually do the raising
        PyErr_SetString(PyExc_ValueError, "Iterable must contain datetime objects");
        // Py_DECREF: since we created a python object "item" we have to
        //    Py_DECREF it so that the gc can deal with it before we go back
        Py_DECREF(item); // clean up the objects
        //go back to python and let the exception be raised
        return NULL;
    }
    // Py_DECREF: since we created a python object "item" we have to Py_DECREF
    //    it since we don't use it
    //    "Your mother does not work here, clean up after yourself"
    Py_DECREF(item);

    // PyArray_ZEROS: make a new numpy array like numpy.zeros()
    // NPY_DOUBLE: the data type of the array
    // 1: the dimension
    // &inval_len: pointer to the number of elements in each dimension
    // 0: 0-C array order, 1-FORTRAN array order (use C)
    outval = PyArray_ZEROS(1, &inval_len, NPY_DOUBLE, 0);
    // PyArray_DATA: return a pointer to the data of the array we created
    // double*: cast it to the type we expect and want to use
    outval_dat = (double*)PyArray_DATA(outval);

    // step thru all the datetimes and  convert them, putting ans in the array
    for (ind=0;ind<inval_len;ind++) {
        // as above get an item from the iterator ival and cast as needed
        item = (PyDateTime_DateTime*)PySequence_GetItem(inval, ind);
        // If this isn't a datetime error
        if (!PyDate_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "Iterable must contain datetime objects");
            // don't forget to clean up python objects that you create
            Py_DECREF(item);
            return NULL;
        }
        // since outval_dat is a double* we can just use it as we want
        outval_dat[ind] = date2num((PyDateTime_DateTime*)item);
        // we are done with item and ready to get the next one, clean up
        Py_DECREF(item);
    }
    /*Giving away our reference to the caller*/
    // cast the outval of type PyObject* to what we want to return
    return (PyArrayObject *)outval;
}

static PyDateTime_DateTime *num2date(double inval){
    npy_int year, month, day, hour, minute, second, usecond;
    npy_double remainder;

    remainder = inval - (npy_int)inval; // get the fraction part
    remainder *= 24;
    hour = (npy_int)(remainder);
    remainder -= hour;
    remainder *= 60;
    minute = (npy_int)(remainder);
    remainder -= minute;
    remainder *= 60;
    second = (npy_int)(remainder);
    remainder -= second;
    remainder *= 1e6;
    usecond = (npy_int)(remainder);

    if (usecond < 10)
        usecond = 0;

    if (usecond>999990)  // compensate for rounding errors
        usecond += (1e6-usecond);

    ord_to_ymd((npy_int)inval, &year, &month, &day);
    return (PyDateTime_DateTime*)PyDateTime_FromDateAndTime(year, month, day, hour, minute, second, usecond );
}


static PyArrayObject *num2date_common(PyObject *self, PyObject *args) {
    PyObject *inval; // what the user passed in
    Py_ssize_t inval_len; // the length of what was passed in
    Py_ssize_t ind;  // for loop index
    PyObject *item;  // item from the user input list
    PyObject *item_out;  // converted item
    PyObject *outval;   // array to return
    PyObject **outval_dat;

    if (!PyArg_ParseTuple(args, "O", &inval))
        return NULL;

    // the input is a number or an array
    if (PyNumber_Check(inval) & !PySequence_Check(inval))
        return (PyArrayObject *)(num2date(PyNumber_AsDouble(inval))); 
    if (!PySequence_Check(inval)) { // is the input a sequence of sorts
        PyErr_SetString(PyExc_ValueError, "Must be a numeric object or iterable of numeric objects");
        return NULL;
    }
    inval_len = PySequence_Length(inval);

    if (!(item_out = PySequence_GetItem(inval, 0)))
	  return NULL;
    // same check as above, is it a datetime object?
    if (!PyNumber_Check(item_out)) {
        PyErr_SetString(PyExc_ValueError, "Iterable must contain numeric objects");
        Py_DECREF(item_out); // clean up the objects
        return NULL;
    }
    //    "Your mother does not work here, clean up after yourself"
    Py_DECREF(item_out);

    outval = PyArray_SimpleNew(1, &inval_len, NPY_OBJECT);
    outval_dat = (PyObject **)PyArray_DATA(outval);
    
    // step thru all the datetimes and  convert them, putting ans in the array
    for (ind=0;ind<inval_len;ind++) {
        // as above get an item from the iterator ival and cast as needed
        if (!(item = PySequence_GetItem(inval, ind))) {
            return NULL;
        }
        // If this isn't a float error
        if (!PyNumber_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "Iterable must contain numeric objects");
            Py_DECREF(item);
            return NULL;
        }
        item_out = (PyObject *)(num2date(PyNumber_AsDouble(item)));

        outval_dat[ind] = item_out;
        Py_DECREF(item);
    }
    /*Giving away our reference to the caller*/
    // cast the outval of type PyObject* to what we want to return
    return (PyArrayObject *)PyArray_FROM_O(outval);
}



// setup h C that the Python knows what to do with it
// PyMethodDef: struct type used to hold int he info, one set per function you
//     want exposed in python  i.e. this is where you set the name
// static: again we don't want Python to have this in the namespace
static PyMethodDef dates_methods[] = {
    // date2num: this is the name as it will show up in python
    // (PyCFunction)date2num_common: cast the c function we created and associate
    //    it in python with the name
    //  METH_VARARGS: it takes in arguments and not kwargs
   { "_date2num", (PyCFunction)date2num_common, METH_VARARGS,
     "Return value is a floating point number (or sequence of floats) which \n"
     "gives the number of days (fraction part represents hours, minutes, seconds)\n"
     "since 0001-01-01 00:00:00 UTC, plus one. The addition of one here is a \n"
     "historical artifact. Also, note that the Gregorian calendar is assumed; \n"
     "this is not universal practice. For details, see the module docstring.\n"},
     {"_num2date", (PyCFunction)num2date_common, METH_VARARGS,
     "Input is a float value which gives the number of days (fraction part \n"
     "represents hours, minutes, seconds) since 0001-01-01 00:00:00 UTC plus one.\n"
     "The addition of one here is a historical artifact. Also, note that the\n"
     "Gregorian calendar is assumed; this is not universal practice. For details,\n"
     "see the module docstring.\n\n"
     "Return value is a datetime instance in timezone tz (default to rcparams TZ value).\n\n"
     "If input is a sequence, a numpy array of datetime objects will be returned.\n"},
    // NULL terminate Python looking at the object
     { NULL, NULL, 0, NULL }
};

// this the what does the exposing to python
// initdate2num: it is normal to name this with a leading "init" and to name the
//    file with a trailing "module"
PyMODINIT_FUNC init_dates(void) {
    // this sets up the module, there are different versions of Py_InitModule3
    //     that take different arguments based on docstring etc
    // date2num: this is the name of the module
    // date2num_methods: the name of the method dict we defined above
    // docstring: in the "3" version you spec the docstring here
    Py_InitModule3("_dates", dates_methods,
                     "Module for computing matplotlib.dates.date2num and num2date a lot faster");
    // this is a required macro for using datetime objects in the module, it goes here
    PyDateTime_IMPORT;
    // this is a required function (macro) for using numpy arrays in the module
    import_array();
}

