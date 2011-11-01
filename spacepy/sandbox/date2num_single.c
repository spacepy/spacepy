/*************************************************************************/
/*
Test compilation:
gcc -DNDEBUG -g -O3 -Wall -Wstrict-prototypes -fPIC -DMAJOR_VERSION=1 -DMINOR_VERISON=0 -I /usr/include/python2.6  -c ctrace2dmodule.c
gcc -shared ctrace2dmodule.o -o ctrace2d.so

Copyright 2010 - 2011  Los Alamos National Security, LLC. */

/*************************************************************************/
#include <Python.h>
//#include <numpy/arrayobject.h>
#include <datetime.h>


#define HOURS_PER_DAY 24
#define MINUTES_PER_DAY 60*HOURS_PER_DAY
#define SECONDS_PER_DAY 60*MINUTES_PER_DAY
#define MUSECONDS_PER_DAY 1e6*SECONDS_PER_DAY
#define SEC_PER_MIN 60
#define SEC_PER_HOUR 3600
#define SEC_PER_DAY SEC_PER_HOUR*24
#define SEC_PER_WEEK SEC_PER_DAY*7

// from python datetimemodule.c
//static PyObject *
//date_toordinal(PyDateTime_Date *self)
//{
//    return PyInt_FromLong(ymd_to_ord(GET_YEAR(self), GET_MONTH(self),
//                                     GET_DAY(self)));
//}

// from python datetimemodule.c
/* year, month, day -> ordinal, considering 01-Jan-0001 as day 1. */
//static int
//ymd_to_ord(int year, int month, int day)
//{
//    return days_before_year(year) + days_before_month(year, month) + day;
//}

static double date2num(PyDateTime_DateTime *inval) {
    int microsecond, second, minute, hour, day, month, year;
    double ord; // the ordinal
    // make sure that the input is of the right type
    if !(PyDate_Check(inval))
        return -999.; // make this a proper error

    microsecond = PyDateTime_DATE_GET_MICROSECOND(inval);
    second = PyDateTime_DATE_GET_SECOND(inval);
    minute = PyDateTime_DATE_GET_MINUTE(inval);
    hour = PyDateTime_DATE_GET_HOUR(inval);
    day = PyDateTime_GET_DAY(inval);
    month = PyDateTime_GET_MONTH(inval);
    year = PyDateTime_GET_YEAR(inval);
    ord = (double)ymd_to_ord(year, month, day); // this is from datetimemodule.c
    // make sure we dont have any int division    
    ord += (double)hour/HOURS_PER_DAY + 
            (double)minute/MINUTES_PER_DAY + 
            (double)second/SECONDS_PER_DAY + 
            (double)microsecond/MUSECONDS_PER_DAY;  
    return (ord);
}


static PyObject *date2num_common(PyObject *self,
				 PyObject *args, PyObject *kwargs) {
//TODO and the macro PyDateTime_IMPORT must be invoked, usually as part of the 
//module initialisation function. The macro puts a pointer to a C structure into 
//a static variable, PyDateTimeAPI
  PyDateTime_IMPORT; // TODO does it matter where thi goes?
  PyDateTime_DateTime *inval;
  double outval;
  static char *kwlist[] = {NULL}; // TODO there are no kwargs so I don't need this right?
    
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!", kwlist, 
                                    &PyDateTime_DateTime, &inval)) 
    return NULL;

//NPY_BEGIN_ALLOW_THREADS
  outval = date2num(inval);
  
  
//NPY_END_ALLOW_THREADS

//  Py_DECREF(fieldy); // TODO what do I need to dec/inc in here?

//  outdims[0] = count;
//  if (!PyArray_Resize(outx, &outshape, 1, NPY_CORDER))
//    return NULL;
//  if (!PyArray_Resize(outy, &outshape, 1, NPY_CORDER))
//    return NULL;
//  /*Giving away our reference to the caller*/
  return Py_BuildValue("d", outval);
}

static PyMethodDef date2num_methods[] = {
   { "date2num", (PyCFunction)date2num, METH_VARARGS | METH_KEYWORDS,
     "PUT DOCS HERE\n"},
   { NULL, NULL, 0, NULL }
};
