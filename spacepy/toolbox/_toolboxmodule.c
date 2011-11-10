/*************************************************************************
# _toolboxmodule
# make some functions in toolbox c extensions to be a lot faster
#
# Brian Larsen
# balarsen@lanl.gov
# Copyright ©2010 - 2011 Los Alamos National Security, LLC.
*************************************************************************/
// the normal include for python extension modules
#include <Python.h>
#include <math.h>
#include <numpy/arrayobject.h>
#include <datetime.h>

#define PyNumber_AsDouble(y) PyFloat_AsDouble(PyNumber_Float(y))
#define TRUE 1
#define FALSE 0

/*Function to call a method, given an object, method name, arguments*/
// Jon Niehof
static PyObject *callmeth(PyObject *obj, const char* methname,
			  PyObject *args, PyObject *kwargs)
{
    PyObject *meth, *retval;

    if (!(meth = PyObject_GetAttrString(obj, methname)))
        return NULL;
    retval = PyEval_CallObjectWithKeywords(meth, args, kwargs);
    Py_DECREF(meth);
    return retval;
}

//=============================================================================
// Make linearly spaced data from start to finish, lso handle datetime objects
//=============================================================================
static PyObject *linspace_tb(PyObject *self, PyObject *args, PyObject *kwargs)
{
//    Grabbed a 15x speedup over the numpy version
//    tb_t = timeit.Timer("tb.linspace(1.5, 10.5, 6)", setup="import spacepy.toolbox as tb")
//    np_t = timeit.Timer("np.linspace(1.5, 10.5, 6)", setup="import numpy as np")
//    print np_t.timeit(10000)/tb_t.timeit(10000)
    double startVal, stopVal, step, num_in;
    PyObject *startDT, *stopDT;
    double *outval_dat;
    Py_ssize_t num, ii;
    short endpoint=TRUE, retstep=FALSE, forcedate=FALSE;
    PyObject *outval, *outtuple, *calltuple;
    PyObject *datetime_module;

    static char *kwlist[] = {"startVal", "stopVal", "num", "endpoint", "retstep", "forcedate", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddd|hhh", kwlist, &startVal, &stopVal,
        &num_in, &endpoint, &retstep, &forcedate)) {
        PyErr_Clear();
        // this was not doubles, test if it is datetimes
        if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOd|hhh", kwlist, &startDT, &stopDT,
            &num_in, &endpoint, &retstep, &forcedate))
            return NULL;
        else {
            if (!PyDateTime_Check(startDT) | !PyDateTime_Check(stopDT)) {
                PyErr_SetString(PyExc_TypeError, "number or datetime objects expected");
                return NULL;
            }
        // need to call the spacepy.time.date2num
        if (!(datetime_module = PyImport_ImportModule("spacepy.time")))
            return NULL;
        calltuple = PyTuple_Pack(1, startDT);
        startVal = PyFloat_AsDouble(callmeth(datetime_module, "date2num", calltuple, NULL));
        Py_DECREF(calltuple);
        calltuple = PyTuple_Pack(1, stopDT);
        stopVal  = PyFloat_AsDouble(callmeth(datetime_module, "date2num", calltuple, NULL));
        Py_DECREF(calltuple);
        forcedate = TRUE;
        }
    }
    // TODO: Right here we ae totally done with Python objects from the outside, if they
    // change we don't care so we can release the GIL

    num = (Py_ssize_t)num_in;
    // if we want 0 (or neg) length return empty array
    if (num<=0) {
        num=0;
        outval = PyArray_SimpleNew(1, &num, NPY_DOUBLE);
        Py_XDECREF(datetime_module); // X since we may not have used datetime_module
        return outval;
    }
    // no neg to PyArray_SimpleNew()
    outval = PyArray_SimpleNew(1, &num, NPY_DOUBLE);
    outval_dat = (double *)PyArray_DATA(outval);
    if (endpoint) {
        if (num==1) {
            outval_dat[0] = startVal;
            if (!forcedate)
                return outval;
            calltuple = PyTuple_Pack(1, PyFloat_FromDouble(startVal));
            outval  = callmeth(datetime_module, "num2date", calltuple, NULL);
            Py_DECREF(calltuple);
            Py_DECREF(datetime_module);
            return outval;
        }
        step = (stopVal-startVal)/((double)(num-1));
        for (ii=0; ii<num; ii++)
            outval_dat[ii] = ii*step + startVal;
        outval_dat[num-1] = stopVal;
    } else {
        step = (stopVal-startVal)/((double)num);
        for (ii=0; ii<num; ii++)
            outval_dat[ii] = ii*step + startVal;
    }
    if (retstep) {
        if (!forcedate)
            outtuple = PyTuple_Pack(2, outval, PyFloat_FromDouble(step));
        else {

            calltuple = PyTuple_Pack(1, outval);
            outval  = callmeth(datetime_module, "num2date", calltuple, NULL);
            Py_DECREF(calltuple);
            outtuple = PyTuple_Pack(2, outval, PyFloat_FromDouble(step));
            Py_DECREF(datetime_module);
        }
        return outtuple;
    }
    if (!forcedate)
        return outval;
    calltuple = PyTuple_Pack(1, outval);
    outval  = callmeth(datetime_module, "num2date", calltuple, NULL);
    Py_DECREF(calltuple);
    Py_DECREF(datetime_module);
    return outval;
}


//=============================================================================
// Go through the input arguments or iterator and return sqrt(x**2 + y**2 + ...)
//=============================================================================
static PyObject *hypot_tb(PyObject *self, PyObject *args)
{
    Py_ssize_t TupleSize = PyTuple_Size(args);
    Py_ssize_t i, seqSize;
    PyObject *temp_p, *temp_seq;
    double tot=0., tmp_d;

    if(!TupleSize) {
        if(!PyErr_Occurred())
            PyErr_SetString(PyExc_TypeError,"hypot expected at least 1 argument, got 0");
        return NULL;
    }
    // if TupleSize is one maybe we got a sequence in
    if (TupleSize == 1) {
    // is it a sequence
        temp_p = PyTuple_GetItem(args,0);
        if(temp_p == NULL)
            return NULL;
        if (PySequence_Check(temp_p)) { // is the input a sequence of sorts
            seqSize =  PySequence_Length(temp_p);
            for (i=0;i<seqSize;i++) {
                temp_seq = PySequence_GetItem(temp_p, i);
                tmp_d = PyNumber_AsDouble(temp_seq);
                tot += (tmp_d*tmp_d);
            }
            return PyFloat_FromDouble(sqrt(tot));
        }
    }

    // If we get here it is a bunch of args
    for (i=0;i<TupleSize;i++) {
        temp_p = PyTuple_GetItem(args,i); // borrowed reference
        if(temp_p == NULL)
            return NULL;

        /* Check if temp_p is numeric */
        if (!(PyNumber_Check(temp_p))) {
            PyErr_SetString(PyExc_TypeError,"Non-numeric argument.");
            return NULL;
        }
        tmp_d = PyNumber_AsDouble(temp_p);
        tot += (tmp_d*tmp_d);
    }
        return PyFloat_FromDouble(sqrt(tot));
}

static PyMethodDef toolbox_methods[] = {
   { "hypot", (PyCFunction)hypot_tb, METH_VARARGS,
    "Compute sqrt(vals[0] **2 + vals[1] **2 ...), ie. n-dimensional hypotenuse\n\n"
    "If the input is a numpy array a c-backend is called for the calculation\n\n"
    "Parameters\n"
    "==========\n"
    "vals : float (arbitary number), or iterable\n"
    "\tarbitary number of float values as arguments or an iterable\n\n"
    "Returns\n"
    "=======\n"
    "out : float\n"
    "\tthe Euclidian distance of the points ot the origin\n\n"
    "Examples\n"
    "========\n"
    ">>> import spacepy.toolbox as tb\n"
    ">>> tb.hypot(3,4)\n"
    "5.0\n"
    ">>> a = [3, 4]\n"
    ">>> tb.hypot(*a)\n"
    "5.0\n"
    ">>> tb.hypot(*range(10))\n"
    "16.88194...\n\n"
    "See Also\n"
    "========\n"
    "math.hypot\n"},
   {"linspace", (PyCFunction)linspace_tb, METH_VARARGS|METH_KEYWORDS,
    "Returns linearly spaced numbers.  Same as numpy.linspace except\n"
    "allows for support of datetime objects\n\n"
    "Parameters\n"
    "==========\n"
    "start : float\n"
    "\tThe starting value of the sequence.\n"
    "stop : float\n"
    "\tThe end value of the sequence, unless `endpoint` is set to False.\n"
    "\tIn that case, the sequence consists of all but the last of ``num + 1``\n"
    "\tevenly spaced samples, so that `stop` is excluded.  Note that the step\n"
    "\tsize changes when `endpoint` is False.\n"
    "num : int (optional)\n"
    "\tNumber of samples to generate. Default is 50.\n"
    "endpoint : bool, optional\n"
    "\tIf True, `stop` is the last sample. Otherwise, it is not included.\n"
    "\tDefault is True.\n"
    "retstep : bool (optional)\n"
    "\tIf True, return (`samples`, `step`), where `step` is the spacing\n"
    "\tbetween samples.\n"
    "forcedate : bool (optional)\n"
    "\tForces linspace to use the date formulation, needed sometimes on 0-d arrays\n\n"
    "Returns\n"
    "=======\n"
    "samples : array\n"
    "\tThere are `num` equally spaced samples in the closed interval\n"
    "\t``[start, stop]`` or the half-open interval ``[start, stop)``\n"
    "\t(depending on whether `endpoint` is True or False).\n"
    "step : float (only if `retstep` is True)\n"
    "\tSize of spacing between samples.\n\n"
    "See Also\n"
    "========\n"
    "toolbox.geomspace\n"
    "toolbox.logspace\n"},
    // NULL terminate Python looking at the object
     { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_toolbox(void) {
    Py_InitModule3("_toolbox", toolbox_methods,
                     "toolbox module");
    PyDateTime_IMPORT;
    // this is a required function (macro) for using numpy arrays in the module
    import_array();
}

