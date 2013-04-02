/*************************************************************************
# _toolboxmodule
# make some functions in toolbox c extensions to be a lot faster
#
# Brian Larsen
# balarsen@lanl.gov
# Copyright 2010 - 2011 Los Alamos National Security, LLC.
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
// Go through the input arguments or iterator and return sqrt(x**2 + y**2 + ...)
//=============================================================================
static PyObject *hypot_tb(PyObject *self, PyObject *args)
{
    Py_ssize_t TupleSize = PyTuple_Size(args);
    Py_ssize_t i, seqSize;
    npy_intp arraySize;
    PyObject *temp_p, *temp_seq;
    PyArrayObject* temp_array;
    double *temp_data, *data_end;
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
        // if the input is an array ravel it then do the same thing
        if (PyArray_Check(temp_p)) {
            temp_p = PyArray_Ravel((PyArrayObject*)temp_p, 0);
            arraySize = PyArray_SIZE((PyArrayObject*)temp_p);

            if(arraySize >= 8000) {
                temp_array = (PyArrayObject*)PyArray_InnerProduct(temp_p, temp_p);
                temp_seq = PyArray_Cast(temp_array, NPY_DOUBLE);
                Py_DECREF(temp_array);
                tot = *((double*)PyArray_GETPTR1(temp_seq, 0));
                Py_DECREF(temp_seq);
            }
            else {
                // this steals a reference to array_type
                // since this steals a reference you can inline it all without a variable (via BAL Numpy list email)
                temp_array = (PyArrayObject*)PyArray_FromAny(temp_p, PyArray_DescrFromType(NPY_DOUBLE), 0, 0,
                                                            NPY_CARRAY_RO, NULL);
                if (temp_array == NULL) return NULL;
                temp_data = (double*)PyArray_DATA(temp_array);
                data_end = temp_data + arraySize;
                for(; temp_data<data_end; temp_data++)
                    tot+=pow(*temp_data, 2); //tested, same speed as x * x
            }
            return PyFloat_FromDouble(sqrt(tot));
        }
        else if (PySequence_Check(temp_p)) { // is the input a sequence of sorts
            seqSize =  PySequence_Length(temp_p);
            for (i=0;i<seqSize;i++) {
                temp_seq = PySequence_GetItem(temp_p, i);  // returns a new reference
                tmp_d = PyNumber_AsDouble(temp_seq);
                Py_XDECREF(temp_seq);
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
    "vals : float (arbitrary number), or iterable\n"
    "\tarbitary number of float values as arguments or an iterable\n\n"
    "Returns\n"
    "=======\n"
    "out : float\n"
    "\tthe Euclidean distance of the points to the origin\n\n"
    "Examples\n"
    "========\n"
    ">>> import spacepy.toolbox as tb\n"
    ">>> tb.hypot(3,4)\n"
    ">>> # 5.0\n"
    ">>> a = [3, 4]\n"
    ">>> tb.hypot(*a)\n"
    ">>> # 5.0\n"
    ">>> tb.hypot(*range(10))\n"
    ">>> # 16.88194...\n"
    ">>> tb.hypot(range(10))\n"
    ">>> # 16.88194...\n\n"
    "See Also\n"
    "========\n"
    "math.hypot\n"},
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

