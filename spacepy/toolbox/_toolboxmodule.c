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

#define PyNumber_AsDouble(y) PyFloat_AsDouble(PyNumber_Float(y))


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
    // NULL terminate Python looking at the object
     { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_toolbox(void) {
    Py_InitModule3("_toolbox", toolbox_methods,
                     "toolbox module");

}

