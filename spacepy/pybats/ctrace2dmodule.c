/*************************************************************************/
/*Fast, accurate tracing routines for tracing a streamline through a 2D
vector field.  The subroutines included here are meant to be called from
python routines.

Test compilation:
gcc -DNDEBUG -g -O3 -Wall -Wstrict-prototypes -fPIC -DMAJOR_VERSION=1 -DMINOR_VERISON=0 -I /usr/include/python2.6 -I /usr/share/pyshared/numpy/core/include/numpy/ -c ctrace2dmodule.c
gcc -shared ctrace2dmodule.o -o ctrace2d.so

Copyright 2010 - 2011  Los Alamos National Security, LLC. */

/*************************************************************************/
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>

/* Bilinear interpolation for x1,y1=0 and x2,y2=1     */
/* Q's are surrounding points such that Q00 = F[0,0], */
/* Q10 = F[1,0], etc.                                 */
static double bilin_reg(double x, double y, double Q00, 
			double Q01, double Q10, double Q11)
{
  double fout;
  //printf("Input = %.3f, %.3f\n",x, y);
  fout = 
    Q00*(1.0-x)*(1.0-y) +
    Q10* x *    (1.0-y) +
    Q01* y *    (1.0-x) +
    Q11* x * y;

//  printf("points = %.5f, %.5f, %.5f, %.5f; ans = %.5f\n",
//	 Q00, Q10, Q01, Q11, fout);

  return fout;
}

/* Check to see if we should break out of an integration */
static int DoBreak(int iloc, int jloc, int iSize, int jSize)
{
  int ibreak = 0;
  //printf("At %i, %i with limits %i, %i\n", iloc, jloc, iSize, jSize);
  if (iloc >= iSize-1 || jloc >= jSize-1)
    ibreak = 1;
  if (iloc <= 0 || jloc <= 0)
    ibreak = 1;

  return ibreak;
}

/* Create unit vectors of field */
static void make_unit(int iSize, int jSize, double *ux, double *uy)
{
  int i, j;
  double magnitude;

  for (i=0; i<iSize; i++)
    for (j=0; j<jSize; j++)
      {
	magnitude = sqrt(pow(*(ux+i*jSize+j),2)+pow(*(uy+i*jSize+j),2));
	*(ux+i*jSize+j) = *(ux+i*jSize+j) / magnitude;
	*(uy+i*jSize+j) = *(uy+i*jSize+j) / magnitude;
      }

  return;
}

/* Simple tracing using Euler's method. */
/* Super fast but not super accurate.   */
static int cEuler(int iSize, int jSize,           /* Grid size and max steps */
		  int maxstep, double ds,         /* Maxsteps and step size  */
		  double xstart, double ystart,   /* Starting locations      */
		  double xGrid[], double yGrid[], /* Actual coord system     */
		  double *ux, double *uy,         /* Field to trace through  */
		  double x[], double y[], int l[])/* x, y of result stream   */
{
  /* Variable declarations */
  int n=0, i, xloc, yloc;
  double dx, dy, fx, fy;

  /* Get starting points in normalized/array coordinates */
  dx = xGrid[1] - xGrid[0];
  dy = yGrid[1] - yGrid[0];
  x[0] = (xstart-xGrid[0]) / dx;
  y[0] = (ystart-yGrid[0]) / dy;
  
  /* Create unit vectors from full vector field */
  make_unit(iSize, jSize, ux, uy);

  /* Perform tracing using Euler's method */
  for(n=0; n<maxstep; n++)
    {
      /* Find surrounding points */
      xloc = floor(x[n]);
      yloc = floor(y[n]);

      /* Break if we leave the domain */
      if (DoBreak(xloc, yloc, iSize, jSize))
	break;

      /* Interpolate unit vectors to current location */
      fx = bilin_reg(x[n]-xloc, y[n]-yloc, *(ux+xloc*jSize+yloc), 
		     *(ux+xloc*jSize+yloc+1), *(ux+(xloc+1)*jSize+yloc), 
		     *(ux+(xloc+1)*jSize+yloc+1));
      fy = bilin_reg(x[n]-xloc, y[n]-yloc, *(uy+xloc*jSize+yloc), 
		     *(uy+xloc*jSize+yloc+1), *(uy+(xloc+1)*jSize+yloc), 
		     *(uy+(xloc+1)*jSize+yloc+1)); 
      
      /* Detect NaNs in function values */
      if (isnan(fx) || isnan(fy) || isinf(fx) || isinf(fy))
	break;

      /* Perform single step */
      x[n+1] = x[n] + ds * fx;
      y[n+1] = y[n] + ds * fy;
    }

  /* Return traced points to original coordinate system. */
  for (i=0; i<=n; i++)
    {
      x[i] = x[i]*dx + xGrid[0];
      y[i] = y[i]*dy + yGrid[0];
    }
  //printf("Used %i points\n",n);
  //printf("Breaking at %.3f, %.3f\n", x[n], y[n]);
  l[0] = n;
  return n;
}

/* Fast and reasonably accurate tracing with */
/* Runge-Kutta 4 method (constant step size) */
static int cRk4(int iSize, int jSize,             /* Grid size and max steps */
		int maxstep, double ds,           /* Maxsteps and step size  */
		double xstart, double ystart,     /* Starting locations      */
		double xGrid[], double yGrid[],   /* Actual coord system     */
		double *ux, double *uy,           /* Field to trace through  */
		double x[], double y[], int l[])  /* x, y of result stream   */
{
  /* Variable declarations */
  int n=0, i, xloc, yloc;
  double dx, dy, xpos, ypos,
    f1x, f2x, f3x, f4x, f1y, f2y, f3y, f4y;

  /* Get starting points in normalized/array coordinates */
  dx = xGrid[1] - xGrid[0];
  dy = yGrid[1] - yGrid[0];
  x[0] = (xstart-xGrid[0]) / dx;
  y[0] = (ystart-yGrid[0]) / dy;
  
  /* Create unit vectors from full vector field */
  make_unit(iSize, jSize, ux, uy);

  /* Perform tracing using RK4 */
  for(n=0; n<maxstep; n++){
    /* See Euler's method for more descriptive comments. */
    /* SUBSTEP #1 */
    xloc = floor(x[n]);
    yloc = floor(y[n]);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    f1x = bilin_reg(x[n]-xloc, y[n]-yloc, *(ux+xloc*jSize+yloc), 
		   *(ux+xloc*jSize+yloc+1), *(ux+(xloc+1)*jSize+yloc), 
		   *(ux+(xloc+1)*jSize+yloc+1));
    f1y = bilin_reg(x[n]-xloc, y[n]-yloc, *(uy+xloc*jSize+yloc), 
		   *(uy+xloc*jSize+yloc+1), *(uy+(xloc+1)*jSize+yloc), 
		   *(uy+(xloc+1)*jSize+yloc+1));
    if (isnan(f1x) || isnan(f1y) || isinf(f1x) || isinf(f1y))
      break;

    /* SUBSTEP #2 */
    xpos = x[n]+f1x*ds/2.0;
    ypos = y[n]+f1y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    f2x = bilin_reg(xpos-xloc, ypos-yloc, *(ux+xloc*jSize+yloc), 
		   *(ux+xloc*jSize+yloc+1), *(ux+(xloc+1)*jSize+yloc), 
		   *(ux+(xloc+1)*jSize+yloc+1));
    f2y = bilin_reg(x[n]-xloc, y[n]-yloc, *(uy+xloc*jSize+yloc), 
		   *(uy+xloc*jSize+yloc+1), *(uy+(xloc+1)*jSize+yloc), 
		   *(uy+(xloc+1)*jSize+yloc+1));
    if (isnan(f2x) || isnan(f2y) || isinf(f2x) || isinf(f2y))
      break;

    /* SUBSTEP #3 */
    xpos = x[n]+f2x*ds/2.0;
    ypos = y[n]+f2y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    f3x = bilin_reg(xpos-xloc, ypos-yloc, *(ux+xloc*jSize+yloc), 
		   *(ux+xloc*jSize+yloc+1), *(ux+(xloc+1)*jSize+yloc), 
		   *(ux+(xloc+1)*jSize+yloc+1));
    f3y = bilin_reg(x[n]-xloc, y[n]-yloc, *(uy+xloc*jSize+yloc), 
		   *(uy+xloc*jSize+yloc+1), *(uy+(xloc+1)*jSize+yloc), 
		   *(uy+(xloc+1)*jSize+yloc+1));
    if (isnan(f3x) || isnan(f3y) || isinf(f3x) || isinf(f3y))
      break;

    /* SUBSTEP #4 */
    xpos = x[n]+f3x*ds;
    ypos = y[n]+f3y*ds;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    f4x = bilin_reg(xpos-xloc, ypos-yloc, *(ux+xloc*jSize+yloc), 
		   *(ux+xloc*jSize+yloc+1), *(ux+(xloc+1)*jSize+yloc), 
		   *(ux+(xloc+1)*jSize+yloc+1));
    f4y = bilin_reg(x[n]-xloc, y[n]-yloc, *(uy+xloc*jSize+yloc), 
		   *(uy+xloc*jSize+yloc+1), *(uy+(xloc+1)*jSize+yloc), 
		   *(uy+(xloc+1)*jSize+yloc+1));
    if (isnan(f4x) || isnan(f4y) || isinf(f4x) || isinf(f4y))
      break;

    /* Peform the full step using all substeps */
    x[n+1] = (x[n] + ds/6.0 * (f1x + f2x*2.0 + f3x*2.0 + f4x));
    y[n+1] = (y[n] + ds/6.0 * (f1y + f2y*2.0 + f3y*2.0 + f4y));

  }

  /* Return traced points to original coordinate system. */
  for (i=0; i<=n; i++)
    {
      x[i] = x[i]*dx + xGrid[0];
      y[i] = y[i]*dy + yGrid[0];
    }
  l[0] = n;
  return n;

}

/*Since the two tracers have identical interfaces,
 *almost all the code is common
 */
static PyObject *ctrace2d_common(PyObject *self, PyObject *args,
				 int (*func)(int, int, int,
					     double, double, double,
					     double*, double*,
					     double*, double*,
					     double*, double*, int*)) {
  PyArrayObject *gridx, *gridy, *fieldx, *fieldy, *outx, *outy;
  /*Data pointers for the above arrays*/
  double *gridxd, *gridyd, *fieldxd, *fieldyd, *outxd, *outyd;
  double xstart, ystart, ds;
  int maxstep, xsize, ysize, count, lout;
  PyArray_Descr *array_type;
  npy_intp outdims[] = {0};
  npy_intp indims[] = {0};
  PyArray_Dims outshape = { outdims, 1 };

  if (!PyArg_ParseTuple(args,
			"idddO!O!O!O!:cEuler",
			&maxstep, &ds, &xstart, &ystart,
			&PyArray_Type, &gridx, &PyArray_Type, &gridy,
			&PyArray_Type, &fieldx, &PyArray_Type, &fieldy))
    return NULL;

  array_type = PyArray_DescrFromType(NPY_DOUBLE);
  /*For all of these, we are throwing away the borrowed ref
   *to the original, and creating a new object with a new ref.
   *So the new ref will be freed, but the borrowed ref is left alone.
   */
  gridx = (PyArrayObject*)PyArray_FromArray(gridx, array_type, NPY_DEFAULT);
  gridy = (PyArrayObject*)PyArray_FromArray(gridy, array_type, NPY_DEFAULT);
  fieldx = (PyArrayObject*)PyArray_FromArray(fieldx, array_type, NPY_DEFAULT);
  fieldy = (PyArrayObject*)PyArray_FromArray(fieldy, array_type, NPY_DEFAULT);
  xsize = (int)PyArray_DIM(gridx, 0);
  ysize = (int)PyArray_DIM(gridy, 0);
  indims[0] = maxstep;
  outx = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, indims, array_type);
  outy = (PyArrayObject *)PyArray_SimpleNewFromDescr(1, indims, array_type);
  Py_DECREF(array_type);

  gridxd = (double*)PyArray_DATA(gridx);
  gridyd = (double*)PyArray_DATA(gridy);
  fieldxd = (double*)PyArray_DATA(fieldx);
  fieldyd = (double*)PyArray_DATA(fieldy);
  outxd = (double*)PyArray_DATA(outx);
  outyd = (double*)PyArray_DATA(outy);

NPY_BEGIN_ALLOW_THREADS
  count = (*func)(xsize, ysize, maxstep, ds, xstart, ystart,
		 gridxd, gridyd, fieldxd, fieldyd, outxd, outyd,
		 &lout);
NPY_END_ALLOW_THREADS

  Py_DECREF(gridx);
  Py_DECREF(gridy);
  Py_DECREF(fieldx);
  Py_DECREF(fieldy);

  outdims[0] = count;
  if (!PyArray_Resize(outx, &outshape, 1, NPY_CORDER))
    return NULL;
  if (!PyArray_Resize(outy, &outshape, 1, NPY_CORDER))
    return NULL;
  /*Giving away our reference to the caller*/
  return Py_BuildValue("NN", outx, outy);
}

static PyObject *ctrace2d_cEuler(PyObject *self, PyObject *args) {
  return ctrace2d_common(self, args, cEuler);
}

static PyObject *ctrace2d_cRk4(PyObject *self, PyObject *args) {
  return ctrace2d_common(self, args, cRk4);
}

static PyMethodDef ctrace2d_methods[] = {
   { "cEuler", (PyCFunction)ctrace2d_cEuler, METH_VARARGS, NULL },
   { "cRk4", (PyCFunction)ctrace2d_cRk4, METH_VARARGS, NULL },
   { NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC initctrace2d(void) {
   Py_InitModule("ctrace2d", ctrace2d_methods);
   import_array();
}
