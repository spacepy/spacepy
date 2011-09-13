/* Test functions for tracing routines
 *
 * Compile and link to test
 * routines in ctrace2dmodule.c
 * e.g.:
 * gcc -I /usr/include/python2.6 -lm ctrace2dtest.c
 *
 * Copyright 2010 - 2011 Los Alamos National Security, LLC
 */


/*Because internal functions in the Python module must use static
 *linkage (not pollute the namespace), this has to be an include
 *rather than linking against another object.
 */
#include "ctrace2dmodule.c"

static int test_arrays(int iSize, int jSize, double xGrid[], double yGrid[],
		       double *ux, double *uy){
  printf("Seems to work.\n");
  return 1;
}

/* A utility for printing test results. */
static void test_check(double result, double answer)
{
  double diff, thresh;

  thresh = 0.00001;
  diff = 100.0*(result - answer) / answer;

  if (diff < thresh)
    printf("TEST PASSED! (Result=%.4f, Answer=%.4f)\n", result, answer);
  else{
    printf("TEST FAILED!\n");
    printf("Result %.20f differs from answer %.20f\n", result, answer);
    printf("Difference of %.3f%% is over threshold of %.5f%%\n", diff, thresh);
  }
}

/* Main simply tests the functionality of all funcs in this file.*/
int main()
{
  /* Declarations */
  double out=0, x=0.1, y=0.2, Q00=3.0, Q01=5.0, Q10=40.0, Q11=60.0;
  double sol1 = 7.460;
  int l[1];

  /* Test bilin_reg */
  printf("Testing bilin_reg\n");
  out = bilin_reg(x, y, Q00, Q01, Q10, Q11);
  printf("TEST 1: ");
  test_check(out, sol1);


  /* Test cEuler 1 */
  int i, j, nx=841, ny=121, maxstep=10000, npoints;
  double xgrid[nx], ygrid[nx], *ux, *uy, xt[maxstep], yt[maxstep], ds=1.0;
  ux = malloc(nx * ny * sizeof(double));
  uy = malloc(nx * ny * sizeof(double));
  
  for (i=0; i<nx; i++)
    {
      xgrid[i] = -10.0+0.25*i;
      ygrid[i] = xgrid[i];
    }

  for (i=0; i<nx; i++)
    for (j=0; j<ny; j++)
      {
	*(ux+i*ny+j) = xgrid[i];
	*(uy+i*ny+j) = -1.0 * ygrid[j];
      }

  npoints = cEuler(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
			 ux, uy, xt, yt, l);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);

  npoints = cRk4(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
		 ux, uy, xt, yt, l);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);  

  /*for (i=0; i<npoints; i++)
    printf("%.3f, %.3f\n", xt[i], yt[i]);*/
  return 0;
}
