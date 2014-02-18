/*************************************************************************/
/*Fast, accurate tracing routines for tracing a streamline through a 2D
vector field.  The subroutines included here are meant to be called from
python routines.

Test compilation/execution:
gcc -DCTRACE2D_TEST ctrace2d.c -lm
./a.out

Copyright 2010 - 2014  Los Alamos National Security, LLC. */

/*************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
/* As set now, code will extrapolate up to 1 point outside of current block. */
int DoBreak(int iloc, int jloc, int iSize, int jSize)
{
  int ibreak = 0;
  //printf("At %i, %i with limits %i, %i\n", iloc, jloc, iSize, jSize);
  if (iloc >= iSize || jloc >= jSize)
    ibreak = 1;
  if (iloc < -1 || jloc < -1)
    ibreak = 1;

  return ibreak;
}

/* Create unit vectors of field */
static void make_unit(int iSize, int jSize, double *ux, double *uy)
{
  int i;
  double magnitude;

  for (i=0; i< iSize * jSize; i++) {
	magnitude = sqrt(pow(ux[i],2) + pow(uy[i],2));
	ux[i] /= magnitude;
	uy[i] /= magnitude;
  }

  return;
}

/* Min and MAX: */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/*Interpolate a value from a grid
 *double x, double y: X and Y of location to interp
 *double *field: Field to interpolate on
 *int xloc, int yloc: position of grid point below x,y
 *int xsize, int ysize: Size of field in X and y
 */
#define grid_interp(x, y, field, xloc, yloc, xsize, ysize) \
  bilin_reg(x-xloc, y-yloc, *(field+yloc*xsize+xloc), \
	    *(field+(yloc+1)*xsize+xloc), *(field+yloc*xsize+xloc+1), \
	    *(field+(yloc+1)*xsize+xloc+1))

/* Simple tracing using Euler's method. */
/* Super fast but not super accurate.   */
int cEuler(int iSize, int jSize,           /* Grid size and max steps */
	   int maxstep, double ds,         /* Maxsteps and step size  */
	   double xstart, double ystart,   /* Starting locations      */
	   double xGrid[], double yGrid[], /* Actual coord system     */
	   double *ux, double *uy,         /* Field to trace through  */
	   double x[], double y[])         /* x, y of result stream   */
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
  for(n=0; n<maxstep-1; n++)
    {
      /* Find surrounding points */
      xloc = floor(x[n]);
      yloc = floor(y[n]);

      /* Break if we leave the domain */
      if (DoBreak(xloc, yloc, iSize, jSize))
      break;
      if(xloc>iSize-2) xloc=iSize-2;
      if(xloc<0)       xloc=0;
      if(yloc>iSize-2) yloc=iSize-2;
      if(yloc<0)       yloc=0;
      /* Interpolate unit vectors to current location */
      fx = grid_interp(x[n], y[n], ux, xloc, yloc, iSize, jSize);
      fy = grid_interp(x[n], y[n], uy, xloc, yloc, iSize, jSize);
      
      /* Detect NaNs in function values */
      //if (isnan(fx) || isnan(fy) || isinf(fx) || isinf(fy))
      //break;

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
  //printf("Used %i of %i points\n",n,maxstep);
  //printf("Breaking at %.3f, %.3f\n", x[n], y[n]);
  return n;
}

/* Fast and reasonably accurate tracing with */
/* Runge-Kutta 4 method (constant step size) */
int cRk4(int iSize, int jSize,             /* Grid size and max steps */
	 int maxstep, double ds,           /* Maxsteps and step size  */
	 double xstart, double ystart,     /* Starting locations      */
	 double xGrid[], double yGrid[],   /* Actual coord system     */
	 double *ux, double *uy,           /* Field to trace through  */
	 double x[], double y[])           /* x, y of result stream   */
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
  for(n=0; n<maxstep-1; n++){
    /* See Euler's method for more descriptive comments. */
    /* SUBSTEP #1 */
    xloc = floor(x[n]);
    yloc = floor(y[n]);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f1x = grid_interp(x[n], y[n], ux, xloc, yloc, iSize, jSize);
    f1y = grid_interp(x[n], y[n], uy, xloc, yloc, iSize, jSize);
    if (isnan(f1x) || isnan(f1y) || isinf(f1x) || isinf(f1y))
      break;

    /* SUBSTEP #2 */
    xpos = x[n]+f1x*ds/2.0;
    ypos = y[n]+f1y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f2x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f2y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);

    if (isnan(f2x) || isnan(f2y) || isinf(f2x) || isinf(f2y))
      break;

    /* SUBSTEP #3 */
    xpos = x[n]+f2x*ds/2.0;
    ypos = y[n]+f2y*ds/2.0;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f3x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f3y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);
    if (isnan(f3x) || isnan(f3y) || isinf(f3x) || isinf(f3y))
      break;

    /* SUBSTEP #4 */
    xpos = x[n]+f3x*ds;
    ypos = y[n]+f3y*ds;
    xloc = floor(xpos);
    yloc = floor(ypos);
    if (DoBreak(xloc, yloc, iSize, jSize))
      break;
    if(xloc>iSize-2) xloc=iSize-2;
    if(xloc<1)       xloc=1;
    if(yloc>iSize-2) yloc=iSize-2;
    if(yloc<1)       yloc=1;
    //xloc=MAX(0, MIN(iSize-2, xloc));
    //yloc=MAX(0, MIN(jSize-2, yloc));
    f4x = grid_interp(xpos, ypos, ux, xloc, yloc, iSize, jSize);
    f4y = grid_interp(xpos, ypos, uy, xloc, yloc, iSize, jSize);
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
  return n;

}

#ifdef CTRACE2D_TEST
static int test_arrays(int iSize, int jSize, double xGrid[], double yGrid[],
		       double *ux, double *uy){
  printf("Seems to work.\n");
  return 1;
}

/* A utility for printing test results. */
void test_check(double result, double answer)
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
			 ux, uy, xt, yt);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);

  npoints = cRk4(nx, ny, maxstep, ds, 1.0, 10.0, xgrid, ygrid,
		 ux, uy, xt, yt);
  printf("Npoints = %i\n", npoints);
  printf("Grid goes from %.2f to %.2f\n", xgrid[0], xgrid[nx-1]);
  printf("Our trace starts at %.2f, %.2f\n", xt[0], yt[0]);
  printf("...and ends at      %.2f, %.2f\n", xt[npoints], yt[npoints]);  

  /*for (i=0; i<npoints; i++)
    printf("%.3f, %.3f\n", xt[i], yt[i]);*/
  return 0;
}
#endif /*CTRACE2D_TEST*/
