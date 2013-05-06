/*
C-code for solving the 1-D diffusion equation:
df/dt = L^2 * d/dL [ D/L^2 (df/dL)]
using a Cranck-Nicholson method.
Originally written by Enrico Camporeale.
Slight edits by the Spacepy team.

How to incorporate the library into a Python script without Spacepy:
1)Build the library:
 gcc -shared -o cn.so cn_solver.c

2)Import ctypes and the library:
 from ctypes import *
 cn=CDLL("cn.so")
 cn.solve_cn.restype=None

How to incorporate the library using Spacepy:
1)Install spacepy as normal.

2)Import the function from the spacepy.lib library.
  from spacepy.lib import solve_cn

3)Check success of compilation/import.
  if spacepy.lib.has_libspacepy then:
      print "Success!"

How to apply the function from Python:
1)Call the routine, using either ctypes or numpy to convert double precision 
arguments to the correct c-type.
  dptr = ctypes.POINTER(ctypes.c_double)
  solve_cn(f.ctypes.data_as(dptr),
           Lgrid.ctypes.data_as(dptr),
	   DLLm.ctypes.data_as(dptr), 
	   DLLm.ctypes.data_as(dptr), 
	   DLLp.ctypes.data_as(dptr), 
	   DLLp.ctypes.data_as(dptr), 
	   ct.c_double(Tdelta), 
	   NL,
	   S.ctypes.data_as(dptr))

The argument list is defined as follows:
f is the distribution function being propagated forward in time; a vector of
length nL.
L is the L grid.
The next four arguments are the temporally and spatially dependent 
diffusion coefficients at present time (old), present time + dt (new), and 
in between each grid point (m for "minus", p for "plus" one half grid spacing).
Each is a vector of the same length as L.
dt is the timestep, a scalar double.
N is the number of gridpoints as a scalar integer.
The f vector is returned at time+dt.
S is a vector of source and loss terms prepared such that the following system:
A*f(t+dt)=B*f(t)+S
...can be solved by simply performing
f(t+dt)=A_inv*(B*f(t)+S)
This means that, following the C-N derivation, S should be dt/2*(S(t)+S(t+dt)).
This must be done before passing the vector to this routine.

caveats: 
This code allows to input the diffusion coefficient D at two times, 
when D is a function of time. One can always put the two equals.
The boundary conditions are Dirichlet type (i.e. f=0)
The routine solves the linear system by performing a LU decomposition, 
and returns the vector f at the new time.

Copyright 2010 Los Alamos National Security, LLC.
*/


#include <math.h>
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#include <stdio.h>

// Cranck-Nicholson solver
void solve_cn(double *f, double *L, double *Dm_old,  double *Dm_new, 
	      double *Dp_old, double *Dp_new, double dt, int N, double *S){
double a[N][3],r[N],al[N];
double betam[N],betap[N],Lp[N],Lm[N],Dllp[N],Dllm[N];
int i,l,k,j;
int indx[N];
double d,dum;
double dL=L[1]-L[0]; // an uniform grid !
int mm=3;

//define gamma=D/L^2
for (i=0;i<N;i++){
  Lp[i]=L[i]+0.5*dL;
  Lm[i]=L[i]-0.5*dL;
  Dllm[i] = 0.5*(Dm_old[i]+Dm_new[i]);
  Dllp[i] = 0.5*(Dp_old[i]+Dp_new[i]);}
  
 for (i=1;i<N-1;i++){
  betam[i] = Dllm[i]*dt / (2*dL*dL*Lm[i]*Lm[i]);
  betap[i] = Dllp[i]*dt / (2*dL*dL*Lp[i]*Lp[i]);}
  
  
// define a tridiagonal matrix with vectors a[..][0] (subdiagonal), 
// a[..][1](diagonal), a[..][2](superdiagonal)
// r is the known term (right-hand side term) plus source...

for (i=1;i<N-1;i++) {
a[i][0]=-betam[i]*L[i]*L[i];
a[i][1]=1 + betam[i]*L[i]*L[i] + betap[i]*L[i]*L[i];
a[i][2]=-betap[i]*L[i]*L[i];
r[i]=betam[i]*L[i]*L[i]*f[i-1]+
  (1 - betam[i]*L[i]*L[i] - betap[i]*L[i]*L[i])*f[i]+betap[i]*L[i]*L[i]*f[i+1] +
  S[i];
}

// dirichlet BC
a[0][1]=1.0;a[0][2]=0.0;r[0]=0.0;
a[N-1][1]=1.0;a[N-1][0]=0.0;r[N-1]=0.0;

// LU decomposition of matrix A

a[0][0]=a[0][1];
a[0][1]=a[0][2];
a[0][2]=0.0;

d=1.0;l=1;
for (k=0;k<N;k++){
    dum=a[k][0];
    i=k;
    if (l<N) l++;
    for (j=k+1;j<l;j++){
      if (fabs(a[j][0]) > fabs(dum)) {
	dum=a[j][0];
	i=j;
      }
    }
    indx[k]=i+1;
    if (dum==0.0) a[k][0]=1.0e-20;
    if (i!=k){
      d=-d;
      for (j=0;j<3;j++) SWAP(a[k][j],a[i][j]);
    }
    for (i=k+1;i<l;i++){
      dum=a[i][0]/a[k][0];
      al[k]=dum;
      for (j=1;j<3;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
      a[i][2]=0.0;
    }
}
// end of LU decomposition; upper triangular is now in a; and lower triangular in al;



// advance timestep  this is the solver: solution overwrites f

l=1;

for (k=0;k<N;k++){
  j=indx[k]-1;
  if (j!=k) SWAP(r[k],r[j]);
  if (l<N) l++;
  for (j=k+1;j<l;j++) r[j]-=al[k]*r[k];
}
l=1;
for (i=N-1;i>=0;i--){
  dum=r[i];
  for (k=1;k<l;k++) dum-=a[i][k]*r[k+i];
  r[i]=dum/a[i][0];//  r[i]=dum/a[i][0];
  if (l<mm) l++;
}
for (i=0;i<N;i++)
f[i]=r[i];
  
    
}
