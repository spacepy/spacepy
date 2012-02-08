
//   gcc -shared -o hypot.so hypot.c

#include <math.h>

double hypot_c(double *inval, long len){
    double ans=0.0;
    long i;
    for (i=0;i<len;i++)
        ans += (inval[i]*inval[i]);
    return (sqrt(ans));
}




