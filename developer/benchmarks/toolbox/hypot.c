
//   gcc -shared -o hypot.so hypot.c

#include <math.h>

double hypot_c(double *inval, long len){
    double *curr;
    double ans=0.0;
    double *max=inval+len;
    for (curr=inval;curr<max;curr++) {
        ans += (*curr * *curr);
    }
    return (sqrt(ans));
}




