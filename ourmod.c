#include <stdio.h>

typedef float (*function)(float x);

float TangentsMethod(function f, function df, float xn, float eps) {
   float x1  = xn - f(xn)/df(xn);
   float x0 = xn;
   while(abs(x0-x1) > eps) {
      x0 = x1;
      x1 = x1 - f(x1)/df(x1);
   }
   return x1;
}
