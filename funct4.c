#include <stdio.h>
#include <math.h>

typedef double (*function)(double x);

double HelleyMethod(function f, function df, function ddf, double xn, double eps) {
   double x1  = xn - (2*f(xn)*df(xn)/(2*pow(df(xn),2)-f(xn)*ddf(xn)));
   double x0 = xn;
   while(abs(x0-x1) > eps) {
      x0 = x1;
      x1 = x1 - (2*f(x1)*df(x1)/(2*pow(df(x1),2)-f(x1)*ddf(x1)));
   }
   return x1;
}
double MyFunction(double x) { return (pow(x, 5) - x - 0.2); }
double MyDerivative(double x) { return (5*pow(x, 4) - 1); }
double My2Derivative(double x) { return (20*pow(x, 3)); }
int main()
{
    double x,xn;
    xn=1.0;
    x = HelleyMethod(MyFunction, MyDerivative, My2Derivative, xn, 0.1);
    printf ("%lf", x);
    return 0;
}
