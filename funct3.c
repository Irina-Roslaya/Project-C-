#include <stdio.h>
#include <math.h>

typedef double (*function)(double x);

double ChebishevMethod(function f, function df, function ddf, double xn, double eps) {
   double x1  = xn - f(xn)/df(xn)-(pow(f(xn),2)*ddf(xn))/(2*pow(df(xn),3));
   double x0 = xn;
   while(abs(x0-x1) > eps) {
      x0 = x1;
      x1 = x1 - f(x1)/df(x1)-(pow(f(x1),2)*ddf(x1))/(2*pow(df(x1),3));
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
    x = ChebishevMethod(MyFunction, MyDerivative, My2Derivative, xn, 0.1);
    printf ("%lf", x);
    return 0;
}

