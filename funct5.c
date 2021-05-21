#include <stdio.h>
#include <math.h>

typedef double (*function)(double x);

double InverseQuadraticInterpolationMethod(function f,double xn,double xnm1, double xnm2, double eps) {
   double x1  = f(xnm1)*f(xn)*xnm2/((f(xnm2)-f(xnm1))*(f(xnm2)-f(xn))) + f(xnm2)*f(xn)*xnm1/((f(xnm1)-f(xnm2))*(f(xnm1)-f(xn)))+f(xnm2)*f(xnm1)*xn/((f(xn)-f(xnm2))*(f(xn)-f(xnm1)));
   double x0 = xn;
   double xm1 = xnm1;
   double xm2 = xnm2;
   while(abs(x0-x1) > eps) {
      xm2 = xm1;
      xm1 = x0;
      x0 = x1;
      x1 = f(xm1)*f(x0)*xm2/((f(xm2)-f(xm1))*(f(xm2)-f(x0))) + f(xnm2)*f(x0)*xm1/((f(xm1)-f(xm2))*(f(xm1)-f(x0)))+f(xm2)*f(xm1)*x0/((f(x0)-f(xm2))*(f(x0)-f(xm1)));
   }
   return x1;
}
double MyFunction(double x) { return (pow(x, 5) - x - 0.2); }
double MyDerivative(double x) { return (5*pow(x, 4) - 1); }
double My2Derivative(double x) { return (20*pow(x, 3)); }
int main()
{
    double x,xn,xnm1,xnm2;
    xn=1.0;
    xnm1 = 1.0005;
    xnm2 = 1.005;
    x = InverseQuadraticInterpolationMethod(MyFunction,xn,xnm1,xnm2, 0.1);
    printf ("%lf", x);
    return 0;
}
