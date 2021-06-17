#include <stdio.h>
#include <math.h>
typedef double (*function)(double x);

double SecantMethod(function f, double x0, double x1, double eps)
{
    int n;
    double t;
    n = 0;
    while(fabs(x1 - x0) > eps)
    {
     t = x1;
     x1 = x1 - (x1 - x0)*f(x1)/(f(x1) - f(x0));
     x0 = t;
     n++;
    }
    printf ("%lf\n", x1);
    return n;
}


double x0 = -1;
double x1 = 1;
double eps = 0.00001;

double MyFunction(double x)
{
    return (sin(x)-x+0.15);
}

int main()
{
    printf ("%lf", SecantMethod(MyFunction, x0, x1, eps));
    return 0;
}
