import methods
import math


def f(x):
    return x**5 - x - 0.2

def df(x):
    return 5*(x**4)-1

def ddf(x):
    return 20*(x**3)

def f1(x):
    return math.sin(x)-x+0.15

def df1(x):
    return -1 + math.cos(x)

def ddf1(x):
    return -math.sin(x)

def f2(x):
    return x** 3 - 18 * x - 83


xn = 1.0
xnm1 = 1.0005
xnm2 = 1.005
eps = 0.1
x0 = -1.00
x1 = 1.00

print("Метод Ньютона",methods.newton(f, df, xn, eps))

print("Метод Чебышева", methods.chebishev(f, df, ddf, xn, eps))

print("Метод Хэлли", methods.helley(f, df, ddf, xn, eps))

print("Метод обратной параболической интерполяции", methods.inverseInterpolation(f, xn, xnm1, xnm2, eps))

print("Для другой функции: ")

print("Метод секущих", methods.secant(f1, x0, x1, eps))

print("Метод Ньютона", methods.newton(f1, df1, xn, eps))

print("Метод Чебышева", methods.chebishev(f1, df1, ddf1, xn, eps))

print("Метод Хэлли", methods.helley(f1, df1, ddf1, xn, eps))

print("Метод обратной параболической интерполяции", methods.inverseInterpolation(f1, xn, xnm1, xnm2, eps))

input()
