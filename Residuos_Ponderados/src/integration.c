#include "integration.h"

double F (Func f1, Func f2, double x)
{
  return (f1(x) * f2(x));
}

double trapezoidalRule (Func f1, Func f2, double a, double b)
{
  int i, n;
  double I, x;
  n = (b-a) / H;
  I =  (F(f1,f2,a) + F(f1,f2,b)) / 2.0;
  for (i = 1; i < n-1; i++)
  {
    x = a + i*H;
    I += F(f1,f2,x);
  }
  return I*H;
}
