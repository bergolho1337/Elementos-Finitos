#include "integration.h"

double F (int i, int j, Func f1, Func f2, double x)
{
  return (f1(i,x) * f2(j,x));
}

double trapezoidalRule (int i, int j, Func f1, Func f2, double a, double b)
{
  int k, n;
  double I, x;
  n = (b-a) / H;
  I =  (F(i,j,f1,f2,a) + F(i,j,f1,f2,b)) / 2.0;
  for (k = 1; k < n-1; k++)
  {
    x = a + k*H;
    I += F(i,j,f1,f2,x);
  }
  return I*H;
}

double integral_polinomium (int i, int j)
{
  return (1/(double)(i+j+1));
}
