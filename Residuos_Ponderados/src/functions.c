#include "functions.h"

double f (double x)
{
  return x;
}

/* ============================================================== */
// Funcoes de ponderacao
double phi_1 (double x)
{
  //return (pow(x,2) - x);
  return (sin(M_PI * x));
}

double phi_2 (double x)
{
  //return (pow(x,3) - pow(x,2));
  return (sin(2 * M_PI * x));
}

double phi_3 (double x)
{
  //return (pow(x,4) - pow(x,3));
  return (sin(3 * M_PI * x));
}

/* ============================================================== */
// Funcoes do operador linear
double L_1 (double x)
{
  //return (pow(x,2) - x - 2);
  return (sin(M_PI * x)*(pow(M_PI,2) + 1));
}

double L_2 (double x)
{
  //return (-6*x + 2 + pow(x,3) - pow(x,2));
  return (sin(2 * M_PI * x)*(4 * pow(M_PI,2) + 1));
}

double L_3 (double x)
{
  //return (-12*pow(x,2) + 6*x + pow(x,4) + pow(x,3));
  return (sin(3 * M_PI * x)*(9 * pow(M_PI,2) + 1));
}
