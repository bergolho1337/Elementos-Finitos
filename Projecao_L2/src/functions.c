#include "functions.h"

// Conjunto de pontos igualmente espacados para montar os polinomios de Lagrange (para o x_i da formula)
double *points;

// Diferenca entre dois pontos consecutivos
double h;

// Numero de pontos
int n;

// >>>>>>>>>>>>>> FUNCAO A SER APROXIMADA <<<<<<<<<<<<<<<<<<<<<<<<<<
double f (int k, double x)
{
  return sin(M_PI*x);
  //return pow(x,3)*(x-1)*(1-2*x);          // Questao 5a
  //return atan(x-0.5) / 0.01;              // Questao 5b
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

double p_k (int k, double x)
{
  // Usando p(x) = x^k
  return pow(x,k);
}

double l_k (int k, double x)
{
  int i;
  double sum = 1;
  for (i = 0; i < n; i++)
  {
    if (i != k)
      sum *= (x - points[i]) / (points[k] - points[i]);
  }
  return sum;
}

double phi_k (int k, double x)
{
  // phi_0
  if (k == 0)
  {
    if (x >= points[0] && x <= points[1])
      return (points[1]-x) / h;
    else
      return 0;
  }
  // phi_N
  else if (k == n)
  {
    if (x >= points[n-2] && x <= points[n-1])
      return (x-points[n-1]) / h;
    else
      return 0;
  }
  // phi_k
  else
  {
    if (x >= points[k] && x <= points[k+1])
      return (points[k+1]-x) / h;
    else if (x >= points[k-1] && x <= points[k])
      return (x-points[k-1]) / h;
    else
      return 0;
  }
}

// Deve ser somente os pontos interiores
void setPoints (int N)
{
  int i;
  double x, dx;
  points = (double*)malloc(sizeof(double)*N);
  dx = 1 / (double)(N-1);
  h = dx;
  n = N;
  for (i = 0; i < N; i++)
  {
    x = i*dx;
    points[i] = x;
  }
  //for (i = 0; i < N; i++)
  //  printf("%e\n",points[i]);
}
