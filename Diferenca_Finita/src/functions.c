#include "functions.h"

// Conjunto de pontos igualmente espacados (para computar a funcao chapeu dos elementos)
double *points;

// Diferenca entre dois pontos consecutivos
double h;

// Numero de pontos
int n;

// >>>>>>>>>>>>>> FUNCAO A SER APROXIMADA <<<<<<<<<<<<<<<<<<<<<<<<<<
double f (int k, double x)
{
  return 0;
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

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

void setPoints (int N, int numeracao)
{
  int i, j;
  double x, dx, aux;
  points = (double*)malloc(sizeof(double)*N);
  dx = 1 / (double)(N-1);
  h = dx;
  n = N;
  // Automatica (numeracao padrao)
  if (numeracao == 1)
  {
    for (i = 0; i < N; i++)
    {
      x = i*dx;
      points[i] = x;
    }
  }
  // Manual
  else
  {
    printf("Digite as coordenadas dos pontos do dominio:\n");
    for (i = 0; i < N; i++)
    {
      printf("Ponto %d: ",i);
      scanf("%lf",&points[i]);
    }
  }
  //for (i = 0; i < N; i++)
  //  printf("%e\n",points[i]);
}

void sortNumeration (double *alfa, int N)
{
  int i, j;
  double aux;
  for (i = 0; i < N; i++)
  {
    for (j = i+1; j < N; j++)
    {
      if (points[i] > points[j])
      {
        aux = points[i];
        points[i] = points[j];
        points[j] = aux;

        aux = alfa[i];
        alfa[i] = alfa[j];
        alfa[j] = aux;
      }
    }
  }
}
