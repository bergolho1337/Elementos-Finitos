#include "plot.h"

double Faprox (double x, Func *f, double *alfa, int N)
{
  int k;
  double eval = 0;
  for (k = 0; k < N; k++)
    eval += alfa[k]*f[0](k,x);
  return eval;
}

double Fanalit (double x, double a, double b, double c)
{
  return ((exp(b*x/a)-1) / (exp(b/a)-1));
}

void writeDataFile (Func *f, double *alfa, double a, double b, double c, double dx, int N)
{
  FILE *file_solution, *file_residue;
  int i, n;
  double x, fanalit, faprox, residue;

  // Se a numeracao foi feita manualmente acertar a posicao dos coeficientes
  sortNumeration(alfa,N);

  file_solution = fopen("solution.dat","w+");
  file_residue = fopen("residue.dat","w+");
  // Escreve a aproximacao
  n = 1 / dx;
  for (i = 0; i <= n; i++)
  {
    x = i*dx;
    fanalit = Fanalit(x,a,b,c);
    faprox = Faprox(x,f,alfa,N);
    residue = fabs(fanalit-faprox);
    // || x || analitica || aproximacao || residuo ||
    fprintf(file_solution,"%e %e %e %e\n",x,fanalit,faprox,residue);
    fprintf(file_residue,"%e %e\n",x,residue);
  }
  fclose(file_solution);
  fclose(file_residue);
}

void printMatrix (char *str, double *A, int N)
{
  int i, j;
  printf("\n%s\n",str);
  for (i = 0; i < N; i++)
  {
    printf("\n");
    for (j = 0; j < N; j++)
    {
      printf("%e ",A[i * N + j]);
    }
  }
  printf("\n");
}

void printVector (char *str, double *b, int N)
{
  int i, j;
  printf("\n%s\n",str);
  for (i = 0; i < N; i++)
    printf("%e\n",b[i]);
  printf("\n");
}
