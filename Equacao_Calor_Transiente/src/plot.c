#include "../include/plot.h"

// Plotar em cada instante de tempo -- Para t = t_k --> u(x,t_k)
void writeDataFile (double *d, double *X, Func *f, double t, int np, int iter)
{
  FILE *file_solution;
  char filename[50];
  int i;
  if (iter % 10 == 0)
  {
    sprintf(filename,"solution%d.dat",iter);
    file_solution = fopen(filename,"w+");
    
    // || Coluna 1 = Xi || Coluna 2 = Solucao analitica || Coluna 3 = Solucao aproximada ||
    for (i = 0; i < np; i++)
      fprintf(file_solution,"%e %e %e\n",X[i],f[2](iter,t,X[i]),d[i]);
    fclose(file_solution);

  }
}

// Plotar em superficie
void writeDataFile2 (double *d, double *X, Func *f, double t, int np, int iter)
{
  FILE *file_solution;
  char filename[50];
  int i;
  file_solution = fopen("solution.dat","a");
  // || Coluna 1 = t || Coluna 2 = Xi || Coluna 3 = Solucao aproximada ||
  for (i = 0; i < np; i++)
    fprintf(file_solution,"%e %e %e\n",t,X[i],d[i]);
    //fprintf(file_solution,"%e %e %e\n",X[i],f[2](iter,t,X[i]),d[i]);
  fclose(file_solution);
}

// Imprime uma matriz NxN
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

// Imprime um vetor Nx1
void printVector (char *str, double *b, int N)
{
  int i;
  printf("\n%s\n",str);
  for (i = 0; i < N; i++)
    printf("%e\n",b[i]);
  printf("\n");
}
