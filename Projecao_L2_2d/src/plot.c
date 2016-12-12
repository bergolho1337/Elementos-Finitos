#include "plot.h"

void writeDataFile (Func f, double *alfa, double *P,  int np)
{
  FILE *file_analit, *file_solution;
  int i, j, n;
  double x, y, dx;
  dx = 1.0 / 500.0;
  file_analit = fopen("analit.dat","w+");
  file_solution = fopen("solution.dat","w+");

  for (i = 0; i < 500; i++)
  {
    x = i*dx;
    for (j = 0; j < 500; j++)
    {
      y = j*dx;
      fprintf(file_analit,"%e %e %e\n",x,y,f(0,x,y));
    }
  }
  fclose(file_analit);

  for (i = 0; i < np; i++)
    fprintf(file_solution,"%e %e %e\n",P[i*2],P[i*2+1],alfa[i]);
  fclose(file_solution);
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
