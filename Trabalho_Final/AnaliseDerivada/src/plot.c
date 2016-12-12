#include "../include/plot.h"

// Plotar a solucao para a funcao e sua derivada 
void writeDataFile (double *d, double *X, Func* f, int np, int typeElem)
{
  FILE *file_solution, *file_derivarive;
  char filename[50];
  double h;
  int i;
  h = X[1] - X[0];

  // Plotar a solucao da funcao primeiro
  sprintf(filename,"solution.dat");
  file_solution = fopen(filename,"w+");
  // || Coluna 1 = Xi || Coluna 2 = Solucao analitica || Coluna 3 = Solucao aproximada ||
  for (i = 0; i < np; i++)
      fprintf(file_solution,"%e %e %e\n",X[i],f[1](0,0.0,X[i]),d[i]);
  fclose(file_solution);

  // PLotar a solucao da derivada
  sprintf(filename,"derivative.dat");
  file_derivarive = fopen("derivative.dat","w+");
  // Elemento linear nao tem a derivada jah calculada
  if (typeElem == 1)
  {
    double du[np-1];
    for (i = 0; i < np-1; i++)
    {
      du[i] = (d[i+1]-d[i])/h;
      fprintf(file_derivarive,"%e %e %e\n",X[i],f[2](0,0.0,X[i]),du[i]);
    }  
  }
  // Elemento de Hermite jah possui derivada calculada
  else
  {
    for (i = 0; i < np; i++)
      //fprintf(file_derivarive,"%e %e %e\n",X[i],f[2](0,0.0,X[i]),d[i+np]);
      fprintf(file_derivarive,"%e %e %e\n",X[i],f[2](0,0.0,X[i]),d[i+np]*(1/h));  // ASSIM DA CERTO !!! Fator de escala !
  }
  fclose(file_derivarive);
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
