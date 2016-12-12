#include "../include/plot.h"

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

void writeVTKFile (double *Vm, double *x, int *map, int np, int ne, int k)
{
  FILE *file;
  int i;
  char filename[50];
  sprintf(filename,"VTK/monoFEM%d.vtk",k);
  file = fopen(filename,"w+");
  fprintf(file,"# vtk DataFile Version 3.0\n");
  fprintf(file,"Monodomain FEM\n");
  fprintf(file,"ASCII\n");
  fprintf(file,"DATASET POLYDATA\n");
  fprintf(file,"POINTS %d float\n",np);
  for (i = 0; i < np; i++)
    fprintf(file,"%e %e %e\n",x[i],0.0,0.0);
  fprintf(file,"LINES %d %d\n",np,ne*3);
  for (i = 0; i < ne; i++)
  {
    fprintf(file,"2 %d %d\n",map[i*2],map[i*2+1]);
  }
    
  fprintf(file,"POINT_DATA %d\n",np);
  fprintf(file,"SCALARS vm float 1\n");
  fprintf(file,"LOOKUP_TABLE default\n");
  for (i = 0; i < np; i++)
    fprintf(file,"%e\n",Vm[i]);
  fclose(file);
}