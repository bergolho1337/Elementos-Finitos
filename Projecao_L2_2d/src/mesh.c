#include "mesh.h"

Mesh* buildMesh ()
{
  printf("[!] Lendo malha ...");
  fflush(stdout);

  int i, j, ret;
  int np, nv, ne;
  Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
  // Ler os pontos
  ret = scanf("%d",&np);
  mesh->np = np;
  mesh->P = (double*)calloc(np*2,sizeof(double));
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < 2; j++)
      ret = scanf("%lf",&mesh->P[i*2+j]);
  }
  // Ler os elementos
  ret = scanf("%d",&ne);
  mesh->nv = nv = ELEM;
  mesh->ne = ne;
  mesh->E = (int*)calloc(ne*nv,sizeof(int));
  for (i = 0; i < ne; i++)
  {
    for (j = 0; j < nv; j++)
      ret = scanf("%d",&mesh->E[i*nv+j]);
  }

  printf(" ok\n");
  return mesh;
}

// Calcula a area de um elemento com 'nv' vertices
double calcArea (int *E, double *P, int k, int nv)
{
  switch (nv)
  {
    // Elemento triangulo
    case 3: {
              double a, b, c, p;
              int v1, v2, v3;
              // Indices dos pontos do triangulo
              v1 = E[k*nv];
              v2 = E[k*nv+1];
              v3 = E[k*nv+2];
              a = sqrt( pow(P[v1*2]-P[v2*2],2) + pow(P[v1*2+1]-P[v2*2+1],2) );
              b = sqrt( pow(P[v1*2]-P[v3*2],2) + pow(P[v1*2+1]-P[v3*2+1],2) );
              c = sqrt( pow(P[v2*2]-P[v3*2],2) + pow(P[v2*2+1]-P[v3*2+1],2) );
              p = (a+b+c)/2;                    // Semi-perimetro
              return sqrt(p*(p-a)*(p-b)*(p-c));
            }
  }
  return -1;
}

// Calcula as funcoes do termo forcante do sistema linear
void calcLoadFunction (Func function, double *f, int *E, double *P, int k, int nv)
{
  switch (nv)
  {
    // Elemento triangulo
    case 3: {
              int v1, v2, v3;
              // Indices dos pontos do triangulo
              v1 = E[k*nv];
              v2 = E[k*nv+1];
              v3 = E[k*nv+2];
              f[0] = function(0,P[v1*2],P[v1*2+1]);
              f[1] = function(0,P[v2*2],P[v2*2+1]);
              f[2] = function(0,P[v3*2],P[v3*2+1]);
            }
  }
}

void printMesh (Mesh *mesh)
{
  int i, j;
  int np, ne, nv;
  np = mesh->np;
  ne = mesh->ne;
  nv = mesh->nv;
  switch (ELEM)
  {
    case 3: {
              printf("\t[*] Usando elementos triangulares\n");
              break;
            }
  }
  printf("\t====== MATRIZ DOS PONTOS ==========\n");
  for (i = 0; i < np; i++)
      printf("\tPonto %d = (%lf,%lf)\n",i+1,mesh->P[i*2],mesh->P[i*2+1]);

  printf("\t====== MATRIZ DE CONECTIVIDADE ====\n");
  for (i = 0; i < ne; i++)
  {
    printf("\tElemento %d = [",i+1);
    for (j = 0; j < nv; j++)
      printf(" %d ",mesh->E[i*2+j]);
    printf("]\n");
  }
}
