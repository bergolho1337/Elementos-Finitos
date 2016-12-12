#include "l2.h"

L2* newL2 (int argc, char *argv[])
{
  L2 *l2 = (L2*)malloc(sizeof(L2));
  l2->mesh = buildMesh();
  l2->functions = f;
  #ifdef DEBUG
  printMesh(l2->mesh);
  #endif

  // Matriz com a solucao analitica
  assembleMassMatrix(l2);
  assemblyLoadVector(l2);
  l2->alfa = solveLinearSystem_CG(l2->A,l2->F,l2->mesh->np);
  writeDataFile(l2->functions,l2->alfa,l2->mesh->P,l2->mesh->np);

  #ifdef DEBUG
  printVector("Solucao do sistema",l2->alfa,l2->mesh->np);
  #endif

  return l2;
}

// Constroi a matriz usando a solucao analitica da integral do elemento
void assembleMassMatrix (L2 *l2)
{
  int k, np, ne, nv;
  int *E;
  int map[l2->mesh->nv];
  double *P;
  double area;
  ne = l2->mesh->ne;
  nv = l2->mesh->nv;
  np = l2->mesh->np;
  P = l2->mesh->P;
  E = l2->mesh->E;
  l2->K = buildLocalMatrix(nv);
  l2->A = (double*)calloc(np*np,sizeof(double));
  // Mapear cada elemento na matriz global
  for (k = 0; k < ne; k++)
  {
    // Calcula a area do elemento
    area = calcArea(E,P,k,nv);
    // Descobre o mapeamento do elemento na matriz global
    checkMap(map,E,k,nv);
    // Mapeia a matriz local para a global
    mapLocalToGlobal_Matrix(l2->A,l2->K,map,area,nv,np);
  }
  #ifdef DEBUG
  printMatrix("Matriz global",l2->A,np);
  #endif
}

// Constroi a matriz local usando a solucao analitica da integral da matriz de rigidez
/*
      | 2 1 1 |
M^K = | 1 2 1 | = \int_{K} [ phi_i phi_j dx ]
      | 1 1 2 |
*/
double* buildLocalMatrix (int nv)
{
  int i, j;
  double *K = (double*)calloc(nv*nv,sizeof(double));
  switch (nv)
  {
    // Triangulo
    case 3: {
              for (i = 0; i < 3; i++)
              {
                for (j = 0; j < 3; j++)
                {
                  if (i == j)
                    K[i*3+j] = 2.0;
                  else
                    K[i*3+j] = 1.0;
                }
              }
              break;
            }
  }
  return K;
}

// Preenche o vetor map com o mapeamento dos nos locais para globais [1,2,3] -> [r,s,t]
void checkMap (int *map, int *E, int k, int nv)
{
  int i;
  for (i = 0; i < nv; i++)
    map[i] = E[k*nv+i];
}

// Mapeia a matriz local do elemento na matriz global do sistema usando a solucao analitica da integral
void mapLocalToGlobal_Matrix (double *A, double *K, int *map, double area, int nv, int np)
{
  int i, j;
  int ig, jg;
  switch (nv)
  {
    case 3:{
              for (i = 0; i < 3; i++)
              {
                ig = map[i];
                for (j = 0; j < 3; j++)
                {
                    jg = map[j];
                    A[ig*np+jg] += (K[i*3+j]*area) / 12.0;
                }
              }
           }
  }
}

void assemblyLoadVector (L2 *l2)
{
  int k, np, ne, nv;
  int *E;
  int map[l2->mesh->nv];
  double *P;
  double area;
  double f[l2->mesh->nv];
  ne = l2->mesh->ne;
  nv = l2->mesh->nv;
  np = l2->mesh->np;
  P = l2->mesh->P;
  E = l2->mesh->E;
  l2->F = (double*)calloc(np,sizeof(double));
  // Mapear cada elemento no vetor global
  for (k = 0; k < ne; k++)
  {
    // Calcula a area do elemento
    area = calcArea(E,P,k,nv);
    // Calcula o valor da funcao 'f'
    calcLoadFunction(l2->functions,f,E,P,k,nv);
    // Descobre o mapeamento do elemento no vetor global
    checkMap(map,E,k,nv);
    // Mapeia o vetor local para a global
    mapLocalToGlobal_Vector(l2->F,f,map,area,nv,np);
  }
  #ifdef DEBUG
  printVector("Vetor global",l2->F,np);
  #endif
}

//
void mapLocalToGlobal_Vector (double *F, double *f, int *map, double area, int nv, int np)
{
  int i;
  int ig;
  switch (nv)
  {
    case 3:{
              for (i = 0; i < 3; i++)
              {
                ig = map[i];
                F[ig] += (f[i]*area) / 3.0;
              }
           }
  }
}
