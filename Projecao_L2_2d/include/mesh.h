#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ELEM 3          // Tipo de elemento a ser usado (3 = Triangulo)

typedef double (*Func) (int k, double x1, double x2);

// Estrutura de malha considerando triangulos
struct Mesh
{
  int np;               // Numero de pontos da malha
  int ne;               // Numero de elementos da malha
  int nv;               // Numero de vertices do elemento
  double *P;            // Matriz dos pontos da malha
  int *E;               // Matriz de conectividade dos elementos
}typedef Mesh;

Mesh* buildMesh ();
void printMesh (Mesh *mesh);
double calcArea (int *E, double *P, int k, int nv);
void calcLoadFunction (Func function, double *f, int *E, double *P, int k, int nv);
