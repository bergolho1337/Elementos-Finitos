#ifndef L2_H_
#define L2_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/plot.h"
#include "../include/functions.h"
#include "../include/linSolver.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
//#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int k, double t, double x);            

// Estrutura do resolvedor da equacao do Calor
struct Calor
{
  int nElem;                                  // Numero de elementos a serem utilizados
  int nPoints;                                // Numero de pontos do dominio
  int *map;                                   // Mapeamento dos elementos locais para os globais
  double dx;                                  // Tamanho da discretizacao no espaco (h)
  double dt;                                  // Tamanho da discretizacao no tempo (k)
  double t_max;                               // Tempo maximo de simulacao
  double x_max;                               // Tamanho maximo do dominio
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double *X;                                  // Vetor com os pontos do dominio
  double *K;                                  // Matriz global do sistema linear
  double *A;                                  // Matriz massa global
  double *B;                                  // Matriz rigidez global
  double *F;                                  // Vetor de carga global
  double *dNew;                               // Vetor dos coeficientes que aproximam a funcao u(x,t) - Instante n
  double *dOld;                               // Vetor dos coeficientes que aproximam a funcao u(x,t) - Instante n-1
}typedef Calor;

Calor* newCalor (int argc, char *argv[]);
Func* buildFunctions ();
int* buildMap (int nElem);
double* buildPoints (int np, double h);
void assembleMatrix (Calor *calor);
double* buildLocalMassMatrix (double h);
double* buildLocalStiffMatrix (double h, double k);
double* buildGlobalMatrixFromLocal (double *local_A, int *map, int np, int ne);
double* buildGlobalMatrix (double *A, double *B, double dt, int np);
void mapLocalToGlobal_Matrix (double *A, double *B, double *K, int *map, double k, int np, int elem);
void assembleLoadVector (Calor *calor);
void buildLocalVector (double *A, double *dOld, double *V, int *map, int elem);
void setBoundaryConditions_Matrix (Calor *calor);
void setBoundaryConditions_Vector (Calor *calor);
double* setInitialConditions (Func *func, double *X, int np);
void solveTransientProblem (Calor *calor);
void printInfoModel (Calor *calor);

#endif
