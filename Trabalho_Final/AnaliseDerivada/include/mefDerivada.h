#ifndef MEFDERIVADA_H_
#define MEFDERIVADA_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/plot.h"
#include "../include/functions.h"
#include "../include/linSolver.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int k, double t, double x);            

// Estrutura do resolvedor da equacao pelo MEF
struct MEFDerivada
{
  int typeElem;                               // Tipo do elemento
  int nElem;                                  // Numero de elementos a serem utilizados
  int nPoints;                                // Numero de pontos do dominio
  int *map;                                   // Mapeamento dos elementos locais para os globais
  double dx;                                  // Tamanho da discretizacao no espaco (h)
  double x_max;                               // Tamanho maximo do dominio
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double *X;                                  // Vetor com os pontos do dominio
  double *K;                                  // Matriz global do sistema linear
  double *A;                                  // Matriz massa global
  double *B;                                  // Matriz rigidez global
  double *F;                                  // Vetor de carga global
  double *d;                                  // Vetor dos coeficientes que aproximam a funcao u(x)
}typedef MEFDerivada;

MEFDerivada* newMEFDerivada (int argc, char *argv[]);
Func* buildFunctions ();
int* buildMap (int nElem, int nPoints, int typeElem);
double* buildPoints (int nPoints, double h);
void allocMemoryVectors (MEFDerivada *mef);
void printInfoModel (MEFDerivada *mef);
void assembleMatrix (MEFDerivada *mef);
double* buildLocalMassMatrix (double h, int typeElem);
double* buildGlobalMatrixFromLocal (double *local_M, int *map, int np, int ne, int typeElem);
void setBoundaryConditions_Matrix (double *K, int np, int typeElem);
void assembleLoadVector (MEFDerivada *mef);
void setBoundaryConditions_Vector (double *F, int np, int typeElem);
/*
double* setInitialConditions (Func *func, double *X, int np);
void solveTransientProblem (Calor *calor);
void printInfoModel (Calor *calor);
*/

#endif
