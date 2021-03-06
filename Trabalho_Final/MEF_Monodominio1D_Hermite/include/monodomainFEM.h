#ifndef MONODOMAINFEM_H_
#define MONODOMAINFEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/fitz.h"
#include "../include/plot.h"
#include "../include/linSolver.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
//#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int elem, double t, double vm, double w);      

/* ============================== CONSTANTES ================================================== */
const double BETA = 0.14;                     // Razao area superficial por volume (cm^-1)
const double Cm = 1.0;                        // Capacitancia da membrana (uF/cm^2)
const double SIGMA = 4000.0;                  // Condutividade citoplasmatica da celula (mS/cm)
/* ============================================================================================ */

// Estrutura do resolvedor da equacao do Monodominio
struct MonodomainFEM
{
  int nElem;                                  // Numero de elementos a serem utilizados
  int nPoints;                                // Numero de pontos do dominio
  int M;                                      // Numero de subintervalos no tempo
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
  double *VNew;                               // Vetor com o valor do potencial transmembranico de cada ponto no tempo n
  double *VOld;                               // Vetor com o valor do potencial transmembranico de cada ponto no tempo n-1
  double *Vstar;                              // Vetor com o valor do potencial transmembranico de cada ponto no tempo intermediario *
  double *wNew;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
  double *wOld;                               // Vetor com o valor da variavel de estado de cada ponto no tempo n
}typedef MonodomainFEM;

/* ================================= FUNCTIONS ======================================================= */
MonodomainFEM* newMonodomainFEM (int argc, char *argv[]);
void freeMonodomain (MonodomainFEM *monoFEM);
int* buildMap (int nElem);
double* buildPoints (int np, double h);
Func* buildFunctions ();
void printInfoModel (MonodomainFEM *monoFEM);
void assembleMatrix (MonodomainFEM *monoFEM);
double* buildLocalMassMatrix (double h);
double* buildLocalStiffMatrix (double h);
double* buildGlobalMatrixFromLocal (double *local_A, int *map, int np, int ne);
double* buildGlobalMatrix (double *A, double *B, double dt, int np);

void solveMonodomain (MonodomainFEM *monoFEM);
void assembleLoadVector (MonodomainFEM *monoFEM);
void solveEDO (MonodomainFEM *monoFEM, double t);

/* =================================================================================================== */

#endif