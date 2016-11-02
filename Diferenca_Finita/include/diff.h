#ifndef DIFF_H_
#define DIFF_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"
#include "linSolver.h"
#include "plot.h"
#include "functions.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
//#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int k, double x);            // Ponteiro para funcao (o 'k' equivale ao grau do polinomio)

// Estrutura do resolvedor de Diferencas Finitas
struct Diff
{
  int N;                                      // Numero de subintervalos da discretizacao
  double h;                                   // Tamanho da discretizacao para o plot
  double Pe;                                  // Valor de Peclet
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double a, b, c;                             // Coeficientes da equacao do problema
  double *A;                                  // Matriz do sistema linear do metodo
  double *F;                                  // Vetor de termos independentes do metodo
  double *u;                                  // Vetor contendo a solucao dos pontos
}typedef Diff;

// Construtor
Diff* newDiff (int argc, char *argv[]);

// Inicializa as referencias para as funcoes utilizadas no metodo
Func* getFunctions ();

// Inicializa os coeficientes da equacao
void setCoeficients (Diff *diff);

// Constroi a matriz global a partir do mapeamento
void makeMatrix (Diff *diff);

// Constroi o vetor de termos independentes
void makeVector (Diff *diff);

// Inicializa as condicoes de contorno do problema
void setBoundaryConditions (Diff *diff);

#endif
