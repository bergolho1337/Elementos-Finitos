#ifndef L2_H_
#define L2_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"
#include "linSolver.h"
#include "plot.h"
#include "functions.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
//#define ANALYZE 1                                     // Flag para realizar a analise do numero de condicionamento (kapla)
//#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int k, double x);            // Ponteiro para funcao (o 'k' equivale ao grau do polinomio)

// Estrutura do resolvedor de Projecao L2
struct L2
{
  int N;                                      // Numero de funcoes de ponderacao a serem utilizadas
  double dx;                                  // Tamanho da discretizacao para o plot
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double *A;                                  // Matriz do sistema linear do metodo
  double *b;                                  // Vetor de termos independentes do metodo
  double *alfa;                               // Vetor dos coeficientes que aproximam a funcao 'P_hf'
  int num_funcao;                             // Numero das funcoes base escolhidas
}typedef L2;

// Construtor
L2* newL2 (int argc, char *argv[]);

// Inicializa as referencias para as funcoes utilizadas no metodo
Func* getFunctions ();

// Constroi a matriz do sistema linear
void makeMatrix_Polinomial (L2 *l2);

// Constroi o vetor de termos independentes
void makeVector_Polinomial (L2 *l2);

// Constroi a matriz do sistema linear usando polinomios de Lagrange
void makeMatrix_Lagrange (L2 *l2);

// Constroi o vetor de termos independentes usando polinomios de Lagrange
void makeVector_Lagrange (L2 *l2);

// Funcao que calcula o valor do numero de condicionamento da matriz A
double analyzeKapla (L2 *l2);

#endif
