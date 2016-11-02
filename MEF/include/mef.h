#ifndef MEF_H_
#define MEF_H_

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

// Estrutura do resolvedor do MEF
struct MEF
{
  int N;                                      // Numero de funcoes de ponderacao a serem utilizadas
  double dx;                                  // Tamanho da discretizacao para o plot
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double a, b, c;                             // Coeficientes da equacao do problema
  int *mapping;                               // Mapeamento dos indices locais para os globais
  double *K;                                  // Matriz local de cada elemento do metodo
  double *A;                                  // Matriz global do sistema linear do metodo
  double *F;                                  // Vetor de termos independentes do metodo
  double *alfa;                               // Vetor dos coeficientes que aproximam a solucao
  int numeracao;                              // Modo como os nos serao numerados. || 1 = Automatica || 2 = Manual ||
}typedef MEF;

// Construtor
MEF* newMEF (int argc, char *argv[]);

// Inicializa as referencias para as funcoes utilizadas no metodo
Func* getFunctions ();

// Inicializa os coeficientes da equacao
void setCoeficients (MEF *mef);

// Inicializa o mapeamento local para o global
void setMapping (MEF *mef);

// Constroi a matriz local de cada elemento
void makeElementMatrix (MEF *mef);

// Constroi a matriz global a partir do mapeamento
void makeGlobalMatrix (MEF *mef);

// Constroi o vetor de termos independentes
void makeVector (MEF *mef);

// Inicializa as condicoes de contorno do problema
void setBoundaryConditions (MEF *mef);

#endif
