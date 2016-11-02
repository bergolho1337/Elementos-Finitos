#ifndef MWR_H_
#define MWR_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integration.h"
#include "linSolver.h"
#include "plot.h"
#include "functions.h"

typedef double (*Func) (double x);            // Ponteiro para funcao

// Estrutura do resolvedor de Residuos Ponderados
struct MWR
{
  int N;                                      // Numero de funcoes de ponderacao a serem utilizadas
  double dx;                                  // Tamanho da discretizacao para o plot
  Func *functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double *A;                                  // Matriz do sistema linear do metodo
  double *b;                                  // Vetor de termos independentes do metodo
  double *alfa;                               // Vetor dos coeficientes que aproximam a funcao 'u'
  int num_metodo;
  double *epsilon;
}typedef MWR;

// Construtor
MWR* newMWR (int argc, char *argv[]);

// Inicializa as referencias para as funcoes de ponderacao do metodo
Func* getFunctions ();

// Constroi a matriz do sistema linear para o metodo de Galerkin
void makeMatrix_Galerkin (MWR *mwr);

// Constroi o vetor de termos independentes para o metodo de Galerkin
void makeVector_Galerkin (MWR *mwr);

// Leitura dos valores dos pontos de coloacao epsilon
void readEpsilon (MWR *mwr);

// Constroi a matriz do sistema linear para o metodo da Colocacao
void makeMatrix_Colocacao (MWR *mwr);

// Constroi o vetor de termos independentes para o metodo de Colocacao
void makeVector_Colocacao (MWR *mwr);

#endif
