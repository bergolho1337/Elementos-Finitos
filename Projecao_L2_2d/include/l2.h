#ifndef L2_H_
#define L2_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mesh.h"
#include "plot.h"
#include "functions.h"

// >>>>>>>>>>> DESCOMENTAR ESSES defines PARA ATIVAR AS FLAGS <<<<<<<<<<<
#define DEBUG 1                                       // Flag para debugacao e imprimir informacoes na tela (matrizes e vetores)
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

typedef double (*Func) (int k, double x1, double x2);            // Ponteiro para funcao (o 'k' equivale ao grau do polinomio)

// Estrutura do resolvedor de Projecao L2
struct L2
{
  Mesh *mesh;                                 // Estrutura da malha
  int N;                                      // Numero de funcoes de ponderacao a serem utilizadas
  double dx;                                  // Tamanho da discretizacao para o plot
  Func functions;                            // Vetor de ponteiros para as funcoes envolvidas no metodo
  double *A;                                  // Matriz global do sistema linear
  double *K;                                  // Matriz local do elemento
  double *F;                                  // Vetor de termos independentes do metodo
  double *alfa;                               // Vetor dos coeficientes que aproximam a funcao 'P_hf'
  int num_funcao;                             // Numero das funcoes base escolhidas
}typedef L2;

L2* newL2 (int argc, char *argv[]);
void assembleMassMatrix (L2 *l2);
double* buildLocalMatrix (int nv);
void checkMap (int *map, int *E, int k, int nv);
void mapLocalToGlobal_Matrix (double *A, double *K, int *map, double area, int nv, int np);
void assemblyLoadVector (L2 *l2);
void mapLocalToGlobal_Vector (double *F, double *f, int *map, double area, int nv, int np);

#endif
