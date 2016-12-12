#ifndef LINSOLVER_H_
#define LINSOLVER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "plot.h"

#define ITER_MAX 1000					// Numero maximo de iteracoes
#define EPSILON 1.0e-08				// Tolerancia para a resposta do sistema linear

double dotProduct (double *a, double *b, int n);
void saxpy (double *x, double *y, double alpha, int n);
void saxpy2 (double *x, double *y, double alpha, int n);
void matvec (double *A, double *b, double *x, int n);
double multiply_row (int row_size, int *Aj, double *Av, double *b);
double norm (double *v, int n);

// Resolve um sistema linear usando o Metodo Gradiente Conjugado
double* solveLinearSystem_CG (double *A, double *b, int N);

// Checar se o sistema foi resolvido corretamente.
double checkSystem (double *A, double *b, double *x, int N);

// Funcao que inverte uma matriz A
double* invertMatrix (double *A, int n);

// Testa se duas matrizes sao inversas uma da outra
void testInverse (double *A, double *Ainv, int N);

// Calcula a norma L_2 de uma matriz (raiz quadrada da soma dos quadrados dos elementos)
double normL21 (double *A, int N);

#endif
