#ifndef LINSOLVER_H_
#define LINSOLVER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Resolve o sistema linear associado pela Decomposicao LU
double* solveLinearSystem_LU (double *A, double *b, int N);

// Checar se o sistema foi resolvido corretamente. (deletar matriz A antes e construir ela de novo --> A = LU)
double checkSystem (double *A, double *b, double *x, int N);

#endif
