#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (int k, double t, double x);            // Ponteiro para funcao

void writeDataFile (double *d, double *X, Func* f, double t, int np, int iter);
void writeDataFile2 (double *d, double *X, Func *f, double t, int np, int iter);
void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);

#endif
