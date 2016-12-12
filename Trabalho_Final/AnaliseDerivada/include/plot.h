#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/mefDerivada.h"

typedef double (*Func) (int k, double t, double x);            // Ponteiro para funcao

void writeDataFile (double *d, double *X, Func* f, int np, int typeElem);
void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);

#endif
