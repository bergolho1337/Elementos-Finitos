#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (int k, double x1, double x2);            // Ponteiro para funcao

void writeDataFile (Func f, double *alfa, double *P, int np);
void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);

#endif
