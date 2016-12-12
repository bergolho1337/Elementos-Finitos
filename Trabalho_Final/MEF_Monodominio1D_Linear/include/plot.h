#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (int elem, double t, double vm, double w);

void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);
void writeVTKFile (double *Vm, double *x, int *map, int np, int ne, int k);

#endif