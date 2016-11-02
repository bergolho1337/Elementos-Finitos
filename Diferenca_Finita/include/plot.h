#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (int k, double x);            // Ponteiro para funcao

void writeDataFile (double *u, double a, double b, double c, double h, int N);
double Fanalit (double x, double a, double b, double c);
double Faprox (double x, Func *f, double *alfa, int N);
void printMatrix (char *str, double *A, int N);
void printVector (char *str, double *b, int N);

#endif
