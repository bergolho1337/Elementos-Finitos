#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double (*Func) (double x);            // Ponteiro para funcao

void writeDataFile (Func *f, double *alfa, double dx, int N);
double Faprox (double x, Func *f, double *alfa, int N);
double Fanalit (double x);

#endif
