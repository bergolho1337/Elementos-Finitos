#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>                                              // <---- M_PI is define here

typedef double (*Func) (int k, double t, double x);            // Ponteiro para funcao

// Funcao do vetor de carga
double f (int k, double t, double x);

// Funcao da condicao inicial
double u0 (int k, double t, double x);

// Funcao da solucao analitica do problema
double analit (int k, double t, double x);

#endif
