#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>                                    // <---- M_PI is define here

typedef double (*Func) (int k, double x1, double x2);            // Ponteiro para funcao

// Funcao a ser aproximada pelo metodo
double f (int k, double x1, double x2);

#endif
