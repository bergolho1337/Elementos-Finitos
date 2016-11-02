#ifndef INTEGRATION_H_
#define INTEGRATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.h"

#define H 1.0e-06                                // Tamanho dos trapezios usados na integracao

typedef double (*Func) (int k, double x);        // Ponteiro para funcao

/* ======================================================================================================== */
// Resolve uma integral no intervalo [a,b] usando a Regra do Trapezio
double trapezoidalRule (int i, int j, Func f1, Func f2, double a, double b);

// Calcula o produto entre as funcoes f1 e f2
double F (int i, int j, Func f1, Func f2, double x);

#endif
