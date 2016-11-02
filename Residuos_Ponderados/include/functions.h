#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>                             // <---- M_PI is define here

typedef double (*Func) (double x);            // Ponteiro para funcao

// Funcao forcante (lado direito do sistema linear)
double f (double x);

// Funcoes de ponderacao
double phi_1 (double x);
double phi_2 (double x);
double phi_3 (double x);

// Funcoes do operador linear
double L_1 (double x);
double L_2 (double x);
double L_3 (double x);

#endif
