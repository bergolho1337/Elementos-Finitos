#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>                                    // <---- M_PI is define here

typedef double (*Func) (int k, double x);            // Ponteiro para funcao

// Funcao a ser aproximada pelo metodo
double f (int k, double x);

// Polinomio simples do tipo: p(x) = x^k
double p_k (int k, double x);

// Polinomio de Lagrange
double l_k (int k, double x);

// Funcoes chapeu
double phi_k (int k, double x);

// Constroi pontos igualmente espacados no dominio
void setPoints (int N);

#endif
