#include "../include/functions.h"

// >>>>>>>>>>>>>>>>>>>>> FUNCAO VETOR FORCA <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
double f (int k, double t, double x)
{
  return 0;                     
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// >>>>>>>>>>>>>>> FUNCAO CONDICAO INICIAL  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
double u0 (int k, double t, double x)
{
  return sin(2*M_PI*x);                     
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

// >>>>>>>>>>>>>>> FUNCAO CONDICAO INICIAL  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
double analit (int k, double t, double x)
{
  return exp(-4*pow(M_PI,2)*t)*sin(2*M_PI*x);
}
// >>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
