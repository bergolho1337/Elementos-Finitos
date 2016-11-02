#include "plot.h"

double Faprox (double x, Func *f, double *alfa, int N)
{
  int k;
  double eval = 0;
  for (k = 0; k < N; k++)
    eval += alfa[k]*f[k](x);
  return eval;
}

double Fanalit (double x)
{
  return x - (sinh(x)/1.17520119364);
}

void writeDataFile (Func *f, double *alfa, double dx, int N)
{
  FILE *file_aprox, *file_analit, *file_residue;
  int i, n;
  double x;

  file_aprox = fopen("f_aprox.dat","w+");
  file_analit = fopen("f.dat","w+");
  file_residue = fopen("residue.dat","w+");

  // Escreve a aproximacao
  n = 1 / dx;
  for (i = 0; i < n; i++)
  {
    x = i*dx;
    fprintf(file_aprox,"%e %e\n",x,Faprox(x,f,alfa,N));
    fprintf(file_analit,"%e %e\n",x,Fanalit(x));
    fprintf(file_residue,"%e %e\n",x,fabs(Faprox(x,f,alfa,N)-Fanalit(x)));
  }
  fclose(file_aprox);
  fclose(file_analit);
  fclose(file_residue);
}
