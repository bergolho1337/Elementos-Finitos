#include <stdio.h>
#include <stdlib.h>
#include "MWR.h"

int main (int argc, char *argv[])
{
  printf("============ METODO DOS RESIDUOS PONDERADOS =====================\n");
  printf("Problema = u'' + u = x,\t\t0 < x < 1\n");
  printf("u(0) = u(1) = 0\n");
  printf("Aproximando por: ~u = sum_{j=1}^N  (alfa_j * phi_j)\n");
  if (argc-1 < 3)
  {
    printf("------------------------------------------------------------------\n");
    printf("Usage:> %s <N> <Dx> <num_metodo>\n",argv[0]);
    printf("<Dx> = Tamanho da discretizacao para o plot\n");
    printf("<N> = Numero de funcoes de ponderacao a serem utilizadas [1,2,3].\n");
    printf("<num_metodo> = Numero do metodo escolhido para resolver o problema\n");
    printf("\t1 - Galerkin\n");
    printf("\t2 - Colocacao\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    MWR *mwr = newMWR(argc,argv);
    return 0;
  }
}
