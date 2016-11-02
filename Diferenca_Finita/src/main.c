#include <stdio.h>
#include <stdlib.h>
#include "diff.h"

int main (int argc, char *argv[])
{
  printf("========= DIFERENCAS FINITAS ----- EXERCICIO DO PECLET ===========\n");
  printf("Problema: Encontrar uma funcao u de tal maneira que satisfaca a equacao.\n");
  printf("-au'' + bu' = 0\n");
  if (argc-1 < 1)
  {
    printf("------------------------------------------------------------------\n");
    printf("Usage:> %s <h>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("<h> = Tamanho da discretizacao da diferenca finita.\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    Diff *diff = newDiff(argc,argv);
    return 0;
  }
}
