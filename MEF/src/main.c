#include <stdio.h>
#include <stdlib.h>
#include "mef.h"

int main (int argc, char *argv[])
{
  printf("========= METODO DOS ELEMENTOS FINITOS ----- UNIDIMENSIONAL 1D ===========\n");
  printf("Problema: Encontrar uma funcao u de tal maneira que satisfaca a equacao.\n");
  printf("-au'' + bu' + cu = f\n");
  if (argc-1 < 3)
  {
    printf("------------------------------------------------------------------\n");
    printf("Usage:> %s <N> <Dx> <numeracao>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("<Dx> = Tamanho da discretizacao para o plot\n");
    printf("<N> = Numero de funcoes de ponderacao a serem utilizadas.\n");
    printf("<numeracao> = Modo de numeracao da malha\n");
    printf("\t1 - Automatica\n");
    printf("\t2 - Manual\n");
    printf("[!] DEFINIR A FUNCAO f QUE SE QUER APROXIMAR NO ARQUIVO \"functions.c\"\n");
    printf("[!] DEBUGACAO EH ATIVADA POR FLAG NO ARQUIVO \"mef.h\".\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    MEF *mef = newMEF(argc,argv);
    return 0;
  }
}
