/*
******************************************
    Autor: Lucas Berg
    Last Update: 23/11/16
******************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include "l2.h"

int main (int argc, char *argv[])
{
  printf("========= METODO DOS ELEMENTOS FINITOS ----- PROJECAO L2 em 2D ===========\n");
  printf("Problema: aproximar uma funcao f por sua projecao L2, P_hf.\n");
  printf("integral{(f - P_hf) phi dA = 0\n");
  printf("--------------------------------------------------------------------------\n");
  if (argc-1 < 0)
  {
    printf("------------------------------------------------------------------\n");
    printf("Usage:> %s < <arquivo_da_malha>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("[!] DEFINIR A FUNCAO f QUE SE QUER APROXIMAR NO ARQUIVO \"functions.h\" E \"functions.c\"\n");
    printf("[!] DEBUGACAO EH ATIVADA POR FLAG NO ARQUIVO \"l2.h\".\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    L2 *l2 = newL2(argc,argv);
    return 0;
  }
}
