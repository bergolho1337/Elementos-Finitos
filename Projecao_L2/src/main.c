#include <stdio.h>
#include <stdlib.h>
#include "l2.h"

int main (int argc, char *argv[])
{
  printf("========= METODO DOS ELEMENTOS FINITOS ----- PROJECAO L2 ===========\n");
  printf("Problema: aproximar uma funcao f por sua projecao L2, P_hf.\n");
  printf("integral{(f - P_hf) phi dx = 0\n");
  if (argc-1 < 3)
  {
    printf("------------------------------------------------------------------\n");
    printf("Usage:> %s <N> <Dx> <num_funcao>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("<Dx> = Tamanho da discretizacao para o plot\n");
    printf("<N> = Numero de funcoes de ponderacao a serem utilizadas.\n");
    printf("<num_funcao> = Numero da familia de funcoes bases escolhida para resolver o problema\n");
    printf("\t1 - Polinomiais simples\n");
    printf("\t2 - Lagrange\n");
    printf("\t3 - Funcoes chapeu\n");
    printf("[!] DEFINIR A FUNCAO f QUE SE QUER APROXIMAR NO ARQUIVO \"functions.h\" E \"functions.c\"\n");
    printf("[!] DEBUGACAO E ANALISE DO NUMERO DE CONDICIONAMENTO SAO ATIVADAS POR FLAGS NO ARQUIVO \"l2.h\".\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    L2 *l2 = newL2(argc,argv);
    return 0;
  }
}
