/*
******************************************
    Autor: Lucas Berg
    Last Update: 09/12/16
******************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include "../include/mefDerivada.h"

int main (int argc, char *argv[])
{
  printf("========= METODO DOS ELEMENTOS FINITOS ----- ANALISE DE DERIVADA ===========\n");
  printf("Problema: resolver a equacao.\n");
  printf("-u_xx = f\n");
  printf("u(0,t) = u(1,t) = 0\n");
  printf("--------------------------------------------------------------------\n");
  if (argc-1 < 3)
  {
    printf("Usage:> %s <nElem> <x_max> <tipoElem>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("<nElem> = Numero de elementos\n");
    printf("<x_max> = Tamanho maximo do dominio\n");
    printf("<tipoElem> = Tipo da funcao do elemento finito\n");
    printf("\t1 = Linear\n");
    printf("\t2 = Quadratico\n");
    printf("\t3 = Hermite\n");
    printf("[!] DEFINIR A FUNCAO f QUE SE QUER APROXIMAR NO ARQUIVO \"functions.c\"\n");
    printf("[!] DEBUGACAO EH ATIVADA POR FLAG NO ARQUIVO \"mefDerivada.h\".\n");
    printf("------------------------------------------------------------------\n");
    printf("Exemplo: ./mefDerivada 5 1.0 1 \n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    MEFDerivada *mefDerivada = newMEFDerivada(argc,argv);
    return 0;
  }
}
