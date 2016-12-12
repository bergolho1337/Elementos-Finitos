/*
******************************************
    Autor: Lucas Berg
    Last Update: 23/11/16
******************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include "../include/calor.h"

int main (int argc, char *argv[])
{
  printf("========= METODO DOS ELEMENTOS FINITOS ----- EQUACAO DO CALOR TRANSIENTE ===========\n");
  printf("Problema: resolver a equacao do calor.\n");
  printf("u_t - u_xx = 0\n");
  printf("u(0,t) = u(1,t) = 0\n");
  printf("u(x,0) = sin(2.pi.x)\n");
  printf("--------------------------------------------------------------------\n");
  if (argc-1 < 4)
  {
    printf("Usage:> %s <nElem> <dt> <t_max> <x_max>\n",argv[0]);
    printf("------------------------------------------------------------------\n");
    printf("<nElem> = Numero de elementos\n");
    printf("<dt> = Tamanho da discretizacao no tempo\n");
    printf("<t_max> = Tempo maximo de simulacao\n");
    printf("<x_max> = Tamanho maximo do dominio\n");
    printf("[!] DEFINIR A FUNCAO f QUE SE QUER APROXIMAR NO ARQUIVO \"functions.c\"\n");
    printf("[!] DEBUGACAO EH ATIVADA POR FLAG NO ARQUIVO \"calor.h\".\n");
    printf("------------------------------------------------------------------\n");
    printf("Exemplo: ./calor 5 0.01 0.1 1.0\n");
    printf("==================================================================\n");
    return 1;
  }
  else
  {
    Calor *calor = newCalor(argc,argv);
    solveTransientProblem(calor);
    return 0;
  }
}
