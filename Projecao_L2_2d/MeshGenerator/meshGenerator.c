/*
******************************************
    Autor: Lucas Berg
    Last Update: 23/11/16
******************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void writePoints (int n, double h);
void writeElements (int n);

int main (int argc, char *argv[])
{
    if (argc-1 < 1)
    {
        printf("=========================== MESH GENERATOR ============================\n");
        printf("   Este programa cria uma malha regular de triangulos sobre o dominio\n");
        printf("   [0,1]x[0,1] a partir do numero de pontos que cada eixo possui.\n");
        printf("-----------------------------------------------------------------------\n");
        printf("Usage:> %s 'num_pontos' > 'arquivo_da_malha' \n",argv[0]);
        printf("-----------------------------------------------------------------------\n");
        printf("'num_pontos' = Numero de pontos no eixo x1 e x2\n");
        printf("'arquivo_da_malha' = Arquivo de saida contendo a malha.\n");
        printf("-----------------------------------------------------------------------\n");
        printf("Exemplo: %s 3 > malha3.in\n",argv[0]);
        printf("Cria a malha de referencia para o Exercicio 1\n");
        printf("=======================================================================\n");
        return -1;
    }
    else
    {
        int n;
        double h;
        
        // Captura a entrada
        n = atoi(argv[1]);
        h = 1.0 / (double)(n-1);

        writePoints(n,h);
        writeElements(n);

        return 0;
    }
}

// Escrever os pontos
void writePoints (int n, double h)
{
    int i, j;
    double x, y;
    printf("%d\n",n*n);
    for (i = 0; i < n; i++)
    {
        x = i*h;
        for (j = 0; j < n; j++)
        {
            y = j*h;
            printf("%e %e\n",x,y);
        }
    }
}

void writeElements (int n)
{
    int i, j;
    int nElem, elem;
    // Calcula o numero de elementos
    nElem = 2*(n-1)*(n-1);
    printf("%d\n",nElem);
    // Criar os elementos triangulos inferiores
    for (elem = 0; elem < n*(n-1); elem++)
    {
        if ((elem+1) % n != 0 || elem == 0)
            printf("%d %d %d\n",elem,elem+1,elem+n+1);       
    }
    // Criar os elementos triangulos superiores
    for (elem = 0; elem < n*(n-1); elem++)
    {
        if ((elem+1) % n != 0 || elem == 0)
            printf("%d %d %d\n",elem,elem+n,elem+n+1);       
    }
}