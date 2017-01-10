#include <cstdio>
#include <cstdlib>
#include "../include/meshGenerator.h"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc-1 < 4)
    {
        printf("Usage:> %s <Xmax> <nElem> <type> <file_name>\n",argv[0]);
        printf("<Xmax> = Tamanho de uma fibra\n");
        printf("<nElem> = Numero de elementos em uma fibra\n");
        printf("<type> = Tipo de malha\n");
        printf("\t1 = Fibra simples\n");
        printf("\t2 = Fibra com uma bifurcacao\n");
        printf("\t3 = Fibra com uma trifircacao\n");
        printf("\t4 = FIbra gerada pelo L-System\n");
        printf("<file_name> = Nome do arquivo de saida\n");
        exit(-1);
    }
    else
    {
        Mesh *mesh = newMesh(argc,argv);
        //printMeshInfo(mesh);
        writeMeshToFile(mesh,argv[4]);
    }
}
