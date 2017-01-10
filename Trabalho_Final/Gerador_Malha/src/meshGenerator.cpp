#include "../include/meshGenerator.h"

Mesh* newMesh (int argc, char *argv[])
{
    Mesh *mesh = NULL;
    int type = atoi(argv[3]);
    switch (type)
    {
        case 1: {
                    mesh = buildMesh_1(argc,argv);
                    break;
                }
        case 2: {
                    mesh = buildMesh_2(argc,argv);
                    break; 
                }
        case 3: {
                    mesh = buildMesh_3(argc,argv);
                    break; 
                }
        case 4: {
                    mesh = buildMesh_4(argc,argv);
                    break; 
                }
    }
    return mesh;
}

// Constroi uma malha como sendo uma fibra simples reta
Mesh* buildMesh_1 (int argc, char *argv[])
{
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    mesh->xMax = atof(argv[1]);
    mesh->nElem = atoi(argv[2]);
    mesh->nPoints = mesh->nElem + 1;
    mesh->h = mesh->xMax / (double)mesh->nElem;
    mesh->points = (Point*)malloc(sizeof(Point)*mesh->nPoints);
    for (int i = 0; i < mesh->nPoints; i++)
    {
        mesh->points[i].x = i*mesh->h;
        mesh->points[i].y = mesh->points[i].z = 0; 
    }
    mesh->elements = (Element*)malloc(sizeof(Element)*mesh->nElem);
    for (int i = 0; i < mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1; 
    }
    return mesh;
}

// Constroi uma malha contendo uma bifurcacao no final
Mesh* buildMesh_2 (int argc, char *argv[])
{
    int i, j;
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    mesh->xMax = atof(argv[1]);
    mesh->nElem = atoi(argv[2]);
    mesh->nPoints = mesh->nElem + 1;
    mesh->h = mesh->xMax / (double)mesh->nElem;
    // Construir os pontos
    mesh->points = (Point*)malloc(sizeof(Point)*(3*mesh->nElem+1));
    // Construir a primeira fibra reta
    for (i = 0; i < mesh->nPoints; i++)
    {
        mesh->points[i].x = i*mesh->h;
        mesh->points[i].y = mesh->points[i].z = 0; 
    }
    // Construir o primeiro ramo rotacionado de 45 em relacao a fibra anterior
    double d_rot[3];
    int id_bif = mesh->nPoints-1; 
    calculateRotation(d_rot,mesh->points[0],mesh->points[1],TETA);
    for (i = mesh->nPoints, j = 1; i < mesh->nPoints*2-1; i++, j++)
    {
       mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
       mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
       mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
    }
    // Construir o segundo ramo rotacionado de -45 em relacao a fibra anterior
    calculateRotation(d_rot,mesh->points[0],mesh->points[1],-TETA);
    for (i = mesh->nPoints*2-1, j = 1; i < mesh->nPoints*3-2; i++, j++)
    {
       mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
       mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
       mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
    }

    // Construir os elementos
    mesh->elements = (Element*)malloc(sizeof(Element)*(3*mesh->nElem));
    // Primeira fibra
    for (i = 0; i < mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }
    // Segunda fibra
    // Bifurcacao
    mesh->elements[mesh->nElem].left = id_bif;
    mesh->elements[mesh->nElem].right = i+1;
    for (i = mesh->nElem+1; i < 2*mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }
    // Terceira fibra
    // Bifurcacao
    mesh->elements[2*mesh->nElem].left = id_bif;
    mesh->elements[2*mesh->nElem].right = i+1;
    for (i = 2*mesh->nElem+1; i < 3*mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }

    // Numero total de elementos e de pontos
    mesh->nElem *= 3;
    mesh->nPoints = mesh->nElem + 1;
    return mesh;
}

// Constroi uma malha contendo uma trifurcacao no final
Mesh* buildMesh_3 (int argc, char *argv[])
{
    int i, j;
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    mesh->xMax = atof(argv[1]);
    mesh->nElem = atoi(argv[2]);
    mesh->nPoints = mesh->nElem + 1;
    mesh->h = mesh->xMax / (double)mesh->nElem;
    // Construir os pontos
    mesh->points = (Point*)malloc(sizeof(Point)*(4*mesh->nElem+1));
    // Construir a primeira fibra reta
    for (i = 0; i < mesh->nPoints; i++)
    {
        mesh->points[i].x = i*mesh->h;
        mesh->points[i].y = mesh->points[i].z = 0; 
    }
    // Construir o primeiro ramo rotacionado de 45 em relacao a fibra anterior
    double d_rot[3];
    int id_bif = mesh->nPoints-1; 
    calculateRotation(d_rot,mesh->points[0],mesh->points[1],TETA);
    for (i = mesh->nPoints, j = 1; i < mesh->nPoints*2-1; i++, j++)
    {
       mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
       mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
       mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
    }
    // Construir o segundo ramo rotacionado de -45 em relacao a fibra anterior
    calculateRotation(d_rot,mesh->points[0],mesh->points[1],-TETA);
    for (i = mesh->nPoints*2-1, j = 1; i < mesh->nPoints*3-2; i++, j++)
    {
       mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
       mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
       mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
    }
    // Construir o terceiro ramo rotacionado de 0 em relacao a fibra anterior
    calculateRotation(d_rot,mesh->points[0],mesh->points[1],0);
    for (i = mesh->nPoints*3-2, j = 1; i < mesh->nPoints*4-3; i++, j++)
    {
       mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
       mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
       mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
    }

    // Construir os elementos
    mesh->elements = (Element*)malloc(sizeof(Element)*(4*mesh->nElem));
    // Primeira fibra
    for (i = 0; i < mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }
    // Segunda fibra
    // Bifurcacao
    mesh->elements[mesh->nElem].left = id_bif;
    mesh->elements[mesh->nElem].right = i+1;
    for (i = mesh->nElem+1; i < 2*mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }
    // Terceira fibra
    // Bifurcacao
    mesh->elements[2*mesh->nElem].left = id_bif;
    mesh->elements[2*mesh->nElem].right = i+1;
    for (i = 2*mesh->nElem+1; i < 3*mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }

    // Quarta fibra
    // Bifurcacao
    mesh->elements[3*mesh->nElem].left = id_bif;
    mesh->elements[3*mesh->nElem].right = i+1;
    for (i = 3*mesh->nElem+1; i < 4*mesh->nElem; i++)
    {
        mesh->elements[i].left = i;
        mesh->elements[i].right = i+1;    
    }

    // Numero total de elementos e de pontos
    mesh->nElem *= 4;
    mesh->nPoints = mesh->nElem + 1;
    return mesh;
}

// Constroi uma malha a partir das regras do L-System
Mesh* buildMesh_4 (int argc, char *argv[])
{
    int i, j, k, cont;
    int id_fiber, iter, id_bif, start_fiber;
    double d_rot[3];
    Mesh *mesh = (Mesh*)malloc(sizeof(Mesh));
    Queue *q = newQueue();

    mesh->xMax = atof(argv[1]);
    mesh->nElem = atoi(argv[2]);
    mesh->nPoints = mesh->nElem + 1;
    mesh->h = mesh->xMax / (double)mesh->nElem;
    // Construir os pontos e os elementos
    int total_points = ((pow(2,MAX_ITER+1)-1)*mesh->nElem) + 1;
    mesh->points = (Point*)malloc(sizeof(Point)*total_points);
    mesh->elements = (Element*)malloc(sizeof(Element)*(total_points-1));
    // Construir a primeira fibra reta
    for (i = 0, k = 0; i < mesh->nPoints; i++)
    {
        if (i < mesh->nPoints-1)
        {
            mesh->elements[k].left = i;
            mesh->elements[k].right = i+1;
            k++;
        }
        mesh->points[i].x = i*mesh->h;
        mesh->points[i].y = mesh->points[i].z = 0; 
    }
    // Enfileirar o ultimo ponto como sendo de crescimento
    Enqueue(&q,i-1);
    // Iteracoes de crescimento do L-System
    id_fiber = 0; iter = 0;
    while (iter < MAX_ITER)
    {
        //printQueue(q);
        cont = q->in_the_queue;
        while (cont > 0)
        {
            id_bif = Dequeue(&q);
            // Crescer 45 graus para cima
            calculateRotation(d_rot,mesh->points[id_bif-1],mesh->points[id_bif],TETA);
            start_fiber = mesh->nPoints + (id_fiber*(mesh->nPoints-1));
            // Definir o elemento na bifurcacao
            mesh->elements[k].left = id_bif;
            mesh->elements[k].right = start_fiber;
            k++;
            for (i = start_fiber, j = 1; i < start_fiber+mesh->nElem; i++, j++)
            {
                if (j < mesh->nElem)
                {
                    mesh->elements[k].left = i;
                    mesh->elements[k].right = i+1;
                    k++;
                }
                mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
                mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
                mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
            }
            // Enfileira o ultimo ponto
            Enqueue(&q,i-1);
            id_fiber++;
            // Crescer 45 graus para baixo
            calculateRotation(d_rot,mesh->points[id_bif-1],mesh->points[id_bif],-TETA);
            start_fiber = mesh->nPoints + (id_fiber*(mesh->nPoints-1));
            // Definir o elemento na bifurcacao
            mesh->elements[k].left = id_bif;
            mesh->elements[k].right = start_fiber;
            k++;
            for (i = start_fiber, j = 1; i < start_fiber+mesh->nElem; i++, j++)
            {
                if (j < mesh->nElem)
                {
                    mesh->elements[k].left = i;
                    mesh->elements[k].right = i+1;
                    k++;
                }
                mesh->points[i].x = mesh->points[id_bif].x + j*d_rot[0]*mesh->h; 
                mesh->points[i].y = mesh->points[id_bif].y + j*d_rot[1]*mesh->h;
                mesh->points[i].z = mesh->points[id_bif].z + j*d_rot[2]*mesh->h;
            }
            // Enfileira o ultimo ponto
            Enqueue(&q,i-1);
            id_fiber++;
            cont--;
        }
        //printQueue(q);
        iter++;
    }
    mesh->nPoints = total_points;
    mesh->nElem = total_points - 1;
    //printf("Numero de pontos = %d\n",mesh->nPoints);
    //printf("Numero de elementos = %d\n",mesh->nElem);
    return mesh;
}

void calculateRotation (double d_rot[], Point p1, Point p2, double ang)
{
    // Calcula o vetor diferenca entre p1 e p2
    double d[3];
    d[0] = p2.x - p1.x;
    d[1] = p2.y - p1.y;
    d[2] = p2.z - p1.z;
    normalizeVector(d);

    // Rotaciona o vetor com angulo THETA em relacao ao eixo z
    d_rot[0] = d[0]*cos(ang) + d[1]*sin(ang);
    d_rot[1] = d[1]*cos(ang) - d[0]*sin(ang);
    d_rot[2] = d[2];

    /*
    printf("Rotacionado %lf\n",ang);
    printf("sin = %lf\n",sin(ang));
    printf("cos = %lf\n",cos(ang));
    printf("x = %lf\n",d_rot[0]);
    printf("y = %lf\n",d_rot[1]);
    printf("z = %lf\n",d_rot[2]);
    */    

}

// Normaliza um vetor
void normalizeVector (double d_rot[])
{
    double norm = sqrt(pow(d_rot[0],2)+pow(d_rot[1],2)+pow(d_rot[2],2));
    for (int i = 0; i < 3; i++) d_rot[i] /= norm;
}

// Escrever a malha em um arquivo
void writeMeshToFile (Mesh *mesh, char *filename)
{
    if (mesh != NULL)
    {
        FILE *file = fopen(filename,"w+");
        fprintf(file,"%d %d %lf\n",mesh->nElem,mesh->nPoints,mesh->h);
        for (int i = 0; i < mesh->nPoints; i++)
            fprintf(file,"%lf %lf %lf\n",mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
        for (int i = 0; i < mesh->nElem; i++)
            fprintf(file,"%d %d\n",mesh->elements[i].left,mesh->elements[i].right);
        fclose(file);
    }   
}

// Imprime informacoes da malha
void printMeshInfo (Mesh *mesh)
{
    if (mesh != NULL)
    {
        printf("=============================== MALHA ================================\n");
        printf("Tamanho da fibra = [0,%lf]\n",mesh->xMax);
        printf("Numero de elementos por fibra = %d\n",mesh->nElem);
        printf("Numero de pontos por fibra = %d\n",mesh->nPoints);
        printf("h = %lf\n",mesh->h);
        printf("----------------------------------------------------------------------\n");
        printf("Pontos:\n");
        for (int i = 0; i < mesh->nPoints; i++)
            printf("%d - %lf %lf %lf\n",i,mesh->points[i].x,mesh->points[i].y,mesh->points[i].z);
        printf("----------------------------------------------------------------------\n");
        printf("Elementos:\n");
        for (int i = 0; i < mesh->nElem; i++)
            printf("%d - %d %d\n",i,mesh->elements[i].left,mesh->elements[i].right);
        printf("======================================================================\n");
    }   
}
