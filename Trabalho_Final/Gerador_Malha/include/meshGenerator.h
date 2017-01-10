#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "../include/queue.h"

struct Mesh;
struct Point;
struct Element;

const double TETA = M_PI / 6.0;     // Angulo de rotacao 
const int MAX_ITER = 3;             // Numero maximo de iteracoes de crescimento

struct Point
{
    double x, y, z;                 // Coordenadas do ponto
}typedef Point;

struct Element
{
    int left;                       // Identificador do ponto a esquerda do elemento
    int right;                      // Identificador do ponto a direita do elemento
}typedef Element;

struct Mesh
{
    int nElem;                      // Numero de elementos
    int nPoints;                    // Numero de pontos
    double xMax;                    // Tamanho da fibra
    double h;                       // Distancia entre dois pontos
    Point *points;                  // Vetor de pontos
    Element *elements;              // Vetor de elementos
}typedef Mesh;

Mesh* newMesh (int argc, char *argv[]);
Mesh* buildMesh_1 (int argc, char *argv[]);
Mesh* buildMesh_2 (int argc, char *argv[]);
Mesh* buildMesh_3 (int argc, char *argv[]);
Mesh* buildMesh_4 (int argc, char *argv[]);
void calculateRotation (double d_rot[], Point p1, Point p2, double ang);
void normalizeVector (double d_rot[]);
void writeMeshToFile (Mesh *mesh, char *filename);
void printMeshInfo (Mesh *mesh);