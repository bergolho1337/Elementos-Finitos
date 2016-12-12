#include "../include/calor.h"

// Construtor do Calor
Calor* newCalor (int argc, char *argv[])
{
  Calor *calor = (Calor*)malloc(sizeof(Calor));
  calor->nElem = atoi(argv[1]);
  calor->dt = atof(argv[2]);
  calor->t_max = atof(argv[3]);
  calor->x_max = atof(argv[4]);
  calor->dx = calor->x_max/(double)calor->nElem;
  calor->nPoints = calor->nElem + 1;
  calor->functions = buildFunctions();
  calor->map = buildMap(calor->nElem);
  calor->X = buildPoints(calor->nPoints,calor->dx);
  calor->dOld = setInitialConditions(calor->functions,calor->X,calor->nPoints);

  calor->dNew = (double*)calloc(calor->nPoints,sizeof(double));
  calor->F = (double*)calloc(calor->nPoints,sizeof(double));
  
  #ifdef DEBUG
  printInfoModel(calor);
  #endif
  
  assembleMatrix(calor);
  
  return calor;
}

// Constroi a matriz usando a solucao analitica da integral do elemento linear
// K = (A+B*dt)
void assembleMatrix (Calor *calor)
{
  printf("[!] Construindo matriz ... ");
  fflush(stdout);
  int np, ne;
  double *local_A, *local_B;
  ne = calor->nElem;
  np = calor->nPoints;

  // Construir as matrizes locais de massa e de rigidez
  local_A = buildLocalMassMatrix(calor->dx);
  local_B = buildLocalStiffMatrix(calor->dx,calor->dt);

  // Construir as matrizes globais A (massa) e B (rigidez)
  calor->A = buildGlobalMatrixFromLocal(local_A,calor->map,np,ne);
  calor->B = buildGlobalMatrixFromLocal(local_B,calor->map,np,ne);

  // Construir a matriz global do sistema linear: K = (A + B*dt)
  calor->K = buildGlobalMatrix(calor->A,calor->B,calor->dt,np);

  printf("ok\n");  
  #ifdef DEBUG
  printMatrix("Global mass matrix A",calor->A,np);
  printMatrix("Global stiff matrix B",calor->B,np);
  printMatrix("Global matrix K",calor->K,np);
  #endif

  // Colocar as condicoes de contorno
  setBoundaryConditions_Matrix(calor);
  LUDecomposition(calor->K,calor->nPoints);

  #ifdef DEBUG
  printMatrix("Global matrix K after LU decomposition",calor->K,np);
  #endif
}

// Constroi a matriz global do sistema linear
double* buildGlobalMatrix (double *A, double *B, double dt, int np)
{
  int i, j;
  double *K = (double*)calloc(np*np,sizeof(double));
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      K[i*np+j] = A[i*np+j] + B[i*np+j]*dt;
  }
  return K;
}

// Constroi uma matriz global a partir da matriz local de cada elemento
double* buildGlobalMatrixFromLocal (double *local_M, int *map, int np, int ne)
{
  int i, j, k;
  int ig, jg;
  double *M = (double*)calloc(np*np,sizeof(double));
  // Para cada elemento
  for (k = 0; k < ne; k++)
  {
    // Percorrer a matriz local
    for (i = 0; i < 2; i++)
    {
      ig = map[k*2+i];
      for (j = 0; j < 2; j++)
      {
        jg = map[k*2+j];
        // Juntar a contribuicao local do elemento na matriz global
        M[ig*np+jg] += local_M[i*2+j];
      }
    }
  }
  return M;
}

// Constroi a matriz de massa local (2x2) de cada elemento
double* buildLocalMassMatrix (double h)
{
  int i, j;
  double *A = (double*)calloc(4,sizeof(double));
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
    {
      if (i == j)
        A[i*2+j] = h / 3.0;
      else
        A[i*2+j] = h / 6.0;
    }
  }
  return A;
}

// Constroi a matriz de rigidez local (2x2) de cada elemento
double* buildLocalStiffMatrix (double h, double k)
{
  int i, j;
  double *B = (double*)calloc(4,sizeof(double));
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
    {
      if (i == j)
        B[i*2+j] = 1 / h;
      else
        B[i*2+j] = -1 / h;
    }
  }
  return B;
}

// Montar o vetor de carga desse problema (A*d^{n-1})
void assembleLoadVector (Calor *calor)
{
  int i, j;
  int np, ne, elem;
  np = calor->nPoints;
  ne = calor->nElem;
  // Iniciliza o vetor global com 0's
  memset(calor->F,0,sizeof(double)*np);
  // Montar o vetor global a partir da multiplicacao A*dOld
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      calor->F[i] += calor->A[i*np+j]*calor->dOld[j];
  }

  #ifdef DEBUG
  printVector("Load vector F",calor->F,np);
  #endif
}


// Constroi o vetor de carga local do elemento
void buildLocalVector (double *A, double *dOld, double *V, int *map, int elem)
{
  int i, j;
  int ig = map[elem*2];
  memset(V,0,sizeof(double)*2);
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 2; j++)
    {
      // O termo 0.0 esta relacionado a integral do termo F^n que eh zero pra esse problema (FAZER GERAL)
      V[i] += A[i*2+j]*dOld[ig+j] + 0.0;
    }
  }
}

// Mapeia o vetor local do elemento para o vetor global
void mapLocalToGlobal_Vector (double *V, double *F, int *map, int elem)
{
  int i, ig;
  ig = map[elem*2];
  for (i = 0; i < 2; i++)
    F[ig+i] += V[i]; 
}

// Colocar as condicoes de contorno na matriz global (so estou considerando Dirichlet)
void setBoundaryConditions_Matrix (Calor *calor)
{
  int np;
  np = calor->nPoints;
  
  // Acertar a solucao no primeiro elemento 
  calor->K[0] = 1.0;
  calor->K[1] = 0.0; 
  calor->K[np] = 0.0;
  //calor->K[np+1] -= calor->A[0] + calor->B[0]*calor->dt;

  // Acertar a solucao no ultimo elemento
  calor->K[(np-1)*np+(np-1)] = 1.0;
  calor->K[(np-1)*np+(np-2)] = 0.0;
  calor->K[(np-2)*np+(np-1)] = 0.0;
  //calor->K[(np-2)*np+(np-2)] -= calor->A[0] + calor->B[0]*calor->dt;

  #ifdef DEBUG
  printMatrix("Global matrix after boundary condition",calor->K,np);
  #endif 

}

// Setar as condicoes de contorno no vetor global (TA GENERICO)
void setBoundaryConditions_Vector (Calor *calor)
{
  int np = calor->nPoints;
  calor->F[0] = 0.0;
  calor->F[np-1] = 0.0;

  #ifdef DEBUG
  printVector("Load vector after boundary conditions",calor->F,np);
  #endif
}

// Constroi as condicoes iniciais
double* setInitialConditions (Func *func, double *X, int np)
{
  int i;
  double *u = (double*)calloc(np,sizeof(double));
  for (i = 0; i < np; i++)
    u[i] = func[1](0,0.0,X[i]);
  
  #ifdef DEBUG
  printVector("Initial condition",u,np);
  #endif // DEBUG
  return u;
}

// Constroi a tabela de mapeamento local-global
int* buildMap (int nElem)
{
  int i, j;
  int *map = (int*)calloc(nElem*2,sizeof(double));
  for (i = 0; i < nElem; i++)
  {
    for (j = 0; j < 2; j++)
      map[i*2+j] = i+j;
  }
  return map;
}

// Constroi o vetor de funcoes que serao utilizados no programa
// 0 = Funcao do vetor de carga
// 1 = Funcao da condicao inicial
// 2 = Funcao da solucao analitica
Func* buildFunctions ()
{
  Func *func = (Func*)malloc(sizeof(Func)*3);
  func[0] = f;
  func[1] = u0;
  func[2] = analit;
  return func;
}

// Constroi o vetor de pontos do dominio igualmente espacados
double* buildPoints (int np, double h)
{
  int i;
  double *x = (double*)calloc(np,sizeof(double));
  for (i = 0; i < np; i++)
    x[i] = i*h;
  return x;
}

// Resolvedor do problema transiente
void solveTransientProblem (Calor *calor)
{
  double t;
  // A matriz global do problema jah se encontra com as condicoes de contorno e LU
  printf("[!] Resolvendo o problema transiente ... ");
  fflush(stdout);
  int i, N;
  N = nearbyint(calor->t_max / calor->dt);
  // Iterar o metodo a cada passo de tempo
  for (i = 0; i < N; i++)
  {
    t = i*calor->dt;
    //writeDataFile(calor->dOld,calor->X,calor->functions,t,calor->nPoints,i);
    writeDataFile2(calor->dOld,calor->X,calor->functions,t,calor->nPoints,i);
    assembleLoadVector(calor);
    setBoundaryConditions_Vector(calor);
    solveLinearSystem_LU(calor->K,calor->F,calor->dNew,calor->nPoints);

    memcpy(calor->dOld,calor->dNew,sizeof(double)*calor->nPoints);
  } 
  printf("ok\n");
}

void printInfoModel (Calor *calor)
{
  int i;
  printf("======================== INFO MODEL ============================\n");
  printf("Number of elements = %d\n",calor->nElem);
  printf("Number of points = %d\n",calor->nPoints);
  printf("dt = %e\n",calor->dt);
  printf("t_max = %e\n",calor->t_max);
  printf("dx = %e\n",calor->dx);
  printf("------------------------- MAPPING ------------------------------\n");
  for (i = 0; i < calor->nElem; i++)
    printf("%d %d\n",calor->map[i*2],calor->map[i*2+1]);
  printf("------------------------- POINTS -------------------------------\n");
  for (i = 0; i < calor->nPoints; i++)
    printf("x%d - %e\n",i,calor->X[i]);
  printf("================================================================\n");
}