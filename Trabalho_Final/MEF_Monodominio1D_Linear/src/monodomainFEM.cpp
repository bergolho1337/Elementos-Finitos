#include "../include/monodomainFEM.h"

// Construtor da estrutura do resolvedor da equacao do monodominio
MonodomainFEM* newMonodomainFEM (int argc, char *argv[])
{
    MonodomainFEM *monoFEM = (MonodomainFEM*)malloc(sizeof(MonodomainFEM));
    monoFEM->nElem = atoi(argv[1]);
    monoFEM->dt = atof(argv[2]);
    monoFEM->t_max = atof(argv[3]);
    monoFEM->x_max = atof(argv[4]);
    monoFEM->nPoints = monoFEM->nElem + 1;
    monoFEM->M = nearbyint(monoFEM->t_max / monoFEM->dt);
    monoFEM->dx = monoFEM->x_max / (double)monoFEM->nElem;
    monoFEM->map = buildMap(monoFEM->nElem);
    monoFEM->X = buildPoints(monoFEM->nPoints,monoFEM->dx);
    monoFEM->functions = buildFunctions();
    
    // Alocar memoria e setar as condicoes iniciais (para o Fitz-Hugh Nagumo eh tudo 0)
    monoFEM->VOld = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->Vstar = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->VNew = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->wOld = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->wNew = (double*)calloc(monoFEM->nPoints,sizeof(double));
    monoFEM->F = (double*)calloc(monoFEM->nPoints,sizeof(double));

    // Construir a matriz global do sistema linear ligado a solucao da EDP
    assembleMatrix(monoFEM);

    #ifdef DEBUG
    printInfoModel(monoFEM);
    #endif

    return monoFEM;
}

// Imprime informacoes sobre o modelo construido
void printInfoModel (MonodomainFEM *monoFEM)
{
    int i;
    printf("======================== INFO MODEL ============================\n");
    printf("Number of elements = %d\n",monoFEM->nElem);
    printf("Number of points = %d\n",monoFEM->nPoints);
    printf("dt = %e\n",monoFEM->dt);
    printf("t_max = %e\n",monoFEM->t_max);
    printf("dx = %e\n",monoFEM->dx);
    printf("------------------------- MAPPING ------------------------------\n");
    for (i = 0; i < monoFEM->nElem; i++)
    printf("%d %d\n",monoFEM->map[i*2],monoFEM->map[i*2+1]);
    printf("------------------------- POINTS -------------------------------\n");
    for (i = 0; i < monoFEM->nPoints; i++)
    printf("x%d - %e\n",i,monoFEM->X[i]);
    printf("================================================================\n");
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

// Constroi o vetor de pontos do dominio igualmente espacados
double* buildPoints (int np, double h)
{
  int i;
  double *x = (double*)calloc(np,sizeof(double));
  for (i = 0; i < np; i++)
    x[i] = i*h;
  return x;
}

// Constroi o vetor de funcoes que serao utilizados no programa
// 0 = Funcao do potencial transmembranico (Vm)
// 1 = Funcao da variavel de estado (w)
Func* buildFunctions ()
{
  Func *func = (Func*)malloc(sizeof(Func)*num_eq);
  func[0] = dvdt__Fitz;
  func[1] = dwdt__Fitz;
  return func;
}

// Constroi a matriz global K do metodo implicito da solucao da EDP pelo MEF
// K = (BETA*Cm*A + SIGMA*dt*B)
void assembleMatrix (MonodomainFEM *monoFEM)
{
  printf("[!] Construindo matriz ... ");
  fflush(stdout);
  int np, ne;
  double *local_A, *local_B;
  ne = monoFEM->nElem;
  np = monoFEM->nPoints;

  // Construir as matrizes locais de massa e de rigidez 
  local_A = buildLocalMassMatrix(monoFEM->dx);
  local_B = buildLocalStiffMatrix(monoFEM->dx);

  // Construir as matrizes globais A (massa) e B (rigidez)
  monoFEM->A = buildGlobalMatrixFromLocal(local_A,monoFEM->map,np,ne);
  monoFEM->B = buildGlobalMatrixFromLocal(local_B,monoFEM->map,np,ne);

  // Construir a matriz global do sistema linear: 
  monoFEM->K = buildGlobalMatrix(monoFEM->A,monoFEM->B,monoFEM->dt,np);

  printf("ok\n");  
  #ifdef DEBUG
  printMatrix("Global mass matrix A",monoFEM->A,np);
  printMatrix("Global stiff matrix B",monoFEM->B,np);
  printMatrix("Global matrix K",monoFEM->K,np);
  #endif

  // !!! As condicoes de contorno saem naturalmente da formulacao variacional !!!
  // Decompor a matriz K em LU
  LUDecomposition(monoFEM->K,monoFEM->nPoints);

  #ifdef DEBUG
  printMatrix("Global matrix K after LU decomposition",monoFEM->K,np);
  #endif

  // Liberar memoria
  free(local_A);
  free(local_B);

}

// Constroi a matriz local de massa de cada elemento (2x2)
// !!! USANDO ELEMENTOS LINEARES !!!
// || phi_1 = 1 - x/h || phi_2 = x/h || --> integral_(0,h) {phi_i . phi_j dx}
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

// Constroi a matriz local de rigidez de cada elemento (2x2)
// !!! USANDO ELEMENTOS LINEARES !!!
// || phi_1 = 1 - x/h || phi_2 = x/h || --> integral_(0,h) {phi_i' . phi_j' dx}
double* buildLocalStiffMatrix (double h)
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

// Constroi a matriz global do sistema linear da EDP
// K = (BETA*Cm*A + SIGMA*dt*B)
double* buildGlobalMatrix (double *A, double *B, double dt, int np)
{
  int i, j;
  double *K = (double*)calloc(np*np,sizeof(double));
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      K[i*np+j] = (BETA*Cm*A[i*np+j]) + (SIGMA*dt*B[i*np+j]);
  }
  return K;
}

// Constroi o vetor de carga do sistema linear relacionado a EDP
// F = (BETA*Cm*A*VOld)
void assembleLoadVector (MonodomainFEM *monoFEM)
{
  int i, j;
  int np;
  np = monoFEM->nPoints;
  // Iniciliza o vetor global com 0's
  memset(monoFEM->F,0,sizeof(double)*np);
  // Montar o vetor global a partir da multiplicacao A*VOld
  for (i = 0; i < np; i++)
  {
    for (j = 0; j < np; j++)
      monoFEM->F[i] += monoFEM->A[i*np+j]*monoFEM->VOld[j];
    // Multiplicar pelo coeficiente BETA*Cm
    monoFEM->F[i] *= BETA*Cm;
  }

  #ifdef DEBUG
  printVector("Load vector F",monoFEM->F,np);
  #endif
}

// Resolver o sistema nao linear de EDOs no tempo atual
void solveEDO (MonodomainFEM *monoFEM, double t)
{
  int np;
  int i, point;
  double f, dt;
  np = monoFEM->nPoints;
  dt = monoFEM->dt;
  // Resolver o sistema de EDO para cada ponto
  for (point = 0; point < np; point++)
  {
    for (i = 0; i < num_eq; i++)
    {
      // V^{n+1} = V^{*} + f*dt
      if (i == 0)
      {
        f = monoFEM->functions[i](point,t,monoFEM->Vstar[point],monoFEM->wOld[point]);
        monoFEM->VNew[point] = monoFEM->Vstar[point] + f*dt;
      }
      // w^{n+1} = w^{n} + f*dt
      else
      {
        f = monoFEM->functions[i](point,t,monoFEM->VOld[point],monoFEM->wOld[point]);
        monoFEM->wNew[point] = monoFEM->wOld[point] + f*dt;
      }     
    }
  } 
}

// Resolve a equacao do monodominio
void solveMonodomain (MonodomainFEM *monoFEM)
{
  double t;
  // A matriz global do problema jah se encontra como LU
  printf("[!] Resolvendo o problema transiente ... ");
  fflush(stdout);
  int i, M;
  M = monoFEM->M;
  // Iterar o metodo a cada passo de tempo
  for (i = 0; i < M; i++)
  {
    t = i*monoFEM->dt;
    // Escrever em arquivo .vtk
    writeVTKFile(monoFEM->VOld,monoFEM->X,monoFEM->map,monoFEM->nPoints,monoFEM->nElem,i);

    // Resolver a EDP (parte difusiva)
    assembleLoadVector(monoFEM);
    //setBoundaryConditions_Vector(monoFEM); // As condicoes de contorno saem naturalmente da formulacao fraca
    solveLinearSystem_LU(monoFEM->K,monoFEM->F,monoFEM->Vstar,monoFEM->nPoints);

    // Resolver as EDOs (parte reativa)
    solveEDO(monoFEM,t);

    // Passa para a proxima iteracao
    memcpy(monoFEM->VOld,monoFEM->VNew,sizeof(double)*monoFEM->nPoints);
    memcpy(monoFEM->wOld,monoFEM->wNew,sizeof(double)*monoFEM->nPoints);
  } 
  printf("ok\n");
}

void freeMonodomain (MonodomainFEM *monoFEM)
{
  printf("[!] Liberando memoria ... ");
  fflush(stdout);
  free(monoFEM->map);
  free(monoFEM->functions);
  free(monoFEM->X);
  free(monoFEM->K);
  free(monoFEM->A);
  free(monoFEM->B);
  free(monoFEM->F);
  free(monoFEM->VNew);
  free(monoFEM->VOld);
  free(monoFEM->Vstar);
  free(monoFEM->wOld);
  free(monoFEM->wNew);
  free(monoFEM);
  printf("ok\n");
}