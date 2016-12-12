#include "../include/mefDerivada.h"

// Construtor do MEFDerivada
MEFDerivada* newMEFDerivada (int argc, char *argv[])
{
  MEFDerivada *mef = (MEFDerivada*)malloc(sizeof(MEFDerivada));
  mef->nElem = atoi(argv[1]);
  mef->x_max = atof(argv[2]);
  mef->typeElem = atoi(argv[3]);
  mef->dx = mef->x_max / (double)mef->nElem;
  mef->nPoints = mef->nElem + 1;
  mef->functions = buildFunctions();
  mef->map = buildMap(mef->nElem,mef->nPoints,mef->typeElem);
  mef->X = buildPoints(mef->nPoints,mef->dx);
  allocMemoryVectors(mef);
  
  #ifdef DEBUG
  printInfoModel(mef);
  #endif
  
  assembleMatrix(mef);
  assembleLoadVector(mef);
  solveLinearSystem(mef->K,mef->F,mef->d,mef->nPoints,mef->typeElem);

  writeDataFile(mef->d,mef->X,mef->functions,mef->nPoints,mef->typeElem);

  return mef;
}

// Constroi o vetor de funcoes que serao utilizados no programa
// 0 = Funcao do vetor de carga
// 1 = Funcao da solucao analitica
// 2 = Funcao da derivada analitica
Func* buildFunctions ()
{
  Func *func = (Func*)malloc(sizeof(Func)*3);
  func[0] = f;
  func[1] = analit;
  func[2] = deriv;
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

// Constroi a tabela de mapeamento local-global
int* buildMap (int nElem, int nPoints, int typeElem)
{
  int i, j;
  int *map = NULL;
  switch (typeElem)
  {
    case 1: {
              printf("[!] Utilizando elemento linear !\n");
              map =  (int*)calloc(nElem*2,sizeof(double));
              for (i = 0; i < nElem; i++)
              {
                for (j = 0; j < 2; j++)
                  map[i*2+j] = i+j;
              }
              break;
            }
    case 2: {
              printf("[!] Utilizando elemento quadratico !\n");
              map =  (int*)calloc(nElem*3,sizeof(double));
              for (i = 0; i < nElem; i++)
              {
                for (j = 0; j < 3; j++)
                  map[i*3+j] = i+j;
              }
              break;
            }
    case 3: {
              printf("[!] Utilizando elemento Hermite !\n");
              map =  (int*)calloc(nElem*4,sizeof(double));
              for (i = 0; i < nElem; i++)
              {
                for (j = 0; j < 4; j++)
                {
                  // Funcao de ponderacao normal
                  if (j < 2)
                    map[i*4+j] = i+j;
                  // Funcao de ponderacao relacionada a derivada
                  else
                    map[i*4+j] = i+j-2+nPoints;
                }
                  
              }
              break;
            }
    default:{
              printf("[-] ERROR! Elemento nao foi cadastrado!\n");
              break;
            }
          
  }
  return map;
}

// Alocar memoria para os vetores
void allocMemoryVectors (MEFDerivada *mef)
{
  int typeElem = mef->typeElem;
  switch (typeElem)
  {
    case 1: {
              mef->d = (double*)calloc(mef->nPoints,sizeof(double));
              mef->F = (double*)calloc(mef->nPoints,sizeof(double));
              break;
            }
    case 2: {
              mef->d = (double*)calloc(mef->nPoints,sizeof(double));
              mef->F = (double*)calloc(mef->nPoints,sizeof(double));
              break;
            }
    // Elementos de Hermite devem possuir 2*N valores, N sao para a funcao normal e outros N para a derivada
    case 3: {
              mef->d = (double*)calloc(mef->nPoints*2,sizeof(double));
              mef->F = (double*)calloc(mef->nPoints*2,sizeof(double));
              break;
            }
  }
}

void printInfoModel (MEFDerivada *mef)
{
  int i;
  printf("======================== INFO MODEL ============================\n");
  printf("Number of elements = %d\n",mef->nElem);
  printf("Number of points = %d\n",mef->nPoints);
  printf("dx = %e\n",mef->dx);
  printf("------------------------- MAPPING ------------------------------\n");
  switch (mef->typeElem)
  {
    case 1: {
              printf("[!] Elementos lineares ...\n");
              for (i = 0; i < mef->nElem; i++)
                printf("%d %d\n",mef->map[i*2],mef->map[i*2+1]);
              break;
            }
    case 2: {
              printf("[!] Elementos quadraticos ... (FALTA MONTAR)\n");
              for (i = 0; i < mef->nElem; i++)
                printf("%d %d %d\n",mef->map[i*3],mef->map[i*3+1],mef->map[i*3+2]);
              exit (-1);
              break;
            }
    case 3: {
              printf("[!] Elementos Hermite ...\n");
              for (i = 0; i < mef->nElem; i++)
                printf("%d %d %d %d\n",mef->map[i*4],mef->map[i*4+1],mef->map[i*4+2],mef->map[i*4+3]);
              break;
            }
  }
  printf("------------------------- POINTS -------------------------------\n");
  for (i = 0; i < mef->nPoints; i++)
    printf("x%d - %e\n",i,mef->X[i]);
  printf("================================================================\n");
}

// Constroi a matriz usando a solucao analitica da integral do elemento associado
void assembleMatrix (MEFDerivada *mef)
{
  printf("[!] Construindo matriz ... ");
  fflush(stdout);
  int np, ne;
  double *local_A;
  ne = mef->nElem;
  np = 0;
  // Acerta o tamanho correto da matriz global dependendo do elemento
  switch (mef->typeElem)
  {
    case 1: np = mef->nPoints;
            break;
    case 2: //FALTA FAZER
            break;
    case 3: np = mef->nPoints*2;
            break;
  }

  // Construir a matriz local de massa
  local_A = buildLocalMassMatrix(mef->dx,mef->typeElem);
  mef->K = buildGlobalMatrixFromLocal(local_A,mef->map,np,ne,mef->typeElem);

  printf("ok\n");  

  // Colocar as condicoes de contorno
  setBoundaryConditions_Matrix(mef->K,np,mef->typeElem);
  LUDecomposition(mef->K,np);

  #ifdef DEBUG
  printMatrix("Global matrix K after LU decomposition",mef->K,np);
  #endif
}

// Constroi a matriz de massa local do elemento
// A_ij = \int phi_i' phi_j' dx
double* buildLocalMassMatrix (double h, int typeElem)
{
  int i, j;
  double *A = NULL;
  switch (typeElem)
  {
    case 1: {
              A = (double*)calloc(4,sizeof(double));
              for (i = 0; i < 2; i++)
              {
                for (j = 0; j < 2; j++)
                {
                  if (i == j)
                    A[i*2+j] = 1 / h;
                  else
                    A[i*2+j] = -1 / h;
                }
              }
              #ifdef DEBUG
              printMatrix("Local matrix",A,2);
              #endif
              break;
            }
    case 2: {
              A = (double*)calloc(9,sizeof(double));
              // FALTA MONTAR
              break;
            }
    case 3: {
              A = (double*)calloc(16,sizeof(double));
              double r = 1.0 / (30.0*h);
              A[0*4] = r*36.0;
              A[1*4+1] = r*36.0;
              A[2*4+2] = r*4.0;
              A[3*4+3] = r*4.0;
              A[0*4+1] = A[1*4] = -36.0*r;
              A[0*4+2] = A[2*4] = r*3.0;
              A[0*4+3] = A[3*4] = r*3.0;
              A[1*4+2] = A[2*4+1] = -3.0*r;
              A[1*4+3] = A[3*4+1] = -3.0*r;
              A[2*4+3] = A[3*4+2] = -1.0*r;

              #ifdef DEBUG
              printMatrix("Local matrix",A,4);
              #endif

              break;
            }
  }
  
  return A;
}

// Constroi uma matriz global a partir da matriz local de cada elemento
double* buildGlobalMatrixFromLocal (double *local_M, int *map, int np, int ne, int typeElem)
{
  int i, j, k;
  int ig, jg;
  double *M = NULL;
  switch (typeElem)
  {
    // Elemento linear
    case 1: {
              M = (double*)calloc(np*np,sizeof(double));
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

              #ifdef DEBUG
              printMatrix("Global matrix K",M,np);
              #endif

              break;
            }
    // Elemento quadratico
    case 2: {
              // FALTA FAZER
              break;
            }
    // Elemento Hermite
    case 3: {
              // Lembrando que a matriz do elemento de Hermite fica (2N)x(2N)
              M = (double*)calloc(np*np,sizeof(double));
              // Para cada elemento
              for (k = 0; k < ne; k++)
              {
                // Percorrer a matriz local
                for (i = 0; i < 4; i++)
                {
                  ig = map[k*4+i];
                  for (j = 0; j < 4; j++)
                  {
                    jg = map[k*4+j];
                    // Juntar a contribuicao local do elemento na matriz global
                    M[ig*np+jg] += local_M[i*4+j];
                  }
                }
              }

              #ifdef DEBUG
              printMatrix("Global matrix K",M,np);
              #endif

              break;
            }
  }
  
  return M;
}

// Colocar as condicoes de contorno na matriz global atraves da penalizacao (so estou considerando Dirichlet)
void setBoundaryConditions_Matrix (double *K, int np, int typeElem)
{
  double epsilon = 1.0e+06;

  K[0] += epsilon;
  switch (typeElem)
  {
    case 1: K[(np-1)*np+(np-1)] += epsilon;
            break;
    case 2: // FALTA FAZER
            break;
    // A condicao de contorno tem que ser imposta na parte da matriz relacionada a funcao u 
    case 3: K[(np/2-1)*np+(np/2-1)] += epsilon;
            break; 
  }

  #ifdef DEBUG
  printMatrix("Global matrix after boundary condition",K,np);
  #endif 

}

// Montar o vetor de carga desse problema
// F = \int 1*phi_i dx
void assembleLoadVector (MEFDerivada *mef)
{
  int i;
  int np = 0;
  switch (mef->typeElem)
  {
    case 1: np = mef->nPoints;
            for (i = 0; i < np; i++)
            {
              // Somente o primeiro e o ultimo
              if (i == 0 || i == np-1)
                mef->F[i] = mef->dx / 2.0;
              else
                mef->F[i] = mef->dx;
            }
            break;
    case 2: // FALTA FAZER
            break;
    case 3: np = mef->nPoints*2;
            // Contribuicao da funcao phi_i
            for (i = 0; i < np/2; i++)
            {
              if (i == 0 || i == np/2-1)
                mef->F[i] = mef->dx / 2.0;
              else
                mef->F[i] = mef->dx;
            }
            // Contribuicao da funcao phi_i'
            for (i = np/2; i < np; i++)
            {
              if (i == np/2)
                mef->F[i] = mef->dx / 12.0;
              else if (i == np-1)
                mef->F[i] = -mef->dx / 12.0;
              else
                mef->F[i] = 0.0;
            }
            break;
  }
  
  #ifdef DEBUG
  printVector("Load vector F",mef->F,np);
  #endif

  setBoundaryConditions_Vector(mef->F,np,mef->typeElem);

}

void setBoundaryConditions_Vector (double *F, int np, int typeElem)
{
  F[0] = 0.0;
  switch (typeElem)
  {
    case 1: F[np-1] = 0.0;
            break;
    case 2: // FALTA FAZER
            break;
    // A condicao de contorno tem que ser imposta na parte do vetor relacionada a funcao u 
    case 3: F[np/2-1] = 0.0;
            break; 
  }

  #ifdef DEBUG
  printVector("Load vector F after boundary conditions",F,np);
  #endif
}

/*

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
*/