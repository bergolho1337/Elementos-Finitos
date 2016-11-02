#include "mef.h"

MEF* newMEF (int argc, char *argv[])
{
  printf("------------------------------------------------------------------\n");
  MEF *mef = (MEF*)malloc(sizeof(MEF));
  mef->N = atoi(argv[1])+1;
  mef->dx = atof(argv[2]);
  mef->numeracao = atoi(argv[3]);
  mef->functions = getFunctions();
  printf(">>>>>>>>>>>> RESOLVENDO USANDO FUNCOES CHAPEU <<<<<<<<<<<<<<<<<<<<<<<<\n");
  setCoeficients(mef);
  setPoints(mef->N,mef->numeracao);
  setMapping(mef);

  makeElementMatrix(mef);
  makeGlobalMatrix(mef);
  makeVector(mef);

  setBoundaryConditions(mef);
  //mef->alfa = solveLinearSystem_CG(mef->A,mef->F,mef->N);
  mef->alfa = LU(mef->A,mef->F,mef->N,mef->N);
  writeDataFile(mef->functions,mef->alfa,mef->a,mef->b,mef->c,mef->dx,mef->N);

  #ifdef DEBUG
  printVector("Solucao do sistema:",mef->alfa,mef->N);
  #endif

  return mef;
}

Func* getFunctions ()
{
  Func *functions = (Func*)malloc(sizeof(Func)*2);
  // Funcao Chapeu
  functions[0] = phi_k;
  functions[1] = f;

  return functions;
}

// >>>>>>> MUDAR OS COEFICIENTES DO PROBLEMA AQUI <<<<<<<<<
void setCoeficients (MEF *mef)
{
  //mef->a = 1.0; mef->b = 0.4; mef->c = 0.0;       // Caso 1
  mef->a = 10.0; mef->b = 1000.0; mef->c = 0.0;      // Caso 2
  printf("[!] Coeficientes ---> || a = %e || b = %e || c = %e ||\n",mef->a,mef->b,mef->c);
}

void setMapping (MEF *mef)
{
  int i, j;
  int num_elem = mef->N - 1;
  int numeracao = mef->numeracao;
  mef->mapping = (int*)malloc(sizeof(int)*num_elem*2);

  // Mapeamento automatico (padrao)
  if (numeracao == 1)
  {
    for (i = 0; i < num_elem; i++)
    {
      for (j = 0; j < 2; j++)
        mef->mapping[i * 2 + j] = i + j;
    }
  }
  // Mapeamento manual
  else
  {
    printf("\nDigite a numeracao dos elementos:\n");
    for (i = 0; i < num_elem; i++)
    {
      printf("Elemento %d = ",i);
      for (j = 0; j < 2; j++)
        scanf("%d",&mef->mapping[i*2+j]);
    }
  }

  #ifdef DEBUG
  printf("Mapeamento dos nos:\n\n");
  for (i = 0; i < num_elem; i++)
  {
    for (j = 0; j < 2; j++)
      printf("%d ",mef->mapping[i*2+j]);
    printf("\n");
  }
  #endif
}

void makeElementMatrix (MEF *mef)
{
    printf("\n[+] Construindo matriz local do sistema ... ");
    fflush(stdout);
    double h;
    h = 1 / (double)(mef->N - 1);

    mef->K = (double*)malloc(sizeof(double)*4);
    mef->K[0] = (mef->a/h) - (mef->b/2.0) + (mef->c*h/3.0);
    mef->K[1] = (-mef->a/h) + (mef->b/2.0) + (mef->c*h/6.0);
    mef->K[2] = (-mef->a/h) - (mef->b/2.0) + (mef->c*h/6.0);
    mef->K[3] = (mef->a/h) + (mef->b/2.0) + (mef->c*h/3.0);

    printf("ok\n");
    #ifdef DEBUG
    printMatrix("Matriz do elemento:",mef->K,2);
    #endif
}

void makeGlobalMatrix (MEF *mef)
{
    printf("\n[+] Construindo matriz global do sistema ... ");
    fflush(stdout);
    int i, j, k;
    int ig, jg;
    int num_elem = mef->N - 1;
    int N = mef->N;
    // Alocar memoria e inicializar todos os elementos da matriz global como zero
    mef->A = (double*)calloc(N*N,sizeof(double));
    // Para cada elemento mapear a matriz local para a global
    for (k = 0; k < num_elem; k++)
    {
      // Percorrer a matriz local
      for (i = 0; i < 2; i++)
      {
        for (j = 0; j < 2; j++)
        {
          // Calcular os indices na matriz global
          ig = mef->mapping[k*2+i];
          jg = mef->mapping[k*2+j];
          mef->A[ig*N+jg] += mef->K[i*2+j];
        }
      }
    }

    printf("ok\n");

    #ifdef DEBUG
    printMatrix("Matriz do global:",mef->A,N);
    #endif
}

void makeVector (MEF *mef)
{
  printf("\n[+] Construindo vetor dos termos independentes ... ");
  fflush(stdout);
  int i;
  int N;
  N = mef->N;
  mef->F = (double*)malloc(sizeof(double)*N);
  for (i = 0; i < N; i++)
    mef->F[i] = trapezoidalRule(i,0,mef->functions[0],mef->functions[1],0,1);

  printf("ok\n");

  #ifdef DEBUG
  printVector("Vetor F",mef->F,N);
  #endif

}

void setBoundaryConditions (MEF *mef)
{
  int i, j;
  int N, numeracao;
  N = mef->N;
  numeracao = mef->numeracao;

  // Numeracao padrao
  if (numeracao == 1)
  {
    // A primeira linha e a ultima contem a condicao de contorno || u_0 = 0 || u_n-1 = 1 ||
    mef->F[0] = 0.0;
    mef->F[N-1] = 1.0;
    // A primeira linha contem somente o valor 1 na coluna do no 0
    for (j = 0; j < N; j++)
    {
      if (j == 0)
        mef->A[j] = 1.0;
      else
        mef->A[j] = 0.0;
    }
    // A ultima linha contem somente o valor 1 na coluna do no N-1
    for (j = 0; j < N; j++)
    {
      if (j == N-1)
        mef->A[(N-1)*N+j] = 1.0;
      else
        mef->A[(N-1)*N+j] = 0.0;
    }
  }
  // Numeracao manual
  else
  {
    int bc_1, bc_2;
    printf("Digite os nos das condicoes de contorno: ");
    scanf("%d %d",&bc_1,&bc_2);

    mef->F[bc_1] = 0.0;
    mef->F[bc_2] = 1.0;

    for (j = 0; j < N; j++)
    {
      if (j == bc_1)
        mef->A[bc_1*N+j] = 1.0;
      else
        mef->A[bc_1*N+j] = 0.0;
    }
    for (j = 0; j < N; j++)
    {
      if (j == bc_2)
        mef->A[bc_2*N+j] = 1.0;
      else
        mef->A[bc_2*N+j] = 0.0;
    }
  }

  #ifdef DEBUG
  printMatrix("Matriz do global apos aplicar as condicoes de contorno:",mef->A,N);
  #endif

  #ifdef DEBUG
  printVector("Vetor dos termos independentes com as condicoes de contorno:",mef->F,N);
  #endif
}
