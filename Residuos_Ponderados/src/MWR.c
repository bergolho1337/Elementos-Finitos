#include "MWR.h"

MWR* newMWR (int argc, char *argv[])
{
  printf("------------------------------------------------------------------\n");
  MWR *mwr = (MWR*)malloc(sizeof(MWR));
  mwr->N = atoi(argv[1]);
  mwr->dx = atof(argv[2]);
  mwr->num_metodo = atoi(argv[3]);
  mwr->functions = getFunctions(mwr->N);
  switch (mwr->num_metodo)
  {
    case 1:{
              printf(">>>>>>>>>>>> RESOLVENDO POR GALERKIN <<<<<<<<<<<<<<<<\n");
              makeMatrix_Galerkin(mwr);
              makeVector_Galerkin(mwr);
              mwr->alfa = solveLinearSystem_LU(mwr->A,mwr->b,mwr->N);
              break;
           }
    case 2:{
             printf(">>>>>>>>>>>> RESOLVENDO POR COLOCACAO <<<<<<<<<<<<<<<<\n");
             readEpsilon(mwr);
             makeMatrix_Colocacao(mwr);
             makeVector_Colocacao(mwr);
             mwr->alfa = solveLinearSystem_LU(mwr->A,mwr->b,mwr->N);
             break;
           }
  }

  writeDataFile(mwr->functions,mwr->alfa,mwr->dx,mwr->N);
  return mwr;
}

Func* getFunctions (int N)
{
  // Para cada phi_j deve haver uma L_j, afinal a matriz do metodo eh quadrada.
  // Mais a funcao do forcante
  int num_functions = N*2 + 1;
  Func *functions = (Func*)malloc(sizeof(Func)*num_functions);
  switch (N)
  {
    case 1: {
              functions[0] = phi_1;
              functions[1] = L_1;
              functions[2] = f;
              break;
            }
    case 2: {
              functions[0] = phi_1;
              functions[1] = phi_2;
              functions[2] = L_1;
              functions[3] = L_2;
              functions[4] = f;
              break;
            }
    case 3: {
              functions[0] = phi_1;
              functions[1] = phi_2;
              functions[2] = phi_3;
              functions[3] = L_1;
              functions[4] = L_2;
              functions[5] = L_3;
              functions[6] = f;
              break;
            }
    default:{
              printf("[!] ERRO! As funcoes nao foram setadas para N = %d\n",N);
              printf("[-] Cancelando operacao ... x(\n");
              exit(1);
            }
  }
  return functions;
}

void makeMatrix_Galerkin (MWR *mwr)
{
  printf("[+] Construindo matriz dos coeficientes ... ");
  int i, j;
  int N;
  mwr->A = (double*)malloc(sizeof(double)*mwr->N*mwr->N);
  N = mwr->N;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      mwr->A[i * N + j] = trapezoidalRule(mwr->functions[i],mwr->functions[N+j],0,1);
    }
  }
  printf("ok\n");
}

void makeVector_Galerkin (MWR *mwr)
{
  printf("[+] Construindo vetor dos termos independentes ... ");
  int i;
  int N, last_func_id;
  N = mwr->N;
  last_func_id = 2*N;
  mwr->b = (double*)malloc(sizeof(double)*N);
  for (i = 0; i < N; i++)
  {
    mwr->b[i] = trapezoidalRule(mwr->functions[i],mwr->functions[last_func_id],0,1);
  }
  printf("ok\n");
}

void readEpsilon (MWR *mwr)
{
  int i;
  int N = mwr->N;
  mwr->epsilon = (double*)malloc(sizeof(double)*N);
  printf("[+] Lendo %d valores para os pontos de colocacao epsilon ...\n",N);
  for (i = 0; i < N; i++)
  {
    printf("Valor %d: ",i+1);
    scanf("%lf",&mwr->epsilon[i]);
  }
}

void makeMatrix_Colocacao (MWR *mwr)
{
  printf("[+] Construindo matriz dos coeficientes ... ");
  int i, j;
  int N;
  N = mwr->N;
  mwr->A = (double*)malloc(sizeof(double)*N*N);
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      mwr->A[i * N + j] = mwr->functions[N+j](mwr->epsilon[i]);
    }
  }
  printf("ok\n");
}

void makeVector_Colocacao (MWR *mwr)
{
  printf("[+] Construindo vetor dos termos independentes ... ");
  int i;
  int N, last_func_id;
  N = mwr->N;
  last_func_id = 2*N;
  mwr->b = (double*)malloc(sizeof(double)*N);
  for (i = 0; i < N; i++)
  {
    mwr->b[i] = mwr->functions[last_func_id](mwr->epsilon[i]);
  }
  printf("ok\n");
}
