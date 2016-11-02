#include "diff.h"

Diff* newDiff (int argc, char *argv[])
{
  printf("------------------------------------------------------------------\n");
  Diff *diff = (Diff*)malloc(sizeof(Diff));
  diff->h = atof(argv[1]);
  diff->N = 1 / diff->h;
  diff->functions = getFunctions();
  printf(">>>>>>>>>>>> RESOLVENDO USANDO DIFERENCAS FINITAS <<<<<<<<<<<<<<<<<<<<<<<<\n");

  setCoeficients(diff);
  makeMatrix(diff);
  makeVector(diff);
  setBoundaryConditions(diff);
  //mef->alfa = solveLinearSystem_CG(mef->A,mef->F,mef->N);
  diff->u = LU(diff->A,diff->F,diff->N,diff->N);
  writeDataFile(diff->u,diff->a,diff->b,diff->c,diff->h,diff->N);

  #ifdef DEBUG
  printVector("Solucao do sistema:",diff->u,diff->N);
  #endif

  return diff;
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
void setCoeficients (Diff *diff)
{
  //mef->a = 1.0; mef->b = 0.4; mef->c = 0.0;            // Caso 1
  diff->a = 10.0; diff->b = 1000.0; diff->c = 0.0;       // Caso 2
  diff->Pe = (diff->h*diff->b) / (2*diff->a);
  printf("[!] Coeficientes ---> || a = %e || b = %e || c = %e ||\n",diff->a,diff->b,diff->c);
  printf("[!] Peclet = %e\n",diff->Pe);
}


void makeMatrix (Diff *diff)
{
    printf("\n[+] Construindo matriz do sistema ... ");
    fflush(stdout);
    int i, j, k;
    int N = diff->N;
    // Alocar memoria e inicializar todos os elementos da matriz como zero
    diff->A = (double*)calloc(N*N,sizeof(double));
    // Preencher somente os elementos fora da primeira e ultima linha
    for (i = 1; i < N-1; i++)
    {
      diff->A[i*N+i] = 2.0;
      diff->A[i*N+i+1] = diff->Pe-1;
      diff->A[i*N+i-1] = -(diff->Pe+1);
    }

    printf("ok\n");

    #ifdef DEBUG
    printMatrix("Matriz do global:",diff->A,N);
    #endif
}

void makeVector (Diff *diff)
{
  printf("\n[+] Construindo vetor dos termos independentes ... ");
  fflush(stdout);
  int i;
  int N;
  N = diff->N;
  diff->F = (double*)calloc(N,sizeof(double));

  printf("ok\n");

  #ifdef DEBUG
  printVector("Vetor F",diff->F,N);
  #endif

}


void setBoundaryConditions (Diff *diff)
{
  int N;
  N = diff->N;

  // A primeira linha e a ultima contem a condicao de contorno || u_0 = 0 || u_n-1 = 1 ||
  diff->F[0] = 0.0;
  diff->F[N-1] = 1.0;

  // A primeira linha contem somente o valor 1 na coluna do no 0
  diff->A[0] = 1.0;

  // A ultima linha contem somente o valor 1 na coluna do no N-1
  diff->A[(N-1)*N+(N-1)] = 1.0;

  #ifdef DEBUG
  printMatrix("Matriz do global apos aplicar as condicoes de contorno:",diff->A,N);
  #endif

  #ifdef DEBUG
  printVector("Vetor dos termos independentes com as condicoes de contorno:",diff->F,N);
  #endif
}
