#include "l2.h"

L2* newL2 (int argc, char *argv[])
{
  printf("------------------------------------------------------------------\n");
  L2 *l2 = (L2*)malloc(sizeof(L2));
  l2->N = atoi(argv[1])+1;
  l2->dx = atof(argv[2]);
  l2->num_funcao = atoi(argv[3]);
  l2->functions = getFunctions(l2->num_funcao);
  switch (l2->num_funcao)
  {
    case 1:  printf(">>>>>>>>>>>> RESOLVENDO USANDO POLINOMIOS SIMPLES <<<<<<<<<<<<<<<<\n");
             makeMatrix_Polinomial(l2);
             makeVector_Polinomial(l2);
             break;

    case 2:  printf(">>>>>>>>>>>> RESOLVENDO USANDO POLINOMIOS DE LAGRANGE <<<<<<<<<<<<<<<<\n");
             setPoints(l2->N);
             makeMatrix_Lagrange(l2);
             makeVector_Lagrange(l2);
             break;
    case 3:  printf(">>>>>>>>>>>> RESOLVENDO USANDO FUNCOES CHAPEU <<<<<<<<<<<<<<<<<<<<<<<<\n");
             setPoints(l2->N);
             makeMatrix_Lagrange(l2);
             makeVector_Lagrange(l2);
             break;
  }
  l2->alfa = solveLinearSystem_CG(l2->A,l2->b,l2->N);
  writeDataFile(l2->functions,l2->alfa,l2->dx,l2->N);

  #ifdef DEBUG
  printVector("Solucao do sistema:",l2->alfa,l2->N);
  #endif

  #ifdef ANALYZE
    printf("-------------------------------------------------------------------------------------------\n");
    printf("[+] Analisando o numero de condicionamento kapla ...\n");
    double kapla;
    kapla = analyzeKapla(l2);
    printf("\nKapla = %e\n",kapla);
  #endif
  return l2;
}

Func* getFunctions (int base_function)
{
  Func *functions = (Func*)malloc(sizeof(Func)*2);
  // Desvia para qual funcao base o usuario escolheu
  switch (base_function)
  {
    // Polinomial simples
    case 1: {
              functions[0] = p_k;
              break;
            }
    // Polinomio de Lagrange
    case 2: {
              functions[0] = l_k;
              break;
            }
    // Funcao Chapeu
    case 3: {
              functions[0] = phi_k;
              break;
            }
  }
  functions[1] = f;
  return functions;
}

void makeMatrix_Polinomial (L2 *l2)
{
  printf("[+] Construindo matriz dos coeficientes ... ");
  fflush(stdout);
  int i, j;
  int N;
  N = l2->N;
  l2->A = (double*)malloc(sizeof(double)*N*N);

  for (i = 0; i < N; i++)
  {
    // Matriz eh triangular superior
    for (j = i; j < N; j++)
    {
      l2->A[i * N + j] = integral_polinomium(i,j);
      // Preenche o elemento simetrico
      if (i != j)
        l2->A[j * N + i] = l2->A[i * N + j];
    }
  }

  #ifdef DEBUG
  printMatrix("Matriz A",l2->A,l2->N);
  #endif

  printf("ok\n");
}

void makeVector_Polinomial (L2 *l2)
{
  printf("[+] Construindo vetor dos termos independentes ... ");
  fflush(stdout);
  int i;
  int N;
  N = l2->N;
  l2->b = (double*)malloc(sizeof(double)*N);
  for (i = 0; i < N; i++)
  {
    l2->b[i] = trapezoidalRule(i,0,l2->functions[0],l2->functions[1],0,1);
  }

  #ifdef DEBUG
  printVector("Vetor b",l2->b,l2->N);
  #endif

  printf("ok\n");
}

void makeMatrix_Lagrange (L2 *l2)
{
  printf("[+] Construindo matriz dos coeficientes ... ");
  fflush(stdout);
  int i, j;
  int N;
  N = l2->N;
  l2->A = (double*)malloc(sizeof(double)*N*N);

  for (i = 0; i < N; i++)
  {
    // Matriz eh triangular superior
    for (j = i; j < N; j++)
    {
      l2->A[i * N + j] = trapezoidalRule(i,j,l2->functions[0],l2->functions[0],0,1);
      // Preenche o elemento simetrico
      if (i != j)
        l2->A[j * N + i] = l2->A[i * N + j];
    }
  }

  #ifdef DEBUG
  printMatrix("Matriz A",l2->A,l2->N);
  #endif

  printf("ok\n");
}

void makeVector_Lagrange (L2 *l2)
{
  printf("[+] Construindo vetor dos termos independentes ... ");
  fflush(stdout);
  int i;
  int N;
  N = l2->N;
  l2->b = (double*)malloc(sizeof(double)*N);
  for (i = 0; i < N; i++)
    l2->b[i] = trapezoidalRule(i,0,l2->functions[0],l2->functions[1],0,1);

  #ifdef DEBUG
  printVector("Vetor b",l2->b,l2->N);
  #endif

  printf("ok\n");
}

double analyzeKapla (L2 *l2)
{
  double *Ainv;
  Ainv = invertMatrix(l2->A,l2->N);

  #ifdef DEBUG
  testInverse(l2->A,Ainv,l2->N);
  #endif

  return normL21(l2->A,l2->N)*normL21(Ainv,l2->N);
}
