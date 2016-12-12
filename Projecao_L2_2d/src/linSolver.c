#include "linSolver.h"

double* solveLinearSystem_CG (double *A, double *b, int N)
{
  printf("[+] Resolvendo o sistema linear associado ... ");
  fflush(stdout);
  int k = 0;
  double alpha, beta, aux;
  double *r, *r_ant, *x, *p, *s;

  // Alocar memoria
  r = (double*)calloc(N,sizeof(double));
  r_ant = (double*)calloc(N,sizeof(double));
  x = (double*)calloc(N,sizeof(double));
  p = (double*)calloc(N,sizeof(double));
  s = (double*)calloc(N,sizeof(double));

  // O primeiro residuo eh o vetor 'b'
	memcpy(r,b,sizeof(double)*N);

	// Faz o primeiro chute da solucao ser o vetor nulo
	memset(x,0,N*sizeof(double));

	// Metodo iterativo
	while (norm(r,N) > EPSILON && k < ITER_MAX)
  {
    k++;
		aux = dotProduct(r,r,N);

    if (k == 1)
      memcpy(p,r,sizeof(double)*N);
    else
    {
      beta = aux / dotProduct(r_ant,r_ant,N);
      saxpy2(p,r,beta,N);
    }
    // Copia o ultimo vetor 'r'
    memcpy(r_ant,r,sizeof(double)*N);

    matvec(A,p,s,N);
    alpha = aux / dotProduct(p,s,N);
    saxpy(x,p,alpha,N);
    saxpy(r,s,-alpha,N);
  }
  printf("ok\n");
	printf("k = %d || Norm = %e\n",k,norm(r,N));

  free(r);
  free(r_ant);
  free(p);
  free(s);
  return x;
}

// Funcao que calcula o produto interno entre dois vetores 'a' e 'b' de tamanho 'n'
double dotProduct (double *a, double *b, int n)
{
  int i;
  double dot = 0.0;
  for (i = 0; i < n; i++)
    dot += a[i]*b[i];
  return dot;
}

// Funcao que calcula SAXPY entre dois vetores 'x' e 'y' e utilizando um escalar 'alpha'
// x_i = x_i + alpha*y_i
void saxpy (double *x, double *y, double alpha, int n)
{
  int i;
  for (i = 0; i < n; i++)
    x[i] += y[i]*alpha;
}

// Funcao que calcula SAXPY entre dois vetores 'x' e 'y' e utilizando um escalar 'alpha'
// x_i = y_i + alpha*x_i
void saxpy2 (double *x, double *y, double alpha, int n)
{
  int i;
  for (i = 0; i < n; i++)
    x[i] = y[i] + alpha*x[i];
}

// Funcao que calcula a multiplicacao matriz-vetor (x = A*b) usando a representacao de matriz esparsa CSR
void matvec (double *A, double *b, double *x, int n)
{
  int i, j;
  double sum;
  for (i = 0; i < n; i++)
  {
    sum = 0;
    for (j = 0; j < n; j++)
      sum += A[i*n+j]*b[j];
    x[i] = sum;
  }
}

// Funcao que calcula a norma euclidiana de um vetor 'v' de tamanho 'n'
double norm (double *v, int n)
{
  int i;
  double norm = 0.0;
  for (i = 0; i < n; i++)
    norm += v[i]*v[i];
  return sqrt(norm);
}

double checkSystem (double *A, double *b, double *x, int N)
{
  int i, j;
  double sum;
  double norm = 0;
  for (i = 0; i < N; i++)
  {
    sum = 0;
    for (j = 0; j < N; j++)
    {
      sum += A[i*N+j]*x[j];
    }
    printf("%e\n",b[i]-sum);
    norm += pow(b[i]-sum,2);
  }
  printf("Norma do sistema = %e\n",sqrt(norm));
  return sqrt(norm);
}

double* invertMatrix (double *A, int n)
{
  int i, j;
  double *b, *x;
  double *Ainv;
  Ainv = (double*)malloc(sizeof(double)*n*n);
  b = (double*)malloc(sizeof(double)*n);
  // Resolver um sistema linear para cada coluna da matriz inversa
  for (i = 0; i < n; i++)
  {
    // Acerta o vetor b
    for (j = 0; j < n; j++)
    {
      if (j == i)
        b[j] = 1;
      else
        b[j] = 0;
    }
    x = solveLinearSystem_CG(A,b,n);

    // Colocar a coluna na matriz
    for (j = 0; j < n; j++)
      Ainv[j * n + i] = x[j];

    free(x);
  }
  return Ainv;
}

void testInverse (double *A, double *Ainv, int N)
{
  printf("\n---------------------------------------------------\n");
  printf("[!] Teste da matriz inversa: A*A^-1 = I\n");
  int i, j, k;
  double sum;
  double I[N*N];
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
    {
      sum = 0;
      for (k = 0; k < N; k++)
        sum += A[i * N + k]*Ainv[k * N + j];
      I[i*N+j] = sum;
    }
  }
  printMatrix("A",A,N);
  printMatrix("A^-1",Ainv,N);
  printMatrix("I",I,N);
}

double normL21 (double *A, int N)
{
  int i, j;
  double norm = 0;
  for (i = 0; i < N; i++)
  {
    for (j = 0; j < N; j++)
      norm += pow(A[i*N+j],2);
  }
  return (sqrt(norm));
}
