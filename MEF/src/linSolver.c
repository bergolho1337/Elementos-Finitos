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

double* LU (double *A, double *b, int m, int n)
{
  int i, j, k, p;
	double *pivot = calloc(n,sizeof(double));
	double Amax, t, aux, r, Mult;
	// 1 PASSO: Transformar a matriz A do problema em duas matrizes triangulares L e U.
	for (i = 0; i < n; i++)
		pivot[i] = i;
	for (j = 0; j < n-1; j++)
	{
		// Escolher pivot
		p = j;
		Amax = abs(A[j*m+j]);
		// Verifica na coluna a ser eliminada qual elemento possui o maior valor absoluto, este elemento será o pivô.
		for (k = j+1; k < n; k++)
		{
			if (abs(A[k*m+j]) > Amax)
			{
				Amax = abs(A[k*m+j]);
				p = k;
			}
		}
		// Se (p != j) então deve-se trocar de linhas
		if (p != j)
		{
			for (k = 0; k < n; k++)
			{
				t = A[j*m+k];
				A[j*m+k] = A[p*m+k];
				A[p*m+k] = t;
			}
			aux = pivot[j];
			pivot[j] = pivot[p];
			pivot[p] = aux;
		}
		if (abs(A[j*m+j]) != 0)
		{
			// Eliminação de Gauss
			r = 1 / A[j*m+j];
			for (i = j+1; i < n; i++)
			{
				Mult = A[i*m+j]*r;
				A[i*m+j] = Mult;
				for (k = j+1; k < n; k++)
					A[i*m+k] = A[i*m+k] - Mult*A[j*m+k];
			}
		}
	}
	// A matriz A agora é L\U. Elementos abaixo da diagonal principal são de L, os da diagonal principal para cima pertencem a U.

	// 2 PASSO: Realizar substituições sucessivas para resolver o sistema triangular inferior: Ly = b
	double *y = calloc(n,sizeof(double));
	double soma;
	k = pivot[0];
	y[0] = b[k];
	// Percorre somente os elementos de L em A.
	for (i = 1; i < n; i++)
	{
		soma = 0;
		for (j = 0; j <= i-1; j++)
		{
			soma += A[i*m+j]*y[j];
		}
		k = pivot[i];
		y[i] = b[k] - soma;
	}
	// 3 PASSO: Realizar substituições retroativas para resolver o sistema triangular superior: Ux = y
	double *x = calloc(n,sizeof(double));
	x[n-1] = y[n-1] / A[(n-1)*m+(n-1)];
	for (i = n-2; i >= 0; i--)
	{
		soma = 0;
		for (j = i+1; j < n; j++)
			soma += A[i*m+j]*x[j];
		x[i] = (y[i] - soma) / A[i*m+i];
	}

	free(pivot);
	free(y);

	return (x);
}
