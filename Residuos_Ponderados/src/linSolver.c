#include "linSolver.h"

double* solveLinearSystem_LU (double *A, double *b, int N)
{
  printf("[+] Resolvendo o sistema linear associado ... ");
  int i, j, k, p;
	double Amax, t, m, r, Mult;
  // Aloca memoria
	double *pivot = calloc(N,sizeof(double));

	// 1 PASSO: Transformar a matriz A do problema em duas matrizes triangulares L e U.
	for (i = 0; i < N; i++)
		pivot[i] = i;
	for (j = 0; j < N-1; j++)
	{
		// Escolher pivot
		p = j;
		Amax = abs(A[j*N+j]);
		// Verifica na coluna a ser eliminada qual elemento possui o maior valor absoluto, este elemento será o pivô.
		for (k = j+1; k < N; k++)
		{
			if (abs(A[k*N+j]) > Amax)
			{
				Amax = abs(A[k*N+j]);
				p = k;
			}
		}
		// Se (p != j) então deve-se trocar de linhas
		if (p != j)
		{
			for (k = 0; k < N; k++)
			{
				t = A[j*N+k];
				A[j*N+k] = A[p*N+k];
				A[p*N+k] = t;
			}
			m = pivot[j];
			pivot[j] = pivot[p];
			pivot[p] = m;
		}
		if (abs(A[j*N+j]) != 0)
		{
			// Eliminação de Gauss
			r = 1 / A[j*N+j];
			for (i = j+1; i < N; i++)
			{
				Mult = A[i*N+j]*r;
				A[i*N+j] = Mult;
				for (k = j+1; k < N; k++)
					A[i*N+k] = A[i*N+k] - Mult*A[j*N+k];
			}
		}
	}
  // 2 PASSO Resolver a matriz da decomposicao
  double *y = calloc(N,sizeof(double));
  double *x = calloc(N,sizeof(double));
	double soma;
	k = pivot[0];
	y[0] = b[k];
	// Realizar substituições sucessivas para resolver o sistema triangular inferior: Ly = b
	for (i = 1; i < N; i++)
	{
		soma = 0;
		for (j = 0; j <= i-1; j++)
			soma += A[i*N+j]*y[j];
		k = pivot[i];
		y[i] = b[k] - soma;
	}
	// Realizar substituições retroativas para resolver o sistema triangular superior: Ux = y
	x[N-1] = y[N-1] / A[(N-1)*N+(N-1)];
	for (i = N-2; i >= 0; i--)
	{
		soma = 0;
		for (j = i+1; j < N; j++)
			soma += A[i*N+j]*x[j];
		x[i] = (y[i] - soma) / A[i*N+i];
	}

  // Libera memoria
	free(y);
  free(pivot);

  printf("ok\n");
  return x;
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
    norm += pow(b[i]-sum,2);
  }
  printf("Norma do sistema = %e\n",sqrt(norm));
  return sqrt(norm);
}
