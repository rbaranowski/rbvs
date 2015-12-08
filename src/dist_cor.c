/* Rafal Baranowski 7 Dec 2015 Public Domain */

#include "dist_cor.h"

SEXP distance_cor_r(SEXP subsamples, SEXP x, SEXP y, SEXP index)
{

	SEXP x_dim, subsamples_dim, cor;

	PROTECT(x_dim = getAttrib(x, R_DimSymbol));
	PROTECT(subsamples_dim = getAttrib(subsamples, R_DimSymbol));

	unsigned int n = INTEGER(x_dim)[0];
	unsigned int p = INTEGER(x_dim)[1];
	unsigned int m = INTEGER(subsamples_dim)[0];
	unsigned int B = INTEGER(subsamples_dim)[1];

	PROTECT(cor = allocMatrix(REALSXP, p, B));

	int *ptr_subsamples = INTEGER(subsamples);
	double *ptr_x = REAL(x);
	double *ptr_y = REAL(y);
	double *ptr_cor = REAL(cor);

	double *ptr_y_dist = Calloc((n * (n + 1)) / 2, double);
	double *ptr_x_dist = Calloc((n * (n + 1)) / 2, double);

	//double *val_index = Calloc(B,double);
	double val_index = 1.0;

	//for(i=0; i < B; i++) val_index[i] = 2*unif_rand();

	//compute the distance matrix for y
	distances(ptr_y, ptr_y_dist, n);

#pragma omp parallel
	{

		double *ptr_A = Calloc((m * (m + 1)) / 2, double);
		double *ptr_B = Calloc((m * (m + 1)) / 2, double);
		int *ptr_rows;
		register unsigned int i, j;

		for (j = 1; j <= p; j++) {

			/// find the distance matrix for x[,j]    
			distances(&ptr_x[IDX(1, j, n)], ptr_x_dist, n);

#pragma omp for
			for (i = 1; i <= B; i++) {

				//select subsample      
				ptr_rows = &ptr_subsamples[IDX(1, i, m)];
				//compute A matrix      
				A_matrix(ptr_y_dist, n, ptr_rows, m, ptr_A,
					 val_index);
				//compute B matrix
				A_matrix(ptr_x_dist, n, ptr_rows, m, ptr_B,
					 val_index);

				// find the correlation
				ptr_cor[IDX(j, i, p)] =
				    dist_cor(ptr_A, ptr_B, m);

			}
		}

		Free(ptr_A);
		Free(ptr_B);
	}

	Free(ptr_y_dist);
	Free(ptr_x_dist);

	//Free(val_index);

	UNPROTECT(3);

	return cor;

}

double dist_cor(double *A, double *B, unsigned int n)
{
	register unsigned int i, j, id;

	double cov_AB = 0.0;
	double var_A = 0.0;
	double var_B = 0.0;

	for (j = 1; j <= n; j++) {

		id = IDXUT(j, j, n);
		cov_AB += A[id] * B[id];
		var_A += A[id] * A[id];
		var_B += B[id] * B[id];

		for (i = j + 1; i <= n; i++) {
			id = IDXUT(i, j, n);

			cov_AB += 2 * A[id] * B[id];
			var_A += 2 * A[id] * A[id];
			var_B += 2 * B[id] * B[id];

		}
	}

	if ((var_A < DBL_EPSILON) || (var_B < DBL_EPSILON))
		return 0.0;
	else
		return sqrt(cov_AB / sqrt(var_A * var_B));
}

/*
SEXP distances_r(SEXP x, SEXP index){
  int n = length(x);
  double val_index = REAL(index)[0];
  
  SEXP dist;
  
  PROTECT(dist = allocVector(REALSXP,(n *(n+1))/2));
  
  distances(REAL(x), REAL(dist), n, val_index);
  
  UNPROTECT(1);
  
  return dist;
  
}*/

void distances(double *x, double *distances, unsigned int n)
{

	unsigned int i, j;

	for (j = 1; j <= n; j++) {
		distances[IDXUT(j, j, n)] = 0.0;
		for (i = j + 1; i <= n; i++) {
			distances[IDXUT(i, j, n)] = fabs(x[i - 1] - x[j - 1]);
		}
	}

}

void A_matrix(double *distances, unsigned int n, int *rows, unsigned int m,
	      double *A, double index)
{

	register unsigned int i, j, id;
	double *A_means = Calloc(m, double);

	double A_mean = 0.0;
	double m_inv = 1.0 / ((double)m);
	if (fabs(index - 1) > DBL_EPSILON) {
		for (j = 1; j <= m; j++) {

			A_means[j - 1] = 0.0;
			id = rows[j - 1];

			for (i = 1; i <= m; i++) {
				A_means[j - 1] +=
				    distances[IDXUT(rows[i - 1], id, n)];
			}

			A_mean += A_means[j - 1];
			A_means[j - 1] *= m_inv;
		}

		A_mean /= (double)(m * m);

		for (j = 1; j <= m; j++) {
			for (i = j; i <= m; i++) {

				A[IDXUT(i, j, m)] =
				    distances[IDXUT
					      (rows[i - 1], rows[j - 1],
					       n)] - A_means[i - 1] -
				    A_means[j - 1] + A_mean;
			}

		}
	} else {
		for (j = 1; j <= m; j++) {

			A_means[j - 1] = 0.0;
			id = rows[j - 1];

			for (i = 1; i <= m; i++) {
				A_means[j - 1] +=
				    distances[IDXUT(rows[i - 1], id, n)];
			}

			A_mean += A_means[j - 1];
			A_means[j - 1] *= m_inv;
		}

		A_mean /= (double)(m * m);

		for (j = 1; j <= m; j++) {
			for (i = j; i <= m; i++) {

				A[IDXUT(i, j, m)] =
				    distances[IDXUT
					      (rows[i - 1], rows[j - 1],
					       n)] - A_means[i - 1] -
				    A_means[j - 1] + A_mean;
			}

		}
	}

	Free(A_means);

}
