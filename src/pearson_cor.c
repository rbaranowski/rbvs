/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "pearson_cor.h"

SEXP pearson_cor_r(SEXP subsamples, SEXP x, SEXP y)
{

	SEXP x_dim, subsamples_dim, cor;

	PROTECT(x_dim = getAttrib(x, R_DimSymbol));
	PROTECT(subsamples_dim = getAttrib(subsamples, R_DimSymbol));

	unsigned int n = INTEGER(x_dim)[0];
	unsigned int p = INTEGER(x_dim)[1];
	unsigned int m = INTEGER(subsamples_dim)[0];
	unsigned int B = INTEGER(subsamples_dim)[1];
	unsigned int k;

	PROTECT(cor = allocMatrix(REALSXP, p, B));

	double *ptr_x = REAL(x);
	double *ptr_y = REAL(y);
	double *ptr_cor = REAL(cor);

	unsigned int *ptr_subsamples = (unsigned int *)INTEGER(subsamples);

#pragma omp parallel for
	for (k = 0; k < B; k++)
		pearson_cor_vector(&ptr_subsamples[k * m], m, ptr_x, n, p,
				   ptr_y, &ptr_cor[k * p]);

	UNPROTECT(3);

	return cor;

}

double pearson_cor(unsigned int *rows, unsigned int m, double *x, double *y,
		   double sum_y, double sd_y)
{

	double sum_x = 0.0;
	double sum_x_sq = 0.0;
	double sd_x = 0.0;
	double sum_xy = 0.0;

	register unsigned int i;
	register unsigned int index;

	for (i = 0; i < m; i++) {

		index = rows[i] - 1;

		sum_x += x[index];
		sum_x_sq += (x[index]) * (x[index]);
		sum_xy += (x[index]) * y[index];

	}

	sd_x = sqrt(m * sum_x_sq - sum_x * sum_x);

	if (sd_x > DBL_EPSILON && sd_y > DBL_EPSILON)
		return fabs((m * sum_xy - (sum_x * sum_y)) / (sd_x * sd_y));
	else
		return 0.0;
}

void pearson_cor_vector(unsigned int *rows, unsigned int m, double *x,
			unsigned int n, unsigned int p, double *y, double *cor)
{

	double sd_y = 0.0;
	double sum_y = 0.0;
	double sum_y_sq = 0.0;

	register unsigned int i, j;

	//find std of y
	for (i = 0; i < m; i++) {
		j = rows[i] - 1;
		sum_y += y[j];
		sum_y_sq += y[j] * y[j];

	}

	sd_y = sqrt(m * sum_y_sq - sum_y * sum_y);

	//find correlations  
	if (sd_y > DBL_EPSILON)
		for (j = 0; j < p; j++)
			cor[j] =
			    pearson_cor(rows, m, &x[j * n], y, sum_y, sd_y);
	else
		for (j = 0; j < p; j++)
			cor[j] = 0.0;

}
