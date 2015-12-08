/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include  "rbvs.h"

void factor_model_row(double *x, unsigned int row, unsigned int n,
		      unsigned int p, unsigned int n_factors, double sigma)
{
	register unsigned int j, k, l;

	double factor = 0.0;

	l = row;
	for (j = 0; j < p; j++) {
		x[l] = 0.0;
		l += n;
	}

	for (k = 0; k < n_factors; k++) {

		factor = rnorm(0.0, sigma);
		l = row;
		for (j = 0; j < p; j++) {
			x[l] += factor * rnorm(0.0, sigma);
			l += n;
		}

	}

	l = row;

	for (j = 0; j < p; j++) {
		x[l] += rnorm(0.0, sigma);
		l += n;
	}

}

SEXP factor_model_r(SEXP n, SEXP p, SEXP n_factors, SEXP sigma)
{

	unsigned int val_n = INTEGER(n)[0];
	unsigned int val_p = INTEGER(p)[0];
	unsigned int val_n_factors = INTEGER(n_factors)[0];
	register unsigned int i;

	double val_sigma = REAL(sigma)[0];

	SEXP x = PROTECT(allocMatrix(REALSXP, val_n, val_p));

	double *ptr_x = REAL(x);

	GetRNGstate();
#pragma omp parallel for
	for (i = 0; i < val_n; i++) {
		factor_model_row(ptr_x, i, val_n, val_p, val_n_factors,
				 val_sigma);
	}

	PutRNGstate();

	UNPROTECT(1);

	return (x);
}
