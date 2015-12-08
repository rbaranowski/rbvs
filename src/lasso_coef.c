/* Rafal Baranowski 7 Dec 2015 Public Domain
  based on the code of Patrick Breheny
 */
#include "lasso_coef.h"

SEXP lasso_coef_gaussian_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero,
			   SEXP alpha, SEXP tol, SEXP max_iterations,
			   SEXP lambda_ratio, SEXP n_lambda)
{

	SEXP x_dim, subsamples_dim, coef;

	PROTECT(x_dim = getAttrib(x, R_DimSymbol));
	PROTECT(subsamples_dim = getAttrib(subsamples, R_DimSymbol));

	register unsigned int n = INTEGER(x_dim)[0];
	register unsigned int p = INTEGER(x_dim)[1];
	register unsigned int m = INTEGER(subsamples_dim)[0];
	register unsigned int B = INTEGER(subsamples_dim)[1];
	register unsigned int k;

	PROTECT(coef = allocMatrix(REALSXP, p, B));

	double *ptr_x = REAL(x);
	double *ptr_y = REAL(y);
	double *ptr_coef = REAL(coef);
	unsigned int *ptr_subsamples = (unsigned int *)INTEGER(subsamples);

	//values 
	unsigned int val_max_nonzero = *INTEGER(max_nonzero);
	double val_alpha = *REAL(alpha);
	double val_tol = *REAL(tol);
	unsigned int val_max_iterations = *INTEGER(max_iterations);
	double val_lambda_ratio = *REAL(lambda_ratio);
	unsigned int val_n_lambda = *INTEGER(n_lambda);

	//allocate matrix for subsampled data  
#pragma omp parallel
	{

		double *ptr_x_subsample = Calloc(p * m, double);
		double *ptr_y_subsample = Calloc(m, double);

#pragma omp for
		for (k = 0; k < B; k++) {
			// select subsample
			select_subsample(&ptr_subsamples[k * m], m, ptr_x, n, p,
					 ptr_y, ptr_x_subsample,
					 ptr_y_subsample);
			//evaluate lasso coefficients
			lasso_coef_vector_gaussian(ptr_x_subsample, m, p,
						   ptr_y_subsample,
						   &ptr_coef[IDX(1, k + 1, p)],
						   val_max_nonzero, val_alpha,
						   val_tol, val_max_iterations,
						   val_lambda_ratio,
						   val_n_lambda);
		}

		Free(ptr_x_subsample);
		Free(ptr_y_subsample);

	}

	UNPROTECT(3);

	return coef;

}

SEXP lasso_coef_binomial_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero,
			   SEXP alpha, SEXP tol, SEXP max_iterations,
			   SEXP lambda_ratio, SEXP n_lambda)
{

	SEXP x_dim, subsamples_dim, coef;

	PROTECT(x_dim = getAttrib(x, R_DimSymbol));
	PROTECT(subsamples_dim = getAttrib(subsamples, R_DimSymbol));

	unsigned int n = INTEGER(x_dim)[0];
	unsigned int p = INTEGER(x_dim)[1];
	unsigned int m = INTEGER(subsamples_dim)[0];
	unsigned int B = INTEGER(subsamples_dim)[1];
	unsigned int k;

	PROTECT(coef = allocMatrix(REALSXP, p, B));

	double *ptr_x = REAL(x);
	double *ptr_y = REAL(y);
	double *ptr_coef = REAL(coef);
	unsigned int *ptr_subsamples = (unsigned int *)INTEGER(subsamples);

	//values 
	unsigned int val_max_nonzero = *INTEGER(max_nonzero);
	double val_alpha = *REAL(alpha);
	double val_tol = *REAL(tol);
	unsigned int val_max_iterations = *INTEGER(max_iterations);
	double val_lambda_ratio = *REAL(lambda_ratio);
	unsigned int val_n_lambda = *INTEGER(n_lambda);

#pragma omp parallel
	{
		//allocate matrix for subsampled data  
		double *ptr_x_subsample = Calloc(p * m, double);
		double *ptr_y_subsample = Calloc(m, double);

#pragma omp for
		for (k = 0; k < B; k++) {
			// select subsample
			select_subsample(&ptr_subsamples[k * m], m, ptr_x, n, p,
					 ptr_y, ptr_x_subsample,
					 ptr_y_subsample);
			//evaluate lasso coefficients
			lasso_coef_vector_binomial(ptr_x_subsample, m, p,
						   ptr_y_subsample,
						   &ptr_coef[k * p],
						   val_max_nonzero, val_alpha,
						   val_tol, val_max_iterations,
						   val_lambda_ratio,
						   val_n_lambda);
		}

		Free(ptr_x_subsample);
		Free(ptr_y_subsample);

	}
	UNPROTECT(3);

	return coef;

}

double lasso_shrinkage(double z, double l1, double l2, double v)
{
	double s = 0.0;
	if (z > 0)
		s = 1.0;
	else if (z < 0)
		s = -1.0;
	if (fabs(z) <= l1)
		return (0.0);
	else
		return (s * (fabs(z) - l1) / (v * (1 + l2)));
}

double lasso_cutoff(double lambda1, double lambda2, double alpha)
{
	return (2 * lambda2 - lambda1) * alpha;
}

int lasso_coef_vector_gaussian(double *x, unsigned int n, unsigned int p,
			       double *y, double *coef,
			       unsigned int max_nonzero, double alpha,
			       double tol, unsigned int max_iterations,
			       double lambda_ratio, unsigned int n_lambda)
{

	// declare variables used in the program
	register unsigned int iterations = 0, i = 0, j = 0, l = 0, n_nonzero =
	    0, active_len = 0, active_strong_len = 0, converged =
	    0, n_violations = 0;
	double cutoff = 0.0;
	double l1, l2, shift;
	double *xj;
	double n_inv = 1.0 / ((double)n);

	// variables requiring allocation
	double *tmp_coef = Calloc(p, double);
	double *z = Calloc(p, double);
	double *residuals = Calloc(n, double);
	unsigned int *active = Calloc(p, unsigned int);
	unsigned int *active_strong = Calloc(p, unsigned int);

	double *lambda = Calloc(n_lambda, double);

	// standardise x and y
	standardise_matrix(x, n, p, 1);
	standardise_matrix(y, n, 1, 0);

	//setting initial values

	for (j = 0; j < p; j++) {

		tmp_coef[j] = 0.0;
		coef[j] = 0.0;
		z[j] = fabs(cross_product(&x[j * n], y, n)) * n_inv;

	}

	for (i = 0; i < n; i++)
		residuals[i] = y[i];

	// Determine lambdas
	select_lambda(z, p, alpha, lambda_ratio, lambda, n_lambda, max_nonzero);

	for (j = 0; j < p; j++)
		z[j] = cross_product(&x[j * n], residuals, n) * n_inv;

	//start the algorithm 

	l = 1;

	while ((l < n_lambda) && (n_nonzero < max_nonzero)) {

		iterations = 0;
		converged = 0;

		l1 = lambda[l] * alpha;
		l2 = lambda[l] * (1 - alpha);

		// Determine the eligible set
		cutoff = lasso_cutoff(lambda[l - 1], lambda[l], alpha);

		for (j = 0; j < p; j++) {
			coef[j] = 0.0;

			if (fabs(z[j]) > cutoff) {
				if (is_in_array
				    (active_strong, j,
				     active_strong_len) == 0) {
					insertion_sort(active_strong, j,
						       active_strong_len);
					active_strong_len++;
				}
			}
		}

		n_violations = -1;

		while ((n_violations != 0) && (iterations < max_iterations)) {

 LOOP:
			while ((converged == 0)
			       && (iterations < max_iterations)) {

				iterations++;
				converged = 1;

				// Solve over the active set
				for (i = 0; i < active_len; i++) {

					j = active[i];
					xj = &x[j * n];
					z[j] =
					    cross_product(xj, residuals,
							  n) * n_inv +
					    tmp_coef[j];

					// Update coef_j                  
					coef[j] =
					    lasso_shrinkage(z[j], l1, l2, 1.0);

					// Update the residuals
					shift = coef[j] - tmp_coef[j];

					if (shift != 0) {

						if ((converged == 1)
						    && (fabs(shift) >
							(tmp_coef[j] * tol)))
							converged = 0;
						update_residuals(residuals, xj,
								 n, -shift);
						tmp_coef[j] = coef[j];

					}

				}

			}

			//check for violations in the strong set

			n_violations = 0;
			converged = 0;

			for (i = 0; i < active_strong_len; i++)
				if (is_in_array
				    (active, active_strong[i],
				     active_len) == 0) {

					j = active_strong[i];
					xj = &x[j * n];
					z[j] =
					    cross_product(xj, residuals,
							  n) * n_inv;

					// Update coef_j            
					coef[j] =
					    lasso_shrinkage(z[j], l1, l2, 1.0);

					if (coef[j] != 0) {

						insertion_sort(active, j,
							       active_len);
						active_len++;

						update_residuals(residuals, xj,
								 n, -coef[j]);
						tmp_coef[j] = coef[j];
						n_violations++;

					}
				}

			if (n_violations > 0)
				goto LOOP;

			//check for violations in the rest

			n_violations = 0;

			for (j = 0; j < p; j++)
				if (is_in_array
				    (active_strong, j,
				     active_strong_len) == 0) {

					xj = &x[j * n];
					z[j] =
					    cross_product(xj, residuals,
							  n) * n_inv;

					// Update coef_j              
					coef[j] =
					    lasso_shrinkage(z[j], l1, l2, 1.0);

					if (coef[j] != 0) {

						insertion_sort(active, j,
							       active_len);
						active_len++;
						insertion_sort(active_strong, j,
							       active_strong_len);
						active_strong_len++;

						update_residuals(residuals, xj,
								 n, -coef[j]);
						tmp_coef[j] = coef[j];
						n_violations++;

					}
				}

		}

		//check the number of non-zero coefficients
		n_nonzero = 0;
		for (j = 0; j < p; j++) {
			tmp_coef[j] = coef[j];
			if (coef[j] != 0)
				n_nonzero++;
		}

		// set previous lambda to current lambda
		l++;

	}

	Free(tmp_coef);
	Free(active);
	Free(active_strong);
	Free(residuals);
	Free(lambda);
	Free(z);

	return LASSO_INFO_OK;

}

int lasso_coef_vector_binomial(double *x, int n, int p, double *y, double *coef,
			       int max_nonzero, double alpha, double tol,
			       int max_iterations, double lambda_ratio,
			       int n_lambda)
{

	// declare variables used in the program
	register unsigned int iterations = 0, i = 0, j = 0, l = 0, k =
	    0, n_nonzero = 0, converged = 0, n_violations = 0;
	register double n_inv = 1.0 / ((double)n);
	double cutoff = 0.0;
	double l1, l2, shift;
	double *xj;
	double tmp = 0, tmp1 = 0, tmp2 = 0;
	double intercept = 0;
	double tmp_intercept = 0;
	double null_dev = 0;
	double dev = 0;
	double pi = 0, xwr = 0, xwx = 0, si = 0, u = 0, v = 0;

	// variables requiring allocation
	double *tmp_coef = Calloc(p, double);
	double *z = Calloc(p, double);
	double *w = Calloc(n, double);
	double *s = Calloc(n, double);
	double *r = Calloc(n, double);
	double *eta = Calloc(n, double);
	double *lambda = Calloc(n_lambda, double);
	unsigned int *active = Calloc(p, unsigned int);
	unsigned int *active_strong = Calloc(p, unsigned int);

	register unsigned int active_len = 0, active_strong_len = 0;

	// standardise x and y
	standardise_matrix(x, n, p, 1);

	//setting initial values

	double ybar = sum(y, n) * n_inv;

	if ((fabs(ybar) < DBL_EPSILON) || (fabs(ybar - 1) < DBL_EPSILON))
		return LASSO_INFO_ZERO_SD;

	intercept = tmp_intercept = log(ybar / (1 - ybar));

	null_dev = 0.0;
	tmp1 = log(ybar);
	tmp2 = log(1 - ybar);

	for (i = 0; i < n; i++) {

		if (y[i] == 1)
			null_dev -= tmp1;
		else
			null_dev -= tmp2;

		s[i] = y[i] - ybar;
		eta[i] = tmp_intercept;

	}

	for (j = 0; j < p; j++) {

		tmp_coef[j] = 0.0;
		coef[j] = 0.0;
		z[j] = fabs(cross_product(&x[j * n], s, n)) / n;

	}

	// Determine lambdas
	select_lambda(z, p, alpha, lambda_ratio, lambda, n_lambda, max_nonzero);

	for (j = 0; j < p; j++)
		z[j] = cross_product(&x[j * n], s, n) / n;

	//start the algorithm 

	l = 1;

	while ((l < n_lambda) && (n_nonzero < max_nonzero)) {

		iterations = 0;
		converged = 0;

		l1 = lambda[l] * alpha;
		l2 = lambda[l] * (1 - alpha);

		// Determine the eligible set
		cutoff = lasso_cutoff(lambda[l - 1], lambda[l], alpha);

		for (j = 0; j < p; j++) {
			coef[j] = 0.0;

			if (fabs(z[j]) > cutoff) {
				if (is_in_array
				    (active_strong, j,
				     active_strong_len) == 0) {
					insertion_sort(active_strong, j,
						       active_strong_len);
					active_strong_len++;

				}
			}
		}

		n_violations = -1;

		while ((n_violations != 0) && (iterations < max_iterations)) {

 LOOP:
			converged = 0;

			while ((converged == 0)
			       && (iterations < max_iterations)) {

				iterations++;
				converged = 1;
				dev = 0.0;

				for (i = 0; i < n; i++) {

					if (eta[i] > 10) {
						pi = 1;
						w[i] = .0001;
					} else if (eta[i] < -10) {
						pi = 0;
						w[i] = .0001;
					} else {
						pi = exp(eta[i]) / (1 +
								    exp(eta
									[i]));
						w[i] = pi * (1 - pi);
					}

					s[i] = y[i] - pi;
					r[i] = s[i] / w[i];

					if (y[i] == 1)
						dev -= log(pi);
					else
						dev -= log(1 - pi);

				}

				if (dev / null_dev < .01)
					return LASSO_INFO_SATURATED_MODEL;

				// Intercept
				xwr = cross_product(w, r, n);
				xwx = sum(w, n);

				tmp = xwr / xwx;
				intercept = tmp + tmp_intercept;

				for (k = 0; k < n; k++) {

					r[k] -= tmp;
					eta[k] += tmp;

				}

				// Covariates - solve over the active set

				for (i = 0; i < active_len; i++) {

					j = active[i];
					xj = &x[j * n];

					// Calculate u, v
					/* xwr = weighted_cross_product(xj, r, w, n);
					   xwx = weighted_sum_squares(xj, w, n); */

					xwr = 0.0;
					xwx = 0.0;

					for (k = 0; k < n; k++) {
						tmp = xj[k] * w[k];
						xwr += tmp * r[k];
						xwx += tmp * xj[k];
					}

					u = xwr / n + (xwx / n) * tmp_coef[j];
					v = xwx / n;

					// Update coef_j            
					coef[j] = lasso_shrinkage(u, l1, l2, v);

					// Update the residuals
					shift = coef[j] - tmp_coef[j];

					if (shift != 0) {

						if ((converged == 1)
						    && (fabs(shift) >
							(tmp_coef[j] * tol)))
							converged = 0;

						for (k = 0; k < n; k++) {
							si = shift * xj[k];
							r[k] -= si;
							eta[k] += si;
						}

						tmp_coef[j] = coef[j];

					}

				}

				tmp_intercept = intercept;

			}

			//check for violations in the strong set

			n_violations = 0;

			for (i = 0; i < active_strong_len; i++)
				if (is_in_array
				    (active, active_strong[i],
				     active_len) == 0) {

					j = active_strong[i];
					xj = &x[j * n];
					z[j] = cross_product(xj, r, n) / n;

					if (fabs(z[j]) > l1) {

						insertion_sort(active, j,
							       active_len);
						active_len++;
						n_violations++;

					}

				}

			if (n_violations > 0)
				goto LOOP;

			//check for violations in the rest

			n_violations = 0;

			for (k = 0; k < p; k++)
				if (is_in_array
				    (active_strong, k,
				     active_strong_len) == 0) {

					j = k;
					xj = &x[j * n];
					z[j] = cross_product(xj, r, n) / n;

					if (fabs(z[j]) > l1) {

						insertion_sort(active, j,
							       active_len);
						active_len++;
						insertion_sort(active_strong, j,
							       active_strong_len);
						active_strong_len++;

						n_violations++;

					}

				}

		}

		//check the number of non-zero coefficients
		n_nonzero = 0;
		for (j = 0; j < p; j++) {
			tmp_coef[j] = coef[j];
			if (coef[j] != 0)
				n_nonzero++;
		}

		// set previous lambda to current lambda
		l++;

	}

	Free(tmp_coef);
	Free(active);
	Free(active_strong);
	Free(lambda);
	Free(z);
	Free(w);
	Free(r);
	Free(s);
	Free(eta);

	return LASSO_INFO_OK;

}
