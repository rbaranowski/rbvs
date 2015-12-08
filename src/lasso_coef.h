/* Rafal Baranowski 7 Dec 2015 Public Domain
 based on the code of Patrick Breheny
 */

#ifndef LASSO_COEF_H
#define LASSO_COEF_H

#include "rbvs.h"
#include "utility_functions.h"
#include "penalised_common.h"


#define LASSO_INFO_OK 0
#define LASSO_INFO_ZERO_SD 100
#define LASSO_INFO_SATURATED_MODEL 101

double lasso_shrinkage(double z, double l1, double l2, double v);
double lasso_cutoff(double lambda1, double lambda2, double alpha);

int lasso_coef_vector_gaussian(double *x, unsigned int n, unsigned int p, double *y, double *coef, unsigned int max_nonzero, double alpha,  double tol, unsigned int max_iterations, double lambda_ratio, unsigned int n_lambda);
int lasso_coef_vector_binomial(double *x, int n, int p, double *y, double *coef, int max_nonzero, double alpha,  double tol, int max_iterations, double lambda_ratio, int n_lambda);
SEXP lasso_coef_gaussian_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero, SEXP alpha, SEXP tol, SEXP max_iterations, SEXP lambda_ratio, SEXP n_lambda);
SEXP lasso_coef_binomial_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero, SEXP alpha, SEXP tol, SEXP max_iterations, SEXP lambda_ratio, SEXP n_lambda);

#endif
