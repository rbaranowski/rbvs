/* Rafal Baranowski 7 Dec 2015 Public Domain
 based on the code of Patrick Breheny
 */
#ifndef MCPLUS_COEF_H
#define MCPLUS_COEF_H

#include "rbvs.h"
#include "utility_functions.h"
#include "penalised_common.h"


#define MCPLUS_INFO_OK 0
#define MCPLUS_INFO_ZERO_SD 100
#define MCPLUS_INFO_SATURATED_MODEL 101

double mcplus_shrinkage(double z, double l1, double l2, double gammma, double v);
double mcplus_cutoff(double lambda1, double lambda2, double alpha, double gamma);

int mcplus_coef_vector_gaussian(double *x, unsigned int n, unsigned int p, double *y, double *coef, unsigned int max_nonzero, double alpha, double gamma, double tol, unsigned int max_iterations, double lambda_ratio, unsigned int n_lambda);
int mcplus_coef_vector_binomial(double *x, int n, int p, double *y, double *coef, int max_nonzero, double alpha, double gamma, double tol, int max_iterations, double lambda_ratio, int n_lambda);
SEXP mcplus_coef_gaussian_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero, SEXP alpha, SEXP gamma, SEXP tol, SEXP max_iterations, SEXP lambda_ratio, SEXP n_lambda);
SEXP mcplus_coef_binomial_r(SEXP subsamples, SEXP x, SEXP y, SEXP max_nonzero, SEXP alpha, SEXP gamma, SEXP tol, SEXP max_iterations, SEXP lambda_ratio, SEXP n_lambda);

#endif
