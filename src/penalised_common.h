/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef PENALISED_COMMON_H
#define PENALISED_COMMON_H

#include "rbvs.h"
#include "utility_functions.h"

void update_residuals(double *residuals, double *x, unsigned int n, double multiplier);
void select_lambda(double *abs_cross_product, unsigned int p, double alpha, double lambda_ratio, double *lambda,  unsigned int n_lambda, unsigned int max_nonzero);

#endif
