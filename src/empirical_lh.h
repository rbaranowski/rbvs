/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef PEARSON_COR_H
#define PEARSON_COR_H

#include "rbvs.h"
#include <Rmath.h>

SEXP empirical_lh_r(SEXP subsamples, SEXP x, SEXP y);

double empirical_lh(unsigned int *rows, unsigned int m, double *x, double *y);

void empirical_lh_vector(unsigned int *rows, unsigned int m, double *x, unsigned int n, unsigned int p, double *y,
                         double *cor);

double lagrangian(double lambda, double *products, unsigned int n);

double lagrangian_root(double *products, unsigned int n);

double
lagrangian_root_rec(double left_lambda, double left_value, double right_lambda, double right_value, double *products,
                    unsigned int n);

double elh(double lambda, const double *products, unsigned int n);

#endif 
