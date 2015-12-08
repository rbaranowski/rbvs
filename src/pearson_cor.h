/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef PEARSON_COR_H
#define PEARSON_COR_H

#include "rbvs.h"
#include <Rmath.h>


SEXP pearson_cor_r(SEXP subsamples, SEXP x, SEXP y);
double pearson_cor(unsigned int *rows, unsigned int m, double *x, double *y, double sum_y, double sd_y);
void pearson_cor_vector(unsigned int *rows,unsigned int m, double *x, unsigned int n, unsigned int p, double *y, double *cor);

#endif 
