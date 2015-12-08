/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef DIST_COR_H
#define DIST_COR_H

#include "rbvs.h"
#include <Rmath.h>

SEXP distance_cor_r(SEXP subsamples, SEXP x, SEXP y, SEXP index);

void distances(double *x, double *distances,  unsigned int n);
void A_matrix(double *distances,  unsigned int n, int *rows,  unsigned int m, double *A,  double index);
double dist_cor(double *A, double *B,  unsigned int n);

#endif

