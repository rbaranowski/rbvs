/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef UTILITY_FUNCTIONS_H
#define UTILITY_FUNCTIONS_H
#include "rbvs.h"

SEXP matrix_projection_in_place_r(SEXP x, SEXP projection, SEXP active, SEXP use_gpu);
SEXP standardise_matrix_r(SEXP x, SEXP scale);

void matrix_projection_in_place(double *x_new, double *x, unsigned int n, unsigned int p, double *projection, unsigned int *active, unsigned int n_active);
void standardise_matrix(double *x, unsigned int n, unsigned int p, int scale);
void select_subsample(unsigned int *subsample, unsigned int m, double *x_src, unsigned int n, unsigned int p, double *y_src, double *x_dst, double *y_dest);
double cross_product(double *x, double *y, unsigned int n);
double weighted_cross_product(double *x, double *y, double *w, unsigned int n);
double sum(double *x, unsigned int n); 
double weighted_sum_squares(double *x, double *w, unsigned int n);
double qselect_double(double *x, unsigned int p, unsigned int k);
int is_in_array(unsigned int *sorted, unsigned int element, unsigned int length);
int unsigned_int_cmp(const void *aa, const void *bb);

#endif
