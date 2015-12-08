/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef RANKINGS_H
#define RANKINGS_H

#include "rbvs.h"

struct omega_with_id_t{
  double omega;
  double x; //randomly assigned to solve ties
  unsigned int id; 
};


int omega_with_id_t_cmp(const void *aa, const void *bb);
struct omega_with_id_t qselect_omega(struct omega_with_id_t *x, unsigned int p, unsigned int k);
void omega_order(double *omega, unsigned int p, unsigned int *index, unsigned int k_max, struct omega_with_id_t *omega_with_id);

SEXP rankings_r(SEXP omega,SEXP k_max);


#endif
