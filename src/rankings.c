/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "rankings.h"

int omega_with_id_t_cmp(const void *aa, const void *bb)
{
	const struct omega_with_id_t *a = (struct omega_with_id_t *)aa;
	const struct omega_with_id_t *b = (struct omega_with_id_t *)bb;

	double tmp = a[0].omega - b[0].omega;

	if (fabs(tmp) > DBL_EPSILON) {
		if (tmp > 0)
			return -1;
		else
			return 1;
	} else {
		tmp = a[0].x - b[0].x;
		if (tmp > 0)
			return -1;
		else if (tmp == 0)
			return 0;
		else
			return 1;
	}
}

struct omega_with_id_t qselect_omega(struct omega_with_id_t *x, unsigned int p,
				     unsigned int k)
{

#define SWAP(a, b) { tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

	struct omega_with_id_t tmp;

	unsigned int st = 0;
	register unsigned int i;

	for (i = 0; i < p - 1; i++) {
		if (omega_with_id_t_cmp(&x[i], &x[p - 1]) > 0)
			continue;
		SWAP(i, st);
		st++;
	}

	SWAP(p - 1, st);

	return k == st ? x[st]
	    : st > k ? qselect_omega(x, st, k)
	    : qselect_omega(&x[st], p - st, k - st);
}

void omega_order(double *omega, unsigned int p, unsigned int *index,
		 unsigned int k_max, struct omega_with_id_t *omega_with_id)
{

	register unsigned int j = 0;

	for (j = 0; j < p; j++) {

		omega_with_id[j].id = j + 1;
		omega_with_id[j].omega = omega[j];
		omega_with_id[j].x = unif_rand();

	}

	///partially sort k_max first elements
	if (k_max < p)
		qselect_omega(omega_with_id, p, k_max);

	qsort(omega_with_id, k_max, sizeof(struct omega_with_id_t),
	      omega_with_id_t_cmp);

	for (j = 0; j < k_max; j++)
		index[j] = omega_with_id[j].id;

}

void rankings(double *omega, unsigned int p, unsigned int B,
	      unsigned int *ranks, unsigned int k_max)
{

	GetRNGstate();

  #pragma omp parallel
	{

		struct omega_with_id_t *omega_with_id = Calloc(p, struct omega_with_id_t);
		unsigned int j;

    #pragma omp for
		for (j = 0; j < B; j++)
			omega_order(&omega[j * p], p, &ranks[j * k_max], k_max,
				    omega_with_id);
    
		Free(omega_with_id);

	}

	PutRNGstate();

}

SEXP rankings_r(SEXP omega, SEXP k_max)
{

	SEXP omega_dim, ranks;

	PROTECT(omega_dim = getAttrib(omega, R_DimSymbol));

	unsigned int p = INTEGER(omega_dim)[0];
	unsigned int B = INTEGER(omega_dim)[1];
	unsigned int val_kmax = INTEGER(k_max)[0];

	if (val_kmax > p)
		val_kmax = p;

	PROTECT(ranks = allocMatrix(INTSXP, val_kmax, B));

	rankings(REAL(omega), p, B, (unsigned int *)INTEGER(ranks), val_kmax);

	UNPROTECT(2);

	return ranks;

}
