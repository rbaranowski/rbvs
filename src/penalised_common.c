/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "penalised_common.h"

void update_residuals(double *residuals, double *x, unsigned int n,
		      double multiplier)
{
	register unsigned int i;
	for (i = 0; i < n; i++)
		residuals[i] += x[i] * multiplier;
}

void select_lambda(double *abs_cross_product, unsigned int p, double alpha,
		   double lambda_ratio, double *lambda, unsigned int n_lambda,
		   unsigned int max_nonzero)
{

	double tmp;
	register unsigned int j;

	tmp = qselect_double(abs_cross_product, p, max_nonzero);

	lambda[0] = qselect_double(abs_cross_product, max_nonzero, 0);

	lambda[1] = (lambda[0] + tmp / alpha) / 2;

	//determine remaining lambdas

	lambda[n_lambda - 1] = lambda[1] * lambda_ratio;

	double log_step =
	    fabs(log(lambda[n_lambda - 1]) - log(lambda[1])) / n_lambda;
	tmp = log(lambda[1]);

	for (j = 1; j < n_lambda; j++) {
		lambda[j] = exp(tmp);
		tmp -= log_step;
	}

}
