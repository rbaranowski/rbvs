/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "utility_functions.h"

SEXP standardise_matrix_r(SEXP x, SEXP scale) {
    SEXP x_dim;

    PROTECT(x_dim = getAttrib(x, R_DimSymbol));

    standardise_matrix(REAL(x), INTEGER(x_dim)[0], INTEGER(x_dim)[1],
                       INTEGER(scale)[0]);

    UNPROTECT(1);

    return (x);
}

void standardise_matrix(double *x, unsigned int n, unsigned int p, int scale) {

    register unsigned int i, j;
    double *xj;

    double mean;
    if (scale) {

        double ss;
        double tmp;

        for (j = 0; j < p; j++) {
            mean = 0;
            xj = &x[j * n];

            for (i = 0; i < n; i++)
                mean += xj[i];

            mean /= n;
            ss = 0;

            for (i = 0; i < n; i++) {
                xj[i] -= mean;
                tmp = xj[i];
                tmp *= tmp;
                ss += tmp;
            }

            ss = sqrt(ss);

            if (ss > DBL_EPSILON)
                for (i = 0; i < n; i++)
                    xj[i] /= ss;

        }
    } else {

        for (j = 0; j < p; j++) {

            mean = 0;
            xj = &x[j * n];

            for (i = 0; i < n; i++)
                mean += xj[i];

            mean /= n;

            for (i = 0; i < n; i++)
                xj[i] -= mean;

        }
    }
}

void select_subsample(unsigned int *subsample, unsigned int m, double *x_src,
                      unsigned int n, unsigned int p, double *y_src,
                      double *x_dst, double *y_dest) {

    register unsigned int i, j, id;

    for (i = 1; i <= m; i++) {

        id = subsample[i - 1];
        y_dest[i - 1] = y_src[id - 1];

        for (j = 1; j <= p; j++)
            x_dst[IDX(i, j, m)] = x_src[IDX(id, j, n)];

    }
}

double cross_product(double *x, double *y, unsigned int n) {
    double val = 0.0;
    register unsigned int i;

    for (i = 0; i < n; i++)
        val += x[i] * y[i];
    return val;
}

double weighted_cross_product(double *x, double *y, double *w, unsigned int n) {
    double val = 0.0;
    register unsigned int i;
    for (i = 0; i < n; i++)
        val += x[i] * y[i] * w[i];
    return val;
}

double sum(double *x, unsigned int n) {
    double val = 0.0;
    register unsigned int i;
    for (i = 0; i < n; i++)
        val += x[i];
    return val;
}

double weighted_sum_squares(double *x, double *w, unsigned int n) {
    double val = 0.0;
    register unsigned int i;
    for (i = 0; i < n; i++)
        val += (x[i] * x[i]) * w[i];
    return val;
}

double qselect_double(double *x, unsigned int p, unsigned int k) {

    #define SWAP(a, b) { tmp = x[a]; x[a] = x[b]; x[b] = tmp; }

    double tmp;
    register unsigned int st = 0, i;

    for (i = 0; i < p - 1; i++) {
        if (x[i] < x[p - 1])
            continue;
        SWAP(i, st);
        st++;
    }

    SWAP(p - 1, st);

    return k == st ? x[st]
                   : st > k ? qselect_double(x, st, k)
                            : qselect_double(&x[st], p - st, k - st);
}

int is_in_array(unsigned int *sorted, unsigned int element, unsigned int length) {

    unsigned int *item =
            (unsigned int *) bsearch(&element, sorted, length,
                                     sizeof(unsigned int), unsigned_int_cmp);
    if (item != NULL)
        return 1;
    else
        return 0;

}

int unsigned_int_cmp(const void *aa, const void *bb) {

    const unsigned int *a = aa, *b = bb;
    return (*a < *b) ? -1 : (*a > *b);

}

unsigned int ceil_div(unsigned int a, unsigned int b) {
    return a > 0 ? 1 + ((a - 1) / b) : 0;
}
