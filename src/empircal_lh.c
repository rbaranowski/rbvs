/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "empirical_lh.h"

SEXP empirical_lh_r(SEXP subsamples, SEXP x, SEXP y) {

    SEXP x_dim, subsamples_dim, elh_mat;

    PROTECT(x_dim = getAttrib(x, R_DimSymbol));
    PROTECT(subsamples_dim = getAttrib(subsamples, R_DimSymbol));

    unsigned int n = (unsigned int) INTEGER(x_dim)[0];
    unsigned int p = (unsigned int) INTEGER(x_dim)[1];
    unsigned int m = (unsigned int) INTEGER(subsamples_dim)[0];
    unsigned int B = (unsigned int) INTEGER(subsamples_dim)[1];
    unsigned int k;

    PROTECT(elh_mat = allocMatrix(REALSXP, p, B));

    double *ptr_x = REAL(x);
    double *ptr_y = REAL(y);
    double *ptr_elh_mat = REAL(elh_mat);

    unsigned int *ptr_subsamples = (unsigned int *) INTEGER(subsamples);

    #pragma omp parallel for
    for (k = 0; k < B; k++)
        empirical_lh_vector(&ptr_subsamples[k * m], m, ptr_x, n, p,
                            ptr_y, &ptr_elh_mat[k * p]);

    UNPROTECT(3);

    return elh_mat;

}

double empirical_lh(unsigned int *rows, unsigned int m, double *x, double *y) {

    unsigned int i, index;

    double *_x = Calloc(m, double);
    double *_y = Calloc(m, double);
    double *products = Calloc(m, double);
    double result;

    double mean_x = 0.0;
    double mean_y = 0.0;
    double std_x = 0.0;

    for (i = 0; i < m; i++) {
        index = rows[i] - 1;
        _x[i] = x[index];
        _y[i] = y[index];

        mean_x += _x[i];
        mean_y += _y[i];
    }

    mean_x /= (double) m;
    mean_y /= (double) m;

    for (i = 0; i < m; i++) {
        _x[i] -= mean_x;
        _y[i] -= mean_y;
        std_x += _x[i] * _x[i];
    }

    std_x = sqrt(std_x / ((double) m - 1.0));

    if (std_x > DBL_EPSILON) {
        for (i = 0; i < m; i++) products[i] = _y[i] * _x[i] / std_x;
        double lambda = lagrangian_root(products, m);

        if ISNAN(lambda) result = R_NegInf;
        else result = elh(lambda, products, m);

    } else {
        result = R_NegInf;
    }

    Free(_x);
    Free(_y);
    Free(products);

    return result;
}

void empirical_lh_vector(unsigned int *rows, unsigned int m, double *x,
                         unsigned int n, unsigned int p, double *y, double *elh_vec) {

    double sd_y = 0.0;
    double sum_y = 0.0;
    double sum_y_sq = 0.0;

    register unsigned int i, j;

    //find std of y
    for (i = 0; i < m; i++) {
        j = rows[i] - 1;
        sum_y += y[j];
        sum_y_sq += y[j] * y[j];

    }

    sd_y = sqrt(m * sum_y_sq - sum_y * sum_y);

    //find correlations
    if (sd_y > DBL_EPSILON)
        for (j = 0; j < p; j++)
            elh_vec[j] =
                    empirical_lh(rows, m, &x[j * n], y);
    else
        for (j = 0; j < p; j++)
            elh_vec[j] = 0.0;

}

double lagrangian(double lambda, double *products, unsigned int n) {
    unsigned int i = 0;
    double result = 0.0;

    for (i = 0; i < n; i++) {
        result += (products[i]) / (1.0 + lambda * products[i]);
    }

    return result;
}

double lagrangian_root(double *products, unsigned int n) {
    unsigned int i = 0;
    double left_lambda = R_NegInf;
    double right_lambda = R_PosInf;
    double left_value, right_value;

    for (i = 0; i < n; i++) {
        if ((products[i] > 0) && (left_lambda < products[i])) {
            left_lambda = products[i];
        } else if ((products[i] < 0) && (right_lambda > products[i])) {
            right_lambda = products[i];
        }
    }
    if (!R_FINITE(left_lambda) && !R_FINITE(right_lambda)) {
        return NA_REAL;
    } else {
        left_lambda = -(1.0 - 1.0 / ((double) n)) / left_lambda;
        right_lambda = -(1.0 - 1.0 / ((double) n)) / right_lambda;

        left_value = lagrangian(left_lambda, products, n);
        right_value = lagrangian(right_lambda, products, n);

        if (fabs(left_value) < sqrt(DBL_EPSILON)) {
            return left_lambda;
        } else if (fabs(right_value) < sqrt(DBL_EPSILON)) {
            return right_lambda;
        } else if (left_value * right_value < 0) {
            return lagrangian_root_rec(left_lambda, left_value, right_lambda, right_value, products, n);
        } else {
            return NA_REAL;
        }

    }

}

double lagrangian_root_rec(double left_lambda, double left_value, double right_lambda, double right_value,
                           double *products, unsigned int n) {

    double lambda = (left_lambda + right_lambda) / 2.0;
    double value = lagrangian(lambda, products, n);

    if (fabs(value) < sqrt(DBL_EPSILON)) {
        return lambda;
    } else if ((right_lambda - left_lambda) < 10.0 * DBL_EPSILON) {
        return NA_REAL;
    } else if ((value * left_value) < 0) {
        return lagrangian_root_rec(left_lambda, left_value, lambda, value, products, n);
    } else if ((value * right_value) < 0) {
        return lagrangian_root_rec(lambda, value, right_lambda, right_value, products, n);
    } else {
        return NA_REAL;
    }

}

double elh(double lambda, const double *products, unsigned int n) {
    unsigned int i = 0;
    double result = 0.0;

    for (i = 0; i < n; i++) {
        result += log(1.0 + lambda * products[i]);
    }

    result *= 2.0;

    return result;

}