/* Rafal Baranowski 7 Dec 2015 Public Domain */
#ifndef RBSS_H
#define RBSS_H

#define IDX(i,j,ld) ((((j)-1) * (ld))+((i)-1))
#define IDXUT(i,j,ld) (i <= j) ? (((((j)-1) * (j))/2)+((i)-1)) : (((((i)-1) * (i))/2)+((j)-1))

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "sort_r.h"

void insertion_sort(unsigned int *sorted, unsigned int element, unsigned int length);

#endif
