/* Rafal Baranowski 7 Dec 2015 Public Domain */
#include "rbvs.h"




/* k is the subset size*/

struct ranks_param {
	unsigned int k;
	unsigned int p;
	unsigned int *ranks;
};



      
void sort_rows(unsigned int *row_id, unsigned int B, unsigned int k, unsigned int p, unsigned int *ranks){
    
  #define COMPARE_ROWS(row_a, row_b)  memcmp(&ranks[row_id[row_a] * p], &ranks[row_id[row_b] * p], (k + 1) * sizeof(unsigned int)) 
  #define COMPARE_ROW_TO_PIVOT(row_a)  memcmp(&ranks[row_id[row_a] * p], &ranks[pivot * p], (k + 1) * sizeof(unsigned int)) 
  #define SWAP_ROWS(row_a, row_b) { tmp = row_id[row_a]; row_id[row_a] = row_id[row_b]; row_id[row_b] = tmp; }
  
  unsigned int i,j, tmp;
  
  /*insertion sort if B<=6*/
  if(B <= 6){
    
    for(i = 1; i < B; i++){
      j = i;
      while(j>=1) if(COMPARE_ROWS(j-1, j) < 0){
        SWAP_ROWS(j-1, j)
        j--;
      }else break;
      
    }
    
  }else if(B>0){
    /*select the pivot*/
    unsigned int pivot = row_id[B-1];

    /*partitioning*/
    i = 0;
    for(j=0; j < B-1; j++){
      if(COMPARE_ROW_TO_PIVOT(j) <= 0){
        SWAP_ROWS(i, j)
        i++;
      }
    }
    
    SWAP_ROWS(i,  B-1)
    
    
    sort_rows(row_id, i, k, p, ranks);
    sort_rows(&row_id[i+1], B-i-1, k, p, ranks);
    
  }
  
  
}

int ranks_rows_cmp(const void *ii, const void *jj, void *arg)
{
	struct ranks_param param = *((struct ranks_param *)arg);
	const unsigned int i = *((int *)ii);
	const unsigned int j = *((int *)jj);
	const unsigned int k = param.k;
	const unsigned int p = param.p;
	const unsigned int *ranks = param.ranks;

	return (memcmp
		(&ranks[i * p], &ranks[j * p], (k + 1) * sizeof(unsigned int)));
}

unsigned int bin_search(unsigned int *sorted, unsigned int element,
			unsigned int low, unsigned int high)
{
	if (high <= low)
		return (element > sorted[low]) ? (low + 1) : low;

	unsigned int mid = (low + high) / 2;

	if (element > sorted[mid])
		return bin_search(sorted, element, mid + 1, high);
	else if (mid == low)
		return (element > sorted[low]) ? (low + 1) : low;
	else
		return bin_search(sorted, element, low, mid - 1);

}

void insertion_sort(unsigned int *sorted, unsigned int element,
		    unsigned int length)
{

	if (length > 0) {

		unsigned int position =
		    bin_search(sorted, element, 0, length - 1);
		memmove(&sorted[position + 1], &sorted[position],
			(length - position) * sizeof(unsigned int));
		sorted[position] = element;

	} else {
		sorted[0] = element;
	}

}

/*function finds susbets which occur the most frequently in the top of the ranking*/

unsigned int k_top_ranked_sets(unsigned int *ranks, unsigned int p,
			       unsigned int B, unsigned int *best_freq,
			       unsigned int *best_subset,
			       unsigned int min_max_freq, unsigned int k_max)
{

	register unsigned int k = 0, i, j;
	register unsigned int current_id, max_id, current_freq, max_freq;
	size_t bytes_to_cmp = 0;

	unsigned int *row_id = Calloc(B, unsigned int);

	for (i = 0; i < B; i++) {
		row_id[i] = i;
	}

	struct ranks_param param = { k, p, ranks };

	for (k = 0; k < k_max; k++) {

		param.k = k;

		for (j = 0; j < B; j++)
			insertion_sort(&ranks[j * p], ranks[j * p + k], k);

		sort_r(row_id, B, sizeof(unsigned int), ranks_rows_cmp, &param);


		current_id = row_id[0];
		max_id = row_id[0];
		current_freq = 0;
		max_freq = 0;

		bytes_to_cmp = (k + 1) * sizeof(unsigned int);

		for (j = 0; j <= (B - 1); j++) {

			if (memcmp
			    (&(ranks[current_id * p]), &(ranks[row_id[j] * p]),
			     bytes_to_cmp)) {
				if (current_freq > max_freq) {
					max_id = current_id;
					max_freq = current_freq;
				}

				current_id = row_id[j];
				current_freq = 1;

			} else
				current_freq++;

		}

		if (current_freq > max_freq) {
			max_freq = current_freq;
			max_id = current_id;
		}

		best_freq[k] = max_freq;

		memcpy(&best_subset[(unsigned int)((k + 1) * k) / 2],
		       &(ranks[max_id * p]), bytes_to_cmp);

		if (max_freq == min_max_freq) {

			Free(row_id);
			return k + 1;
		}

	}

	Free(row_id);

	return k;

}

SEXP k_top_ranked_sets_r(SEXP ranks, SEXP k_max, SEXP min_max_freq, SEXP active)
{

	SEXP ranks_dim;
	int no_protect = 0;

	PROTECT(ranks_dim = getAttrib(ranks, R_DimSymbol));

	no_protect++;

	unsigned int p = INTEGER(ranks_dim)[0];
	unsigned int B = INTEGER(ranks_dim)[1];

	unsigned int val_k_max = INTEGER(k_max)[0];
	unsigned int val_min_max_freq = INTEGER(min_max_freq)[0];
	unsigned int *best_subset_tmp =
	    Calloc((val_k_max * (val_k_max + 1)) / 2, unsigned int);
	unsigned int *best_freq_tmp = Calloc(val_k_max, unsigned int);
	unsigned int *ranks_tmp = Calloc(p * B, unsigned int);

	memcpy(ranks_tmp, INTEGER(ranks), p * B * sizeof(unsigned int));

	unsigned int val_k =
	    k_top_ranked_sets(ranks_tmp, p, B, best_freq_tmp, best_subset_tmp,
			      val_min_max_freq, val_k_max);

	SEXP subsets;
	PROTECT(subsets = allocVector(VECSXP, val_k));
	SEXP freq;
	PROTECT(freq = allocVector(INTSXP, val_k));
	no_protect += 2;
	SEXP subset;
	unsigned int *ptr_subset;
	unsigned int *ptr_freq = (unsigned int *)INTEGER(freq);

	register unsigned int i = 0, l, j;
	register unsigned int len_active = length(active);
	unsigned int *ptr_active = (unsigned int *)INTEGER(active);

	if (len_active > 0) {

		for (l = 0; l < len_active; l++) {
			for (j = 0; j < (val_k * (val_k + 1)) / 2; j++)
				if (best_subset_tmp[j] >= ptr_active[l])
					best_subset_tmp[j]++;
		}

	}

	i = 0;

	for (j = 0; j < val_k; j++) {

		subset = SET_VECTOR_ELT(subsets, j, allocVector(INTSXP, j + 1));
		ptr_subset = (unsigned int *)INTEGER(subset);

		for (l = 0; l < j + 1; l++) {

			ptr_subset[l] = best_subset_tmp[i];
			i++;
		}

	}

	for (l = 0; l < val_k; l++)
		ptr_freq[l] = best_freq_tmp[l];

	SEXP res;
	PROTECT(res = allocVector(VECSXP, 2));
	SEXP res_names;
	PROTECT(res_names = allocVector(STRSXP, 2));
	char *names[2] = { "frequencies", "subsets" };
	no_protect += 2;

	SET_VECTOR_ELT(res, 0, freq);
	SET_VECTOR_ELT(res, 1, subsets);

	for (j = 0; j < 2; j++)
		SET_STRING_ELT(res_names, j, mkChar(names[j]));
	setAttrib(res, R_NamesSymbol, res_names);

	UNPROTECT(no_protect);

	Free(best_subset_tmp);
	Free(best_freq_tmp);
	Free(ranks_tmp);

	return (res);

}
