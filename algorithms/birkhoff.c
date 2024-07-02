/***************************************************************************************************
 *                                                                                                 *
 *  birkhoff.c - implements classic birkhoff decomposition algorithm                               *
 *  Copyright (C) 2024  SteezyE                                                                    *
 *                                                                                                 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of     *
 *  the GNU Affero General Public License as published by the Free Software Foundation, either     *
 *  version 3 of the License, or (at your option) any later version.                               *
 *                                                                                                 *
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;      *
 *  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.      *
 *  See the GNU Affero General Public License for more details.                                    *
 *                                                                                                 *
 *  You should have received a copy of the GNU Affero General Public License along with this       *
 *  program. If not, see <https://www.gnu.org/licenses/>.                                          *
 *                                                                                                 *
 ***************************************************************************************************/

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "cscio.h"
#include "extern/matchmaker/matchmaker.h"

#define EPS 1e-3
#define ZERO_EPS 1e-8
#define ONE_EPS 1e-4

#define ALGO_ID 5
#define CHEAP_ID 1
#define RELABEL_P 1.0

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

#define DEBUG 0
#define DEBUG_SMALL 0
#define BENCHMARK 0

void decomposition(int *col_ptrs, int *col_ids, double *col_vals, int n, int *k, int **_matchings,
                   double **_factors, int *perm_count)
{
#if DEBUG
    printf("-- main loop begin --\n");
#endif
    int lb = 0;
    int * row_count = (int *) calloc(n, sizeof(int));
    int i, j, h;
    for(i=0; i<(*k); ++i)
        row_count[col_ids[i]]++;
    for(i=0; i<n; ++i)
        lb = MAX(lb, MAX(col_ptrs[i+1]-col_ptrs[i], row_count[i]));
    int max_elements = lb;
    // allocate memory for matchings and factors (lowerbound of columns / rows)
    int * matchings = *_matchings = (int *) malloc(max_elements * n * sizeof(int));
    double * factors = *_factors = (double *) malloc(max_elements * sizeof(double));
    int * col_to_row = row_count;
    int * row_to_col = (int *) malloc(n * sizeof(int));
    int cardm = 0;
    double fac_sum = 0.0;
    *perm_count = 0;
    while(1)
    {
        // get matching
        cardm = 0;
        for(i=0; i<n; ++i)
            row_to_col[i] = -1;
        matching(col_ptrs, col_ids, col_to_row, row_to_col, n, n, ALGO_ID, CHEAP_ID, RELABEL_P);
        for(i=0; i<n; ++i)
            cardm += row_to_col[i] > -1;
        if(cardm < n)
        {
#if DEBUG || DEBUG_SMALL
            printf("Found matching was not perfect, leaving main loop.\n");
#endif
            break;
        }
        // realloc twice the size if necessary
        if((*perm_count) == max_elements)
        {
            max_elements <<= 1;
            matchings = *_matchings = (int *) realloc(matchings, max_elements * n * sizeof(int));
            factors = *_factors = (double *) realloc(factors, max_elements * sizeof(double));
        }
        // match array mapping columns to rows will be reused (to update residual)
        int * update_pos = col_to_row;
        // compute factor and store matching
        double mini = 1.0;
        for(i=0; i<n; ++i)
        {
            int a = row_to_col[i];
            matchings[(*perm_count)*n+i] = a;
            for(j=col_ptrs[a]; j<col_ptrs[a+1]; ++j)
            {
                if(col_ids[j] == i)
                {
                    update_pos[i] = j;
                    mini = MIN(mini, col_vals[j]);
                    break;
                }
            }
        }
        // store factor
        factors[(*perm_count)++] = mini;
        fac_sum += mini;
#if DEBUG_SMALL
        printf("%d %lg %f\n", *perm_count, mini, fac_sum);
#endif
        // update values in csc
        int compact_sz = 0, compact_min = *k;
        for(i=0; i<n; ++i)
        {
            col_vals[update_pos[i]] -= mini;
            if(col_vals[update_pos[i]] <= ZERO_EPS)
            {
                compact_min = MIN(compact_min, update_pos[i]);
                update_pos[compact_sz++] = row_to_col[i];
            }
        }
        // check end condition
        double sos = 0.0;
        for(i=0; i<*k; ++i)
            sos += col_vals[i] * col_vals[i];
        if(sqrt(sos) <= EPS) break;
        // move zero elements to end and decrease size (maintaining relative positions)
        i = compact_min; h = 1;
        while(h <= compact_sz)
        {
            while(i != (*k)-1 && col_vals[i] > ZERO_EPS)
                i++;
            for(j=i; j<(*k)-h; j++)
            {
                col_vals[j] = col_vals[j+1];
                col_ids[j] = col_ids[j+1];
            }
            h++;
        }
        (*k) -= compact_sz;
        // update sizes in col_ptrs array
        for(i=n; i>=1; --i)
            col_ptrs[i] -= col_ptrs[i-1];
        for(i=0; i<compact_sz; ++i)
            col_ptrs[update_pos[i]+1] -= 1;
        for(i=1; i<n+1; ++i)
            col_ptrs[i] += col_ptrs[i-1];
#if DEBUG
        printf("perm_number: %d; factor: %f; factor_sum: %f;\n", *perm_count, mini, fac_sum);
        print_csc(col_ptrs, col_ids, col_vals, n, *k);
        printf("matching: ");
        for(i=0; i<n; ++i)
            printf("%d ", matchings[((*perm_count)-1)*n+i]);
        printf("\n");
#endif
    }
    free(col_to_row); free(row_to_col);
#if DEBUG
    printf("-- main loop end --\n\n");
#endif
}

int main()
{
    int n, m, k, perm_count, *col_ids, *col_ptrs, *matchings;
    double *col_vals, *factors;
    read_mm(&n, &m, &k, &col_ptrs, &col_ids, &col_vals);
#if DEBUG
    printf("-- initialization begin --\n");
    print_csc(col_ptrs, col_ids, col_vals, n, k);
    printf("-- initialization end --\n\n");
#endif
    is_bistochastic(col_ptrs, col_ids, col_vals, n, ONE_EPS);
#if BENCHMARK
    clock_t tic = clock();
#endif
    decomposition(col_ptrs, col_ids, col_vals, n, &k, &matchings, &factors, &perm_count);
#if BENCHMARK
    clock_t toc = clock();
    printf("It took %fs to find the decomposition.\n", (double) (toc - tic) / CLOCKS_PER_SEC);
#endif
#if DEBUG || DEBUG_SMALL || BENCHMARK
    double sos = 0.0;
    int i;
    for(i=0; i<k; ++i)
        sos += col_vals[i] * col_vals[i];
    printf("Frobenius norm of %f was achieved.\n", sqrt(sos));
#endif
    print_matches(matchings, factors, n, perm_count);
    free(col_ids); free(col_ptrs); free(col_vals); free(matchings); free(factors);
    return 0;
}
