/***************************************************************************************************
 *                                                                                                 *
 *  kpack.c - implements k-pack birkhoff decomposition algorithm                                   *
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
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "extern/perfect-matching-enumeration/uno.h"
#include "extern/matchmaker/matchmaker.h"
#include "bottleneck_matching.h"
#include "cscio.h"

#define FOV 64
#define EPS 1e-3
#define ZERO_EPS 1e-8
#define ONE_EPS 1e-4

#define MAX(X, Y) ((X) > (Y) ? (X) : (Y))
#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

#define DEBUG 0
#define DEBUG_SMALL 0
#define BENCHMARK 0

enum q_element_state
{
    EMPTY, GEN1, GEN2
};

void thresh_graph(int *col_ptrs, int *col_ids, double *col_vals, int *col_ptrs2, int *col_ids2,
                  double v, int n)
{
    // creat copy of csc without vals and only elements >= threshold v
    int i, j, h;
    for(i=n; i>=1; --i)
        col_ptrs2[i] = col_ptrs[i] - col_ptrs[i-1];
    for(j=0, h=0; j<n; ++j)
    {
        for(i=col_ptrs[j]; i<col_ptrs[j+1]; ++i)
        {
            if(col_vals[i] >= v)
                col_ids2[h++] = col_ids[i];
            else col_ptrs2[j+1]--;
        }
    }
    for(i=1; i<n+1; ++i)    
        col_ptrs2[i] += col_ptrs2[i-1];
}

void update_residual(int *col_ptrs, int *col_ids, double *col_vals, int n, int *k, int *row_to_col,
                     double v)
{
    int * update_pos = (int *) malloc(n * sizeof(int));
    int i, j, h;
    for(i=0; i<n; ++i)
    {
        int a = row_to_col[i];
        for(j=col_ptrs[a]; j<col_ptrs[a+1]; ++j)
        {
            if(col_ids[j] == i)
            {
                col_vals[j] -= v;
                update_pos[i] = j;
                break;
            }
        }
    }
    int compact_sz = 0, compact_min = *k;
    for(i=0; i<n; ++i)
    {
        if(col_vals[update_pos[i]] <= ZERO_EPS)
        {
            compact_min = MIN(compact_min, update_pos[i]);
            update_pos[compact_sz++] = row_to_col[i];
        }
    }
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
    free(update_pos);
}

double frobenius_norm(double *col_vals, int k)
{
    int i;
    double sos = 0.0;
    for(i=0; i<k; ++i)
        sos += col_vals[i] * col_vals[i];
    return sqrt(sos);
}

void decomposition(int *col_ptrs, int *col_ids, double *col_vals, int n, int *k,
                   int **_matchings, double **_factors, int *perm_count)
{
    int * qstate = (int *) calloc(FOV, sizeof(int));
    int * res_col_ptrs = (int *) calloc(FOV * (n+1), sizeof(int));
    int * res_col_ids = (int *) malloc(FOV * (*k) * sizeof(int));
    double * res_col_vals = (double *) malloc(FOV * (*k) * sizeof(double));
    int * res_cur_size = (int *) malloc(FOV * sizeof(int));
    int * col_ids2 = (int *) malloc((*k) * sizeof(int));
    int * col_ptrs2 = (int *) calloc(n+1, sizeof(int));
    int i, j, h, lb = 0, max_perm_count;
    for(i=0; i<FOV; ++i)
        res_cur_size[i] = *k;
    int * row_count = (int *) calloc(n, sizeof(int));
    for(i=0; i<(*k); ++i)
        row_count[col_ids[i]]++;
    for(i=0; i<n; ++i)
        lb = MAX(lb, MAX(col_ptrs[i+1]-col_ptrs[i], row_count[i]));
    max_perm_count = lb;
    int * inactive = NULL;
    int * best_col_to_row = row_count;
    int * best_row_to_col = (int *) malloc(n * sizeof(int));
    int * matchings = (int *) malloc(lb * FOV * n * sizeof(int));
    double * factors = (double *) malloc(lb * FOV * sizeof(double));
    // add initial element to queue
    qstate[0] = GEN1;
    for(i=0; i<n+1; ++i)
        res_col_ptrs[i] = col_ptrs[i];
    for(i=0; i<(*k); ++i)
    {
        res_col_ids[i]  = col_ids[i];
        res_col_vals[i] = col_vals[i];
    }
    int qsize = 1;
    *perm_count = 0;
    bool solution_found = 0;
    while(!solution_found)
    {
        int min_free_index = 0;
        double max_factor = 0.0;
        int req_lb = FOV / qsize;
        int rest = FOV % qsize;
        // realloc twice the size if necessary
        if((*perm_count) == max_perm_count)
        {
            max_perm_count <<= 1;
            matchings = (int *) realloc(matchings, max_perm_count * FOV * n * sizeof(int));
            factors = (double *) realloc(factors, max_perm_count * FOV * sizeof(double));
        }
        for(i=0; i<FOV; ++i)
        {
            if(qstate[i] != GEN1)
                continue;
            int req_limit = req_lb + (i < rest);
            int * ptrs = res_col_ptrs+i*(n+1);
            int * ids = res_col_ids+i*(*k);
            double * vals = res_col_vals+i*(*k);
            int * k_ = res_cur_size+i;
            double bottleneck_v = bottleneck_matching(ptrs, ids, vals, n, *k_, best_col_to_row,
                                                      best_row_to_col);
            max_factor = MAX(max_factor, bottleneck_v);
            if(bottleneck_v <= ZERO_EPS)
            {
#if DEBUG || DEBUG_SMALL
                printf("Values of q-element %d have degenerated to the point of no "
                                       "improvement.\n", i+1);
#endif
                qstate[i] = EMPTY;
                qsize--;
                continue;
            }
            if(req_limit == 1)
            {
                // update q-element with already found bottleneck matching
                // residual update
                update_residual(ptrs, ids, vals, n, k_, best_row_to_col, bottleneck_v);
                // add factor
                factors[(*perm_count)*FOV+i] = bottleneck_v;
                // add matching
                memcpy(matchings+(*perm_count)*FOV*n+i*n, best_row_to_col, n * sizeof(int));
                // set q-state
                qstate[i] = GEN2;
                // check end condition
                double frob = frobenius_norm(vals, *k_);
                if(frob <= EPS)
                {
#if DEBUG || DEBUG_SMALL || BENCHMARK
                    printf("Frobenius norm of %f was achieved.\n", frob);
#endif
                    solution_found = 1;
                    // allocate solution storage
                    int p = (*perm_count) + 1;
                    int * match = *_matchings = (int *) malloc(p * n * sizeof(int));
                    double * fac = *_factors = (double *) malloc(p * sizeof(double));
                    // copy solution to storage
                    for(j=0; j<p; ++j)
                    {
                        fac[j] = factors[j*FOV+i];
                        memcpy(match+j*n, matchings+j*FOV*n+i*n, n*sizeof(int));
                    }
                    break;
                }
            }
            else // req_limit > 1
            {
                // compute up to req_limit-1 additional matchings
                // build thresholded bipartite graph
                thresh_graph(ptrs, ids, vals, col_ptrs2, col_ids2, bottleneck_v, n);
                // run uno to get matchings
                if(!inactive)
                    inactive = (int *) calloc(n*n, sizeof(int));
                else
                    memset(inactive, 0, n*n*sizeof(int));
                int * req_matchings = (int *) malloc(req_limit * n * sizeof(int));
                int matching_count = 1;
                memcpy(req_matchings, best_row_to_col, n * sizeof(int));
                uno(col_ids2, col_ptrs2, best_col_to_row, best_row_to_col, inactive, n, req_limit,
                    req_matchings, &matching_count);
                // find matching_count - 1 empty slots
                int * element_pos = (int *) malloc(matching_count * sizeof(int));
                element_pos[0] = i;
                for(j=1; min_free_index<FOV && j<matching_count; ++min_free_index)
                    if(qstate[min_free_index] == EMPTY)
                        element_pos[j++] = min_free_index;
                // copy base to matching_count - 1 empty slots
                // copy residual, matchings and factors
                for(j=1; j<matching_count; ++j)
                {
                    int t = element_pos[j];
                    memcpy(res_col_ptrs+t*(n+1), res_col_ptrs+i*(n+1), (n+1) * sizeof(int));
                    memcpy(res_col_ids+t*(*k), res_col_ids+i*(*k), (*k) * sizeof(int));
                    memcpy(res_col_vals+t*(*k), res_col_vals+i*(*k), (*k) * sizeof(double));
                    res_cur_size[t] = res_cur_size[i];
                    for(h=0; h<(*perm_count); ++h)
                    {
                        memcpy(matchings+h*FOV*n+t*n, matchings+h*FOV*n+i*n,
                                                       n * sizeof(int));
                        factors[h*FOV+t] = factors[h*FOV+i];
                    }
                    qsize++;
                }
                // insert new matchings, bottle factor, update residual and qstate
                for(j=0; j<matching_count; ++j)
                {
                    int t = element_pos[j];
                    int pos = *perm_count;
                    memcpy(matchings+pos*FOV*n+t*n, req_matchings+j*n, n * sizeof(int));
                    factors[pos*FOV+t] = bottleneck_v;
                    update_residual(res_col_ptrs+t*(n+1), res_col_ids+t*(*k),
                                                        res_col_vals+t*(*k), n, res_cur_size+t,
                                                        req_matchings+j*n, bottleneck_v);
                    // check for end condition
                    double frob = frobenius_norm(res_col_vals+t*(*k), res_cur_size[t]);
                    if(frob <= EPS)
                    {
#if DEBUG || DEBUG_SMALL || BENCHMARK
                        printf("Frobenius norm of %f was achieved.\n", frob);
#endif
                        solution_found = 1;
                        // allocate solution storage
                        int p = (*perm_count) + 1;
                        int * match = *_matchings = (int *) malloc(p * n * sizeof(int));
                        double * fac = *_factors = (double *) malloc(p * sizeof(double));
                        // copy solution to storage
                        for(h=0; h<p; ++h)
                        {
                            fac[h] = factors[h*FOV+t];
                            memcpy(match+h*n, matchings+h*FOV*n+t*n, n * sizeof(int));
                        }
                        break;
                    }
                    qstate[t] = GEN2;
                }
                free(req_matchings);
                free(element_pos);
            }
            if(solution_found)
                break;
        }
        // filtering: use max factor to set non max elements to inactive (qsize-x)
        for(i=0; !solution_found && i<FOV; ++i)
        {
            int pos = (*perm_count);
            if(qstate[i] == EMPTY)
                continue;
            if(factors[pos*FOV+i] < max_factor-ZERO_EPS)
            {
                qstate[i] = EMPTY;
                qsize--;
            }
            else
                qstate[i] = GEN1;
        }
        if(!solution_found && !qsize)
        {
            *_matchings = NULL;
            *_factors = NULL;
            *perm_count = 0;
            break;
        }
        (*perm_count)++;
    }
    free(qstate); free(res_col_ptrs); free(res_col_ids); free(res_col_vals); free(res_cur_size);
    free(matchings); free(factors); free(best_col_to_row); free(best_row_to_col); free(col_ids2);
    free(col_ptrs2); free(inactive);
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
    print_matches(matchings, factors, n, perm_count);
    free(col_ids); free(col_ptrs); free(col_vals); free(matchings); free(factors);
    return 0;
}
