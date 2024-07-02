/***************************************************************************************************
 *                                                                                                 *
 *  cscio.c - implements input, output and verification of compressed sparse column (csc) data     *
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

#include "cscio.h"

void is_bistochastic(int *col_ptrs, int *col_ids, double *col_vals, int n, double one_eps)
{
    int i, j;
    double * row_sums = (double *) calloc(n, sizeof(double));
    for(j=0; j<n; ++j)
    {
        double sum = 0.0;
        for(i=col_ptrs[j]; i<col_ptrs[j+1]; ++i)
        {
            if(col_vals[i] < 0.0)
            {
                printf("Only non-negative matrix entries are supported.\n");
                free(row_sums);
                exit(5);
            }
            sum += col_vals[i]; row_sums[col_ids[i]] += col_vals[i];
        }
        if(sum < 1.0 - one_eps || sum > 1.0 + one_eps)
        {
            printf("Column sum was %lg, but has to be 1.0 +- %lg.\n", sum, one_eps);
            printf("Consider adjusting ONE_EPS tolerance.\n");
            free(row_sums);
            exit(6);
        }
    }
    for(i=0; i<n; ++i)
    {
        if(row_sums[i] < 1.0 - one_eps || row_sums[i] > 1.0 + one_eps)
        {
            printf("Row sum was %lg, but has to be 1.0 +- %lg.\n", row_sums[i], one_eps);
            printf("Consider adjusting ONE_EPS tolerance.\n");
            free(row_sums);
            exit(7);
        }
    }
    free(row_sums);
}

void print_csc_to_mm(int *col_ptrs, int *col_ids, double *col_vals, int n, int k)
{
    printf("%%MatrixMarket matrix coordinate real general\n%d %d %d\n", n, n, k);
    int i, j;
    for(i=0,j=1; i<k; ++i)
    {
        if(i==col_ptrs[j]) j++;
        printf("%d %d %lg\n", col_ids[i]+1, j, col_vals[i]);
    }
}

void print_matches(int * matchings, double * factors, int n, int perm_count)
{
    printf("%d %d\n", n, perm_count);
    int i, j;
    for(i=0; i<perm_count; ++i) 
    {
        printf("\n%lg\n", factors[i]);
        for(j=0; j<n; ++j)
            printf("%d ", matchings[i*n+j]);
        printf("\n");
    }
}

void print_matchings(int * matchings, int matching_count, int n)
{
    int i, j;
    for(i=0; i<matching_count; ++i)
    {
        printf("matching %d: ", i+1);
        for(j=0; j<n; ++j)
            printf("(%d,%d) ", j, matchings[i*n+j]);
        printf("\n");
    }
}

void print_csc(int *col_ptrs, int *col_ids, double *col_vals, int n, int k)
{
    int i;
    printf("col_vals: ");
    for(i=0; i<k; ++i)
        printf("%f ", col_vals[i]);
    printf("\ncol_ids: ");
    for(i=0; i<k; ++i)
        printf("%d ", col_ids[i]);
    printf("\ncol_ptrs: ");
    for(i=0; i<n+1; ++i)
        printf("%d ", col_ptrs[i]);
    printf("\n");
}

void read_mm(int *n, int *m, int *k, int **col_ptrs, int **col_ids, double **col_vals)
{
    MM_typecode matcode;
    if(mm_read_banner(stdin, &matcode))
    {
        printf("Matrix Market banner couldn't been read from stdin.\n");
        exit(1);
    }
    if(!mm_is_matrix(matcode) || !mm_is_sparse(matcode) || !mm_is_real(matcode))
    {
        printf("Only sparse matrices with float values are supported.\n");
        exit(2);
    }
    if(mm_read_mtx_crd_size(stdin, n, m, k))
    {
        printf("Matrix size parsing failed.\n");
        exit(3);
    }
    if(*n != *m)
    {
        printf("Only square matrices are supported.\n");
        exit(4);
    }
    int * I = (int *) malloc((*k) * sizeof(int));
    int * J = (int *) malloc((*k) * sizeof(int));
    double * v = (double *) malloc((*k) * sizeof(double));
    int i;
    for(i=0; i<*k; ++i)
    {
        scanf("%d %d %lg\n", &I[i], &J[i], &v[i]);
        I[i]--; J[i]--;
    }
    int * w = (int *) calloc((*n)+1, sizeof(int));
    int * ci = *col_ids = (int *) malloc((*k) * sizeof(int));
    int * cp = *col_ptrs = (int *) malloc(((*n)+1) * sizeof(int));
    double * cv = *col_vals = (double *) malloc((*k) * sizeof(double));
    for(i=0; i<*k; ++i)
        w[J[i]]++;
    cp[0] = 0;
    for(i=1; i<=*n; ++i)
    {
        cp[i] = w[i-1];
        w[i] = w[i-1] + w[i];
    }
    for(i=(*k)-1; i>=0; --i)
    {
        int a = I[i], b = J[i];
        ci[--w[b]] = a;
        cv[w[b]] = v[i];
    }
    free(I); free(J); free(v); free(w);
}

void print_cycle(int * cycle, int cycle_len)
{
    printf("cycle: ");
    int i;
    for(i=0; i<cycle_len; ++i)
        printf("%d ", cycle[i]);
    printf("\n");
}

void print_scc(int * scc_ptrs, int * scc_nodes, int scc_count)
{
    int i, j;
    for(i=0; i<scc_count; ++i)
    {
        printf("cc %d: ", i+1);
        for(j=scc_ptrs[i]; j<scc_ptrs[i+1]; ++j)
            printf("%d ", scc_nodes[j]);
        printf("\n");
    }
}
