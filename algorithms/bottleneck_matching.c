/***************************************************************************************************
 *                                                                                                 *
 *  bottleneck_matching.h - implements bottleneck matching algorithm (based on binary search)      *
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
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "extern/quadsort/quadsort.h"
#include "extern/matchmaker/matchmaker.h"

#define ALGO_ID 5
#define CHEAP_ID 1
#define RELABEL_P 1.0

#define MIN(X, Y) ((X) < (Y) ? (X) : (Y))

double bottleneck_matching(int *col_ptrs, int *col_ids, double *col_vals, int n, int k,
                           int * best_col_to_row, int * best_row_to_col)
{
    double * sorted_vals = (double *) malloc(k * sizeof(double));
    int * col_ids2 = (int *) malloc(k * sizeof(int));
    int * col_ptrs2 = (int *) malloc((n+1) * sizeof(int));
    int * col_to_row = (int *) malloc(n * sizeof(int));
    int * row_to_col = (int *) malloc(n * sizeof(int));
    col_ptrs2[0] = 0;
    memcpy(sorted_vals, col_vals, k * sizeof(double));
    quadsort_prim(sorted_vals, k, sizeof(double));
    double best_factor = 0.0, mini, v;
    int l = 0, r = k-1, mid, cardm, i, j, h;
    while(l <= r)
    {
        mid = l + (r-l) / 2;
        v = sorted_vals[mid];
        // skip current value <= best value
        if(v <= best_factor)
        {
            l = mid + 1;
            continue;
        }
        // create copy of csc without vals and only elements >= threshold
        for(i=n; i>=1; --i)
            col_ptrs2[i] = col_ptrs[i] - col_ptrs[i-1];
        bool flag = 0;
        for(j=0, h=0; j<n; ++j)
        {
            int m = h;
            for(i=col_ptrs[j]; i<col_ptrs[j+1]; ++i)
            {
                if(col_vals[i] >= v)
                    col_ids2[h++] = col_ids[i];
                else col_ptrs2[j+1]--;
            }
            // break early if there can't be a perfect matching
            if(m == h)
            {
                r = mid - 1;
                flag = 1;
                break;
            }
        }
        if(flag)
            continue;
        for(i=1; i<n+1; ++i)
            col_ptrs2[i] += col_ptrs2[i-1];
        // find matching
        cardm = 0;
        for(i=0; i<n; ++i)
            row_to_col[i] = -1;
        matching(col_ptrs2, col_ids2, col_to_row, row_to_col, n, n, ALGO_ID, CHEAP_ID, RELABEL_P);
        for(i=0; i<n; ++i)
            cardm += row_to_col[i] > -1;
        if(cardm < n)
        {
            r = mid - 1;
            continue;
        }
        // else if cardm == n
        mini = 1.0;
        for(i=0; i<n; ++i)
        {
            int a = row_to_col[i];
            best_row_to_col[i] = a;
            best_col_to_row[i] = col_to_row[i];
            for(j=col_ptrs[a]; j<col_ptrs[a+1]; ++j)
            {
                if(col_ids[j] == i)
                {
                    mini = MIN(mini, col_vals[j]);
                    break;
                }
            }
        }
        // update best with factor of current matching
        best_factor = mini;
        l = mid + 1;
    }
    free(sorted_vals); free(col_ptrs2); free(col_ids2); free(col_to_row); free(row_to_col);
    return best_factor;
}
