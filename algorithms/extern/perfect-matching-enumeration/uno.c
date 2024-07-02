/***************************************************************************************************
 *                                                                                                 *
 *  uno.c - enumerates perfect matchings in a bipartite graph                                      *
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
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

#include "../matchmaker/matchmaker.h"

#define MIN(X,Y) ( (X) > (Y) ? (Y) : (X) )

void scc_iter(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
              int n, int v, int * i, int * st, int * st_size, int * index, int * lowlink,
              bool * onstack, int * scc_ptrs, int * scc_nodes, int * scc_count, int * scc_pos)
{
    index[v] = *i;
    lowlink[v] = (*i)++;
    st[(*st_size)++] = v;
    onstack[v] = 1;
    int j, w, t;
    bool flag = v < n;
    if(!flag)
        t = col_to_row[v-n];
    for(j=(flag ? 0 : col_ptrs[v-n]); j<(flag ? 1 : col_ptrs[v+1-n]); ++j)
    {
        w = flag ? row_to_col[v]+n : col_ids[j];
        if(!flag && (w == t || inactive[(v-n)*n+w]))
            continue;
        if(index[w] == -1)
        {
            scc_iter(col_ids, col_ptrs, col_to_row, row_to_col, inactive, n,
                     w, i, st, st_size, index, lowlink, onstack,
                     scc_ptrs, scc_nodes, scc_count, scc_pos);
            lowlink[v] = MIN(lowlink[v], lowlink[w]);
        }
        else if(onstack[w])
            lowlink[v] = MIN(lowlink[v], index[w]);
    }
    if(lowlink[v] == index[v])
    {
        do
        {
            w = st[--(*st_size)];
            onstack[w] = 0; 
            scc_nodes[(*scc_pos)++] = w;
        } while(w != v);
        scc_ptrs[++(*scc_count)] = *scc_pos;
    }
}

// csc: edges from column to rows
// row_to_col: edges from row to column
// col_to_row: edges from column to row
void scc(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive, int n,
         int ** _scc_ptrs, int ** _scc_nodes, int *scc_count)
{
    int * scc_nodes = *_scc_nodes = (int *) malloc(2 * n * sizeof(int));
    int * scc_ptrs = *_scc_ptrs = (int *) malloc((2 * n + 1) * sizeof(int));
    scc_ptrs[0] = 0;
    *scc_count = 0;

    int i = 0, v, st_size = 0, scc_pos = 0;
    int * st = (int *) malloc(2 * n * sizeof(int));
    int * index = (int *) malloc(2 * n * sizeof(int));
    for(v=0; v<2*n; ++v)
        index[v] = -1;
    int * lowlink = (int *) malloc(2 * n * sizeof(int));
    bool * onstack = (bool *) calloc(2 * n, sizeof(bool));

    for(v=0; v<2*n; ++v)
    {
        if(index[v] == -1)
            scc_iter(col_ids, col_ptrs, col_to_row, row_to_col, inactive, n,
                                 v, &i, st, &st_size, index, lowlink, onstack,
                                 scc_ptrs, scc_nodes, scc_count, &scc_pos);
    }
    free(st); free(index); free(lowlink); free(onstack);
}

bool dfs(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
         bool * visited, int * st, bool * onstack, int * st_size, int n, int v)
{
    visited[v] = 1;
    st[(*st_size)++] = v;
    onstack[v] = 1;
    int j, w, t;
    bool flag = v < n;
    if(!flag)
        t = col_to_row[v-n];
    for(j=(flag ? 0 : col_ptrs[v-n]); j<(flag ? 1 : col_ptrs[v+1-n]); ++j)
    {
        w = flag ? row_to_col[v]+n : col_ids[j];
        if(!flag && (w == t || inactive[(v-n)*n+w]))
            continue;
        if(!visited[w])
            if(dfs(col_ids, col_ptrs, col_to_row, row_to_col, inactive, visited, st, onstack,
                   st_size, n, w))
                return 1;
        if(onstack[w])
        {
            st[(*st_size)++] = w;
            return 1;
        }
    }
    (*st_size)--;
    onstack[v] = 0;
    return 0;
}

// cycle is returned "backwards"
void find_cycle(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
                bool * visited, int * cycle, int * cycle_len, int n, int start)
{
    int * st = (int *) malloc((2 * n + 1) * sizeof(int));
    bool * onstack = (bool *) calloc(2 * n, sizeof(bool));
    int st_size = 0;

    dfs(col_ids, col_ptrs, col_to_row, row_to_col, inactive, visited, st, onstack, &st_size, n,
        start);

    *cycle_len=0;
    start = st[--st_size];
    do cycle[(*cycle_len)++] = st[--st_size];
    while(st[st_size] != start);
    int last = (*cycle_len)-1;
    if(cycle[last] < n)
    {
        int t = cycle[0];
        int i;
        for(i=0; i<last; ++i)
            cycle[i] = cycle[i+1];
        cycle[last] = t;
    }

    free(st); free(onstack);
}

void uno(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive, int n,
         int r, int * matchings, int * matching_count)
{
    if(*matching_count >= r)
        return;
    int i, j;
    int *scc_ptrs, *scc_nodes, scc_count;
    scc(col_ids, col_ptrs, col_to_row, row_to_col, inactive, n, &scc_ptrs, &scc_nodes, &scc_count);
    if(scc_count == 2*n)
    {
        free(scc_ptrs); free(scc_nodes);
        return;
    }
    int * new_row_to_col = (int *) malloc(n * sizeof(int));
    int * new_col_to_row = (int *) malloc(n * sizeof(int));
    memcpy(new_row_to_col, row_to_col, n * sizeof(int));
    memcpy(new_col_to_row, col_to_row, n * sizeof(int));
    int min_cycle_size = INT_MAX, min_cycle_start;
    // marks node for deletion
    bool * del = (bool *) calloc(2*n, sizeof(bool));
    for(i=0; i<scc_count; ++i)
    {
        int low = scc_ptrs[i];
        int scc_size = scc_ptrs[i+1] - low;
        if(scc_size > 1)
        {
            if(scc_size < min_cycle_size)
            {
                min_cycle_size = scc_size;
                min_cycle_start = low;
            }
        }
        else del[scc_nodes[low]] = 1;
    }
    while(scc_nodes[min_cycle_start] < n)
        min_cycle_start++;
    int start = scc_nodes[min_cycle_start];
    // begin dfs from start, return cycle when encountering a node once again
    // (don't visit del nodes)
    int * cycle = (int *) malloc((2 * n + 1) * sizeof(int));
    int cycle_len;
    find_cycle(col_ids, col_ptrs, col_to_row, row_to_col, inactive, del, cycle, &cycle_len, n,
               start);
    for(j=0; j<cycle_len; ++j)
    {
        if(j % 2 == 0)
        {
            int t = (cycle_len<<1)-j-1;
            int v_1 = cycle[t % cycle_len] - n;
            int v_2 = cycle[(t-1) % cycle_len];
            new_col_to_row[v_1] = v_2;
            new_row_to_col[v_2] = v_1;
        }
    }
    for(j=0; j<n; ++j)
        matchings[(*matching_count)*n+j] = new_row_to_col[j];
    (*matching_count)++;
    if(*matching_count >= r)
    {
        free(cycle); free(new_row_to_col); free(new_col_to_row); free(scc_ptrs);
                free(scc_nodes); free(del);
        return;
    }
    // set edges inactive for subproblem 1
    int e_1 = cycle[cycle_len-1]-n;
    int e_2 = cycle[cycle_len-2];
    for(i=0; i<n; ++i)
    {
        for(j=col_ptrs[i]; j<col_ptrs[i+1]; ++j)
        {
            // delete e, both endpoints and all adjacent edges
            if(i == e_1 || col_ids[j] == e_2)
                inactive[i * n + col_ids[j]]++;
        }
    }
    // call uno for subproblem 1
    uno(col_ids, col_ptrs, new_col_to_row, new_row_to_col, inactive, n, r, matchings,
            matching_count);
    // reset edges that have been set inactive
    for(i=0; i<n; ++i)
    {
        for(j=col_ptrs[i]; j<col_ptrs[i+1]; ++j)
        {
            // reincarnate edges that have been just deleted
            if(i == e_1 || col_ids[j] == e_2)
                inactive[i * n + col_ids[j]]--;
        }
    }
    // set edges inactive for subproblem 2
    inactive[e_1 * n + e_2]++;
    // call uno for subproblem 2
    uno(col_ids, col_ptrs, col_to_row, row_to_col, inactive, n, r, matchings,
        matching_count);
    // reset nodes that have been set inactive
    inactive[e_1 * n + e_2]--;
    free(cycle); free(new_row_to_col); free(new_col_to_row); free(scc_ptrs); free(scc_nodes);
    free(del);
}
