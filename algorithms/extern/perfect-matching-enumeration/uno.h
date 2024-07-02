/***************************************************************************************************
 *                                                                                                 *
 *  uno.h - provides function declarations for uno.c                                               *
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

#ifndef __UNO_H
#define __UNO_H

void scc_iter(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
              int n, int v, int * i, int * st, int * st_size, int * index, int * lowlink,
              bool * onstack, int * scc_ptrs, int * scc_nodes, int * scc_count, int * scc_pos);

void scc(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive, int n,
         int ** _scc_ptrs, int ** _scc_nodes, int *scc_count);

bool dfs(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
         bool * visited, int * st, bool * onstack, int * st_size, int n, int v);

void find_cycle(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive,
                bool * visited, int * cycle, int * cycle_len, int n, int start);

void uno(int * col_ids, int * col_ptrs, int * col_to_row, int * row_to_col, int * inactive, int n,
         int r, int * matchings, int * matching_count);

#endif
