/***************************************************************************************************
 *                                                                                                 *
 *  cscio.h - provides function declarations and includes for cscio.c                              *
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

#ifndef __CSCIO_H
#define __CSCIO_H

#include <stdio.h>
#include <stdlib.h>

#include "extern/mmio/mmio.h"

void is_bistochastic(int *col_ptrs, int *col_ids, double *col_vals, int n, double one_eps);
void read_mm(int *n, int *m, int *k, int **col_ptrs, int **col_ids, double **col_vals);
void print_csc_to_mm(int *col_ptrs, int *col_ids, double *col_vals, int n, int k);
void print_matches(int * matchings, double * factors, int n, int perm_count);
void print_csc(int *col_ptrs, int *col_ids, double *col_vals, int n, int k);
void print_matchings(int *matchings, int matching_count, int n);
void print_scc(int * scc_ptrs, int * scc_nodes, int scc_count);
void print_cycle(int * cycle, int cycle_len);

#endif
