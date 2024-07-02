/***************************************************************************************************
 *                                                                                                 *
 *  bottleneck_matching.h - provides function declaration for bottleneck_matching.c                *
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

#ifndef __BOTTLENECK_H
#define __BOTTLENECK_H

double bottleneck_matching(int *col_ptrs, int *col_ids, double *col_vals, int n, int k,
                           int * best_col_to_row, int * best_row_to_col);

#endif
