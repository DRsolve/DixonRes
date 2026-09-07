/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_COMPONENT_TESTS_H
#define DRSOLVE_COMPONENT_TESTS_H

void test_unified_mpoly_det(void);
void test_polynomial_solver(void);
void test_rational_polynomial_solver(void);
int test_dixon_complexity(void);
void test_iterative_elimination(void);
void test_iterative_elimination_str(void);
void test_iterative_elimination_str2(void);
void test_fq_nmod_correctness(void);
void test_fq_nmod_benchmarks(void);
int test_sparse_interpolation(void);

#endif
