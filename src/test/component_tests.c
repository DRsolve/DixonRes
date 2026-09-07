/* SPDX-License-Identifier: GPL-2.0-or-later */
#include <stdio.h>
#include <string.h>

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

int main(int argc, char **argv)
{
    if (argc != 2 || strcmp(argv[1], "--list") == 0) {
        printf("Usage: %s <test>\n", argv[0]);
        puts("Tests: determinant, roots, complexity, rational-solver, polynomial-solver,");
        puts("       ideal, ideal-string, resultant-ideal, sparse, roots-benchmark");
        puts("The sparse interpolation test and roots benchmark can take a long time.");
        return argc == 2 ? 0 : 1;
    }
    if (strcmp(argv[1], "determinant") == 0) test_unified_mpoly_det();
    else if (strcmp(argv[1], "roots") == 0) test_fq_nmod_correctness();
    else if (strcmp(argv[1], "complexity") == 0) return test_dixon_complexity();
    else if (strcmp(argv[1], "rational-solver") == 0) test_rational_polynomial_solver();
    else if (strcmp(argv[1], "polynomial-solver") == 0) test_polynomial_solver();
    else if (strcmp(argv[1], "ideal") == 0) test_iterative_elimination();
    else if (strcmp(argv[1], "ideal-string") == 0) test_iterative_elimination_str();
    else if (strcmp(argv[1], "resultant-ideal") == 0) test_iterative_elimination_str2();
    else if (strcmp(argv[1], "sparse") == 0) return test_sparse_interpolation();
    else if (strcmp(argv[1], "roots-benchmark") == 0) test_fq_nmod_benchmarks();
    else {
        fprintf(stderr, "Unknown test: %s\n", argv[1]);
        return 1;
    }
    return 0;
}
