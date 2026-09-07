/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "dixon_complexity.h"

int test_dixon_complexity(void) {
    // Test data
    long a1[] = {1000, 1000, 1000, 1001, 1002, 1003};
    int len = sizeof(a1) / sizeof(a1[0]);
    double omega = DIXON_OMEGA;
    
    printf("Dixon Complexity Results (Hessenberg Method):\n");
    for (int n = 5; n < 10; n++) {
        double complexity = dixon_complexity(a1, len, n, omega);
        printf("n=%d: %.6f\n", n, complexity);
    }
    
    // Test dixon_size function
    printf("\nTesting dixon_size with Hessenberg method:\n");
    fmpz_t test_result;
    fmpz_init(test_result);
    dixon_size(test_result, a1, len, 1);
    fmpz_clear(test_result);
    
    return 0;
}
