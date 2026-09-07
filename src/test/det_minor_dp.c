/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "unified_mpoly_det.h"
#include "fq_mpoly_mat_det.h"
#include <flint/nmod_mat.h>

/* Observe root scheduling without relying on noisy wall-time assertions.
 * Tagged first-row operands identify the expensive final multiplications. */
typedef ulong scheduling_poly_t[1];
typedef struct { int parallel, root_products, wrong_team, base_calls; } scheduling_ctx_t;
static void scheduling_base(scheduling_poly_t result, scheduling_poly_t **matrix,
                            slong size, scheduling_ctx_t *ctx)
{
    nmod_mat_t numeric;
    nmod_mat_init(numeric, size, size, 101);
    for (slong i = 0; i < size; i++)
        for (slong j = 0; j < size; j++)
            nmod_mat_entry(numeric, i, j) = matrix[i][j][0] % 101;
    result[0] = nmod_mat_det(numeric);
    nmod_mat_clear(numeric);
    ctx->base_calls++;
}
static ulong scheduling_mul(ulong a, ulong b, scheduling_ctx_t *ctx)
{
    if (a >= 10000) {
#ifdef _OPENMP
        #pragma omp atomic update
#endif
        ctx->root_products++;
#ifdef _OPENMP
        if ((omp_get_num_threads() > 1) != ctx->parallel) {
            #pragma omp atomic update
            ctx->wrong_team++;
        }
#endif
    }
    return ((a % 101) * b) % 101;
}
#define DET_DP_NAME scheduling_dp
#define DET_DP_MINOR scheduling_minor
#define DET_DP_BASE scheduling_base
#define DET_DP_POLY scheduling_poly_t
#define DET_DP_CTX scheduling_ctx_t *
#define DET_DP_LABEL "scheduling-test"
#define DET_DP_INIT(p, c) ((p)[0] = 0, 1)
#define DET_DP_CLEAR(p, c) ((void) (p))
#define DET_DP_SET(p, q, c) ((p)[0] = (q)[0])
#define DET_DP_ZERO(p, c) ((p)[0] = 0)
#define DET_DP_IS_ZERO(p, c) ((p)[0] == 0)
#define DET_DP_MUL(p, a, b, c) ((p)[0] = scheduling_mul((a)[0], (b)[0], c))
#define DET_DP_ADD(p, a, b, c) ((p)[0] = ((a)[0] + (b)[0]) % 101)
#define DET_DP_SUB(p, a, b, c) ((p)[0] = ((a)[0] + 101 - (b)[0]) % 101)
#define DET_DP_SWAP(p, q, c) do { ulong t = (p)[0]; (p)[0] = (q)[0]; (q)[0] = t; } while (0)
#include "determinant/det_minor_dp.h"

static void check_root_scheduling(void)
{
    for (slong n = 4; n <= 7; n++) {
        scheduling_poly_t entries[7][7], *rows[7], result = {0};
        nmod_mat_t numeric;
        nmod_mat_init(numeric, n, n, 101);
        for (slong i = 0; i < n; i++) {
            rows[i] = entries[i];
            for (slong j = 0; j < n; j++) {
                ulong value = 1;
                for (slong e = 0; e < (i == 0 ? n - 1 : i - 1); e++)
                    value = (value * (j + 1)) % 101;
                if (i == 0) value += 10000;
                entries[i][j][0] = value;
                nmod_mat_entry(numeric, i, j) = value % 101;
            }
        }
        ulong expected = nmod_mat_det(numeric);
        for (int parallel = 0; parallel <= 1; parallel++) {
            scheduling_ctx_t ctx = {parallel, 0, 0, 0};
            if (n == 6 && scheduling_dp(result, rows, n, &ctx, parallel, 34)) abort();
            if (!scheduling_dp(result, rows, n, &ctx, parallel, n == 6 ? 35 : 1024) ||
                result[0] != expected || ctx.root_products != n || ctx.wrong_team) {
                fprintf(stderr, "root scheduling/reduction regression: n=%ld parallel=%d\n", n, parallel);
                abort();
            }
        }
        if (n == 6) {
            scheduling_ctx_t ctx = {0, 0, 0, 0};
            scheduling_poly_t saved[7][7];
            for (slong i = 0; i < n; i++)
                memcpy(saved[i], entries[i], (size_t) n * sizeof(scheduling_poly_t));
            /* 34 entries cannot fit n=6, but each n=5 child fits. Do not
             * silently revert to the uncached base algorithm for this case. */
            scheduling_minor(result, rows, n, &ctx, 0, 34);
            if (result[0] != expected || ctx.base_calls != 0 || ctx.root_products != n) abort();
            for (slong i = 0; i < n; i++)
                if (memcmp(saved[i], entries[i], (size_t) n * sizeof(scheduling_poly_t))) abort();
            scheduling_minor(result, rows, n, &ctx, 0, 0);
            if (result[0] != expected || ctx.base_calls != 1) abort();
        }
        nmod_mat_clear(numeric);
    }
}

static void require_equal(unified_mpoly_t actual, unified_mpoly_t expected,
                          slong n, slong limit, int parallel)
{
    if (!unified_mpoly_equal(actual, expected)) {
        fprintf(stderr, "minor DP mismatch: field=%d n=%ld limit=%ld parallel=%d\n",
                actual->field_id, n, limit, parallel);
        abort();
    }
}

/* Exercise the public prime-field route, which bypasses unified_mpoly. */
static void check_direct_nmod(unified_mpoly_t **matrix, slong n,
                              unified_mpoly_t expected, const fq_nmod_ctx_t fq)
{
    unified_mpoly_ctx_t ctx = expected->ctx_ptr;
    fq_mvpoly_t **input = malloc((size_t) n * sizeof(*input));
    fq_mvpoly_t result;
    nmod_mpoly_t actual;
    nmod_mpoly_init(actual, GET_NMOD_CTX(ctx));
    for (slong i = 0; i < n; i++) {
        input[i] = malloc((size_t) n * sizeof(**input));
        for (slong j = 0; j < n; j++)
            nmod_mpoly_to_fq_mvpoly(&input[i][j], GET_NMOD_POLY(matrix[i][j]),
                                    2, 0, GET_NMOD_CTX(ctx), fq);
    }
    compute_fq_det_unified_interface(&result, input, n);
    fq_mvpoly_to_nmod_mpoly(actual, &result, GET_NMOD_CTX(ctx));
    if (!nmod_mpoly_equal(actual, GET_NMOD_POLY(expected), GET_NMOD_CTX(ctx))) abort();
    fq_mvpoly_clear(&result);
    nmod_mpoly_clear(actual, GET_NMOD_CTX(ctx));
    for (slong i = 0; i < n; i++) {
        for (slong j = 0; j < n; j++) fq_mvpoly_clear(&input[i][j]);
        free(input[i]);
    }
    free(input);
}

static void check_field(ulong prime, slong degree, ulong zech_limit)
{
    fq_nmod_ctx_t fq;
    field_ctx_t field;
    flint_rand_t state;
    fq_nmod_ctx_init_ui(fq, prime, degree, "a");
    field_ctx_init_enhanced(&field, fq, zech_limit);
    unified_mpoly_ctx_t ctx = unified_mpoly_ctx_init(2, ORD_LEX, &field);
    flint_rand_init(state);
    flint_rand_set_seed(state, 123, 456);
    for (slong n = 0; n <= 7; n++) {
        unified_mpoly_t **matrix = unified_mpoly_mat_init(n, n, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx);
        unified_mpoly_t actual = unified_mpoly_init(ctx);
        for (int variant = 0; variant < 4; variant++) {
            for (slong i = 0; i < n; i++) {
                for (slong j = 0; j < n; j++) {
                    unified_mpoly_randtest(matrix[i][j], state, 2, 2);
                    if ((variant == 1 && (i + j) % 3 == 0) ||
                        (variant == 2 && i > j)) unified_mpoly_zero(matrix[i][j]);
                }
            }
            if (variant == 3 && n > 1)
                for (slong j = 0; j < n; j++) unified_mpoly_set(matrix[1][j], matrix[0][j]);
            compute_unified_mpoly_det_recursive(expected, matrix, n, ctx);
            /* n=6 needs exactly 35 entries for two adjacent layers. */
            const slong limits[] = {0, 1, 34, 35, 1024};
            for (size_t l = 0; l < sizeof(limits) / sizeof(limits[0]); l++) {
                g_dixon_det_cache_limit = limits[l];
                for (int parallel = 0; parallel <= 1; parallel++) {
                    compute_unified_mpoly_det(actual, matrix, n, ctx, parallel);
                    require_equal(actual, expected, n, limits[l], parallel);
                }
                if (field.field_id == FIELD_ID_NMOD && n >= 4)
                    check_direct_nmod(matrix, n, expected, fq);
            }
        }
        unified_mpoly_clear(expected);
        unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(matrix, n, n);
    }

    if (field.field_id == FIELD_ID_NMOD) {
        /* Larger layers, checked against independent numeric elimination. */
        slong n = 11;
        nmod_mat_t numeric;
        nmod_mat_init(numeric, n, n, prime);
        nmod_mat_randtest(numeric, state);
        unified_mpoly_t **matrix = unified_mpoly_mat_init(n, n, ctx);
        unified_mpoly_t expected = unified_mpoly_init(ctx);
        unified_mpoly_t actual = unified_mpoly_init(ctx);
        for (slong i = 0; i < n; i++)
            for (slong j = 0; j < n; j++)
                nmod_mpoly_set_ui(GET_NMOD_POLY(matrix[i][j]),
                                  nmod_mat_entry(numeric, i, j), GET_NMOD_CTX(ctx));
        nmod_mpoly_set_ui(GET_NMOD_POLY(expected), nmod_mat_det(numeric), GET_NMOD_CTX(ctx));
        for (slong limit = 923; limit <= 924; limit++) {
            g_dixon_det_cache_limit = limit;
            for (int parallel = 0; parallel <= 1; parallel++) {
                compute_unified_mpoly_det(actual, matrix, n, ctx, parallel);
                require_equal(actual, expected, n, limit, parallel);
            }
        }
        unified_mpoly_clear(expected);
        unified_mpoly_clear(actual);
        unified_mpoly_mat_clear(matrix, n, n);
        nmod_mat_clear(numeric);
    }
    printf("minor DP tests passed: field=%d\n", field.field_id);
    flint_rand_clear(state);
    unified_mpoly_ctx_clear(ctx);
    field_ctx_clear(&field);
    fq_nmod_ctx_clear(fq);
}

int main(void)
{
    g_dixon_verbose_level = 0;
    g_field_equation_reduction = 0;
#ifdef _OPENMP
    omp_set_dynamic(0);
    omp_set_num_threads(4);
    omp_set_max_active_levels(1);
#endif
    check_root_scheduling();
    check_field(101, 1, 0);
    check_field(2, 3, 1024);
    check_field(3, 2, 0);
    check_field(2, 1, 0);
    return 0;
}
