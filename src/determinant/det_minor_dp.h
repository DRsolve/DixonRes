/* SPDX-License-Identifier: GPL-2.0-or-later */
/* Internal template, instantiated for unified_mpoly and nmod_mpoly.
 * No include guard: each inclusion defines DET_DP_* operations below.
 * A state is the last k rows and a k-subset of columns, in colex order.
 * Only adjacent layers are retained; the root uses one cofactor per column
 * so that its largest multiplications can also run in parallel.
 * limit counts polynomial values,
 * excluding arithmetic temporaries and matrix views.
 */
static int DET_DP_NAME(DET_DP_POLY result, DET_DP_POLY **matrix,
                       slong size, DET_DP_CTX ctx, int use_parallel, slong limit)
{
    ulong choose[FLINT_BITS][FLINT_BITS] = {{0}};
    DET_DP_POLY *previous = NULL, *current = NULL;
    slong previous_count = 0, current_count = 0;
    slong peak = 0;
    int ok = 0;

    if (size <= 3 || size >= FLINT_BITS || limit <= 0) return 0;

    /* Saturation avoids overflow even when the requested matrix is far too
     * large. All ranks used after this preflight fit within the entry limit. */
    for (slong n = 0; n <= size; n++) {
        choose[n][0] = 1;
        for (slong k = 1; k <= n; k++) {
            ulong a = choose[n - 1][k - 1], b = choose[n - 1][k];
            choose[n][k] = (a > (ulong) limit || b > (ulong) limit - a)
                            ? (ulong) limit + 1 : a + b;
        }
    }
    for (slong k = 1; k <= size; k++) {
        ulong count = k == size ? (ulong) size : choose[size][k];
        ulong prev = k == 1 ? 0 : choose[size][k - 1];
        if (count > (ulong) limit || prev > (ulong) limit - count)
            return 0;
        if (count > (size_t) -1 / sizeof(DET_DP_POLY)) return 0;
        if ((slong) (count + prev) > peak) peak = (slong) (count + prev);
    }

    /* Keep cheap sparse expansions on the demand-driven path. The product
     * of row nonzero counts bounds the number of recursive branches; stop
     * counting once it already exceeds the full DP multiplication count. */
    {
        double dp_work = 0, recursive_work = 0, branches = 1;
        for (slong k = 2; k <= size; k++) dp_work += k * (double) choose[size][k];
        for (slong row = 0; row < size - 3; row++) {
            slong nonzero = 0;
            for (slong col = 0; col < size; col++)
                if (!DET_DP_IS_ZERO(matrix[row][col], ctx)) nonzero++;
            branches *= FLINT_MIN(nonzero, size - row);
            recursive_work += branches;
            if (recursive_work >= dp_work) break;
        }
        if (recursive_work + 12 * branches < dp_work) return 0;
    }

    if (g_dixon_verbose_level >= 2)
        printf("  determinant layered DP (%s): size=%ld, peak entries=%ld, limit=%ld, parallel=%s\n",
               DET_DP_LABEL, size, peak, limit, use_parallel ? "yes" : "no");

    for (slong k = 1; k <= size; k++) {
        slong count = k == size ? size : (slong) choose[size][k];
        int failed = 0;
        current = malloc((size_t) count * sizeof(DET_DP_POLY));
        if (current == NULL) goto cleanup;
        for (current_count = 0; current_count < count; current_count++) {
            if (!DET_DP_INIT(current[current_count], ctx)) goto cleanup;
        }

        if (k == 1) {
            for (slong i = 0; i < count; i++)
                DET_DP_SET(current[i], matrix[size - 1][i], ctx);
        } else {
            /* One writer per output, immutable previous layer, implicit
             * barrier before freeing it. No hash locks or duplicate work. */
#ifdef _OPENMP
            #pragma omp parallel if(use_parallel && count > 1 && !omp_in_parallel()) num_threads(FLINT_MIN(count, omp_get_max_threads())) reduction(|:failed)
#endif
            {
                DET_DP_POLY product, sum;
                int product_ok = DET_DP_INIT(product, ctx);
                int sum_ok = DET_DP_INIT(sum, ctx);
                if (!product_ok || !sum_ok) failed = 1;
#ifdef _OPENMP
                #pragma omp for schedule(dynamic, 1)
#endif
                for (slong index = 0; index < count; index++) {
                    slong cols[FLINT_BITS];
                    ulong prefix[FLINT_BITS], suffix[FLINT_BITS];
                    ulong rank = (ulong) index;
                    slong col = size - 1;
                    if (!product_ok || !sum_ok) continue;

                    /* There is only one root state. Parallelize its expansion
                     * terms instead of serializing all n large products in
                     * one worker. The colex rank of the (n-1)-subset missing
                     * column index is n-1-index. For n>=4, these n outputs and
                     * n inputs fit below the middle-layer entry peak. */
                    if (k == size) {
                        if (DET_DP_IS_ZERO(matrix[0][index], ctx) ||
                            DET_DP_IS_ZERO(previous[size - 1 - index], ctx)) continue;
                        DET_DP_MUL(current[index], matrix[0][index],
                                   previous[size - 1 - index], ctx);
                        if (index & 1) {
                            DET_DP_ZERO(product, ctx);
                            DET_DP_SUB(current[index], product, current[index], ctx);
                        }
                        continue;
                    }

                    /* Unrank sum_i C(cols[i], i+1) in O(n) time. */
                    for (slong j = k; j > 0; j--) {
                        while (choose[col][j] > rank) col--;
                        cols[j - 1] = col;
                        rank -= choose[col][j];
                        col--;
                    }
                    prefix[0] = 0;
                    for (slong j = 0; j < k; j++)
                        prefix[j + 1] = prefix[j] + choose[cols[j]][j + 1];
                    suffix[k] = 0;
                    for (slong j = k - 1; j > 0; j--)
                        suffix[j] = suffix[j + 1] + choose[cols[j]][j];

                    DET_DP_ZERO(current[index], ctx);
                    for (slong j = 0; j < k; j++) {
                        ulong child = prefix[j] + suffix[j + 1];
                        if (DET_DP_IS_ZERO(matrix[size - k][cols[j]], ctx) ||
                            DET_DP_IS_ZERO(previous[child], ctx)) continue;
                        DET_DP_MUL(product, matrix[size - k][cols[j]], previous[child], ctx);
                        if (j & 1)
                            DET_DP_SUB(sum, current[index], product, ctx);
                        else
                            DET_DP_ADD(sum, current[index], product, ctx);
                        DET_DP_SWAP(current[index], sum, ctx);
                    }
                }
                /* The worksharing barrier above publishes all root terms.
                 * A balanced reduction keeps large sums parallel and avoids
                 * repeatedly merging one term into an ever-growing prefix. */
                if (k == size) {
                    for (slong stride = 1; stride < count; stride *= 2) {
#ifdef _OPENMP
                        #pragma omp for schedule(static)
#endif
                        for (slong left = 0; left < count; left += 2 * stride) {
                            if (!product_ok || !sum_ok || left + stride >= count) continue;
                            DET_DP_ADD(sum, current[left], current[left + stride], ctx);
                            DET_DP_SWAP(current[left], sum, ctx);
                        }
                    }
                }
                if (product_ok) DET_DP_CLEAR(product, ctx);
                if (sum_ok) DET_DP_CLEAR(sum, ctx);
            }
        }
        if (failed) goto cleanup;
        for (slong i = 0; i < previous_count; i++) DET_DP_CLEAR(previous[i], ctx);
        free(previous);
        previous = current;
        previous_count = current_count;
        current = NULL;
        current_count = 0;
    }
    DET_DP_SWAP(result, previous[0], ctx);
    ok = 1;

cleanup:
    for (slong i = 0; i < current_count; i++) DET_DP_CLEAR(current[i], ctx);
    for (slong i = 0; i < previous_count; i++) DET_DP_CLEAR(previous[i], ctx);
    free(current);
    free(previous);
    return ok;
}

/* Method 0 backend. If a complete layer does not fit, expand one row and
 * retry DP on each child. Siblings execute sequentially, so their DP storage
 * never multiplies the entry budget; each child may still use layer threads.
 * The submatrix is a shallow, read-only view of the input polynomials. */
#ifdef DET_DP_MINOR
static void DET_DP_MINOR(DET_DP_POLY result, DET_DP_POLY **matrix,
                         slong size, DET_DP_CTX ctx, int use_parallel, slong limit)
{
    DET_DP_POLY **rows = NULL, *entries = NULL;
    DET_DP_POLY accum, child, product, sum;
    int accum_ok, child_ok, product_ok, sum_ok;
    size_t width;

    if (size <= 3 || limit <= 0) {
        DET_DP_BASE(result, matrix, size, ctx);
        return;
    }
    if (DET_DP_NAME(result, matrix, size, ctx, use_parallel, limit)) return;

    width = (size_t) (size - 1);
    if (width > (size_t) -1 / sizeof(*rows) ||
        width > (size_t) -1 / sizeof(*entries) / width) {
        DET_DP_BASE(result, matrix, size, ctx);
        return;
    }
    rows = malloc(width * sizeof(*rows));
    entries = malloc(width * width * sizeof(*entries));
    if (rows == NULL || entries == NULL) {
        free(rows);
        free(entries);
        DET_DP_BASE(result, matrix, size, ctx);
        return;
    }
    for (size_t i = 0; i < width; i++) rows[i] = entries + i * width;
    accum_ok = DET_DP_INIT(accum, ctx);
    child_ok = DET_DP_INIT(child, ctx);
    product_ok = DET_DP_INIT(product, ctx);
    sum_ok = DET_DP_INIT(sum, ctx);
    if (accum_ok && child_ok && product_ok && sum_ok) {
        DET_DP_ZERO(accum, ctx);
        for (slong col = 0; col < size; col++) {
            if (DET_DP_IS_ZERO(matrix[0][col], ctx)) continue;
            for (slong i = 1; i < size; i++) {
                slong dst = 0;
                for (slong j = 0; j < size; j++) {
                    if (j == col) continue;
                    memcpy(&rows[i - 1][dst++], &matrix[i][j], sizeof(*entries));
                }
            }
            DET_DP_MINOR(child, rows, size - 1, ctx, use_parallel, limit);
            if (DET_DP_IS_ZERO(child, ctx)) continue;
            DET_DP_MUL(product, matrix[0][col], child, ctx);
            if (col & 1)
                DET_DP_SUB(sum, accum, product, ctx);
            else
                DET_DP_ADD(sum, accum, product, ctx);
            DET_DP_SWAP(accum, sum, ctx);
        }
        DET_DP_SWAP(result, accum, ctx);
    } else {
        DET_DP_BASE(result, matrix, size, ctx);
    }
    if (accum_ok) DET_DP_CLEAR(accum, ctx);
    if (child_ok) DET_DP_CLEAR(child, ctx);
    if (product_ok) DET_DP_CLEAR(product, ctx);
    if (sum_ok) DET_DP_CLEAR(sum, ctx);
    free(entries);
    free(rows);
}
#endif

#undef DET_DP_MINOR
#undef DET_DP_BASE
#undef DET_DP_NAME
#undef DET_DP_POLY
#undef DET_DP_CTX
#undef DET_DP_LABEL
#undef DET_DP_INIT
#undef DET_DP_CLEAR
#undef DET_DP_SET
#undef DET_DP_ZERO
#undef DET_DP_IS_ZERO
#undef DET_DP_MUL
#undef DET_DP_ADD
#undef DET_DP_SUB
#undef DET_DP_SWAP
