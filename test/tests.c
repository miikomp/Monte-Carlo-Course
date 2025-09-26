#include "header.h"
#include <inttypes.h>

#define N 10000000
#define MEAN_TOL 0.002
#define VAR_TOL 0.002
#define PROD_TOL 0.002
#define BIT_TOL  0.002
#define TRIG_TOL 0.00001
#define CDF_TOL  0.00001

int test_xoshiro256ss_uniformity();
int test_xoshiro256ss_independence();
int test_splitmix64_distribution();
int test_initTrigTables();
int test_rotl_bit_rotation();
int test_xoshiro256ss_reproducibility();
int test_splitmix64_known_values();
int test_randd_bounds();

/**
 * @brief Main function for the testing environment. Compiled and run with "make test"
 * 
 * @return int 0 if all tests passed
 */
int main(void) {
    printf("Running tests...\n");

    int i = 0;

    /* Test pRNGs */
    
    i += test_xoshiro256ss_uniformity();
    i += test_xoshiro256ss_independence();
    i += test_splitmix64_distribution();
    i += test_rotl_bit_rotation();
    i += test_xoshiro256ss_reproducibility();
    i += test_splitmix64_known_values();
    i += test_randd_bounds();

    /* Test initializers */
    
    i += test_initTrigTables();

    /* Final result */

    if (i != 0) {
        printf("\n%d test(s) failed!\n", i);
    }
    else {
        printf("\nAll tests passed!\n");
    }

    return i;
}

int test_xoshiro256ss_uniformity() {
    xoshiro256ss_state state;
    xoshiro256ss_seed(&state, 123456789);
    double sum = 0.0, sumsq = 0.0, min = 1.0, max = 0.0;
    for (int i = 0; i < N; ++i) {
        double r = randd(&state);
        sum += r;
        sumsq += r * r;
        if (r < min) min = r;
        if (r > max) max = r;
    }
    double mean = sum / N;
    double var = (sumsq / N) - (mean * mean);
    printf("\nxoshiro256ss_double() uniformity test:\n");
    printf("  Mean: %.6f (expected 0.5)\n", mean);
    printf("  Variance: %.6f (expected ~0.0833)\n", var);
    printf("  Min: %.6f, Max: %.6f (expected [0,1))\n", min, max);
    int fail = 0;
    if (fabs(mean - 0.5) > MEAN_TOL) fail = 1;
    if (fabs(var - 1.0/12.0) > VAR_TOL) fail = 1;
    if (min < 0.0 || max >= 1.0) fail = 1;
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_xoshiro256ss_independence() {
    xoshiro256ss_state state;
    xoshiro256ss_seed(&state, 987654321);
    double prev = randd(&state);
    double sum = 0.0;
    for (int i = 1; i < N; ++i) {
        double curr = randd(&state);
        sum += prev * curr;
        prev = curr;
    }
    double corr = sum / (N - 1);
    printf("\nxoshiro256ss_double() independence test:\n");
    printf("  Mean product of consecutive values: %.6f (expected ~0.25)\n", corr);
    int fail = fabs(corr - 0.25) > PROD_TOL;
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_splitmix64_distribution() {
    uint64_t state = 123456789;
    double sum = 0.0, sumsq = 0.0;
    for (int i = 0; i < N; ++i) {
        uint64_t r = splitmix64(&state);
        double d = (double)(r >> 11) * (1.0 / 9007199254740992.0);
        sum += d;
        sumsq += d * d;
    }
    double mean = sum / N;
    double var = (sumsq / N) - (mean * mean);
    printf("\nsplitmix64() uniformity test:\n");
    printf("  Mean: %.6f (expected 0.5)\n", mean);
    printf("  Variance: %.6f (expected ~0.0833)\n", var);
    int fail = 0;
    if (fabs(mean - 0.5) > MEAN_TOL) fail = 1;
    if (fabs(var - 1.0/12.0) > VAR_TOL) fail = 1;
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_initTrigTables() {
    initTrigTables();
    int fail = 0;
    for (int i = 0; i < 1000; ++i) {
        double angle = (M_PI * i) / (TRIG_LOOKUP_TABLE_SIZE - 1);
        double s = sin_table[i];
        double c = cos_table[i];
        double t = tan_table[i];
        if (fabs(s - sin(angle)) > TRIG_TOL) fail = 1;
        if (fabs(c - cos(angle)) > TRIG_TOL) fail = 1;
        if (fabs(t - tan(angle)) > TRIG_TOL) fail = 1;
    }
    printf("\ninitTrigTables() test:\n");
    printf("  First 1000 sin/cos/tan values checked against math.h\n");
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_rotl_bit_rotation() {
    printf("\nrotl() bit rotation test:\n");
    int fail = 0;
    uint64_t v1 = rotl(UINT64_C(1), 1);
    uint64_t v2 = rotl(UINT64_C(0x8000000000000000), 1);
    uint64_t v3 = rotl(UINT64_C(0x0123456789ABCDEF), 16);
    printf("  rotl(1,1) -> 0x%016" PRIx64 " (expected 0x0000000000000002)\n", v1);
    printf("  rotl(0x8000...0,1) -> 0x%016" PRIx64 " (expected 0x0000000000000001)\n", v2);
    printf("  rotl(0x0123456789ABCDEF,16) -> 0x%016" PRIx64 "\n", v3);
    if (v1 != UINT64_C(0x2)) fail = 1;
    if (v2 != UINT64_C(0x1)) fail = 1;
    if (v3 != UINT64_C(0x456789ABCDEF0123)) fail = 1;
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_xoshiro256ss_reproducibility() {
    printf("\nxoshiro256ss() reproducibility test:\n");
    xoshiro256ss_state a, b;
    xoshiro256ss_seed(&a, 0xDEADBEEFCAFEBABEULL);
    xoshiro256ss_seed(&b, 0xDEADBEEFCAFEBABEULL);

    int fail = 0;
    for (int i = 0; i < 8; ++i) {
        uint64_t ra = xoshiro256ss(&a);
        uint64_t rb = xoshiro256ss(&b);
        if (ra != rb) fail = 1;
        printf("  draw[%d] = 0x%016" PRIx64 "\n", i, ra);
    }

    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_splitmix64_known_values() {
    printf("\nsplitmix64() known value test:\n");
    uint64_t state = 0;
    uint64_t first = splitmix64(&state);
    uint64_t second = splitmix64(&state);
    printf("  first  = 0x%016" PRIx64 "\n", first);
    printf("  second = 0x%016" PRIx64 "\n", second);
    int fail = 0;
    if (first != UINT64_C(0xE220A8397B1DCDAF)) fail = 1;
    if (second != UINT64_C(0x6E789E6AA1B965F4)) fail = 1;
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}

int test_randd_bounds() {
    printf("\nrandd() bounds test:\n");
    xoshiro256ss_state state;
    xoshiro256ss_seed(&state, 0x123456789ABCDEF0ULL);
    int fail = 0;
    double min = 1.0, max = 0.0;
    const int samples = 1000000;
    for (int i = 0; i < samples; ++i) {
        double r = randd(&state);
        if (r < min) min = r;
        if (r > max) max = r;
        if (r < 0.0 || r >= 1.0) fail = 1;
    }
    printf("  min = %.12f, max = %.12f\n", min, max);
    printf("  Result: %s\n", fail ? "Fail" : "Success");
    return fail;
}
