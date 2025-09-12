#include "header.h"

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