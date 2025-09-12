#include "header.h"

double normalCDF(double x);

double ksTest(const double *data, long size);

int compareDoubles(const void *a, const void *b);

void summarizePiResultsArray(const double *results) {
    const long N = GLOB.n_generations;

    /* Calculate mean of results */
    
    double sum = 0.0L, sumsq = 0.0L;
    for (long i = 0; i < N; i++)
    {
        sum += results[i];
        sumsq += results[i] * results[i];
    }

    const double mean = sum / (double)N;

    /* Variance and sample standard deviation */

    const double ss = sumsq - (sum * sum) / (double)N;
    const double var = ss / (double)(N - 1);
    const double sd_s = sqrt(var);

    /* Standard error of mean (SEM) and 95% CI for the mean */

    const double sem  = sd_s / sqrt((double)N);
    const double lo   = mean - 1.96 * sem;
    const double hi   = mean + 1.96 * sem;

    fprintf(stdout, "\n\nMean value = %.6lf +/- %.6lf [%.6lf, %.6lf]\n\n", mean, sem, lo, hi);

    /* Figure-of-Merit from relative error */

    const double t0 = GLOB.t0;
    const double t1 = omp_get_wtime();
    const double rerr = sem / mean;                               
    const double FOM = 1e-9 * (1.0L / (rerr * rerr * (t1 - t0)));     // in billions

    fprintf(stdout, "FOM = %.2lf\n\n", FOM);

    /* Normalise the results by sample standard deviation */

    double *norm = calloc(N, sizeof(double));
    if (!norm) {
        fprintf(stderr, "[ERROR] Memory allocation error.\n");
        exit(EXIT_FAILURE);
    }
    for (long i = 0; i < N; i++) {
        norm[i] = (results[i] - mean) / sd_s;
    }

    /* K-S test against standard normal */

    const double d_max = ksTest(norm, N);

    /* Output the K-S test result */

    const double ks_crit = 1.36 / sqrt((double)N);
    const double lf_crit = 0.886 / sqrt((double)N);

    fprintf(stdout, "Kolmogorov-Smirnov test for normality:\n");
    fprintf(stdout, "          K-S D: %.6f\n", d_max);
    fprintf(stdout, "  Standard crit: %.6f: %s\n",
            ks_crit, (d_max < ks_crit ? "TRUE" : "FALSE"));
    fprintf(stdout, "Lilliefors crit: %.6f: %s\n\n",
            lf_crit, (d_max < lf_crit ? "TRUE" : "FALSE"));

    free(norm);

    return;
}

double normalCDF(double x) {
    return 0.5 * (1.0 + erf(x / sqrt(2.0)));
}

double ksTest(const double *data, long size) {
    if (size <= 0) {
        fprintf(stderr, "[ERROR] Invalid size for \"ksTest\".\n");
        exit(EXIT_FAILURE);
    }
    /* Sort the data */

    double *sorted = malloc(size * sizeof(double));
    if (!sorted) {
        fprintf(stderr, "[ERROR] Memory allocation error in ksTest.\n");
        exit(EXIT_FAILURE);
    }
    memcpy(sorted, data, size * sizeof(double));
    qsort(sorted, size, sizeof(double), compareDoubles);

    /* Compute the empirical CDF and compare to the theoretical CDF */

    double d_max = 0.0;
    for (long i = 0; i < size; i++) {
        double empirical_cdf = (double)(i + 1) / size;
        double theoretical_cdf = normalCDF(sorted[i]);
        double d_plus = fabs(empirical_cdf - theoretical_cdf);
        double d_minus = fabs(theoretical_cdf - (double)i / size);
        double d = fmax(d_plus, d_minus);
        if (d > d_max) {
            d_max = d;
        }
    }

    free(sorted);
    return d_max;
}

int compareDoubles(const void *a, const void *b) {
    double diff = *(double *)a - *(double *)b;
    return (diff > 0) - (diff < 0);
}