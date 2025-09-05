#include "header.h"

void summarizeResultsArray(const double *results) {
    const long N = GLOB.n_outer;

    /* Calculate mean of results */
    
    double sum = 0.0L, sumsq = 0.0L;
    for (long i = 0; i < GLOB.n_outer; i++)
    {
        sum += results[i];
        sumsq += results[i] * results[i];
    }

    double mean = sum / (double)N;

    /* Standard deviation */

    double bracket = sumsq - (sum * sum) / (double)N;

    double sd = sqrtl(bracket / ((double)N * (N - 1)));

    /* 95% confidence interval 1.96 * sd */

    double lo = mean - 1.96 * sd;
    double hi = mean + 1.96 * sd;

    fprintf(stdout, "\n\nMean value = %.6lf +- %.6lf [%.6lf, %.6lf]\n\n", mean, sd, lo, hi);

    /* Figure-of-Merit from relative error */

    double t0 = GLOB.t0;
    double t1 = omp_get_wtime();

    double rerr = sd / mean;

    /* FOM is scaled to billions */

    double FOM = 1e-9 * (1.0L / (rerr * rerr * (t1 - t0)));

    fprintf(stdout, "FOM = %.2lf\n", FOM);
}