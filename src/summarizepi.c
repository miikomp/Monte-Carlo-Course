#include "header.h"

void summarizeResultsArray(const long double *results) {
    const long N = GLOB.n_outer;

    /* Calculate mean of results */
    
    long double sum = 0.0L, sumsq = 0.0L;
    for (long i = 0; i < GLOB.n_outer; i++)
    {
        sum += results[i];
        sumsq += results[i] * results[i];
    }

    long double mean = sum / (long double)N;

    /* Standard deviation */

    long double bracket = sumsq - (sum * sum) / (long double)N;

    long double sd = sqrtl(bracket / ((long double)N * (N - 1)));

    /* 95% confidence interval 1.96 * sd */

    long double lo = mean - 1.96 * sd;
    long double hi = mean + 1.96 * sd;

    fprintf(stdout, "\nmean = %.6Lf +- %.6Lf [%.6Lf, %.6Lf]\n\n", mean, sd, lo, hi);
}