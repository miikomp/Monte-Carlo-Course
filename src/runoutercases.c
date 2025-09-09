#include "header.h"

void updateProgressBar(long current, long total);

int runCirclePi(uint64_t sm, uint64_t *seeds) {

    /* Run the specified number of outer iterations and approximate pi */
    /* Outer loop is single threaded and discpatched inner loops are multihreaded */

    fprintf(stdout, "Approximating pi uisng the 1/4 circle method with %ld outer, and %ld inner iterations...\n\n", GLOB.n_outer, GLOB.n_inner);

    /* Allocate for results array */

    double *results = calloc(GLOB.n_outer, sizeof(double));
    if (!results)
    {
        fprintf(stderr, "Memory allocation error.\n");
        free(seeds);
        return EXIT_FAILURE;
    }
    
    /* Run outer loop */

    for (long m = 0; m < GLOB.n_outer; m++)
    {

        /* Update progress bar */

        updateProgressBar(m+1, GLOB.n_outer);


        /* Derive seeds for run */

        for (long i = 0; i < GLOB.n_threads; i++) 
        {
            seeds[i] = splitmix64(&sm);
            if (!seeds[i]) 
                seeds[i] = 0x9E3779B97F4A7C15ULL;
        }

        /* Run inner loop */

        PiResult res = {0,0};
        
        if (runCirclePiInner(&res, seeds) != 0)
        {
            fprintf(stderr, "Failed at outer iteration %ld.\n", m);
            free(seeds); 
            free(results);
            return EXIT_FAILURE;
        }

        /* Calculate the estimate for pi */
        if (res.n_hits == 0) 
            results[m] = 0.0;
        else 
            results[m] = 4.0L * (double)res.n_hits / (double)res.n_tot;
    }

    /* Process results */

    summarizeResultsArray(results);

    /* Free memory */

    free(results);

    return 0;
}

int runBuffonsPi(uint64_t sm, uint64_t *seeds) {
    /* Run the specified number of outer iterations for Buffon's Needle */

    fprintf(stdout, "Approximating pi using Buffon's Needle with %ld outer, and %ld inner iterations...\n\n", GLOB.n_outer, GLOB.n_inner);

    /* Allocate results array */

    double *results = calloc(GLOB.n_outer, sizeof(double));
    if (!results) {
        fprintf(stderr, "Memory allocation error.\n");
        free(seeds);
        return EXIT_FAILURE;
    }

    /* Outer loop */

    for (long m = 0; m < GLOB.n_outer; m++) 
    {

        /* Update progress bar */

        updateProgressBar(m + 1, GLOB.n_outer);

        /* Derive seeds for this run */

        for (long i = 0; i < GLOB.n_threads; i++) 
        {
            seeds[i] = splitmix64(&sm);
            if (!seeds[i]) 
                seeds[i] = 0x9E3779B97F4A7C15ULL;
        }

        /* Run inner loop */

        PiResult res = {0, 0};
        if (runBuffonsPiInner(&res, seeds) != 0) 
        {
            fprintf(stderr, "Failed at outer iteration %ld.\n", m);
            free(seeds);
            free(results);
            return EXIT_FAILURE;
        }

        /* Calculate the estimate for pi */

        if (res.n_hits == 0) 
            results[m] = 0.0;
        else 
            results[m] = (2.0 * GLOB.needle_length * res.n_tot) / (res.n_hits * GLOB.line_spacing);
    }

    /* Process results */

    summarizeResultsArray(results);

    /* Free memory */

    free(results);

    return 0;
}

void updateProgressBar(long current, long total) {
    long barWidth = PRG_BAR_WIDTH;
    long progress = (long)((100.0 * current) / total);
    long pos = (current * barWidth) / total;

    char bar[barWidth + 50];
    long idx = 0;

    /* Start the progress bar */
    idx += sprintf(&bar[idx], "\r[");

    // Fill the progress bar
    for (long i = 0; i < barWidth; i++) {
        if (i < pos)
            bar[idx++] = '=';
        else if (i == pos)
            bar[idx++] = '>';
        else
            bar[idx++] = ' ';
    }

    /* Close and put percentage */
    sprintf(&bar[idx], "] %ld%%", progress);

    /* Flush */
    fprintf(stdout, "%s", bar);
    fflush(stdout);
}