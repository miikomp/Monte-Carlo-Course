#include "header.h"

/**
 * @brief Updates the progress bar in the console.
 * 
 * @param current Current progress value
 * @param total Total progress value
 */
void updateProgressBar(long current, long total);

/**
 * @brief Runs the quarter circle approximation for pi
 * 
 * @param out ptr to result struct
 * @param seeds array of thread private seeds
 * @return int 0 on success 1 on failure
 */
int runCirclePiInner(PiResult *out, xoshiro256ss_state *state);

/**
 * @brief Runs the specified number of outer iterations for Buffon's Needle
 * 
 * @param out ptr to result struct
 * @param seeds array of thread private seeds
 * @return int 0 on success 1 on failure
 */
int runBuffonsPiInner(PiResult *out, xoshiro256ss_state *state);

int runCirclePi(uint64_t sm) {

    /* Allocate memory for thread-private seeds */

    xoshiro256ss_state *states = calloc(GLOB.n_threads, sizeof(xoshiro256ss_state));
    if (!states)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    /* Run the specified number of outer iterations for the quarter circle approximation of pi */
    
    fprintf(stdout, "Approximating pi using the 1/4 circle method with %ld outer, and %ld inner iterations...\n\n", GLOB.n_generations, GLOB.n_particles);

    /* Allocate for results array */

    double *results = calloc(GLOB.n_generations, sizeof(double));
    if (!results)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        free(states);
        exit(EXIT_FAILURE);
    }

    /* Seed thread PRNG states */

    for (long i = 0; i < GLOB.n_threads; i++)
        xoshiro256ss_seed(&states[i], sm + i * 0x9E3779B97F4A7C15ULL);
    
    /* Run outer loop */

    for (long m = 0; m < GLOB.n_generations; m++)
    {
        /* Update progress bar */

        updateProgressBar(m+1, GLOB.n_generations);

        /* Run inner loop */

        PiResult res = {0,0};
        if (runCirclePiInner(&res, states) != 0)
        {
            fprintf(stderr, "[ERROR] Failed at outer iteration %ld.\n", m);
            free(states); 
            free(results);
            exit(EXIT_FAILURE);
        }

        /* Calculate the estimate for pi */

        if (res.n_hits == 0) 
            results[m] = 0.0;
        else 
            results[m] = 4.0L * (double)res.n_hits / (double)res.n_tot;
    }

    /* Process results */

    summarizePiResultsArray(results);

    /* Free memory */

    free(results);
    free(states);

    return 0;
}

int runBuffonsPi(uint64_t sm) {

    /* Allocate memory for thread-private seeds */
    
    xoshiro256ss_state *states = calloc(GLOB.n_threads, sizeof(xoshiro256ss_state));
    if (!states)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    /* Run the specified number of outer iterations for Buffon's Needle */

    fprintf(stdout, "Approximating pi using Buffon's Needle with %ld outer, and %ld inner iterations...\n\n", GLOB.n_generations, GLOB.n_particles);

    /* Allocate results array */

    double *results = calloc(GLOB.n_generations, sizeof(double));
    if (!results) {
        fprintf(stderr, "[ERROR Memory allocation failed.\n");
        free(states);
        exit(EXIT_FAILURE);
    }

    /* Seed thread PRNG states */

    for (long i = 0; i < GLOB.n_threads; i++)
        xoshiro256ss_seed(&states[i], sm + i * 0x9E3779B97F4A7C15ULL);

    /* Outer loop */

    for (long m = 0; m < GLOB.n_generations; m++) 
    {
        /* Update progress bar */

        updateProgressBar(m + 1, GLOB.n_generations);

        /* Run inner loop */

        PiResult res = {0, 0};
        if (runBuffonsPiInner(&res, states) != 0) 
        {
            fprintf(stderr, "[WARNING] Failed at outer iteration %ld.\n", m);
            free(states);
            free(results);
            exit(EXIT_FAILURE);
        }
        /* Calculate the estimate for pi */

        if (res.n_hits == 0) 
            results[m] = 0.0;
        else 
            results[m] = (2.0 * GLOB.needle_length * res.n_tot) / (res.n_hits * GLOB.line_spacing);
    }

    /* Process results */

    summarizePiResultsArray(results);

    /* Free memory */

    free(results);
    free(states);

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

int runCirclePiInner(PiResult *out, xoshiro256ss_state *states)
{
    long n_tot = 0, n_hits = 0;

    #pragma omp parallel default(none) \
            shared(GLOB, states) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        xoshiro256ss_state *state = &states[id];

        #pragma omp for schedule(dynamic, 32000)
        for (long m = 0; m < GLOB.n_particles; ++m) {
            
            /* Generate a random point */

            double x = randd(state);
            double y = randd(state);

            n_tot++;

            /* Check if it falls inside the quarter circle */

            if (x*x + y*y < 1.0L) 
                n_hits++;
        }
    }

    /* Combine */

    out->n_tot = n_tot;
    out->n_hits = n_hits;

    return 0;
}

int runBuffonsPiInner(PiResult *out, xoshiro256ss_state *states) {
    long n_tot = 0, n_hits = 0;

    /* Precompute constant values */

    const double half_line_spacing = GLOB.line_spacing / 2.0;

    #pragma omp parallel default(none) \
            shared(GLOB, states, sin_table, stderr, half_line_spacing) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        xoshiro256ss_state *state = &states[id];

        #pragma omp for schedule(auto) 
        for (long m = 0; m < GLOB.n_particles; ++m)
        {
            /* Generate random center position and angle */

            double center_dist = randd(state) * half_line_spacing;
            double angle = randd(state) * M_PI;

            n_tot++;

            /* Look-up sin */

            int index = (int)((angle / M_PI) * (TRIG_LOOKUP_TABLE_SIZE - 1));

            /* Check if the needle crosses a line */

            double projection = sin_table[index];

            if (projection >= center_dist) {
                n_hits++;
            }
        }
    }

    /* Combine */

    out->n_tot = n_tot;
    out->n_hits = n_hits;

    return 0;
}