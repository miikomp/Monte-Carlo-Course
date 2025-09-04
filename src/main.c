#include "header.h"

int main(int argc, char **argv) {
    uint64_t np;

    /* ########################################################################################## */

    /* Check for valid usage of executable */

    if (argc < 2) 
    {
        fprintf(stderr, "[ERROR] Usage: %s [OPTS] filename\n", argv[0]);
        return EXIT_FAILURE;
    }

    /* Filename is the last argument */

    GLOB.fname = argv[argc - 1];

    /* Try to parse commandline arguments */
    
    for (long i = 1; i < argc - 1; i++) 
    {
        if (!strcmp(argv[i], "-omp")) 
        {
            if (i + 1 >= argc - 1) 
            {
                fprintf(stderr, "[ERROR] Number of threads not given for \"%s\"\n", argv[i]);
                return EXIT_FAILURE;
            }

            /* Put number of threads */

            GLOB.n_threads = (int)fmin((strtol(argv[++i], NULL, 10)), 16);
        }
        
    }

    /* ########################################################################################## */

    /* Read input file */
    
    np = readInput();

    if (GLOB.n_outer < 1 || GLOB.n_inner < 1)
    {
        fprintf(stdout, "Number of iterations not given!\n");
        return EXIT_FAILURE;
    }

    GLOB.n_kwargs = np;
    fprintf(stdout, "DONE. %ld keyword arguments succesfully parsed\n\n", np);

    /* Get random seed if not provided */
    
    if (!GLOB.seed) 
    {
        GLOB.seed = (uint64_t)time(NULL);
        fprintf(stdout, "Using random seed %lu...\n\n", GLOB.seed);
    }

    /* ########################################################################################## */

    /* Set desired number of threads and get actual */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();

    /* Disable dynamic teaming */

    omp_set_dynamic(0);

    /* Thread-private seeds */
    
    uint64_t *seeds = calloc(nt, sizeof(uint64_t));
    if (!seeds)
    {
        fprintf(stderr, "Memory allocation error\n");
        return EXIT_FAILURE;
    }

    /* ########################################################################################## */

    /* Start timer */

    double t0 = omp_get_wtime();

    /* --- Main loop --- */

    uint64_t sm = ((uint64_t)GLOB.seed << 32) ^ UINT64_C(0x94D049BB133111EB);
    
    switch (GLOB.mode) {
    case RUNMODE_CIRCLE_PI: 
    {
        fprintf(stdout, "Approximating pi with %ld outer, and %ld inner iterations...\n\n", GLOB.n_outer, GLOB.n_inner);
        long double *results = calloc(GLOB.n_outer, sizeof(long double));
        if (!results)
        {
            fprintf(stderr, "Memory allocation error\n");
            free(seeds);
            return EXIT_FAILURE;
        }
        
        for (long m = 0; m < GLOB.n_outer; m++)
        {
            /* Derive seeds for run */

            for (long i = 0; i < nt; i++) 
            {
                seeds[i] = splitmix64(&sm);
                if (!seeds[i]) 
                    seeds[i] = 0x9E3779B97F4A7C15ULL;
            }

            PiResult res = {0,0};
            if (runCirclePi(&res, seeds) != 0)
            {
                fprintf(stderr, "Failed at outer iteration %ld\n", m);
                free(seeds); 
                free(results);
                return EXIT_FAILURE;
            }

            results[m] = 4.0L * (long double)res.n_hits / (long double)res.n_tot;
        }

        /* Process results */

        summarizeResultsArray(results);

        break;
    }
    default: 
    {
        fprintf(stderr, "Mode %ld not implemented\n", GLOB.mode);

        break;
    }
    }

    /* ########################################################################################## */

    /* Stop timer */

    double t1 = omp_get_wtime();
    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Runtime: %.4lfs\n", t1 - t0);
    fprintf(stdout, "------------------------\n");

    /* Prepare for termination */

    free(seeds);

    /* Exit */

    return EXIT_SUCCESS;
}