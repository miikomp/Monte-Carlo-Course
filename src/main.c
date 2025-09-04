#include "header.h"

int main(int argc, char **argv) {
    long n, m;
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

            GLOB.n_threads = (int)strtol(argv[++i], NULL, 10);
        }
        
    }

    /* ########################################################################################## */

    /* Read input file */
    
    np = readInput();
    GLOB.n_kwargs = np;
    fprintf(stdout, "DONE. %ld keyword arguments succesfully parsed\n", np);

    /* Get random seed if not provided */
    
    if (!GLOB.seed) 
    {
        GLOB.seed = (uint64_t)time(NULL);
        fprintf(stdout, "Using random seed %lu...\n", GLOB.seed);
    }

    srand((uint32_t)GLOB.seed);

    /* ########################################################################################## */

    /* Set desired number of threads and get actual */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();

    /* Allocate memory for thread private data */

    Tallies *partials = calloc(nt, sizeof(Tallies));
    if (!partials)
    {
        fprintf(stderr, "Memory allocation error\n");
        return EXIT_FAILURE;
    }

    /* Create thread-private seeds */
    
    uint64_t *seeds = calloc(nt, sizeof(uint64_t));
    if (!seeds)
    {
        fprintf(stderr, "Memory allocation error\n");
        free(partials);
        return EXIT_FAILURE;
    }

    for (long i = 0; i < nt; i++)
        seeds[i] = splitmix64(GLOB.seed);
    
    /* ########################################################################################## */

    /* Start timer */
    double t0 = omp_get_wtime();

    /* --- Main loop --- */

    #pragma omp parallel default(none) \
        private (n, m) \
        shared  ( partials, GLOB, seeds)
    {
        int id = omp_get_thread_num();

        Tallies local = initTallies();

        /* Get thread-local seed */

        uint64_t rs = seeds[id];

        /* Loop over outer and inner iterations */

        #pragma omp for schedule(static) collapse(2)
        for (n = 0; n < GLOB.n_outer; n++)
        {
            for (m = 0; m < GLOB.n_inner; m++)
            {
                local.n_tot += 2;
                local.n_hits += 1;
                local.dis += randd(&rs);
            }
        }

        /* Write cumulated thread-private data */

        partials[id].n_tot  += local.n_tot;
        partials[id].n_hits += local.n_hits;
        partials[id].dis    += local.dis;

    }

    /* Merge thread-private data */

    Tallies total = initTallies();
    for (int t = 0; t < nt; t++)
    {
        total.n_tot     += partials[t].n_tot;
        total.n_hits    += partials[t].n_hits;
        total.dis       += partials[t].dis;
    }

    /* ########################################################################################## */
    
    /* Process results */

    fprintf(stdout, "tot=%lld dis=%lf\n", total.n_tot, total.dis);

    /* Stop timer */

    double t1 = omp_get_wtime();
    fprintf(stdout, "Runtime: %.4lf s\n", t1 - t0);

    /* ########################################################################################## */

    /* Prepare for termination */

    free(seeds);
    free(partials);

    /* Exit */

    return EXIT_SUCCESS;
}