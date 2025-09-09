#include "header.h"

int main(int argc, char **argv) {
    /* ########################################################################################## */

    /* Check for valid usage of executable */

    if (argc < 2) 
    {
        fprintf(stderr, "[ERROR] Usage: %s [OPTS] filename.\n", argv[0]);
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
                fprintf(stderr, "[ERROR] Number of threads not given for \"%s\".\n", argv[i]);
                return EXIT_FAILURE;
            }

            /* Put number of threads */

            GLOB.n_threads = (int)fmin((strtol(argv[++i], NULL, 10)), 16);
        }
        
    }

    /* ########################################################################################## */

    /* Read input file */
    
    long np = readInput();
    GLOB.n_kwargs = np;

    /* Process input data */

    if (processInput() != 0) 
        return EXIT_FAILURE;

    fprintf(stdout, "[NOTE] %ld keyword arguments succesfully parsed.\nDONE. \n\n", np);

    /* ########################################################################################## */

    /* Set desired number of threads and get actual provided */

    omp_set_num_threads(GLOB.n_threads);
    int nt = omp_get_max_threads();
    GLOB.n_threads = nt;

    /* Disable dynamic teaming to not mess up thread-private seeding */

    omp_set_dynamic(0);

    /* Allocate memory for thread-private seeds */
    
    xoshiro256ss_state *states = calloc(nt, sizeof(xoshiro256ss_state));;
    if (!states)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        return EXIT_FAILURE;
    }

    /* Initialize look-up tables */

    initTrigTables();

    /* ########################################################################################## */

    /* Start timer */

    GLOB.t0 = omp_get_wtime();

    /* --- Main loop --- */

    /* Get scrambler seed from global seed */

    uint64_t sm = (GLOB.seed << 32) ^ UINT64_C(0x94D049BB133111EB);
    
    /* Dispatch case to correct sub-routine */

    switch (GLOB.mode) {
    case RUNMODE_CIRCLE_PI: 
    {
        if (runCirclePi(sm, states) != 0)
        {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            free(states);
            return EXIT_FAILURE;
        }

        break;
    }
    case RUNMODE_BUFFONS_PI: 
    {
        if (runBuffonsPi(sm, states) != 0) {
            fprintf(stderr, "[ERROR] Computation failed.\n");
            free(states);
            return EXIT_FAILURE;
        }

        break;
    }
    default: 
    {
        fprintf(stderr, "[ERROR] Mode %ld not implemented.\n", GLOB.mode);
        free(states);
        return EXIT_FAILURE;
    }
    }

    /* ########################################################################################## */

    /* Stop timer */

    GLOB.t1 = omp_get_wtime();

    /* Print runtime information */

    fprintf(stdout, "\n------------------------\n");
    fprintf(stdout, "  Runtime: %.4lfs\n", GLOB.t1 - GLOB.t0);
    fprintf(stdout, "------------------------\n\n");

    /* Prepare for termination */

    free(states);

    /* Exit */

    return EXIT_SUCCESS;
}