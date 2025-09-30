#include "header.h"

int processInput() {

    /* Check for valid number of iterations */
    
    if (GLOB.n_generations < 1 || GLOB.n_particles < 1) {
        fprintf(stderr, "[ERROR] Number of iterations must be > 0.\n");
        return EXIT_FAILURE;
    }

    /* Get and put random seed if not provided */
    
    if (!GLOB.seed) 
    {
        GLOB.seed = (uint64_t)time(NULL);
        if (VERBOSITY >= 1)
            fprintf(stdout, "Using random seed %lu.\n", GLOB.seed);
    }
    else 
    {
        if (VERBOSITY >= 1)
            fprintf(stdout, "Using provided seed %lu.\n", GLOB.seed);
    }

    /* Initialize RNG state for global use */
    
    xoshiro256ss_seed(&GLOB.rng_state, GLOB.seed);

    /* ########################################################################################## */
    /* Check that requirements for running a transport simulation are set */
    if (GLOB.mode == RUNMODE_CRITICALITY) 
    {
        /* Check that atleast one material is defined */

        if (DATA.n_mats < 1) 
        {
            fprintf(stderr, "[ERROR] No materials provided for transport mode.\n");
            return EXIT_FAILURE;
        }

        /* Check that cross section library path is provided */

        if (GLOB.xslibpath[0] == '\0') 
        {
            fprintf(stderr, "[ERROR] No cross section library path provided for transport mode.\n");
            return EXIT_FAILURE;
        }
    }

    /* Return success */
    return EXIT_SUCCESS;
}
