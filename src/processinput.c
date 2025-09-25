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

    /* Check that requirements for running a transport simulation are set */
    if (GLOB.mode == RUNMODE_TRANSPORT) 
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

    /* Check Buffon's needle related parameters */

    if (GLOB.mode != RUNMODE_BUFFONS_PI) 
    {
        if (GLOB.needle_length > 0.0 || GLOB.line_spacing > 0.0) 
        {
            fprintf(stderr, "[WARNING] Buffon's needle parameters provided but mode is not set to %d. Ignoring.\n", RUNMODE_BUFFONS_PI);
        }
    }
    else 
    {
        /* Check for valid Buffon's needle parameters */

        if (GLOB.needle_length <= 0.0) 
        {
            fprintf(stderr, "[ERROR] Invalid needle length %.3f for Buffon's needle mode.\n", GLOB.needle_length);
            return EXIT_FAILURE;
        }

        if (GLOB.line_spacing <= 0.0) 
        {
            fprintf(stderr, "[ERROR] Invalid line spacing %.3f for Buffon's needle mode.\n", GLOB.line_spacing);
            return EXIT_FAILURE;
        }

        if (GLOB.needle_length > GLOB.line_spacing) 
        {
            fprintf(stderr, "[ERROR] Needle length %.3f cannot be larger than line spacing %.3f.\n", GLOB.needle_length, GLOB.line_spacing);
            return EXIT_FAILURE;
        }
    }

    /* Return success */
    return EXIT_SUCCESS;
}
