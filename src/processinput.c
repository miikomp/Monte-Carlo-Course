#include "header.h"

int processInput() {

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
    /* Check that requirements for running any transport simulation are set */
    if (GLOB.mode == RUNMODE_CRITICALITY || GLOB.mode == RUNMODE_EXTERNAL_SOURCE) 
    {
        /* Check for valid number of iterations */

        if (!(GLOB.n_generations > 0 || GLOB.n_cycles > 0) || GLOB.n_particles < 1) {
            fprintf(stderr, "[ERROR] Number of iterations must be > 0.\n");
            return EXIT_FAILURE;
        }

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

    /* Check that parameters needed for succesfully running an external source simulation are set */

    if (GLOB.mode == RUNMODE_EXTERNAL_SOURCE)
    {
       /* 
       Check if a cut-off has been set and set default to avoid infinite loops.
       Neutron energy and collision count cutoffs are not considered here since
       the relate more to limiting what is scored than to what is simulated.
       If neither generation nor time cutoffs are set, set generation cutoff to 30.
       */

    if (!(GLOB.time_cutoff < LONG_MAX) && !(GLOB.generation_cutoff < LONG_MAX)) 
        GLOB.generation_cutoff = 30l;

    /* Unless rewritten by user, increase neutron buffer size in external source simulation */
    if (GLOB.nbuf_factor <= 1.0)
        GLOB.nbuf_factor = 10.0;
    }

    /* Return success */
    return EXIT_SUCCESS;
}
