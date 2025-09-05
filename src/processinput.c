#include "header.h"

int processInput() {

    /* Check for valid number of iterations */

    if (GLOB.n_outer < 1 || GLOB.n_inner < 1) {
        fprintf(stderr, "[ERROR] Number of iterations must be > 0.\n");
        return 1;
    }

    /* Get and put random seed if not provided */
    
    if (!GLOB.seed) 
    {
        GLOB.seed = (uint64_t)time(NULL);
        fprintf(stdout, "Using random seed %lu...\n\n", GLOB.seed);
    }

    return 0;
}