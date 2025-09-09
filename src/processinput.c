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
        fprintf(stdout, "[NOTE] Using random seed %lu.\n", GLOB.seed);
    }
    else 
    {
        fprintf(stdout, "[NOTE] Using provided seed %lu.\n", GLOB.seed);
    }

    /* Check if runmode specific parameters */
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
            return 1;
        }

        if (GLOB.line_spacing <= 0.0) 
        {
            fprintf(stderr, "[ERROR] Invalid line spacing %.3f for Buffon's needle mode.\n", GLOB.line_spacing);
            return 1;
        }

        if (GLOB.needle_length > GLOB.line_spacing) 
        {
            fprintf(stderr, "[ERROR] Needle length %.3f cannot be larger than line spacing %.3f.\n", GLOB.needle_length, GLOB.line_spacing);
            return 1;
        }
    }

    /* Return success */
    return 0;
}