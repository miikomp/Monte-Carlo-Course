#include "header.h"

int runExternalSourceSimulation(void) {

    /* Check if a cut-off has been set and set default to avoid infinite loops.
       Neutron energy and collision count cutoffs are not considered here since
       the relate more to limiting what is scored than to what is simulated.
       If neither generation nor time cutoffs are set, set generation cutoff to 30.
    */

    if (!(GLOB.time_cutoff < LONG_MAX) && !(GLOB.generation_cutoff < LONG_MAX)) 
        GLOB.generation_cutoff = 30l;
    return EXIT_FAILURE;
}