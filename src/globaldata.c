#include "header.h"

int VERBOSITY = 0;

/* Initialize global data structs */

runInfo GLOB = {
    .fname          = NULL,
    .outfname       = NULL,
    .errfname       = NULL,
    .xslibpath      = {0},
    .n_kwargs       = 0l,
    .mode           = 0l,
    .seed           = 0u,
    .t0             = 0.0,
    .t1             = 0.0,
    .n_threads      = 1,

    .n_generations  = 0l,
    .n_particles    = 0l,

    .energy_cutoff  = 0.0,

    .needle_length  = 0.0,
    .line_spacing   = 0.0,
};

runData DATA = { 
    .generation = 1ul,
    .n_mats = 0ul,
    .mats   = NULL,
    .n_bank = 0ul,
    .bank_cap = 0ul,
    .bank   = NULL,
    .src = NULL,
    .src_type = 0,
    .detector = NULL,
    .detector_type = 0
};

ResultsData RES = {
    .n_generations = 0l,
    .avg_scores = {0.0, 0.0, 0u, 0u}
};

double sin_table[TRIG_LOOKUP_TABLE_SIZE];
double cos_table[TRIG_LOOKUP_TABLE_SIZE];
double tan_table[TRIG_LOOKUP_TABLE_SIZE];

void initTrigTables() {
    for (int i = 0; i < TRIG_LOOKUP_TABLE_SIZE; i++) 
    {
        double angle = (M_PI * i) / (TRIG_LOOKUP_TABLE_SIZE - 1);
        if (GLOB.mode == RUNMODE_BUFFONS_PI) 
        {   
            /* Scale sin by needle length / 2 for Buffon's needle problem */

            sin_table[i] = (GLOB.needle_length / 2.0) * sin(angle);
            cos_table[i] = cos(angle);
            tan_table[i] = tan(angle);
            continue;
        }
        else 
        {
            /* Just put raw sin, cos and tan for other runmodes */

            sin_table[i] = sin(angle);
            cos_table[i] = cos(angle);
            tan_table[i] = tan(angle);
            continue;
        }
    }
}
