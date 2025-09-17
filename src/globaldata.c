#include "header.h"

/* Initialize global data structs */

runInfo GLOB = {
    .fname          = NULL,
    .outfname       = NULL,
    .errfname       = NULL,
    .xslibpath      = {0},
    .n_threads      = 1,
    .n_generations  = 0l,
    .n_particles    = 0l,
    .n_kwargs       = 0l,
    .mode           = 0l,
    .seed           = 0u,
    .t0             = 0.0,
    .t1             = 0.0,
    .needle_length  = 0.0,
    .line_spacing   = 0.0,
};

runData DATA = { 
    .n_mats = 0ul,
    .mats   = NULL,
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