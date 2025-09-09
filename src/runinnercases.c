#include "header.h"

int runCirclePiInner(PiResult *out, xoshiro256ss_state *states)
{
    long n_tot = 0, n_hits = 0;

    #pragma omp parallel default(none) \
            shared(GLOB, states) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        xoshiro256ss_state *state = &states[id];

        #pragma omp for schedule(dynamic, 32000)
        for (long m = 0; m < GLOB.n_inner; ++m) {
            
            /* Generate a random point */

            double x = randd(state);
            double y = randd(state);

            n_tot++;

            /* Check if it falls inside the quarter circle */

            if (x*x + y*y < 1.0L) 
                n_hits++;
        }
    }

    /* Combine */

    out->n_tot = n_tot;
    out->n_hits = n_hits;

    return 0;
}

int runBuffonsPiInner(PiResult *out, xoshiro256ss_state *states) {
    long n_tot = 0, n_hits = 0;

    /* Precompute constant values */

    const double half_line_spacing = GLOB.line_spacing / 2.0;

    #pragma omp parallel default(none) \
            shared(GLOB, states, sin_table, stderr, half_line_spacing) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        xoshiro256ss_state *state = &states[id];

        #pragma omp for schedule(auto) 
        for (long m = 0; m < GLOB.n_inner; ++m)
        {
            /* Generate random center position and angle */

            double center_dist = randd(state) * half_line_spacing;
            double angle = randd(state) * M_PI;

            n_tot++;

            /* Look-up sin */

            int index = (int)((angle / M_PI) * (TRIG_LOOKUP_TABLE_SIZE - 1));

            /* Check if the needle crosses a line */

            double projection = sin_table[index];

            if (projection >= center_dist) {
                n_hits++;
            }
        }
    }

    /* Combine */

    out->n_tot = n_tot;
    out->n_hits = n_hits;

    return 0;
}