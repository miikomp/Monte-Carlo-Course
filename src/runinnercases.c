#include "header.h"

int runCirclePiInner(PiResult *out, uint64_t *seeds)
{
    long n_tot = 0, n_hits = 0;

    #pragma omp parallel default(none) \
            shared(GLOB, seeds) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        uint64_t rs = seeds[id];

        #pragma omp for schedule(static)
        for (long m = 0; m < GLOB.n_inner; ++m) {
            
            /* Generate a random point */

            uint64_t rand64 = xorShift64(&rs);

            /* Split into two random doubles */
            
            double x = (double)(rand64 >> 32) * INV_INT32_MAX;
            double y = (double)(rand64 & 0xFFFFFFFF) * INV_INT32_MAX;

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

int runBuffonsPiInner(PiResult *out, uint64_t *seeds) {
    long n_tot = 0, n_hits = 0;

    /* Precompute constant values */

    const double half_line_spacing = GLOB.line_spacing / 2.0;

    #pragma omp parallel default(none) \
            shared(GLOB, seeds, sin_table, stderr, half_line_spacing) \
            reduction(+:n_tot, n_hits)
    {
        int id = omp_get_thread_num();
        uint64_t rs = seeds[id];

        #pragma omp for schedule(static)
        for (long m = 0; m < GLOB.n_inner; ++m)
        {
            /* Generate random center position and angle */

            uint64_t rand64 = xorShift64(&rs);

            /* Split into two random doubles */

            double center_dist = (double)(rand64 >> 32) * INV_INT32_MAX * half_line_spacing;
            double angle = (double)(rand64 & 0xFFFFFFFF) * INV_INT32_MAX * M_PI;;

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