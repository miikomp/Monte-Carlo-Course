#include "header.h"

int runCirclePi(PiResult *out, uint64_t *seeds)
{
    PiResult total = {0,0};
    #pragma omp parallel default(none) \
            shared(total, GLOB, seeds)
    {
        int id = omp_get_thread_num();
        uint64_t rs = seeds[id];
        long long n_tot = 0, n_hits = 0;

        #pragma omp for schedule(static)
        for (long m = 0; m < GLOB.n_inner; ++m) {
            
            /* Generate a random point */

            long double x = randd(&rs);
            long double y = randd(&rs);

            n_tot++;

            /* Check if it falls inside the quarter circle */

            if (x*x + y*y < 1.0L) 
            n_hits++;
        }
        

        /* Combine once per thread */

        #pragma omp atomic
        total.n_tot += n_tot;

        #pragma omp atomic
        total.n_hits += n_hits;
    }

    *out = total;
    return 0;
}