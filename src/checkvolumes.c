#include "header.h"

int checkVolumes() {
    long n_points = GLOB.n_points;
    
    if (n_points <= 0)
        return EXIT_SUCCESS;

    fprintf(stdout, "\nChecking material volumes by sampling %ld random points...\n", n_points);

    double t0 = omp_get_wtime();

    /* Allocate thread local buffers */

    double **vol_bufs = (double **)calloc((size_t)GLOB.n_threads, sizeof(double*));
    if (!vol_bufs)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    long *hit_counts = (long *)calloc((size_t)GLOB.n_threads, sizeof(long));
    if (!hit_counts)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        free(vol_bufs);
        return EXIT_FAILURE;
    }

    for (int t = 0; t < GLOB.n_threads; t++)
    {
        vol_bufs[t] = (double *)calloc(DATA.n_mats, sizeof(double));
        if (!vol_bufs[t])
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            for (int j = 0; j < t; ++j)
                free(vol_bufs[j]);
            free(vol_bufs);
            free(hit_counts);
            return EXIT_FAILURE;
        }
    }

    #pragma omp parallel default(none) \
            shared(DATA, vol_bufs, hit_counts, n_points)
    {
        /* Get outer boundaries to sample within */

        double xmin = DATA.x_min;
        double xmax = DATA.x_max;
        double ymin = DATA.y_min;
        double ymax = DATA.y_max;
        double zmin = DATA.z_min;
        double zmax = DATA.z_max;
        
        /* Get thread-local buffer */

        int tid = omp_get_thread_num();

        double *volumes = vol_bufs[tid];

        /* Get random state for thread */

        uint64_t seed = (uint64_t)time(NULL) ^ ((uint64_t)tid << 32);
        xoshiro256ss_state state;
        xoshiro256ss_seed(&state, seed);

        /* Loop over the specified number of points */
        
        #pragma omp for schedule(dynamic)
        for (long i = 0; i < n_points; i++)
        {
            /* Get random point */

            double x = xmin + randd(&state) * (xmax - xmin);
            double y = ymin + randd(&state) * (ymax - ymin);
            double z = zmin + randd(&state) * (zmax - zmin);

            /* Get current cell */
            cellSearchRes res = cellSearch(x, y, z, 0.0, 0.0, 0.0);
            long c = res.cell_idx;
            int err = res.err;
            if (c < 0 || err != CELL_ERR_OK)
                continue;

            Cell *cell = &DATA.cells[c];

            /* Get material idx */

            int m = cell->mat_idx;

            if (m >= 0)
            {
                volumes[m] += 1.0;
                hit_counts[tid]++;
            }
        }
    } // End of parallel region

    /* Combine buffers */

    double *volumes = (double *)calloc(DATA.n_mats, sizeof(double));
    long total_hits = 0;

    for (int t = 0; t < GLOB.n_threads; t++)
    {
        double *buffer = vol_bufs[t];

        for (size_t m = 0; m < DATA.n_mats; m++)
            volumes[m] += buffer[m];
        free(buffer);
        total_hits += hit_counts[t];
    }

    free(vol_bufs);
    free(hit_counts);

    /* Calculate total volume */

    double dimx = DATA.x_max - DATA.x_min;
    double dimy = DATA.y_max - DATA.y_min;
    double dimz = DATA.z_max - DATA.z_min;
    double vol0;

    if (dimx == 0.0)
        vol0 = dimy * dimz;
    else if (dimy == 0.0)
        vol0 = dimx * dimz;
    else if (dimz == 0.0)
        vol0 = dimx * dimy;
    else
        vol0 = dimx * dimy * dimz;

    if (total_hits == 0)
    {
        fprintf(stderr, "[ERROR] Volume sampling produced zero valid hits.\n");
        free(volumes);
        return EXIT_FAILURE;
    }

    if (VERBOSITY >= 1)
        fprintf(stdout, "Sampling efficiency: %.2f%%\n", 100.0 * (double)total_hits / (double)n_points);
    fprintf(stdout, "\n");

    /* Divide by total samples, scale to bounding volume and put to material */

    for (size_t m = 0; m < DATA.n_mats; m++)
    {
        Material *mat = &DATA.mats[m];
        double hits = volumes[m];
        double p = hits / (double)n_points;
        mat->vol = vol0 * p;

        double rel_unc = (p > 0.0) ? sqrt(fmax(0.0, (1.0 - p) / (p * (double)n_points))) : 0.0;
        fprintf(stdout, "%10s: %.4E  +/- %.3E\n", mat->name, mat->vol, rel_unc);
    }

    free(volumes);

    double ttot = omp_get_wtime() - t0;

    fprintf(stdout, "\nTime elapsed: %.3f s\n", ttot);
    fprintf(stdout, "DONE.\n");

    return EXIT_SUCCESS;
}
