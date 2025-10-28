#include "header.h"

int checkVolumes2()
{
    long n_lines = GLOB.n_lines;
    if (n_lines <= 0)
        return EXIT_SUCCESS;

    fprintf(stdout, "\nChecking material volumes by sampling %ld random lines...\n", n_lines);

    double t0 = omp_get_wtime();

    /* Allocate thread local arrays */

    double **length_bufs = (double **)calloc((size_t)GLOB.n_threads, sizeof(double *));
    if (!length_bufs)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    double *tot_lengths = (double *)calloc((size_t)GLOB.n_threads, sizeof(double));
    if (!tot_lengths)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        free(length_bufs);
        return EXIT_FAILURE;
    }

    for (int t = 0; t < GLOB.n_threads; ++t)
    {
        length_bufs[t] = (double *)calloc(DATA.n_mats, sizeof(double));
        if (!length_bufs[t])
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            for (int j = 0; j < t; ++j)
                free(length_bufs[j]);
            free(length_bufs);
            free(tot_lengths);
            return EXIT_FAILURE;
        }
    }

    double xmin = DATA.x_min;
    double xmax = DATA.x_max;
    double ymin = DATA.y_min;
    double ymax = DATA.y_max;
    double zmin = DATA.z_min;
    double zmax = DATA.z_max;

    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;

    bool planar = (fabs(dz) < 1e-12);
    double step_eps = 1e-6;

    #pragma omp parallel default(none) \
            shared(length_bufs, tot_lengths, n_lines, xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz, planar, step_eps, DATA, stdout)
    {
        int tid = omp_get_thread_num();
        double *lengths = length_bufs[tid];
        double total = 0.0;

        /* Initialise a random state */

        xoshiro256ss_state state;
        uint64_t seed = (uint64_t)time(NULL) ^ ((uint64_t)tid << 32);
        xoshiro256ss_seed(&state, seed);

        /* Allocate a material-wise line length arrays */

        double *line_lengths = (double *)malloc(DATA.n_mats * sizeof(double));
        if (!line_lengths)
            exit(EXIT_FAILURE);
        
        #pragma omp for schedule(static)
        for (long i = 0; i < n_lines; i++)
        {
            /* Sample a random point in XYZ or XY within the defined geometry*/

            double x, y, z, f;
            Surface *bb_surf = &DATA.surfs[DATA.outside_surf_idx];
            
            do 
            {
                x = xmin + randd(&state) * dx;
                y = ymin + randd(&state) * dy;
                z = planar ? 0.0 : (zmin + randd(&state) * dz);

                f = surfaceTest(bb_surf->type, bb_surf->params, bb_surf->n_params, x, y, z);
            } 
            while (f >= 0.0);

            /* Clear scores */

            memset(line_lengths, 0, DATA.n_mats * sizeof(double));
            double total_len = 0.0;

            double dirs[2][3];

            /* Sample a random direction in 3D or 2D */

            if (planar)
            {
                double phi = 2.0 * M_PI * randd(&state);
                double u = cos(phi);
                double v = sin(phi);
                dirs[0][0] =  u;
                dirs[0][1] =  v; 
                dirs[0][2] = 0.0;
                dirs[1][0] = -u; 
                dirs[1][1] = -v; 
                dirs[1][2] = 0.0;
            }
            else
            {
                double u, v, w;
                sampleIsotropicDirection(&state, &u, &v, &w);
                dirs[0][0] =  u; 
                dirs[0][1] =  v; 
                dirs[0][2] =  w;
                dirs[1][0] = -u; 
                dirs[1][1] = -v; 
                dirs[1][2] = -w;
            }

            /* Track ray forwards and backwards to geometry boundaries */

            for (int k = 0; k < 2; k++)
            {
                double ru = dirs[k][0];
                double rv = dirs[k][1];
                double rw = dirs[k][2];

                double rx = x;
                double ry = y;
                double rz = z;

                double dir_len = 0.0;

                while (1)
                {
                    int err;
                    long cell_idx = cellSearch(rx, ry, rz, &err, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
                    if (cell_idx < 0 || err != CELL_ERR_OK || cell_idx >= (long)DATA.n_cells)
                        break;

                    double d = distanceToNearestBoundary(rx, ry, rz, ru, rv, rw);

                    if (!(d > 0.0 && isfinite(d)))
                        break;

                    Cell *cell = &DATA.cells[cell_idx];
                    if (cell->mat_idx >= 0)
                    {
                        line_lengths[cell->mat_idx] += d;
                        dir_len += d;
                    }
                    
                    rx += ru * d;
                    ry += rv * d;
                    rz += rw * d;

                    rx += ru * step_eps;
                    ry += rv * step_eps;
                    rz += rw * step_eps;
                }

                total_len += dir_len;
            }

            /* Normalise material-wise lengths */

            if (total_len > 0.0)
            {
                double inv = 1.0 / total_len;
                for (size_t m = 0; m < DATA.n_mats; ++m)
                    lengths[m] += line_lengths[m] * inv;
                total += 1.0;
            }
        }

        tot_lengths[tid] = total;
        free(line_lengths);
    }

    /* Combine buffers */

    double *lengths = (double *)calloc(DATA.n_mats, sizeof(double));
    double total_length = 0.0;

    for (int t = 0; t < GLOB.n_threads; ++t)
    {
        double *buffer = length_bufs[t];
        for (size_t m = 0; m < DATA.n_mats; ++m)
            lengths[m] += buffer[m];
        free(buffer);
        total_length += tot_lengths[t];
    }

    free(length_bufs);
    free(tot_lengths);

    /* Get bounding surface volume */

    double vol0 = DATA.tot_vol;

    if (vol0 < 0.0)
    {
        fprintf(stderr, "[ERROR] Invalid bounding box volume.\n");
        return EXIT_FAILURE;
    }

    if (total_length == 0.0)
    {
        fprintf(stderr, "[ERROR] Volume sampling produced zero valid tracks.\n");
        free(lengths);
        return EXIT_FAILURE;
    }

    /* Calculate and put results */

    for (size_t m = 0; m < DATA.n_mats; ++m)
    {
        Material *mat = &DATA.mats[m];
        double p = lengths[m] / total_length;
        mat->vol = vol0 * p;

        double rel_unc = (p > 0.0) ? sqrt(fmax(0.0, (1.0 - p) / (p * total_length))) : 0.0;
        fprintf(stdout, "%10s: %.4E  +/- %.3E\n", mat->name, mat->vol, rel_unc);
    }

    free(lengths);

    double ttot = omp_get_wtime() - t0;
    fprintf(stdout, "\nTime elapsed: %.3f s\n", ttot);
    fprintf(stdout, "DONE.\n");

    return EXIT_SUCCESS;
}
