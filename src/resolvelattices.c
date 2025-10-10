#include "header.h"

int resolveLattices() {
    fprintf(stdout, "\nProcessing lattices...\n");

    /* First pass: ensure each lattice has a corresponding universe slot */
    
    for (size_t l = 0; l < DATA.n_lats; ++l)
    {
        Lattice *lat = &DATA.lats[l];

        if (lat->uni_idx >= 0)
            continue;

        for (size_t u = 0; u < DATA.n_unis; ++u)
        {
            if (!strcmp(DATA.unis[u].name, lat->name))
            {
                fprintf(stderr, "[ERROR] Universe name '%s' clashes with lattice definition.\n",
                        lat->name);
                return EXIT_FAILURE;
            }
        }

        Universe *tmp = (Universe*)realloc(DATA.unis, (DATA.n_unis + 1) * sizeof(Universe));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed while adding lattice universe '%s'.\n",
                    lat->name);
            return EXIT_FAILURE;
        }
        DATA.unis = tmp;

        Universe *uni = &DATA.unis[DATA.n_unis];
        memset(uni, 0, sizeof(*uni));
        snprintf(uni->name, sizeof(uni->name), "%s", lat->name);
        uni->type = UNI_LATTICE;
        uni->lat_idx = (int)l;
        uni->t_idx = -1;
        uni->n_cells = 0;
        uni->cell_idxs = NULL;

        lat->uni_idx = (int)DATA.n_unis;
        DATA.n_unis++;
    }

    /* Second pass: resolve fills */
    for (size_t l = 0; l < DATA.n_lats; ++l)
    {
        Lattice *lat = &DATA.lats[l];

        switch (lat->type)
        {
            case LAT_SQUARE_INFINITE:
            {
                const char *fill_name = lat->uni_names;
                int fill_idx = -1;
                for (size_t u = 0; u < DATA.n_unis; ++u)
                {
                    if (!strcmp(fill_name, DATA.unis[u].name))
                    {
                        fill_idx = (int)u;
                        break;
                    }
                }

                if (fill_idx < 0)
                {
                    fprintf(stderr, "[ERROR] Lattice '%s' references unknown universe '%s'.\n",
                            lat->name, fill_name);
                    return EXIT_FAILURE;
                }

                lat->uni_idxs = (long*)calloc(1, sizeof(long));
                if (!lat->uni_idxs)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed while resolving lattice '%s'.\n",
                            lat->name);
                    return EXIT_FAILURE;
                }
                lat->uni_idxs[0] = fill_idx;

                break;
            }
            case LAT_SQUARE_FINITE:
            {
                /* Resolve universe indeces */

                lat->uni_idxs = (long*)calloc(lat->n_unis, sizeof(long));
                memset(lat->uni_idxs, -1, sizeof(long) * lat->n_unis);

                for (size_t n = 0; n < lat->n_unis; n++)
                {
                    char* uname = lat->uni_names + MAX_STR_LEN * n;

                    if (!strcmp(uname, lat->name))
                    {
                        fprintf(stderr, "[ERROR] Lattice \"%s\" references itself.\n", lat->name);
                        return EXIT_FAILURE;
                    }

                    for (size_t u = 0; u < DATA.n_unis; u++)
                    {
                        if (!strcmp(uname, DATA.unis[u].name))
                            lat->uni_idxs[n] = u;
                    }
                    if (lat->uni_idxs[n] < 0)
                    {
                        fprintf(stderr, "[ERROR] Universe \"%s\" for lattice \"%s\" not found.\n", uname, lat->name);
                        return EXIT_FAILURE;
                    }
                }

                break;
            }
            default:
                fprintf(stderr, "[ERROR] Unsupported lattice type %d for lattice '%s'.\n",
                        (int)lat->type, lat->name);
                return EXIT_FAILURE;
        }

        Universe *uni = &DATA.unis[lat->uni_idx];
        uni->type = UNI_LATTICE;
        uni->lat_idx = (int)l;
        uni->n_cells = 0;
        free(uni->cell_idxs);
        uni->cell_idxs = NULL;

        if (VERBOSITY >= 1)
        {
            fprintf(stdout,
                    "  Lattice %zu: %s (type=%d, origin=(%.3f, %.3f, %.3f), pitch=(%.3f, %.3f, %.3f), size=(%ld %ld %ld))\n",
                    l, lat->name, (int)lat->type,
                    lat->x0, lat->y0, lat->z0,
                    lat->dx, lat->dy, lat->dz,
                    lat->nx, lat->ny, lat->nz);
        }
        if (VERBOSITY >= 2 && lat->nz == 1)
        {
            /* Print 2D lattice ofr debugging */
 
            for (long i = 0; i < lat->nx; i++)
            {
                fprintf(stdout, "  ");

                for (long j = 0; j < lat->ny; j++)
                {
                    fprintf(stdout, "%ld ", lat->uni_idxs[j * lat->nx + i]);
                }
                fprintf(stdout, "\n");
            }
        }
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}
