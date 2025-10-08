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
            if (!strcmp(DATA.unis[u].name, lat->uni_name))
            {
                fprintf(stderr, "[ERROR] Universe name '%s' clashes with lattice definition.\n",
                        lat->uni_name);
                return EXIT_FAILURE;
            }
        }

        Universe *tmp = (Universe*)realloc(DATA.unis, (DATA.n_unis + 1) * sizeof(Universe));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed while adding lattice universe '%s'.\n",
                    lat->uni_name);
            return EXIT_FAILURE;
        }
        DATA.unis = tmp;

        Universe *uni = &DATA.unis[DATA.n_unis];
        memset(uni, 0, sizeof(*uni));
        snprintf(uni->name, sizeof(uni->name), "%s", lat->uni_name);
        uni->type = UNI_LATTICE;
        uni->lat_idx = (int)l;
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
                if (lat->n_unis != 1)
                {
                    fprintf(stderr, "[ERROR] Infinite square lattice '%s' must define exactly one fill universe.\n",
                            lat->uni_name);
                    return EXIT_FAILURE;
                }

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
                            lat->uni_name, fill_name);
                    return EXIT_FAILURE;
                }

                lat->uni_idxs = (int*)calloc(1, sizeof(int));
                if (!lat->uni_idxs)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed while resolving lattice '%s'.\n",
                            lat->uni_name);
                    return EXIT_FAILURE;
                }
                lat->uni_idxs[0] = fill_idx;

                lat->nx = lat->ny = lat->nz = -1;
                lat->dz = 0.0;

                break;
            }
            default:
                fprintf(stderr, "[ERROR] Unsupported lattice type %d for lattice '%s'.\n",
                        (int)lat->type, lat->uni_name);
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
            const char *fill_name = (lat->n_unis > 0) ? lat->uni_names : "";
            fprintf(stdout,
                    "  Lattice %zu: %s (type=%d, origin=(%.3f, %.3f, %.3f), pitch=(%.3f, %.3f, %.3f), fill=%s)\n",
                    l, lat->uni_name, (int)lat->type,
                    lat->x0, lat->y0, lat->z0,
                    lat->dx, lat->dy, lat->dz,
                    fill_name);
        }
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}
