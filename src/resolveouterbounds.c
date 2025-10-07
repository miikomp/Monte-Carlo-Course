#include "header.h"

int resolveOuterBounds() {

    fprintf(stdout, "\nCalculating outer boundaries...\n");

    /* First find the outside cell */

    int i = -1;

    for (size_t c = 0; c < DATA.n_cells; c++)
    {
        Cell *cell = &DATA.cells[c];
        if (cell->mat_idx == -1)
        {
            if (i != -1)
            {
                fprintf(stderr, "[ERROR] More than one cell defined as outside cell (material \"outside\").\n");
                return EXIT_FAILURE;
            }
            i = (int)c;
            break;
        }
    }

    /* Get the outside cell */

    Cell c0 = DATA.cells[i];

    /* Then find the surface that defines the outer bounds */
    /* Currently it is required that the outside cell confined by a single surface */
    int j = -1;
    for (size_t s = 0; s < c0.n_surfs; s++)
    {
        if (c0.side[s] > 0)
        {
            if (j != -1)
            {
                fprintf(stderr, "[ERROR] Outside cell must be bounded by a single surface.\n");
                return EXIT_FAILURE;
            }
            j = c0.surf_idxs[s];
        }
    }

    /* Get the bounding surface */

    Surface *s0 = &DATA.surfs[j];

    /* Check for valid surface type */

    if (s0->type != SURF_CYLX && s0->type != SURF_CYLY && s0->type != SURF_CYLZ &&
        s0->type != SURF_SPH && s0->type != SURF_CUBE && s0->type != SURF_CUBOID)
    {
        fprintf(stderr, "[ERROR] Unsupported surface type for outer bounds: %d.\n", s0->type);
        return EXIT_FAILURE;
    }

    /* Set the outer bounds */
    size_t n_params = s0->n_params;
    double *params = s0->params;
    SurfaceTypes type = s0->type;

    switch (type)
    {
    case SURF_CYLX:
    {
        DATA.y_min = params[0] - params[2];
        DATA.y_max = params[0] + params[2];
        DATA.z_min = params[1] - params[2];
        DATA.z_max = params[1] + params[2];

        if (n_params == 5)
        {
            DATA.x_min = params[3];
            DATA.x_max = params[4];
        }
        else
        {
            DATA.x_min = 0.0;
            DATA.x_max = 0.0;
        }

        break;
    }  
    
    case SURF_CYLY:
    {
        DATA.x_min = params[0] - params[2];
        DATA.x_max = params[0] + params[2];
        DATA.z_min = params[1] - params[2];
        DATA.z_max = params[1] + params[2];

        if (n_params == 5)
        {
            DATA.y_min = params[3];
            DATA.y_max = params[4];
        }
        else
        {
            DATA.y_min = 0.0;
            DATA.y_max = 0.0;
        }

        break;
    }

    case SURF_CYLZ:
    {
        DATA.x_min = params[0] - params[2];
        DATA.x_max = params[0] + params[2];
        DATA.y_min = params[1] - params[2];
        DATA.y_max = params[1] + params[2];

        if (n_params == 5)
        {
            DATA.z_min = params[3];
            DATA.z_max = params[4];
        }
        else
        {
            DATA.z_min = 0.0;
            DATA.z_max = 0.0;
        }

        break;
    }
    default:
    {
        fprintf(stderr, "[ERROR] Surface type %d not implemented.\n", type);
        return EXIT_FAILURE;
    }
    }

    fprintf(stdout, "X: [%lf, %lf] Y: [%lf %lf] Z: [%lf %lf]\n", DATA.x_min, DATA.x_max, DATA.y_min, DATA.y_max, DATA.z_min, DATA.z_max);

    fprintf(stdout, "DONE.\n");

    return EXIT_SUCCESS;
}