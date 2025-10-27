#include "header.h"

int resolveOuterBounds() {

    const double EPS = 1e-10;

    fprintf(stdout, "\nCalculating outer boundaries...\n");

    /* First find the outside cell in the root universe */

    int i = -1;

    for (size_t c = 0; c < DATA.unis[0].n_cells; c++)
    {
        int cell_idx = DATA.unis[0].cell_idxs[c];
        Cell *cell = &DATA.cells[cell_idx];
        if (!strcmp(cell->mat_name, "outside"))
        {
            if (i != -1)
            {
                fprintf(stderr, "[ERROR] More than one cell defined as outside cell (material \"outside\").\n");
                return EXIT_FAILURE;
            }
            i = cell_idx;
            break;
        }
    }

    /* Get the outside cell */

    if (i == -1)
    {
        fprintf(stderr, "[ERROR] Outside cell (material \"outside\") not found in root universe.\n");
        return EXIT_FAILURE;
    }

    Cell c0 = DATA.cells[i];

    /* Then find the surface that defines the outer bounds */
    /* Currently it is required that the outside cell confined by a single surface */
    int j = -1;
    for (size_t s = 0; s < c0.n_surfs; s++)
    {
        if (c0.sides[s] > 0)
        {
            if (j != -1)
            {
                fprintf(stderr, "[ERROR] Outside cell must be bounded by a single surface.\n");
                return EXIT_FAILURE;
            }
            j = c0.surf_idxs[s];
        }
    }

    if (j < 0 || j > (int)DATA.n_surf)
    {
        fprintf(stderr, "[ERROR] Outer bounding surface not found.\n");
        return EXIT_FAILURE;
    }

    /* Get the bounding surface and save its index */

    Surface *s0 = &DATA.surfs[j];
    DATA.outside_surf_idx = j;

    /* Set the outer bounds */
    size_t n_params = s0->n_params;
    double *params = s0->params;
    SurfaceTypes type = s0->type;

    switch (type)
    {
    case SURF_CYLX:
    {
        double r = params[2];
        DATA.y_min = params[0] - r;
        DATA.y_max = params[0] + r;
        DATA.z_min = params[1] - r;
        DATA.z_max = params[1] + r;

        if (n_params == 5)
        {
            DATA.x_min = params[3];
            DATA.x_max = params[4];

            DATA.tot_vol = M_PI * r * r * (DATA.x_max - DATA.x_min);
        }
        else
        {
            DATA.x_min = 0.0;
            DATA.x_max = 0.0;

            DATA.tot_vol = M_PI * r * r;
        }

        break;
    }  
    
    case SURF_CYLY:
    {
        double r = params[2];
        DATA.x_min = params[0] - r;
        DATA.x_max = params[0] + r;
        DATA.z_min = params[1] - r;
        DATA.z_max = params[1] + r;

        if (n_params == 5)
        {
            DATA.y_min = params[3];
            DATA.y_max = params[4];

            DATA.tot_vol = M_PI * r * r * (DATA.y_max - DATA.y_min);
        }
        else
        {
            DATA.y_min = 0.0;
            DATA.y_max = 0.0;

            DATA.tot_vol = M_PI * r * r;
        }

        break;
    }

    case SURF_CYLZ:
    {
        double r = params[2];
        DATA.x_min = params[0] - r;
        DATA.x_max = params[0] + r;
        DATA.y_min = params[1] - r;
        DATA.y_max = params[1] + r;

        if (n_params == 5)
        {
            DATA.z_min = params[3];
            DATA.z_max = params[4];

            DATA.tot_vol = M_PI * r * r * (DATA.z_max - DATA.z_min);
        }
        else
        {
            DATA.z_min = 0.0;
            DATA.z_max = 0.0;

            DATA.tot_vol = M_PI * r * r;
        }

        break;
    }

    case SURF_HEXX:
    {
        double d = params[2];
        DATA.x_min = params[0] - d;
        DATA.x_max = params[0] + d;
        DATA.y_min = (params[1] - d) / (0.5 * SQRT3);
        DATA.y_max = (params[1] + d) / (0.5 * SQRT3);

        if (n_params == 5)
        {
            DATA.z_min = params[3];
            DATA.z_max = params[4];

            DATA.tot_vol = (2 * SQRT3 * d * d) * (DATA.z_max - DATA.z_min);
        }
        else
        {
            DATA.z_min = 0.0;
            DATA.z_max = 0.0;

            DATA.tot_vol = (2 * SQRT3 * d * d);
        }

        break;
    }

    case SURF_HEXY:
    {
        double d  = params[2];
        DATA.x_min = (params[0] - d) / (0.5 * SQRT3);
        DATA.x_max = (params[0] + d) / (0.5 * SQRT3);
        DATA.y_min = params[1] - d;
        DATA.y_max = params[1] + d;

        if (n_params == 5)
        {
            DATA.z_min = params[3];
            DATA.z_max = params[4];

            DATA.tot_vol = (2 * SQRT3 * d * d) * (DATA.z_max - DATA.z_min);
        }
        else
        {
            DATA.z_min = 0.0;
            DATA.z_max = 0.0;

            DATA.tot_vol = (2 * SQRT3 * d * d);
        }

        break;
    }

    case SURF_SQR:
    {
        double a = params[2];
        DATA.x_min = params[0] - a;
        DATA.x_max = params[0] + a;
        DATA.y_min = params[1] - a;
        DATA.y_max = params[1] + a;
        if (n_params == 5)
        {
            DATA.z_min = params[3];
            DATA.z_max = params[4];

            DATA.tot_vol = (DATA.x_max - DATA.x_min) * (DATA.y_max - DATA.y_min) * (DATA.z_max - DATA.z_min);
        }
        else
        {
            DATA.z_min = 0.0;
            DATA.z_max = 0.0;

            DATA.tot_vol = (DATA.x_max - DATA.x_min) * (DATA.y_max - DATA.y_min);
        }

        break;
    }
    case SURF_SPH:
    {
        double r = params[3];
        DATA.x_min = -r;
        DATA.x_max =  r;
        DATA.y_min = -r;
        DATA.y_max =  r;
        DATA.z_min = -r;
        DATA.z_max =  r;

        DATA.tot_vol = 4 * M_PI * r * r * r / 3;

        break;
    }
    case SURF_CUBE:
    {
        double a = params[3];
        DATA.x_min = -a;
        DATA.x_max =  a;
        DATA.y_min = -a;
        DATA.y_max =  a;
        DATA.z_min = -a;
        DATA.z_max =  a;

        DATA.tot_vol = (DATA.x_max - DATA.x_min) * (DATA.y_max - DATA.y_min) * (DATA.z_max - DATA.z_min);

        break;
    }
    case SURF_CUBOID:
    {
        DATA.x_min = params[0];
        DATA.x_max = params[1];
        DATA.y_min = params[2];
        DATA.y_max = params[3];
        DATA.z_min = params[4];
        DATA.z_max = params[5];

        DATA.tot_vol = (DATA.x_max - DATA.x_min) * (DATA.y_max - DATA.y_min) * (DATA.z_max - DATA.z_min);

        break;
    }
    default:
    {
        fprintf(stderr, "[ERROR] Surface type %d not implemented for outer boundaries.\n", type);
        return EXIT_FAILURE;
    }
    }

    /* Check for legal geometry */
    
    if (((fabs(DATA.x_min) < EPS) && (fabs(DATA.x_max) < EPS)) || 
        ((fabs(DATA.y_min) < EPS) && (fabs(DATA.y_max) < EPS)))
    {
        fprintf(stderr, "[ERROR] Universe must be finite in the XY plane.\n");
        return EXIT_FAILURE;
    }

    fprintf(stdout, "  X: [%8.4lf, %8.4lf] cm\n  Y: [%8.4lf, %8.4lf] cm\n  Z: [%8.4lf, %8.4lf] cm\nVol: %.4lf\n", DATA.x_min, DATA.x_max, DATA.y_min, DATA.y_max, DATA.z_min, DATA.z_max, DATA.tot_vol);

    fprintf(stdout, "DONE.\n");

    return EXIT_SUCCESS;
}
