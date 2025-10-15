#include "header.h"

enum {
    PLOT_YZ = 1,
    PLOT_XZ,
    PLOT_XY
};

enum {
    BOUNDS_NONE,
    BOUNDS_MATERIAL,
    BOUNDS_CELL,
};

long plotGeometry() {

    if (GLOB.noplot)
        return EXIT_SUCCESS;

    fprintf(stdout, "\nPlotting geometry...\n");

    /* Loop over geometryPlotters */

    for (size_t g = 0; g < DATA.n_gpls; g++)
    {
        GeometryPlotter *gpl = &DATA.gpls[g];

        /* Get plot boundaries or use global bounds */

        double min1, max1, min2, max2;

        switch (gpl->axis)
        {
            case PLOT_YZ:
            {
                min1 = (gpl->ymin != 0) ? gpl->ymin : DATA.y_min;
                max1 = (gpl->ymax != 0) ? gpl->ymax : DATA.y_max;
                min2 = (gpl->zmin != 0) ? gpl->zmin : DATA.z_min;
                max2 = (gpl->zmax != 0) ? gpl->zmax : DATA.z_max;

                break;
            }
            case PLOT_XZ:
            {
                min1 = (gpl->xmin != 0) ? gpl->xmin : DATA.x_min;
                max1 = (gpl->xmax != 0) ? gpl->xmax : DATA.x_max;
                min2 = (gpl->zmin != 0) ? gpl->zmin : DATA.z_min;
                max2 = (gpl->zmax != 0) ? gpl->zmax : DATA.z_max;

                break;
            }
            case PLOT_XY:
            {
                min1 = (gpl->xmin != 0) ? gpl->xmin : DATA.x_min;
                max1 = (gpl->xmax != 0) ? gpl->xmax : DATA.x_max;
                min2 = (gpl->ymin != 0) ? gpl->ymin : DATA.y_min;
                max2 = (gpl->ymax != 0) ? gpl->ymax : DATA.y_max;

                break;
            }
            default:
            {
                fprintf(stderr, "WTF?\n");
                return EXIT_FAILURE;
            }
        }
        
        /* Calculate correct pixel counts to preserve aspect ratio */

        long pix = gpl->pixx;
        double dimx = max1 - min1;
        double dimy = max2 - min2;

        if ((dimx) > (dimy))
        {
            gpl->pixx = (dimx > 0.0 && dimy > 0.0) ? pix * (dimx / dimy) : pix;
            gpl->pixy = pix;
        }
        else
        {
            gpl->pixx = pix;
            gpl->pixy = (dimx > 0.0 && dimy > 0.0) ? pix * (dimy / dimx) : pix;
        }

        /* Create image file (.ppm format) */

        char filename[MAX_STR_LEN];
        snprintf(filename, sizeof(filename), "%s_geom%zu.ppm", GLOB.inputfname, g + 1);

        FILE *imgf = fopen(filename, "w");
        if (!imgf)
        {
            fprintf(stderr, "[ERROR] Could not open file \"%s\" for writing", filename);
            return EXIT_FAILURE;
        }

        fprintf(imgf, "P3\n%ld %ld\n255\n", gpl->pixx, gpl->pixy);

        /* Plot image */

        double di = (dimx > 0.0) ? dimx / (double)gpl->pixx : 0.0;
        double dj = (dimy > 0.0) ? dimy / (double)gpl->pixy : 0.0;

        double x, y, z;
        x = y = z = 0.0;

        /* Put constant coordinate */

        if (gpl->axis == PLOT_YZ)
            x = gpl->pos;
        else if (gpl->axis == PLOT_XZ)
            y = gpl->pos;
        else
            z = gpl->pos;


        /* Allocate buffer for mat or cell indeces of a single row to use for boundary checking */
        
        long *prev_row = (long*)calloc(gpl->pixx, sizeof(long));
        memset(prev_row, -1, sizeof(long) * gpl->pixx);

        /* Cell boundary plotting requires tracking lattice boundaries as well */
        long *prev_row2 = (long*)calloc(gpl->pixx, sizeof(long));
        memset(prev_row2, -1, sizeof(long) * gpl->pixx);
        
        long last_mat = -1, last_cell = -1, last_lat_elem = -1;

        for (long j = 0; j < gpl->pixy; j++)
        {
            /* Advance next row */

            if (gpl->axis == PLOT_YZ)
                z = max2 - (j + 0.5) * dj;
            else if (gpl->axis == PLOT_XZ)
                z = max2 - (j + 0.5) * dj;
            else
                y = max2 - (j + 0.5) * dj;

            for (long i = 0; i < gpl->pixx; i++)
            {
                /* Advance next pixel */

                if (gpl->axis == PLOT_YZ)
                    y = max1 - (i + 0.5) * di;
                else if (gpl->axis == PLOT_XZ)
                    x = max1 - (i + 0.5) * di;
                else
                    x = max1 - (i + 0.5) * di;

                /* Get material colour at position */

                int r = 0, g = 0, b = 0, err = CELL_ERR_OK;
                long mat_idx = getMaterialAtPosition(x, y, z, &err);

                if (err == CELL_ERR_OVERLAP)
                {
                    r = 255;
                    g = 0;
                    b = 0;
                }
                else if (err == CELL_ERR_UNDEFINED)
                {
                    r = 0;
                    g = 255;
                    b = 0;
                }
                else
                {
                    if (mat_idx >= 0)
                    {
                        Material *mat = &DATA.mats[mat_idx];
                        r = (int)mat->rgb[0];
                        g = (int)mat->rgb[1];
                        b = (int)mat->rgb[2];
                    }
                }

                /* Check if boundary crossed and draw black pixel */

                switch (gpl->bounds)
                {
                    case BOUNDS_MATERIAL:
                    {
                        if (j == 0)
                        {
                            prev_row[i] = mat_idx;
                            if (i == 0)
                                last_mat = mat_idx;
                        }
                        else if (prev_row[i] != mat_idx || last_mat != mat_idx)
                        {
                            prev_row[i] = mat_idx;
                            last_mat = mat_idx;
                            r = 0;
                            g = 0;
                            b = 0;
                        }
                        break;
                    }
                    case BOUNDS_CELL:
                    {
                        cellSearchRes res = cellSearch(x, y, z, 0.0, 0.0, 0.0);
                        long cell_idx = res.cell_idx;
                        long lat_elem_idx = res.lattice_eidx;

                        if (j == 0)
                        {
                            prev_row[i] = cell_idx;
                            prev_row2[i] = lat_elem_idx;
                            if (i == 0)
                            {
                                last_cell = cell_idx;
                                last_lat_elem = lat_elem_idx;
                            }
                        }
                        else if ((prev_row[i] != cell_idx || last_cell != cell_idx) ||
                                 (prev_row2[i] >= 0 && prev_row2[i] != lat_elem_idx) ||
                                 (last_lat_elem >= 0 && last_lat_elem != lat_elem_idx))
                        {
                            r = 0;
                            g = 0;
                            b = 0;
                        }

                        last_lat_elem = lat_elem_idx;
                        last_cell = cell_idx;
                        prev_row[i] = cell_idx;
                        prev_row2[i] = lat_elem_idx;
                        break;
                    }
                }

                fprintf(imgf, "%d %d %d ", r, g, b);
            }
            fprintf(imgf, "\n");
        }

        fprintf(stdout, "  %5.1lf%%\n", (100 * (g + 1.0) / (double)DATA.n_gpls));
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}