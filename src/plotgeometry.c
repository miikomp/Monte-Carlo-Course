#include "header.h"

enum {
    PLOT_YZ = 1,
    PLOT_XZ,
    PLOT_XY
};

enum {
    BOUNDS_NONE,
    BOUNDS_MATERIAL,
    BOUNDS_CELL
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
                min1 = isnan(gpl->ymin) ? DATA.y_min : gpl->ymin;
                max1 = isnan(gpl->ymax) ? DATA.y_max : gpl->ymax;
                min2 = isnan(gpl->zmin) ? DATA.z_min : gpl->zmin;
                max2 = isnan(gpl->zmax) ? DATA.z_max : gpl->zmax;

                break;
            }
            case PLOT_XZ:
            {
                min1 = isnan(gpl->xmin) ? DATA.x_min : gpl->xmin;
                max1 = isnan(gpl->xmax) ? DATA.x_max : gpl->xmax;
                min2 = isnan(gpl->zmin) ? DATA.z_min : gpl->zmin;
                max2 = isnan(gpl->zmax) ? DATA.z_max : gpl->zmax;

                break;
            }
            case PLOT_XY:
            {
                min1 = isnan(gpl->xmin) ? DATA.x_min : gpl->xmin;
                max1 = isnan(gpl->xmax) ? DATA.x_max : gpl->xmax;
                min2 = isnan(gpl->ymin) ? DATA.y_min : gpl->ymin;
                max2 = isnan(gpl->ymax) ? DATA.y_max : gpl->ymax;

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

        /* Create image file (.png format) */

        gdImagePtr img = gdImageCreateTrueColor(gpl->pixx, gpl->pixy);
        if (!img)
        {
            fprintf(stderr, "[ERROR] Could not create image file.\n");
            return EXIT_FAILURE;
        }

        /* Create colour palette, colour for each material + special colours for overlap, undefined and outside */
        
        int *palette = calloc(DATA.n_mats + 3, sizeof(int));
        if (!palette)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            return EXIT_FAILURE;
        }

        for (size_t m = 0; m < DATA.n_mats; m++)
        {
            Material *mat = &DATA.mats[m];
            palette[m] = gdImageColorAllocate(img, mat->rgb[0], mat->rgb[1], mat->rgb[2]);
        }

        /* Add the special colours after the material colours */
        size_t idx = DATA.n_mats;

        size_t overlap_colour_idx = idx;
        palette[idx++] = gdImageColorAllocate(img, 255, 0, 0);

        size_t undefined_colour_idx = idx;
        palette[idx++] = gdImageColorAllocate(img, 0, 255, 0);
        
        size_t outside_colour_idx = idx;
        palette[idx++] = gdImageColorAllocate(img, 0, 0, 0);

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

                int err = CELL_ERR_OK;
                long mat_idx = getMaterialAtPosition(x, y, z, &err);
                int colour;

                if (err == CELL_ERR_OVERLAP)
                    colour = palette[overlap_colour_idx];
                else if (err == CELL_ERR_UNDEFINED)
                    colour = palette[undefined_colour_idx];
                else
                {
                    if (mat_idx >= 0 && mat_idx < (long)DATA.n_mats)
                        colour = palette[mat_idx];
                    else
                        colour = palette[outside_colour_idx];
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
                            colour = palette[outside_colour_idx];
                        }
                        break;
                    }
                    case BOUNDS_CELL:
                    {
                        long lat_idx = -1;
                        long cell_idx = cellSearch(x, y, z, NULL, NULL, NULL, NULL, &lat_idx, NULL, NULL, NULL, NULL);
                        if (j == 0)
                        {
                            prev_row[i] = cell_idx;
                            prev_row2[i] = lat_idx;
                            if (i == 0)
                            {
                                last_cell = cell_idx;
                                last_lat_elem = lat_idx;
                            }
                        }
                        else if ((prev_row[i] != cell_idx || last_cell != cell_idx) ||
                                 (prev_row2[i] >= 0 && prev_row2[i] != lat_idx) ||
                                 (last_lat_elem >= 0 && last_lat_elem != lat_idx))
                        {
                            colour = palette[outside_colour_idx];
                        }

                        last_lat_elem = lat_idx;
                        last_cell = cell_idx;
                        prev_row[i] = cell_idx;
                        prev_row2[i] = lat_idx;
                        break;
                    }
                }
                gdImageSetPixel(img, i, j, colour);
            }
        }

        /* Save to image file*/

        char filename[MAX_STR_LEN];
        snprintf(filename, sizeof(filename), "%s_geom%zu.png", GLOB.inputfname, g + 1);
        FILE *fp = fopen(filename, "wb");
        if (!fp)
        {
            fprintf(stderr, "[ERROR] Could not open file \"%s\" for editing.\n", filename);
            return EXIT_FAILURE;
        }

        gdImagePng(img, fp);
        fclose(fp);
        gdImageDestroy(img);

        fprintf(stdout, "  %5.1lf%%\n", (100 * (g + 1.0) / (double)DATA.n_gpls));
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}