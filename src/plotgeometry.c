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

    if (GLOB.trackplotmode && DATA.tracks != NULL)
        fprintf(stdout, "\nPlotting tracks...\n");
    else
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

        /* If doing tracks, draw them onto the image now */
        
        if (GLOB.trackplotmode && DATA.tracks != NULL && DATA.track_counts != NULL)
        {
            const size_t points_per_track = (size_t)(MAX_COLLISION_BINS + 1u);
            const int track_colour = gdTrueColor(255, 255, 255); // White
            const int start_colour = gdTrueColor(0, 255, 0); // Green
            const int end_colour   = gdTrueColor(255, 0, 0); // Red
            const int mark_size  = (pix > 1500) ? 12 : 4;

            if (dimx > 0.0 && dimy > 0.0 && di > 0.0 && dj > 0.0)
            {
                for (size_t t = 0; t < (size_t)GLOB.n_tracks; t++)
                {
                    size_t count = DATA.track_counts[t];
                    if (count < 2)
                        continue;

                    double *points = DATA.tracks + t * points_per_track * 3u;

                    /* --- Draw track segments as lines --- */

                    for (size_t p = 0; p < count - 1; p++)
                    {
                        double *p1 = points + p * 3u;
                        double *p2 = p1 + 3u;

                        double coord1_start = 0.0, coord2_start = 0.0;
                        double coord1_end = 0.0, coord2_end = 0.0;

                        switch (gpl->axis)
                        {
                            case PLOT_YZ:
                            {
                                coord1_start = p1[1];
                                coord2_start = p1[2];
                                coord1_end   = p2[1];
                                coord2_end   = p2[2];
                                break;
                            }
                            case PLOT_XZ:
                            {
                                coord1_start = p1[0];
                                coord2_start = p1[2];
                                coord1_end   = p2[0];
                                coord2_end   = p2[2];
                                break;
                            }
                            case PLOT_XY:
                            default:
                            {
                                coord1_start = p1[0];
                                coord2_start = p1[1];
                                coord1_end   = p2[0];
                                coord2_end   = p2[1];
                                break;
                            }
                        }

                        double fx1 = (max1 - coord1_start) / di - 0.5;
                        double fy1 = (max2 - coord2_start) / dj - 0.5;
                        double fx2 = (max1 - coord1_end) / di - 0.5;
                        double fy2 = (max2 - coord2_end) / dj - 0.5;

                        if (!isfinite(fx1) || !isfinite(fy1) ||
                            !isfinite(fx2) || !isfinite(fy2))
                            continue;

                        int x1 = (int)lround(fx1);
                        int y1 = (int)lround(fy1);
                        int x2 = (int)lround(fx2);
                        int y2 = (int)lround(fy2);

                        if (x1 < 0) 
                            x1 = 0;
                        if (x1 >= gpl->pixx) 
                            x1 = gpl->pixx - 1;
                        if (x2 < 0) 
                            x2 = 0;
                        if (x2 >= gpl->pixx) 
                            x2 = gpl->pixx - 1;

                        if (y1 < 0) 
                            y1 = 0;
                        if (y1 >= gpl->pixy) 
                            y1 = gpl->pixy - 1;
                        if (y2 < 0) 
                            y2 = 0;
                        if (y2 >= gpl->pixy) 
                            y2 = gpl->pixy - 1;

                        gdImageLine(img, x1, y1, x2, y2, track_colour);
                    }

                    /* --- Draw points at the start and end of the tracks --- */

                    double *p_start = points;
                    double *p_end   = points + (count - 1u) * 3u;

                    double start_c1 = 0.0, start_c2 = 0.0;
                    double end_c1   = 0.0, end_c2   = 0.0;

                    switch (gpl->axis)
                    {
                        case PLOT_YZ:
                        {
                            start_c1 = p_start[1]; 
                            start_c2 = p_start[2];
                            end_c1   = p_end[1];   
                            end_c2   = p_end[2];
                            break;
                        }
                        case PLOT_XZ:
                        {
                            start_c1 = p_start[0]; 
                            start_c2 = p_start[2];
                            end_c1   = p_end[0];   
                            end_c2   = p_end[2];
                            break;
                        }
                        case PLOT_XY:
                        default:
                        {
                            start_c1 = p_start[0]; 
                            start_c2 = p_start[1];
                            end_c1   = p_end[0];   
                            end_c2   = p_end[1];
                            break;
                        }
                    }

                    double fx = (max1 - start_c1) / di - 0.5;
                    double fy = (max2 - start_c2) / dj - 0.5;

                    if (isfinite(fx) && isfinite(fy))
                    {
                        int sx = (int)lround(fx);
                        int sy = (int)lround(fy);
                        if (sx < 0) 
                            sx = 0;
                        else if (sx >= gpl->pixx) 
                            sx = gpl->pixx - 1;
                        if (sy < 0) 
                            sy = 0;
                        else if (sy >= gpl->pixy) 
                            sy = gpl->pixy - 1;

                        gdImageFilledEllipse(img, sx, sy, mark_size, mark_size, start_colour);
                    }

                    fx = (max1 - end_c1) / di - 0.5;
                    fy = (max2 - end_c2) / dj - 0.5;

                    if (isfinite(fx) && isfinite(fy))
                    {
                        int ex = (int)lround(fx);
                        int ey = (int)lround(fy);
                        if (ex < 0) 
                            ex = 0;
                        else if (ex >= gpl->pixx) 
                            ex = gpl->pixx - 1;
                        if (ey < 0) 
                            ey = 0;
                        else if (ey >= gpl->pixy) 
                            ey = gpl->pixy - 1;
                        
                        gdImageFilledEllipse(img, ex, ey, mark_size, mark_size, end_colour);
                    }

                }
            }
        }

        /* Save image to file*/

        char filename[MAX_STR_LEN];
        if (GLOB.trackplotmode && DATA.tracks != NULL)
            snprintf(filename, sizeof(filename), "%s_tracks%zu.png", GLOB.inputfname, g + 1);
        else
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

    if (GLOB.trackplotmode && DATA.tracks != NULL)
    {
        free(DATA.tracks);
        free(DATA.track_counts);
        exit(EXIT_SUCCESS);
    }
    else
        return EXIT_SUCCESS;
}
