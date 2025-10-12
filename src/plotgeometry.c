#include "header.h"

static const int DRAW_CELL_BOUNDARIES = 1;

void getPixelColor(double x, double y, double z, int *r, int *g, int *b, long *cell_idx_out)
{
    int err = CELL_ERR_OK;
    long cell_idx = cellSearch(x, y, z, &err, NULL, NULL, NULL);
    if (cell_idx_out)
        *cell_idx_out = cell_idx;

    *r = 0;
    *g = 0;
    *b = 0;

    if (err == CELL_ERR_OVERLAP)
    {
        /* Overlaps highlighted in red */
        *r = 255;
        return;
    }
    if (err == CELL_ERR_UNDEFINED)
    {
        /* Undefined regions highlighted in green */
        *g = 255;
        return;
    }

    if (cell_idx >= 0 && err == CELL_ERR_OK)
    {
        Cell *cell = &DATA.cells[cell_idx];
        if (cell->mat_idx >= 0)
        {
            Material *mat = &DATA.mats[cell->mat_idx];
            *r = (int)mat->rgb[0];
            *g = (int)mat->rgb[1];
            *b = (int)mat->rgb[2];
        }
    }
}

int writeXYSlice(const char *filename, size_t width, size_t height, double z_slice)
{
    double dimx = DATA.x_max - DATA.x_min;
    double dimy = DATA.y_max - DATA.y_min;

    if (dimx < 0.0 || dimy < 0.0)
    {
        fprintf(stderr, "[ERROR] Invalid XY bounds: dimx=%g dimy=%g.\n", dimx, dimy);
        return 1;
    }

    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("[ERROR] Failed to open XY slice file");
        return 1;
    }

    fprintf(fp, "P3\n%zu %zu\n255\n", width, height);

    double dx = (dimx > 0.0) ? dimx / (double)width : 0.0;
    double dy = (dimy > 0.0) ? dimy / (double)height : 0.0;

    long *prev_row = NULL;
    unsigned char *prev_row_valid = NULL;
    const long sentinel = LONG_MIN;

    if (DRAW_CELL_BOUNDARIES && width > 0)
    {
        prev_row = (long*)malloc(sizeof(long) * width);
        prev_row_valid = (unsigned char*)malloc(sizeof(unsigned char) * width);
        if (!prev_row || !prev_row_valid)
        {
            fprintf(stderr, "[WARN] Failed to allocate memory for boundary plotting. Disabling boundaries.\n");
            free(prev_row);
            free(prev_row_valid);
            prev_row = NULL;
            prev_row_valid = NULL;
        }
        else
        {
            for (size_t i = 0; i < width; ++i)
            {
                prev_row[i] = sentinel;
                prev_row_valid[i] = 0;
            }
        }
    }

    for (size_t j = 0; j < height; ++j)
    {
        double y = (dimy > 0.0) ? DATA.y_max - (j + 0.5) * dy : DATA.y_min;
        long prev_cell = sentinel;
        int prev_cell_valid = 0;

        for (size_t i = 0; i < width; ++i)
        {
            double x = (dimx > 0.0) ? DATA.x_min + (i + 0.5) * dx : DATA.x_min;

            int r, g, b;
            long cell_idx = -1;
            getPixelColor(x, y, z_slice, &r, &g, &b, &cell_idx);

            long normalized_idx = (cell_idx >= 0) ? cell_idx : sentinel;
            int boundary = 0;

            if (DRAW_CELL_BOUNDARIES && prev_row && prev_row_valid)
            {
                if (prev_cell_valid && normalized_idx != prev_cell)
                    boundary = 1;
                if (prev_row_valid[i] && normalized_idx != prev_row[i])
                    boundary = 1;
            }

            if (boundary)
                r = g = b = 0;

            fprintf(fp, "%d %d %d ", r, g, b);

            prev_cell = normalized_idx;
            prev_cell_valid = (cell_idx >= 0);
            if (prev_row && prev_row_valid)
            {
                prev_row[i] = normalized_idx;
                prev_row_valid[i] = (cell_idx >= 0);
            }
        }
        fputc('\n', fp);
    }

    free(prev_row);
    free(prev_row_valid);

    fclose(fp);
    fprintf(stdout, "Wrote XY slice image to %s (%zux%zu).\n", filename, width, height);
    return 0;
}

int writeXZSlice(const char *filename, size_t width, size_t height, double y_slice)
{
    double dimx = DATA.x_max - DATA.x_min;
    double dimz = DATA.z_max - DATA.z_min;

    if (dimx < 0.0 || dimz < 0.0)
    {
        fprintf(stderr, "[ERROR] Invalid XZ bounds: dimx=%g dimz=%g.\n", dimx, dimz);
        return 1;
    }

    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("[ERROR] Failed to open XZ slice file");
        return 1;
    }

    fprintf(fp, "P3\n%zu %zu\n255\n", width, height);

    double dx = (dimx > 0.0) ? dimx / (double)width : 0.0;
    double dz = (dimz > 0.0) ? dimz / (double)height : 0.0;

    long *prev_row = NULL;
    unsigned char *prev_row_valid = NULL;
    const long sentinel = LONG_MIN;

    if (DRAW_CELL_BOUNDARIES && width > 0)
    {
        prev_row = (long*)malloc(sizeof(long) * width);
        prev_row_valid = (unsigned char*)malloc(sizeof(unsigned char) * width);
        if (!prev_row || !prev_row_valid)
        {
            fprintf(stderr, "[WARN] Failed to allocate memory for boundary plotting. Disabling boundaries.\n");
            free(prev_row);
            free(prev_row_valid);
            prev_row = NULL;
            prev_row_valid = NULL;
        }
        else
        {
            for (size_t i = 0; i < width; ++i)
            {
                prev_row[i] = sentinel;
                prev_row_valid[i] = 0;
            }
        }
    }

    for (size_t j = 0; j < height; ++j)
    {
        double z = (dimz > 0.0) ? DATA.z_max - (j + 0.5) * dz : DATA.z_min;
        long prev_cell = sentinel;
        int prev_cell_valid = 0;

        for (size_t i = 0; i < width; ++i)
        {
            double x = (dimx > 0.0) ? DATA.x_min + (i + 0.5) * dx : DATA.x_min;

            int r, g, b;
            long cell_idx = -1;
            getPixelColor(x, y_slice, z, &r, &g, &b, &cell_idx);

            long normalized_idx = (cell_idx >= 0) ? cell_idx : sentinel;
            int boundary = 0;

            if (DRAW_CELL_BOUNDARIES && prev_row && prev_row_valid)
            {
                if (prev_cell_valid && normalized_idx != prev_cell)
                    boundary = 1;
                if (prev_row_valid[i] && normalized_idx != prev_row[i])
                    boundary = 1;
            }

            if (boundary)
                r = g = b = 0;

            fprintf(fp, "%d %d %d ", r, g, b);

            prev_cell = normalized_idx;
            prev_cell_valid = (cell_idx >= 0);
            if (prev_row && prev_row_valid)
            {
                prev_row[i] = normalized_idx;
                prev_row_valid[i] = (cell_idx >= 0);
            }
        }
        fputc('\n', fp);
    }

    free(prev_row);
    free(prev_row_valid);

    fclose(fp);
    fprintf(stdout, "Wrote XZ slice image to %s (%zux%zu).\n", filename, width, height);
    return 0;
}

long plotGeometry() {
    fprintf(stdout, "\nPlotting geometry...\n");

    const size_t pixx = 1000;
    const size_t pixy = 1000;

    int status = 0;
    status |= writeXYSlice("geometry_xy.ppm", pixx, pixy, 0.0);
    status |= writeXZSlice("geometry_xz.ppm", pixx, pixy, 0.0);

    if (status == 0)
        fprintf(stdout, "DONE.\n");
    else
        fprintf(stderr, "[ERROR] Geometry plotting encountered issues.\n");

    return status;
}
