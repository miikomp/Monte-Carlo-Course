#include "header.h"

static void getPixelColor(double x, double y, double z, int *r, int *g, int *b)
{
    int err = CELL_ERR_OK;
    int cell_idx = cellSearch(x, y, z, &err);

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

static int writeXYSlice(const char *filename, size_t width, size_t height, double z_slice)
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

    for (size_t j = 0; j < height; ++j)
    {
        double y = (dimy > 0.0) ? DATA.y_max - (j + 0.5) * dy : DATA.y_min;
        for (size_t i = 0; i < width; ++i)
        {
            double x = (dimx > 0.0) ? DATA.x_min + (i + 0.5) * dx : DATA.x_min;

            int r, g, b;
            getPixelColor(x, y, z_slice, &r, &g, &b);
            fprintf(fp, "%d %d %d ", r, g, b);
        }
        fputc('\n', fp);
    }

    fclose(fp);
    fprintf(stdout, "Wrote XY slice image to %s (%zux%zu).\n", filename, width, height);
    return 0;
}

static int writeXZSlice(const char *filename, size_t width, size_t height, double y_slice)
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

    for (size_t j = 0; j < height; ++j)
    {
        double z = (dimz > 0.0) ? DATA.z_max - (j + 0.5) * dz : DATA.z_min;
        for (size_t i = 0; i < width; ++i)
        {
            double x = (dimx > 0.0) ? DATA.x_min + (i + 0.5) * dx : DATA.x_min;

            int r, g, b;
            getPixelColor(x, y_slice, z, &r, &g, &b);
            fprintf(fp, "%d %d %d ", r, g, b);
        }
        fputc('\n', fp);
    }

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
