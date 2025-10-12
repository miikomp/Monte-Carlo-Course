#include "header.h"

double distanceToNearestBoundary(double x, double y, double z, double u, double v, double w) {

    /* Find current cell */

    int err;
    double lx, ly, lz;
    long cell_idx = cellSearch(x, y, z, &err, &lx, &ly, &lz);
    Cell *cell;

    if (cell_idx >= 0 && err == CELL_ERR_OK)
        cell = &DATA.cells[cell_idx];
    else
        return INFINITY;

    Universe *uni = &DATA.unis[cell->uni_idx];

    /* Check if cell is inside a lattice element filling universe */

    long parent_lat_idx = uni->parent_lat_idx;
    double dmax = INFINITY;
    if (parent_lat_idx >= 0 && parent_lat_idx < (long)DATA.n_lats)
    {
        Lattice *lat = &DATA.lats[parent_lat_idx];
        long bb_idx = lat->bb_surf_idx;
        if (bb_idx >= 0)
        {
            Surface *S = &DATA.surfs[bb_idx];
            dmax = surfaceDistance(S->type, S->params, S->n_params, lx, ly, lz, u, v, w);
        }
        else
        {
            fprintf(stderr, "[ERROR] Lattice \"%s\" without element boundaries.\n", lat->name);
            exit(EXIT_FAILURE);
        }
    }

    /* Compute distance to system bounding surface */

    Surface *S = &DATA.surfs[DATA.outside_surf_idx];
    double tx = x;
    double ty = y;
    double tz = z;
    double tu = u;
    double tv = v;
    double tw = w;

    if (S->t_idx >= 0)
        applyTransformation(&DATA.transforms[S->t_idx], &tx, &ty, &tz, &tu, &tv, &tw);

    double d0 = surfaceDistance(S->type, S->params, S->n_params, tx, ty, tz, tu, tv, tw);

    if (dmax > d0)
        dmax = d0;

    /* Find nearest surface in cell */

    double d = INFINITY;

    for (size_t s = 0; s < cell->n_surfs; s++)
    {
        size_t idx = cell->surf_idxs[s];
        Surface *S = &DATA.surfs[idx];

        double tx = lx;
        double ty = ly;
        double tz = lz;
        double tu = u;
        double tv = v;
        double tw = w;

        if (S->t_idx >= 0)
            applyTransformation(&DATA.transforms[S->t_idx], &tx, &ty, &tz, &tu, &tv, &tw);

        double d0 = surfaceDistance(S->type, S->params, S->n_params, tx, ty, tz, tu, tv, tw);
        if (d0 < d)
            d = d0;
    }

    if (d > dmax)
        d = dmax;

    return d;
}

/* Cell search should return lattice index!*/