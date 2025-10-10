#include "header.h"

double distanceToNearestBoundary(double x, double y, double z, double u, double v, double w) {

    /* Find current cell */

    int err;
    long cell_idx = cellSearch(x, y, z, &err);
    Cell *cell;

    if (cell_idx >= 0 && err == CELL_ERR_OK)
        cell = &DATA.cells[cell_idx];
    else
        return INFINITY;

    /* Find nearest surface in cell */

    double d = INFINITY;

    for (size_t s = 0; s < cell->n_surfs; s++)
    {
        size_t idx = cell->surf_idxs[s];
        Surface *S = &DATA.surfs[idx];

        double tx = x;
        double ty = y;
        double tz = z;
        double tu = u;
        double tv = v;
        double tw = w;

        if (S->t_idx >= 0)
            applyTransformation(&DATA.transforms[S->t_idx], &tx, &ty, &tz, &tu, &tv, &tw);

        double d0 = surfaceDistance(S->type, S->params, S->n_params, tx, ty, tz, tu, tv, tw);
        if (d0 < d)
            d = d0;
    }

    return d;
}