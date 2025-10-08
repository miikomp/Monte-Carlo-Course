#include "header.h"

int locateCellInUniverse(size_t uni_idx, double x, double y, double z, int *err, double *min_abs_out);

int cellSearch(double x, double y, double z, int *err) {

    if (err)
        *err = CELL_ERR_OK;

    /* Start at the root universe 0 */

    size_t universe_idx = 0;
    size_t safety_counter = 0;

    while (true)
    {
        if (universe_idx >= DATA.n_unis)
        {
            if (err)
                *err = CELL_ERR_UNDEFINED;
            return -1;
        }

        Universe *current_uni = &DATA.unis[universe_idx];

        /* If lattice, determine which universe to go to next */

        if (current_uni->type == UNI_LATTICE)
        {
            if (current_uni->lat_idx < 0 || (size_t)current_uni->lat_idx >= DATA.n_lats)
            {
                if (err)
                    *err = CELL_ERR_UNDEFINED;
                return -1;
            }

            Lattice *lat = &DATA.lats[current_uni->lat_idx];

            if (!lat->uni_idxs || lat->n_unis == 0)
            {
                if (err)
                    *err = CELL_ERR_UNDEFINED;
                return -1;
            }

            if (lat->dx == 0.0 || lat->dy == 0.0)
            {
                if (err)
                    *err = CELL_ERR_UNDEFINED;
                return -1;
            }

            switch (lat->type)
            {
                case LAT_SQUARE_INFINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;
                    long ix = (long)floor(rel_x / lat->dx + 0.5);
                    long iy = (long)floor(rel_y / lat->dy + 0.5);

                    x -= (lat->x0 + ix * lat->dx);
                    y -= (lat->y0 + iy * lat->dy);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                default:
                {
                    if (err)
                        *err = CELL_ERR_UNDEFINED;
                    return -1;
                }
            }

            /* Forward loop to move into the lattice universe */

            continue;
        }

        int local_err = CELL_ERR_OK;
        double min_abs = INFINITY;
        int cell_idx = locateCellInUniverse(universe_idx, x, y, z, &local_err, &min_abs);

        if (local_err == CELL_ERR_OVERLAP && err && *err != CELL_ERR_OVERLAP)
            *err = CELL_ERR_OVERLAP;

        if (cell_idx < 0)
        {
            if (err && *err == CELL_ERR_OK)
                *err = local_err;
            return -1;
        }

        Cell *cell = &DATA.cells[cell_idx];

        /* If material filled, return cell index */

        if (!cell->unifilled)
            return cell_idx;

        /* if filled with universe get the filling universe */

        if (cell->filluni_idx < 0 || (size_t)cell->filluni_idx >= DATA.n_unis)
        {
            if (err)
                *err = CELL_ERR_UNDEFINED;
            return -1;
        }

        /* Next universe (material or further fills) */
        universe_idx = (size_t)cell->filluni_idx;

        if (++safety_counter > DATA.n_unis * 16UL)
        {
            if (err && *err == CELL_ERR_OK)
                *err = CELL_ERR_OVERLAP;
            return cell_idx;
        }
    }
}


/**
 * @brief Locate the cell in the given universe that contains the point (x, y, z).
 * 
 * @param uni_idx index of universe in DATA.unis
 * @param x x-coordinate of the point
 * @param y y-coordinate of the point
 * @param z z-coordinate of the point
 * @param err pointer to an error code
 * @param min_abs_out pointer to store the minimum absolute distance
 * @return int index of the cell containing the point, or -1 if not found
 */
int locateCellInUniverse(size_t uni_idx, double x, double y, double z,
                                int *err, double *min_abs_out)
{
    const double EPS = 1e-12;

    if (err)
        *err = CELL_ERR_OK;
    if (min_abs_out)
        *min_abs_out = INFINITY;

    if (uni_idx >= DATA.n_unis)
    {
        if (err)
            *err = CELL_ERR_UNDEFINED;
        return -1;
    }

    Universe *uni = &DATA.unis[uni_idx];
    int best_cell = -1;
    double best_min_abs = INFINITY;

    for (size_t i = 0; i < uni->n_cells; ++i)
    {
        int cell_idx = uni->cell_idxs[i];
        Cell *cell = &DATA.cells[cell_idx];

        double min_abs = INFINITY;
        size_t s;
        for (s = 0; s < cell->n_surfs; ++s)
        {
            int side = cell->sides[s];
            Surface *surf = &DATA.surfs[cell->surf_idxs[s]];
            double res = surfaceTest(surf->type, surf->params, surf->n_params, x, y, z);
            double abs_res = fabs(res);
            if (abs_res < min_abs)
                min_abs = abs_res;

            if (side * res < -EPS)
                break;
        }

        if (s == cell->n_surfs)
        {
            if (best_cell < 0)
            {
                best_cell = cell_idx;
                best_min_abs = min_abs;
            }
            else
            {
                if (!(best_min_abs <= EPS && min_abs <= EPS))
                {
                    if (err && *err != CELL_ERR_OVERLAP)
                        *err = CELL_ERR_OVERLAP;
                }
                else if (min_abs < best_min_abs)
                {
                    best_cell = cell_idx;
                    best_min_abs = min_abs;
                }
            }
        }
    }

    if (best_cell < 0)
    {
        if (err && *err == CELL_ERR_OK)
            *err = CELL_ERR_UNDEFINED;
    }
    else if (min_abs_out)
    {
        *min_abs_out = best_min_abs;
    }

    return best_cell;
}
