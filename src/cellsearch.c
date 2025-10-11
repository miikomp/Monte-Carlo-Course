#include "header.h"

static void roundHexAxial(double q, double r, long *q_round, long *r_round)
{
    double s = -q - r;

    long rq = (long)lround(q);
    long rr = (long)lround(r);
    long rs = (long)lround(s);

    double dq = fabs((double)rq - q);
    double dr = fabs((double)rr - r);
    double ds = fabs((double)rs - s);

    if (dq > dr && dq > ds)
        rq = -rr - rs;
    else if (dr > ds)
        rr = -rq - rs;
    else
        rs = -rq - rr;

    if (q_round)
        *q_round = rq;
    if (r_round)
        *r_round = rr;
}

long locateCellInUniverse(size_t uni_idx, double x, double y, double z, int *err, double *min_abs_out);

long cellSearch(double x, double y, double z, int *err) {

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

        /* If lattice, determine which sub-universe to go to next */

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

                    long ix = floor(rel_x / lat->dx + 0.5);
                    long iy = floor(rel_y / lat->dy + 0.5);

                    x -= (lat->x0 + ix * lat->dx);
                    y -= (lat->y0 + iy * lat->dy);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_HEXX_INFINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double pitch = lat->dx;
                    double hex_size = pitch / SQRT3;

                    double q = (SQRT3 / 3.0 * rel_x - rel_y / 3.0) / hex_size;
                    double r = (2.0 / 3.0 * rel_y) / hex_size;

                    long iq = 0;
                    long ir = 0;
                    roundHexAxial(q, r, &iq, &ir);

                    double center_x = hex_size * SQRT3 * ((double)iq + 0.5 * (double)ir);
                    double center_y = hex_size * 1.5 * (double)ir;

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_HEXY_INFINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double pitch = lat->dy;
                    double hex_size = pitch / SQRT3;

                    double q = (2.0 / 3.0 * rel_x) / hex_size;
                    double r = (-rel_x / 3.0 + SQRT3 / 3.0 * rel_y) / hex_size;

                    long iq = 0;
                    long ir = 0;
                    roundHexAxial(q, r, &iq, &ir);

                    double center_x = hex_size * 1.5 * (double)iq;
                    double center_y = pitch * ((double)ir + 0.5 * (double)iq);

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_SQUARE_FINITE:
                {
                    double ix_f = (x - lat->x0) / lat->dx + 0.5 * ((double)lat->nx - 1.0);
                    double iy_f = (y - lat->y0) / lat->dy + 0.5 * ((double)lat->ny - 1.0);

                    long ix = (long)floor(ix_f + 0.5);
                    long iy = (long)floor(iy_f + 0.5);

                    if (ix < 0 || ix >= lat->nx || iy < 0 || iy >= lat->ny)
                    {
                        if (err)
                            *err = CELL_ERR_UNDEFINED;
                        return -1;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double cell_x0 = lat->x0 + ((double)ix - 0.5 * ((double)lat->nx - 1.0)) * lat->dx;
                    double cell_y0 = lat->y0 + ((double)iy - 0.5 * ((double)lat->ny - 1.0)) * lat->dy;

                    x -= cell_x0;
                    y -= cell_y0;
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];

                    break;
                }
                case LAT_HEXX_FINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double pitch = lat->dx;
                    double hex_size = pitch / SQRT3;

                    double q = (SQRT3 / 3.0 * rel_x - rel_y / 3.0) / hex_size;
                    double r = (2.0 / 3.0 * rel_y) / hex_size;

                    long iq = 0;
                    long ir = 0;
                    roundHexAxial(q, r, &iq, &ir);

                    double q_offset = 0.5 * ((double)lat->nx - 1.0);
                    double r_offset = 0.5 * ((double)lat->ny - 1.0);

                    long ix = (long)floor(((double)iq + q_offset) + 0.5);
                    long iy = (long)floor(((double)ir + r_offset) + 0.5);

                    if (ix < 0 || ix >= lat->nx || iy < 0 || iy >= lat->ny)
                    {
                        if (err)
                            *err = CELL_ERR_UNDEFINED;
                        return -1;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double center_x = hex_size * SQRT3 * ((double)iq + 0.5 * (double)ir);
                    double center_y = hex_size * 1.5 * (double)ir;

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];
                    break;
                }
                case LAT_HEXY_FINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double pitch = lat->dy;
                    double hex_size = pitch / SQRT3;

                    double q = (2.0 / 3.0 * rel_x) / hex_size;
                    double r = (-rel_x / 3.0 + SQRT3 / 3.0 * rel_y) / hex_size;

                    long iq = 0;
                    long ir = 0;
                    roundHexAxial(q, r, &iq, &ir);

                    double q_offset = 0.5 * ((double)lat->nx - 1.0);
                    double r_offset = 0.5 * ((double)lat->ny - 1.0);

                    long ix = (long)floor(((double)iq + q_offset) + 0.5);
                    long iy = (long)floor(((double)ir + r_offset) + 0.5);

                    if (ix < 0 || ix >= lat->nx || iy < 0 || iy >= lat->ny)
                    {
                        if (err)
                            *err = CELL_ERR_UNDEFINED;
                        return -1;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double center_x = hex_size * 1.5 * (double)iq;
                    double center_y = pitch * ((double)ir + 0.5 * (double)iq);

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];
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
 * @return long index of the cell containing the point, or -1 if not found
 */
long locateCellInUniverse(size_t uni_idx, double x, double y, double z,
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

    if (uni->t_idx >= 0)
        applyTransformation(&DATA.transforms[uni->t_idx], &x, &y, &z, NULL, NULL, NULL);

    long best_cell = -1;
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
            double tx = x;
            double ty = y;
            double tz = z;

            if (surf->t_idx >= 0)
                applyTransformation(&DATA.transforms[surf->t_idx], &tx, &ty, &tz, NULL, NULL, NULL);
            
            double res = surfaceTest(surf->type, surf->params, surf->n_params, tx, ty, tz);
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
