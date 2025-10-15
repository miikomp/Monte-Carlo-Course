#include "header.h"
void roundHexAxial(double q, double r, long *q_round, long *r_round);

long locateCellInUniverse(size_t uni_idx, double x, double y, double z, int *err, double *min_abs_out);

cellSearchRes cellSearch(double x, double y, double z, double u, double v, double w) {

    cellSearchRes res = {0};
    res.err = CELL_ERR_OK;
    res.lattice_eidx = -1;

    /* Start at the root universe 0 */

    size_t universe_idx = 0;

    while (true)
    {
        if (universe_idx >= DATA.n_unis)
        {
            
            res.err = CELL_ERR_UNDEFINED;
            return res;
        }

        Universe *current_uni = &DATA.unis[universe_idx];

        /* If lattice, determine which sub-universe to go to next */

        if (current_uni->type == UNI_LATTICE)
        {
            if (current_uni->lat_idx < 0 || (size_t)current_uni->lat_idx >= DATA.n_lats)
            {
                
                res.err = CELL_ERR_UNDEFINED;
                return res;
            }

            Lattice *lat = &DATA.lats[current_uni->lat_idx];

            if (!lat->uni_idxs || lat->n_unis == 0)
            {
                res.err = CELL_ERR_UNDEFINED;
                return res;
            }

            if (lat->pitch == 0.0)
            {
                res.err = CELL_ERR_UNDEFINED;
                return res;
            }

            switch (lat->type)
            {
                case LAT_SQUARE_INFINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    long ix = floor(rel_x / lat->pitch + 0.5);
                    long iy = floor(rel_y / lat->pitch + 0.5);

                    x -= (lat->x0 + ix * lat->pitch);
                    y -= (lat->y0 + iy * lat->pitch);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_HEXX_INFINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double hex_size = lat->pitch / SQRT3;

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

                    double hex_size = lat->pitch / SQRT3;

                    double q = (2.0 / 3.0 * rel_x) / hex_size;
                    double r = (-rel_x / 3.0 + SQRT3 / 3.0 * rel_y) / hex_size;

                    long iq = 0;
                    long ir = 0;
                    roundHexAxial(q, r, &iq, &ir);

                    double center_x = hex_size * 1.5 * (double)iq;
                    double center_y = lat->pitch * ((double)ir + 0.5 * (double)iq);

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_TRI_INFINITE:
                {
                    const double s = lat->pitch;

                    const double inv_s = 1.0 / s;
                    const double inv_s_sqrt3 = inv_s / SQRT3;
                    const double two_inv_s_sqrt3 = 2.0 * inv_s / SQRT3;
                    const double eps = 1e-12;

                    double origin_x = lat->x0 - 0.5 * s;
                    double origin_y = lat->y0 - (SQRT3 / 6.0) * s;

                    double rel_x = x - origin_x;
                    double rel_y = y - origin_y;

                    double i_coord = rel_x * inv_s - rel_y * inv_s_sqrt3;
                    double j_coord = rel_y * two_inv_s_sqrt3;

                    long i_base = (long)floor(i_coord);
                    long j_base = (long)floor(j_coord);

                    double fi = i_coord - (double)i_base;
                    double fj = j_coord - (double)j_base;

                    double vx = origin_x + (double)i_base * s + (double)j_base * (0.5 * s);
                    double vy = origin_y + (double)j_base * (0.5 * SQRT3 * s);

                    double cx, cy;
                    bool invert = false;

                    if (fi + fj <= 1.0 + eps)
                    {
                        cx = vx + 0.5 * s;
                        cy = vy + (SQRT3 / 6.0) * s;
                    }
                    else
                    {
                        cx = vx + s;
                        cy = vy + (SQRT3 / 3.0) * s;
                        invert = true;
                    }

                    x -= cx;
                    y -= cy;
                    z -= lat->z0;

                    if (invert)
                    {
                        x = -x;
                        y = -y;
                        u = -u;
                        v = -v;
                    }

                    universe_idx = (size_t)lat->uni_idxs[0];
                    break;
                }
                case LAT_SQUARE_FINITE:
                {
                    double ix_f = (x - lat->x0) / lat->pitch + 0.5 * ((double)lat->nx - 1.0);
                    double iy_f = (y - lat->y0) / lat->pitch + 0.5 * ((double)lat->ny - 1.0);

                    long ix = (long)floor(ix_f + 0.5);
                    long iy = (long)floor(iy_f + 0.5);

                    if (ix < 0 || ix >= lat->nx || iy < 0 || iy >= lat->ny)
                    {
                        res.err = CELL_ERR_UNDEFINED;
                        return res;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double cell_x0 = lat->x0 + ((double)ix - 0.5 * ((double)lat->nx - 1.0)) * lat->pitch;
                    double cell_y0 = lat->y0 + ((double)iy - 0.5 * ((double)lat->ny - 1.0)) * lat->pitch;

                    x -= cell_x0;
                    y -= cell_y0;
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];
                    res.lattice_eidx = idx;

                    break;
                }
                case LAT_HEXX_FINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double hex_size = lat->pitch / SQRT3;

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
                        res.err = CELL_ERR_UNDEFINED;
                        return res;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double center_x = hex_size * SQRT3 * ((double)iq + 0.5 * (double)ir);
                    double center_y = hex_size * 1.5 * (double)ir;

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];
                    res.lattice_eidx = idx;

                    break;
                }
                case LAT_HEXY_FINITE:
                {
                    double rel_x = x - lat->x0;
                    double rel_y = y - lat->y0;

                    double hex_size = lat->pitch / SQRT3;

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
                        res.err = CELL_ERR_UNDEFINED;
                        return res;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    double center_x = hex_size * 1.5 * (double)iq;
                    double center_y = lat->pitch * ((double)ir + 0.5 * (double)iq);

                    x -= (lat->x0 + center_x);
                    y -= (lat->y0 + center_y);
                    z -= lat->z0;

                    universe_idx = (size_t)lat->uni_idxs[idx];
                    res.lattice_eidx = idx;

                    break;
                }
                case LAT_TRI_FINITE:
                {
                    const double s = lat->pitch;
                    if (s <= 0.0)
                    {
                        res.err = CELL_ERR_UNDEFINED;
                        return res;
                    }

                    const double inv_s = 1.0 / s;
                    const double inv_s_sqrt3 = inv_s / SQRT3;
                    const double two_inv_s_sqrt3 = 2.0 * inv_s / SQRT3;
                    const double eps = 1e-12;

                    double origin_x = lat->x0 - 0.5 * s;
                    double origin_y = lat->y0 - (SQRT3 / 6.0) * s;

                    double rel_x = x - origin_x;
                    double rel_y = y - origin_y;

                    double i_coord = rel_x * inv_s - rel_y * inv_s_sqrt3;
                    double j_coord = rel_y * two_inv_s_sqrt3;

                    long i_base = (long)floor(i_coord);
                    long j_base = (long)floor(j_coord);

                    double fi = i_coord - (double)i_base;
                    double fj = j_coord - (double)j_base;

                    double vx = origin_x + (double)i_base * s + (double)j_base * (0.5 * s);
                    double vy = origin_y + (double)j_base * (0.5 * SQRT3 * s);

                    double cx, cy;
                    bool invert = false;

                    if (fi + fj <= 1.0 + eps)
                    {
                        cx = vx + 0.5 * s;
                        cy = vy + (SQRT3 / 6.0) * s;
                    }
                    else
                    {
                        cx = vx + s;
                        cy = vy + (SQRT3 / 3.0) * s;
                        invert = true;
                    }

                    long col = 2 * i_base + (invert ? 1 : 0);
                    long row = j_base;

                    double col_f = (double)col + 0.5 * ((double)lat->nx - 1.0);
                    double row_f = (double)row + 0.5 * ((double)lat->ny - 1.0);

                    long ix = (long)floor(col_f + 0.5);
                    long iy = (long)floor(row_f + 0.5);

                    if (ix < 0 || ix >= lat->nx || iy < 0 || iy >= lat->ny)
                    {
                        res.err = CELL_ERR_UNDEFINED;
                        return res;
                    }

                    size_t idx = (size_t)ix + (size_t)iy * (size_t)lat->nx;

                    x -= cx;
                    y -= cy;
                    z -= lat->z0;

                    if (invert)
                    {
                        x = -x;
                        y = -y;
                        u = -u;
                        v = -v;
                    }

                    universe_idx = (size_t)lat->uni_idxs[idx];
                    res.lattice_eidx = (long)idx;

                    break;
                }
                default:
                {
                    res.err = CELL_ERR_UNDEFINED;
                    return res;
                }
            }

            /* Forward loop to move into the lattice universe */

            continue;
        }

        int err = CELL_ERR_OK;
        double min_abs = INFINITY;
        int cell_idx = locateCellInUniverse(universe_idx, x, y, z, &err, &min_abs);

        if (err == CELL_ERR_OVERLAP && res.err != CELL_ERR_OVERLAP)
            res.err = CELL_ERR_OVERLAP;

        if (cell_idx < 0)
        {
            if (res.err == CELL_ERR_OK)
                res.err = err;
            return res;
        }

        Cell *cell = &DATA.cells[cell_idx];

        /* If material filled, put local coordinates, universe index and return cell index */

        if (!cell->unifilled)
        {
            res.lx = x;
            res.ly = y;
            res.lz = z;
            res.lu = u;
            res.lv = v;
            res.lw = w;
            res.cell_idx = cell_idx;
            return res;
        }
        /* if filled with universe get the filling universe */

        if (cell->filluni_idx < 0 || (size_t)cell->filluni_idx >= DATA.n_unis)
        {
            res.err = CELL_ERR_UNDEFINED;
            return res;
        }

        /* Next universe (material or further fills) */
        universe_idx = (size_t)cell->filluni_idx;
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

void roundHexAxial(double q, double r, long *q_round, long *r_round)
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
