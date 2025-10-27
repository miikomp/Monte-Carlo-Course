#include "header.h"

void applyBoundaryConditions(double *x, double *y, double *z, double *u, double *v, double *w)
{
    if (!x || !y || !z || DATA.outside_surf_idx < 0 || DATA.outside_surf_idx >= (long)DATA.n_surf)
        return;

    Surface *S = &DATA.surfs[DATA.outside_surf_idx];

    if (surfaceTest(S->type, S->params, S->n_params, *x, *y, *z) <= 0.0)
        return;

    const double xmin = DATA.x_min;
    const double xmax = DATA.x_max;
    const double ymin = DATA.y_min;
    const double ymax = DATA.y_max;
    const double zmin = DATA.z_min;
    const double zmax = DATA.z_max;

    switch (DATA.boundary_coef)
    {
        case BC_BLACK:
            return;

        case BC_REFLECTIVE:
        {
            const SurfaceTypes stype = S->type;

            if (stype == SURF_SQR)
            {
                const double EPS = 1.0e-12;
                const double span_x = xmax - xmin;
                const double span_y = ymax - ymin;
                const double span_z = zmax - zmin;

                if (span_x > 0.0)
                {
                    const double span2 = 2.0 * span_x;
                    double offset = fmod(*x - xmin, span2);
                    if (offset < 0.0)
                        offset += span2;

                    bool reflected = false;
                    if (offset > span_x)
                    {
                        offset = span2 - offset;
                        reflected = true;
                    }

                    *x = xmin + offset;
                    if (reflected && u)
                        *u = -*u;

                    if (*x < xmin + EPS)
                        *x = xmin + EPS;
                    else if (*x > xmax - EPS)
                        *x = xmax - EPS;
                }

                if (span_y > 0.0)
                {
                    const double span2 = 2.0 * span_y;
                    double offset = fmod(*y - ymin, span2);
                    if (offset < 0.0)
                        offset += span2;

                    bool reflected = false;
                    if (offset > span_y)
                    {
                        offset = span2 - offset;
                        reflected = true;
                    }

                    *y = ymin + offset;
                    if (reflected && v)
                        *v = -*v;

                    if (*y < ymin + EPS)
                        *y = ymin + EPS;
                    else if (*y > ymax - EPS)
                        *y = ymax - EPS;
                }

                if (span_z > 0.0)
                {
                    const double span2 = 2.0 * span_z;
                    double offset = fmod(*z - zmin, span2);
                    if (offset < 0.0)
                        offset += span2;

                    bool reflected = false;
                    if (offset > span_z)
                    {
                        offset = span2 - offset;
                        reflected = true;
                    }

                    *z = zmin + offset;
                    if (reflected && w)
                        *w = -*w;

                    if (*z < zmin + EPS)
                        *z = zmin + EPS;
                    else if (*z > zmax - EPS)
                        *z = zmax - EPS;
                }
            }

            return;
        }

        case BC_PERIODIC:
        {
            const SurfaceTypes stype = S->type;

            if (stype == SURF_SQR)
            {
                const double span_x = xmax - xmin;
                const double span_y = ymax - ymin;
                const double span_z = zmax - zmin;

                if (span_x > 0.0)
                {
                    double delta = fmod(*x - xmin, span_x);
                    if (delta < 0.0)
                        delta += span_x;
                    *x = xmin + delta;
                }

                if (span_y > 0.0)
                {
                    double delta = fmod(*y - ymin, span_y);
                    if (delta < 0.0)
                        delta += span_y;
                    *y = ymin + delta;
                }

                if (span_z > 0.0)
                {
                    double delta = fmod(*z - zmin, span_z);
                    if (delta < 0.0)
                        delta += span_z;
                    *z = zmin + delta;
                }

                return;
            }

            if (stype == SURF_HEXX)
            {
                const double cx = S->params[0];
                const double cy = S->params[1];
                const double a = S->params[2];

                if (a > 0.0)
                {
                    double dx = *x - cx;
                    double dy = *y - cy;

                    for (int iter = 0; iter < 12; ++iter)
                    {
                        const double planes[6] = {
                            dx - a,
                            -dx - a,
                            dx + SQRT3 * dy - 2.0 * a,
                            -dx - SQRT3 * dy - 2.0 * a,
                            dx - SQRT3 * dy - 2.0 * a,
                            -dx + SQRT3 * dy - 2.0 * a
                        };

                        int max_idx = -1;
                        double max_val = 0.0;

                        for (int p = 0; p < 6; ++p)
                        {
                            if (planes[p] > max_val)
                            {
                                max_val = planes[p];
                                max_idx = p;
                            }
                        }

                        if (max_idx < 0)
                            break;

                        switch (max_idx)
                        {
                            case 0:
                                dx -= 2.0 * a;
                                break;
                            case 1:
                                dx += 2.0 * a;
                                break;
                            case 2:
                                dx -= a;
                                dy -= SQRT3 * a;
                                break;
                            case 3:
                                dx += a;
                                dy += SQRT3 * a;
                                break;
                            case 4:
                                dx -= a;
                                dy += SQRT3 * a;
                                break;
                            case 5:
                                dx += a;
                                dy -= SQRT3 * a;
                                break;
                        }
                    }

                    *x = cx + dx;
                    *y = cy + dy;
                }

                const double span_z = zmax - zmin;
                if (span_z > 0.0)
                {
                    double delta = fmod(*z - zmin, span_z);
                    if (delta < 0.0)
                        delta += span_z;
                    *z = zmin + delta;
                }

                return;
            }

            if (stype == SURF_HEXY)
            {
                const double cx = S->params[0];
                const double cy = S->params[1];
                const double a = S->params[2];

                if (a > 0.0)
                {
                    double dx = *x - cx;
                    double dy = *y - cy;

                    for (int iter = 0; iter < 12; ++iter)
                    {
                        const double planes[6] = {
                            dy - a,
                            -dy - a,
                            SQRT3 * dx + dy - 2.0 * a,
                            -SQRT3 * dx - dy - 2.0 * a,
                            SQRT3 * dx - dy - 2.0 * a,
                            -SQRT3 * dx + dy - 2.0 * a
                        };

                        int max_idx = -1;
                        double max_val = 0.0;

                        for (int p = 0; p < 6; ++p)
                        {
                            if (planes[p] > max_val)
                            {
                                max_val = planes[p];
                                max_idx = p;
                            }
                        }

                        if (max_idx < 0)
                            break;

                        switch (max_idx)
                        {
                            case 0:
                                dy -= 2.0 * a;
                                break;
                            case 1:
                                dy += 2.0 * a;
                                break;
                            case 2:
                                dx -= SQRT3 * a;
                                dy -= a;
                                break;
                            case 3:
                                dx += SQRT3 * a;
                                dy += a;
                                break;
                            case 4:
                                dx -= SQRT3 * a;
                                dy += a;
                                break;
                            case 5:
                                dx += SQRT3 * a;
                                dy -= a;
                                break;
                        }
                    }

                    *x = cx + dx;
                    *y = cy + dy;
                }

                const double span_z = zmax - zmin;
                if (span_z > 0.0)
                {
                    double delta = fmod(*z - zmin, span_z);
                    if (delta < 0.0)
                        delta += span_z;
                    *z = zmin + delta;
                }

                return;
            }

            return;
        }
    }
}
