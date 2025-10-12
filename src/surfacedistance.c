#include "header.h"

double surfaceDistance(SurfaceTypes type, double* params, size_t n_params,
                       double x, double y, double z,
                       double u, double v, double w)
{
    if (u == 0.0 && v == 0.0 && w == 0.0)
        return INFINITY;

    const double EPS = 1e-12;

    switch (type)
    {
        /* YZ-plane, normal axis along X-axis*/
        case SURF_PLANEX:
        {
            double x0 = params[0];
            if (fabs(u) < EPS)
                return INFINITY;

            double d = (x0 - x) / u;
            return (d > EPS) ? d : INFINITY;
        }

        /* XZ-plane, normal axis along Y-axis */
        case SURF_PLANEY:
        {
            double y0 = params[0];
            if (fabs(v) < EPS)
                return INFINITY;

            double d = (y0 - y) / v;
            return (d > EPS) ? d : INFINITY;
        }

        /* XY-plane, normal axis along Z-axis */
        case SURF_PLANEZ:
        {
            double z0 = params[0];
            if (fabs(w) < EPS)
                return INFINITY;

            double d = (z0 - z) / w;
            return (d > EPS) ? d : INFINITY;
        }

        /* Arbitrary plane */
        case SURF_PLANE:
        {
            double G, H, I, J;

            if (n_params == 4)
            {
                G = params[0];
                H = params[1];
                I = params[2];
                J = params[3];
            }
            else
            {
                double x1 = params[0];
                double y1 = params[1];
                double z1 = params[2];
                double x2 = params[3];
                double y2 = params[4];
                double z2 = params[5];
                double x3 = params[6];
                double y3 = params[7];
                double z3 = params[8];

                double v1x = x2 - x1;
                double v1y = y2 - y1;
                double v1z = z2 - z1;
                double v2x = x3 - x1;
                double v2y = y3 - y1;
                double v2z = z3 - z1;

                G = v1y * v2z - v1z * v2y;
                H = v1z * v2x - v1x * v2z;
                I = v1x * v2y - v1y * v2x;
                J = -(G * x1 + H * y1 + I * z1);
            }

            double denom = G * u + H * v + I * w;
            if (fabs(denom) < EPS)
                return INFINITY;

            double numer = G * x + H * y + I * z + J;
            double t = -numer / denom;

            return (t > EPS) ? t : INFINITY;
        }

        /* Sphere */
        case SURF_SPH:
        {
            double x0 = params[0];
            double y0 = params[1];
            double z0 = params[2];
            double R = params[3];
            double R2 = R * R;

            double dx = x - x0;
            double dy = y - y0;
            double dz = z - z0;

            double A = u * u + v * v + w * w;
            double B = dx * u + dy * v + dz * w;
            double C = dx * dx + dy * dy + dz * dz - R2;

            if (A < EPS)
                return INFINITY;

            double disc = B * B - A * C;
            if (disc < -EPS)
                return INFINITY;

            if (disc < 0.0)
                disc = 0.0;

            double sqrt_disc = sqrt(disc);
            double t1 = (-B - sqrt_disc) / A;
            double t2 = (-B + sqrt_disc) / A;

            if (t1 > EPS)
                return t1;
            else if (t2 > EPS)
                return t2;
            else
                return INFINITY;
        }

        /* Cylinder along the X -axis */
        case SURF_CYLX:
        {
            double y0 = params[0];
            double z0 = params[1];
            double R = params[2];
            double R2 = R * R;

            double dy = y - y0;
            double dz = z - z0;

            double A = v * v + w * w;
            double B = dy * v + dz * w;
            double C = dy * dy + dz * dz - R2;

            double d = INFINITY;
            int truncated = (n_params == 5);
            double lower = truncated ? params[3] : 0.0;
            double upper = truncated ? params[4] : 0.0;

            if (A > EPS)
            {
                double disc = B * B - A * C;
                if (disc >= -EPS)
                {
                    if (disc < 0.0)
                        disc = 0.0;

                    double sqrt_disc = sqrt(disc);
                    double t1 = (-B - sqrt_disc) / A;
                    double t2 = (-B + sqrt_disc) / A;

                    if (t1 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t1 < d)
                                d = t1;
                        }
                        else
                        {
                            double x_hit = x + u * t1;
                            if (x_hit >= lower - EPS && x_hit <= upper + EPS && t1 < d)
                                d = t1;
                        }
                    }

                    if (t2 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t2 < d)
                                d = t2;
                        }
                        else
                        {
                            double x_hit = x + u * t2;
                            if (x_hit >= lower - EPS && x_hit <= upper + EPS && t2 < d)
                                d = t2;
                        }
                    }
                }
            }

            if (truncated)
            {
                if (fabs(u) >= EPS)
                {
                    double t_cap = (lower - x) / u;
                    if (t_cap > EPS)
                    {
                        double y_cap = y + v * t_cap;
                        double z_cap = z + w * t_cap;
                        double dyy = y_cap - y0;
                        double dzz = z_cap - z0;
                        if (dyy * dyy + dzz * dzz <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }

                    t_cap = (upper - x) / u;
                    if (t_cap > EPS)
                    {
                        double y_cap = y + v * t_cap;
                        double z_cap = z + w * t_cap;
                        double dyy = y_cap - y0;
                        double dzz = z_cap - z0;
                        if (dyy * dyy + dzz * dzz <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }
                }
            }

            return d;
        }

        /* Cylinder along the Y -axis */
        case SURF_CYLY:
        {
            double x0 = params[0];
            double z0 = params[1];
            double R = params[2];
            double R2 = R * R;

            double dx = x - x0;
            double dz = z - z0;

            double A = u * u + w * w;
            double B = dx * u + dz * w;
            double C = dx * dx + dz * dz - R2;

            double d = INFINITY;
            int truncated = (n_params == 5);
            double lower = truncated ? params[3] : 0.0;
            double upper = truncated ? params[4] : 0.0;

            if (A > EPS)
            {
                double disc = B * B - A * C;
                if (disc >= -EPS)
                {
                    if (disc < 0.0)
                        disc = 0.0;

                    double sqrt_disc = sqrt(disc);
                    double t1 = (-B - sqrt_disc) / A;
                    double t2 = (-B + sqrt_disc) / A;

                    if (t1 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t1 < d)
                                d = t1;
                        }
                        else
                        {
                            double y_hit = y + v * t1;
                            if (y_hit >= lower - EPS && y_hit <= upper + EPS && t1 < d)
                                d = t1;
                        }
                    }

                    if (t2 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t2 < d)
                                d = t2;
                        }
                        else
                        {
                            double y_hit = y + v * t2;
                            if (y_hit >= lower - EPS && y_hit <= upper + EPS && t2 < d)
                                d = t2;
                        }
                    }
                }
            }

            if (truncated)
            {
                if (fabs(v) >= EPS)
                {
                    double t_cap = (lower - y) / v;
                    if (t_cap > EPS)
                    {
                        double x_cap = x + u * t_cap;
                        double z_cap = z + w * t_cap;
                        double dxx = x_cap - x0;
                        double dzz = z_cap - z0;
                        if (dxx * dxx + dzz * dzz <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }

                    t_cap = (upper - y) / v;
                    if (t_cap > EPS)
                    {
                        double x_cap = x + u * t_cap;
                        double z_cap = z + w * t_cap;
                        double dxx = x_cap - x0;
                        double dzz = z_cap - z0;
                        if (dxx * dxx + dzz * dzz <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }
                }
            }

            return d;
        }

        /* Cylinder along the Z -axis */
        case SURF_CYLZ:
        {
            double x0 = params[0];
            double y0 = params[1];
            double R = params[2];
            double R2 = R * R;

            double dx = x - x0;
            double dy = y - y0;

            double A = u * u + v * v;
            double B = dx * u + dy * v;
            double C = dx * dx + dy * dy - R2;

            double d = INFINITY;
            int truncated = (n_params == 5);
            double lower = truncated ? params[3] : 0.0;
            double upper = truncated ? params[4] : 0.0;

            if (A > EPS)
            {
                double disc = B * B - A * C;
                if (disc >= -EPS)
                {
                    if (disc < 0.0)
                        disc = 0.0;

                    double sqrt_disc = sqrt(disc);
                    double t1 = (-B - sqrt_disc) / A;
                    double t2 = (-B + sqrt_disc) / A;

                    if (t1 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t1 < d)
                                d = t1;
                        }
                        else
                        {
                            double z_hit = z + w * t1;
                            if (z_hit >= lower - EPS && z_hit <= upper + EPS && t1 < d)
                                d = t1;
                        }
                    }

                    if (t2 > EPS)
                    {
                        if (!truncated)
                        {
                            if (t2 < d)
                                d = t2;
                        }
                        else
                        {
                            double z_hit = z + w * t2;
                            if (z_hit >= lower - EPS && z_hit <= upper + EPS && t2 < d)
                                d = t2;
                        }
                    }
                }
            }

            if (truncated)
            {
                if (fabs(w) >= EPS)
                {
                    double t_cap = (lower - z) / w;
                    if (t_cap > EPS)
                    {
                        double x_cap = x + u * t_cap;
                        double y_cap = y + v * t_cap;
                        double dxx = x_cap - x0;
                        double dyy = y_cap - y0;
                        if (dxx * dxx + dyy * dyy <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }

                    t_cap = (upper - z) / w;
                    if (t_cap > EPS)
                    {
                        double x_cap = x + u * t_cap;
                        double y_cap = y + v * t_cap;
                        double dxx = x_cap - x0;
                        double dyy = y_cap - y0;
                        if (dxx * dxx + dyy * dyy <= R2 + EPS && t_cap < d)
                            d = t_cap;
                    }
                }
            }

            return d;
        }

        /* Square prism along the Z-axis */
        case SURF_SQR:
        {
            double x0 = params[0];
            double y0 = params[1];
            double h = params[2];

            double xmin = x0 - h;
            double xmax = x0 + h;
            double ymin = y0 - h;
            double ymax = y0 + h;

            int truncated = (n_params == 5);
            double zmin = truncated ? params[3] : -INFINITY;
            double zmax = truncated ? params[4] : INFINITY;

            double d = INFINITY;

            if (fabs(u) >= EPS)
            {
                double x_plane = (u > 0.0) ? xmax : xmin;
                double tx = (x_plane - x) / u;
                if (tx > EPS)
                {
                    double y_hit = y + v * tx;
                    double z_hit = z + w * tx;
                    if (y_hit >= ymin - EPS && y_hit <= ymax + EPS &&
                        z_hit >= zmin - EPS && z_hit <= zmax + EPS &&
                        tx < d)
                        d = tx;
                }
            }
            else if (x < xmin - EPS || x > xmax + EPS)
                return INFINITY;

            if (fabs(v) >= EPS)
            {
                double y_plane = (v > 0.0) ? ymax : ymin;
                double ty = (y_plane - y) / v;
                if (ty > EPS)
                {
                    double x_hit = x + u * ty;
                    double z_hit = z + w * ty;
                    if (x_hit >= xmin - EPS && x_hit <= xmax + EPS &&
                        z_hit >= zmin - EPS && z_hit <= zmax + EPS &&
                        ty < d)
                        d = ty;
                }
            }
            else if (y < ymin - EPS || y > ymax + EPS)
                return INFINITY;

            if (truncated)
            {
                if (fabs(w) >= EPS)
                {
                    double z_plane = (w > 0.0) ? zmax : zmin;
                    double tz = (z_plane - z) / w;
                    if (tz > EPS)
                    {
                        double x_hit = x + u * tz;
                        double y_hit = y + v * tz;
                        if (x_hit >= xmin - EPS && x_hit <= xmax + EPS &&
                            y_hit >= ymin - EPS && y_hit <= ymax + EPS &&
                            tz < d)
                            d = tz;
                    }
                }
                else if (z < zmin - EPS || z > zmax + EPS)
                    return INFINITY;
            }

            return d;
        }

        /* Hexagonal prism X-type along the Z-axis */
        case SURF_HEXX:
        {
            double x0 = params[0];
            double y0 = params[1];
            double a = params[2];
            int truncated = (n_params == 5);
            double zmin = truncated ? params[3] : -INFINITY;
            double zmax = truncated ? params[4] : INFINITY;

            const double normals[6][2] = {
                { 1.0,        0.0},
                {-1.0,        0.0},
                { 1.0,  SQRT3},
                {-1.0, -SQRT3},
                { 1.0, -SQRT3},
                {-1.0,  SQRT3}
            };

            double offsets[6];
            offsets[0] = -(x0 + a);
            offsets[1] =  (x0 - a);
            offsets[2] = -(x0 + SQRT3 * y0 + 2.0 * a);
            offsets[3] =  (x0 + SQRT3 * y0 - 2.0 * a);
            offsets[4] = -(x0 - SQRT3 * y0 + 2.0 * a);
            offsets[5] =  (x0 - SQRT3 * y0 - 2.0 * a);

            double d = INFINITY;

            for (int i = 0; i < 6; ++i)
            {
                double nx = normals[i][0];
                double ny = normals[i][1];
                double denom = nx * u + ny * v;
                if (fabs(denom) < EPS)
                    continue;

                double numer = nx * x + ny * y + offsets[i];
                double t = -numer / denom;
                if (t <= EPS || t >= d)
                    continue;

                double hx = x + u * t;
                double hy = y + v * t;
                double hz = z + w * t;

                double dx = hx - x0;
                double dy = hy - y0;

                double side = fabs(dx) - a;
                double diag1 = fabs(dx + SQRT3 * dy) - 2.0 * a;
                double diag2 = fabs(dx - SQRT3 * dy) - 2.0 * a;

                if (side <= EPS && diag1 <= EPS && diag2 <= EPS &&
                    hz >= zmin - EPS && hz <= zmax + EPS)
                    d = t;
            }

            if (truncated)
            {
                if (fabs(w) >= EPS)
                {
                    double z_plane = (w > 0.0) ? zmax : zmin;
                    double tz = (z_plane - z) / w;
                    if (tz > EPS && tz < d)
                    {
                        double hx = x + u * tz;
                        double hy = y + v * tz;

                        double dx = hx - x0;
                        double dy = hy - y0;

                        double side = fabs(dx) - a;
                        double diag1 = fabs(dx + SQRT3 * dy) - 2.0 * a;
                        double diag2 = fabs(dx - SQRT3 * dy) - 2.0 * a;

                        if (side <= EPS && diag1 <= EPS && diag2 <= EPS)
                            d = tz;
                    }
                }
                else if (z < zmin - EPS || z > zmax + EPS)
                    return INFINITY;
            }

            return d;
        }

        /* Hexagonal prism Y-type along the Z-axis */
        case SURF_HEXY:
        {
            double x0 = params[0];
            double y0 = params[1];
            double a = params[2];
            int truncated = (n_params == 5);
            double zmin = truncated ? params[3] : -INFINITY;
            double zmax = truncated ? params[4] : INFINITY;

            const double normals[6][2] = {
                {0.0,         1.0},
                {0.0,        -1.0},
                {SQRT3,       1.0},
                {-SQRT3,     -1.0},
                {SQRT3,      -1.0},
                {-SQRT3,      1.0}
            };

            double offsets[6];
            offsets[0] = -(y0 + a);
            offsets[1] =  (y0 - a);
            offsets[2] = -(SQRT3 * x0 + y0 + 2.0 * a);
            offsets[3] =  (SQRT3 * x0 + y0 - 2.0 * a);
            offsets[4] = -(SQRT3 * x0 - y0 + 2.0 * a);
            offsets[5] =  (SQRT3 * x0 - y0 - 2.0 * a);

            double d = INFINITY;

            for (int i = 0; i < 6; ++i)
            {
                double nx = normals[i][0];
                double ny = normals[i][1];
                double denom = nx * u + ny * v;
                if (fabs(denom) < EPS)
                    continue;

                double numer = nx * x + ny * y + offsets[i];
                double t = -numer / denom;
                if (t <= EPS || t >= d)
                    continue;

                double hx = x + u * t;
                double hy = y + v * t;
                double hz = z + w * t;

                double dx = hx - x0;
                double dy = hy - y0;

                double side = fabs(dy) - a;
                double diag1 = fabs(SQRT3 * dx + dy) - 2.0 * a;
                double diag2 = fabs(SQRT3 * dx - dy) - 2.0 * a;

                if (side <= EPS && diag1 <= EPS && diag2 <= EPS &&
                    hz >= zmin - EPS && hz <= zmax + EPS)
                    d = t;
            }

            if (truncated)
            {
                if (fabs(w) >= EPS)
                {
                    double z_plane = (w > 0.0) ? zmax : zmin;
                    double tz = (z_plane - z) / w;
                    if (tz > EPS && tz < d)
                    {
                        double hx = x + u * tz;
                        double hy = y + v * tz;

                        double dx = hx - x0;
                        double dy = hy - y0;

                        double side = fabs(dy) - a;
                        double diag1 = fabs(SQRT3 * dx + dy) - 2.0 * a;
                        double diag2 = fabs(SQRT3 * dx - dy) - 2.0 * a;

                        if (side <= EPS && diag1 <= EPS && diag2 <= EPS)
                            d = tz;
                    }
                }
                else if (z < zmin - EPS || z > zmax + EPS)
                    return INFINITY;
            }

            return d;
        }

        /* Cube */
        case SURF_CUBE:
        {
            double x0 = params[0];
            double y0 = params[1];
            double z0 = params[2];
            double h  = params[3];

            double xmin = x0 - h;
            double xmax = x0 + h;
            double ymin = y0 - h;
            double ymax = y0 + h;
            double zmin = z0 - h;
            double zmax = z0 + h;

            double t_enter = -INFINITY;
            double t_exit  = INFINITY;

            if (fabs(u) < EPS)
            {
                if (x < xmin - EPS || x > xmax + EPS)
                    return INFINITY;
            }
            else
            {
                double t1 = (xmin - x) / u;
                double t2 = (xmax - x) / u;
                if (t1 > t2)
                {
                    double tmp = t1;
                    t1 = t2;
                    t2 = tmp;
                }
                if (t1 > t_enter)
                    t_enter = t1;
                if (t2 < t_exit)
                    t_exit = t2;
            }

            if (fabs(v) < EPS)
            {
                if (y < ymin - EPS || y > ymax + EPS)
                    return INFINITY;
            }
            else
            {
                double t1 = (ymin - y) / v;
                double t2 = (ymax - y) / v;
                if (t1 > t2)
                {
                    double tmp = t1;
                    t1 = t2;
                    t2 = tmp;
                }
                if (t1 > t_enter)
                    t_enter = t1;
                if (t2 < t_exit)
                    t_exit = t2;
            }

            if (fabs(w) < EPS)
            {
                if (z < zmin - EPS || z > zmax + EPS)
                    return INFINITY;
            }
            else
            {
                double t1 = (zmin - z) / w;
                double t2 = (zmax - z) / w;
                if (t1 > t2)
                {
                    double tmp = t1;
                    t1 = t2;
                    t2 = tmp;
                }
                if (t1 > t_enter)
                    t_enter = t1;
                if (t2 < t_exit)
                    t_exit = t2;
            }

            if (t_enter > t_exit)
                return INFINITY;

            if (t_exit <= EPS)
                return INFINITY;

            double t_hit = (t_enter > EPS) ? t_enter : t_exit;

            return (t_hit > EPS && t_hit < INFINITY) ? t_hit : INFINITY;
        }

        /* Cuboid between (x1, x2) & (y1 y2) & (z1 z2) */
        case SURF_CUBOID:
        {
            double xmin = params[0];
            double xmax = params[1];
            double ymin = params[2];
            double ymax = params[3];
            double zmin = params[4];
            double zmax = params[5];

            double d = INFINITY;

            if (fabs(u) >= EPS)
            {
                double x_plane = (u > 0.0) ? xmax : xmin;
                double tx = (x_plane - x) / u;
                if (tx > EPS)
                {
                    double y_hit = y + v * tx;
                    double z_hit = z + w * tx;
                    if (y_hit >= ymin - EPS && y_hit <= ymax + EPS &&
                        z_hit >= zmin - EPS && z_hit <= zmax + EPS &&
                        tx < d)
                        d = tx;
                }
            }
            else if (x < xmin - EPS || x > xmax + EPS)
                return INFINITY;

            if (fabs(v) >= EPS)
            {
                double y_plane = (v > 0.0) ? ymax : ymin;
                double ty = (y_plane - y) / v;
                if (ty > EPS)
                {
                    double x_hit = x + u * ty;
                    double z_hit = z + w * ty;
                    if (x_hit >= xmin - EPS && x_hit <= xmax + EPS &&
                        z_hit >= zmin - EPS && z_hit <= zmax + EPS &&
                        ty < d)
                        d = ty;
                }
            }
            else if (y < ymin - EPS || y > ymax + EPS)
                return INFINITY;

            if (fabs(w) >= EPS)
            {
                double z_plane = (w > 0.0) ? zmax : zmin;
                double tz = (z_plane - z) / w;
                if (tz > EPS)
                {
                    double x_hit = x + u * tz;
                    double y_hit = y + v * tz;
                    if (x_hit >= xmin - EPS && x_hit <= xmax + EPS &&
                        y_hit >= ymin - EPS && y_hit <= ymax + EPS &&
                        tz < d)
                            d = tz;
                }
            }
            else if (z < zmin - EPS || z > zmax + EPS)
                return INFINITY;

            return d;
        }
        default:
        {
            fprintf(stderr, "[ERROR] Surface type %d not implemented.\n", type);
            return INFINITY;
        }
    }
}
