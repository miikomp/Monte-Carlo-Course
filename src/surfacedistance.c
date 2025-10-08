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

        default:
        {
            fprintf(stderr, "[ERROR] Surface type %d not implemented.\n", type);
            return INFINITY;
        }
    }
}
