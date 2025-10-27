# include "header.h"

double surfaceTest(SurfaceTypes type, double *params, size_t n_params, double x, double y, double z) {

    /* Note: validity of params or the number of them is not checked here, it must be done prior! */
    /* Return < 0 if inside, > 0 if outside, 0 if on surface */

    switch (type)
    {   
        /* YZ-plane, normal axis along X-axis*/
        case SURF_PLANEX:
            return x - params[0];

        /* XZ-plane, normal axis along Y-axis */
        case SURF_PLANEY:
            return y - params[0];

        /* XY-plane, normal axis along Z-axis */
        case SURF_PLANEZ:
            return z - params[0];

        /* General plane, either parametric or by three points*/
        case SURF_PLANE:
        {
            double G, H, I, J;

            if (n_params == 4)
            {
                /* Parametric form */
            
                G = params[0];
                H = params[1];
                I = params[2];
                J = params[3];
            }
            else
            {
                /* Defined by three points */

                double x1 = params[0];
                double y1 = params[1];
                double z1 = params[2];
                double x2 = params[3];
                double y2 = params[4];
                double z2 = params[5];
                double x3 = params[6];
                double y3 = params[7];
                double z3 = params[8];

                G = y2*z3 - y3*z2 - y1*(z3 - z2) + z1*(y3 - y2);
                H = z2*x3 - x3*x2 - z1*(x3 - x2) + x1*(z3 - z2);
                I = x2*y3 - x3*y2 - x1*(y3 - y2) + y1*(x3 - x2);
                J = -x1*(y2*z3 - y3*z2) + y1*(x2*z3 + x3*z2) - z1*(x2*y3 - x3*y2);
            }
            
            return G*x + H*y + I*z + J;
        } 

        /* Sphere */
        case SURF_SPH:
        {
            double dx = x - params[0];
            double dy = y - params[1];
            double dz = z - params[2];
            double r2 = dx*dx + dy*dy + dz*dz;
            double R2 = params[3] * params[3];
            return r2 - R2;
        }

        /* Cylinder along the X -axis */
        case SURF_CYLX:
        {
            /* First check if the point is inside the infinite cylinder */

            double dy = y - params[0];
            double dz = z - params[1];
            double r2 = dy*dy + dz*dz;
            double R2 = params[2] * params[2];

            /* If cylinder is truncated check point against heights */

            double res = r2 - R2;

            if (n_params == 5) 
            {
                double x1 = params[3];
                double x2 = params[4];
                double lower = x1 - x;
                double upper = x - x2;
                res = fmax(res, fmax(lower, upper));
            }

            return res;

        }

        /* Cylinder along the Y-axis */
        case SURF_CYLY:
        {
            /* First check if the point is inside the infinite cylinder */

            double dx = x - params[0];
            double dz = z - params[1];
            double r2 = dx*dx + dz*dz;
            double R2 = params[2] * params[2];

            /* If cylinder is truncated check point against heights */

            double res = r2 - R2;

            if (n_params == 5) 
            {
                double y1 = params[3];
                double y2 = params[4];
                double lower = y1 - y;
                double upper = y - y2;
                res = fmax(res, fmax(lower, upper));
            }

            return res;
        }

        /* Cylinder along the Z-axis */
        case SURF_CYLZ:
        {

            /* First check if the point is inside the infinite cylinder */

            double dx = x - params[0];
            double dy = y - params[1];
            double r2 = dx*dx + dy*dy;
            double R2 = params[2] * params[2];

            /* If cylinder is truncated check point against heights */

            double res = r2 - R2;

            if (n_params == 5) 
            {
                double z1 = params[3];
                double z2 = params[4];
                double lower = z1 - z;
                double upper = z - z2;
                res = fmax(res, fmax(lower, upper));
            }

            return res;
        }

        /* Square prism parallel to Z-axis*/
        case SURF_SQR:
        {
            /* First check if the point is inside an infinite square prism */

            double dx = x - params[0];
            double dy = y - params[1];
            double h = params[2];

            double a = fmax(fabs(dx), fabs(dy)) - h;

            /* If truncated check against truncating heights */

            if (n_params == 5)
            {
                double z1 = params[3];
                double z2 = params[4];

                double b = z1 - z;
                double c = z - z2;

                return fmax(a, fmax(b, c));
            }

            return a;
        }

        /* Equilateral triangular prism along the Z-axis */
        case SURF_TRI:
        {
            double x0 = params[0];
            double y0 = params[1];
            double d = params[2];
            double dx = x - x0;
            double dy = y - y0;

            double base = -dy - 0.5 * d;
            double left = dy - SQRT3 * dx - d;
            double right = dy + SQRT3 * dx - d;

            double res = fmax(base, fmax(left, right));

            if (n_params == 5)
            {
                double zmin = params[3];
                double zmax = params[4];
                double lower = zmin - z;
                double upper = z - zmax;
                res = fmax(res, fmax(lower, upper));
            }

            return res;
        }

        /* Hexagonal prism X-type */
        case SURF_HEXX:
        {
            double dx = x - params[0];
            double dy = y - params[1];
            double a = params[2];

            double side = fabs(dx) - a;
            double diag1 = fabs(dx + SQRT3 * dy) - 2.0 * a;
            double diag2 = fabs(dx - SQRT3 * dy) - 2.0 * a;

            double res = fmax(side, fmax(diag1, diag2));

            if (n_params == 5)
            {
                double z1 = params[3];
                double z2 = params[4];
                double lower = z1 - z;
                double upper = z - z2;
                res = fmax(res, fmax(lower, upper));
            }

            return res;
        }

        /* Hexagonal prism X-type */
        case SURF_HEXY:
        {
            double dx = x - params[0];
            double dy = y - params[1];
            double a = params[2];

            double side = fabs(dy) - a;
            double diag1 = fabs(SQRT3 * dx + dy) - 2.0 * a;
            double diag2 = fabs(SQRT3 * dx - dy) - 2.0 * a;

            double res = fmax(side, fmax(diag1, diag2));

            if (n_params == 5)
            {
                double z1 = params[3];
                double z2 = params[4];
                double lower = z1 - z;
                double upper = z - z2;
                res = fmax(res, fmax(lower, upper));
            }

            return res;
        }

        /* Semi-infinite cone parallel to z-axis */
        case SURF_CONE:
        {
            double x0 = params[0];
            double y0 = params[1];
            double z_section = params[2];
            double r = params[3];
            double h = params[4];

            if (h <= 0.0)
                return INFINITY;

            double z_apex = z_section + h;
            double slope = r / h;

            if (z > z_apex)
                return z - z_apex;

            double dx = x - x0;
            double dy = y - y0;
            double rho = sqrt(dx * dx + dy * dy);
            double radius = slope * (z_apex - z);

            return rho - radius;
        }

        /* Cube */
        case SURF_CUBE:
        {
            double h = params[3];
            double dx = fabs(x - params[0]) - h;
            double dy = fabs(y - params[1]) - h;
            double dz = fabs(z - params[2]) - h;

            return fmax(dx, fmax(dy, dz));
        }

        /* Cuboid */
        case SURF_CUBOID:
        {
            double xmin = params[0];
            double xmax = params[1];
            double ymin = params[2];
            double ymax = params[3];
            double zmin = params[4];
            double zmax = params[5];

            double fx = fmax(xmin - x, x - xmax);
            double fy = fmax(ymin - y, y - ymax);
            double fz = fmax(zmin - z, z - zmax);

            return fmax(fx, fmax(fy, fz));
        }

        /* Elliptical torus with major radius perpendicular to Z-axis */
        case SURF_TORUS:
        {
            double dx = x - params[0];
            double dy = y - params[1];
            double dz = z - params[2];
            double R = params[3];
            double a = params[4];
            double b = params[5];

            if (a <= 0.0 || b <= 0.0)
                return INFINITY;

            double rho = sqrt(dx * dx + dy * dy);
            double radial = rho - R;
            double term1 = (radial * radial) / (a * a);
            double term2 = (dz * dz) / (b * b);
            return term1 + term2 - 1.0;
        }

        case SURF_INF:
        {
            return -INFINITY;
        }
        default:
        {
            fprintf(stderr, "[ERROR] Surface type %d not implemented.\n", type);
            return INFINITY;
        }    
    }
}
