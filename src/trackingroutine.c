#include "header.h"

#define SURFACE_TRACKING 1
#define DELTA_TRACKING 2

double trackingRoutine(Neutron *n)
{
    long method = DELTA_TRACKING;
    switch (method)
    {
        case SURFACE_TRACKING:
        {
            double dt = 0.0;

            /* Loop surface tracking loop until next collision is reached */

            while(1)
            {
                /* Do boundary conditions */

                applyBoundaryConditions(&n->x, &n->y, &n->z, &n->u, &n->v, &n->w);

                /* Sample distance to next collision */

                double d = sampleDistanceToCollision(n);

                if (d < 0.0)
                {
                    /* Neutron is outside the geometry */

                    n->status = NEUTRON_DEAD_LEAKAGE;

                    return (dt > 0.0) ? -dt : -1.0;
                }

                /* Get distance to nearest boundary */

                double d0 = distanceToNearestBoundary(n->x, n->y, n->z, n->u, n->v, n->w);

                /* Check if collision is beyond the boundary crossing */

                if (isfinite(d0) && d0 < d)
                {
                    /* Move over boundary */

                    double step = d0 + STEP_INTPL;

                    n->x += n->u * step;
                    n->y += n->v * step;
                    n->z += n->w * step;

                    /* Add to total distance */

                    dt += step;

                    /* Continue loop */

                    continue;
                }

                /* Move neutron to collision site */

                n->x += d * n->u;
                n->y += d * n->v;
                n->z += d * n->w;

                /* Add to total distance */

                dt += d;

                return dt;
            }
        }
        case DELTA_TRACKING:
        {
            double dt = 0.0;

            /* Loop delta-tracking loop until collision is accepted */

            while(1)
            {
                /* Get majorant cross-section */

                double sigma_m = getMajorantXS(n->E);
                if (sigma_m <= 0.0)
                {
                    n->status = NEUTRON_DEAD_LEAKAGE;
                    return (dt > 0.0) ? -dt : -1.0;
                }

                /* Sample distance to collision using majorant cross section */

                double xi = randd(&n->state);
                const double eps = 1.0e-12;
                xi = fmin(fmax(xi, eps), 1.0 - eps);
                double d = -log(1.0 - xi) / sigma_m;   

                /* Move neutron to collision site */

                n->x += d * n->u;
                n->y += d * n->v;
                n->z += d * n->w;  

                /* Add to total distance */

                dt += d;

                /* Do boundary conditions */

                applyBoundaryConditions(&n->x, &n->y, &n->z, &n->u, &n->v, &n->w);

                /* Get current material */

                int err;
                if (((n->mat_idx = getMaterialAtPosition(n->x, n->y, n->z, &err)) < 0) || err != CELL_ERR_OK)
                {
                    n->status = NEUTRON_DEAD_LEAKAGE;
                    return (dt > 0.0) ? -dt : -1.0;
                }

                /* Get total macroscopic cross section at current material */

                double sigma_t = getTotalMacroscopicXS(n->E, &DATA.mats[n->mat_idx]);
                if (sigma_t <= 0.0)
                {
                    n->status = NEUTRON_DEAD_LEAKAGE;
                    return -1.0;
                }
                
                /* Compute rejection probability */

                double P = sigma_t / sigma_m;
                if (P > 1.0)
                    P = 1.0;
                else if (P < 0.0)
                    P = 0.0;

                xi = randd(&n->state);
                xi = fmin(fmax(xi, eps), 1.0 - eps);

                /* Break loop if accepted, else experience a virtual collision and continue */

                if (xi < P)
                    return dt;
            }
        }
        default:
        {
            fprintf(stderr, "[ERROR]: Invalid tracking mode? WTF...\n");
            exit(EXIT_FAILURE);
        }
    }

    return -1.0;
}
