#include "header.h"

#define SURFACE_TRACKING 1
#define DELTA_TRACKING 2

double trackingRoutine(Neutron *n, uint64_t *dt_count, uint64_t *dt_vcount, uint64_t *count)
{
    /* Counter for total travelled distance */

    double dt = 0.0;

    /* Reset counters */

    if (dt_count)
        *dt_count = 0;
    if (dt_vcount)
        *dt_vcount = 0;
    if (count)
        *count = 0;

    while(1)
    {
        /* Apply boundary conditions */
        
        applyBoundaryConditions(&n->x, &n->y, &n->z, &n->u, &n->v, &n->w);

        /* Get current material */

        int err;
        if (((n->mat_idx = getMaterialAtPosition(n->x, n->y, n->z, &err)) < 0) || err != CELL_ERR_OK)
        {
            n->status = NEUTRON_DEAD_LEAKAGE;
            return -1.0;
        }

        /* Get majorant cross-section for the incident neutron energy */

        double sigma_m = getMajorantXS(n->E);
        if (sigma_m <= 0.0)
        {
            n->status = NEUTRON_DEAD_LEAKAGE;
            return -1.0;
        }

        /* Get total macroscopic cross section at current material */

        double sigma_t = getTotalMacroscopicXS(n->E, &DATA.mats[n->mat_idx]);
        if (sigma_t < 0.0)
        {
            n->status = NEUTRON_DEAD_LEAKAGE;
            return -1.0;
        }

        /* Compute delta-tracking rejection probability */

        double P = sigma_t / sigma_m;
        if (P > 1.0)
            P = 1.0;
        else if (P < 0.0)
            P = 0.0;

        if (count)
            (*count)++;

        /* --- Use either delta- or surface-tracking --- */

        if (P >= 1.0 - DT_THRESHOLD)
        {
            /* Delta-tracking */

            if (dt_count)
                (*dt_count)++;

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
            long new_mat;
            if (((new_mat = getMaterialAtPosition(n->x, n->y, n->z, &err)) < 0) || err != CELL_ERR_OK)
            {
                n->status = NEUTRON_DEAD_LEAKAGE;
                return -1.0;
            }

            /* Check if material boundary crossed */

            if (new_mat != n->mat_idx)
            {
                n->mat_idx = new_mat;

                /* Get new total macroscopic cross section at current material */

                sigma_t = getTotalMacroscopicXS(n->E, &DATA.mats[n->mat_idx]);
                if (sigma_t <= 0.0)
                {
                    n->status = NEUTRON_DEAD_LEAKAGE;
                    return -1.0;
                }
                
                /* Compute new rejection probability */

                P = sigma_t / sigma_m;
                if (P > 1.0)
                    P = 1.0;
                else if (P < 0.0)
                    P = 0.0;
            } 

            /* Check if collision is accepted */

            xi = randd(&n->state);
            xi = fmin(fmax(xi, eps), 1.0 - eps);

            /* Break loop if accepted, else experience a virtual collision and continue */

            if (xi < P)
                return dt;

            if (dt_vcount)
                (*dt_vcount)++;

            continue;
        }
        else
        {
            /* Surface tracking */

            /* Do boundary conditions */

            applyBoundaryConditions(&n->x, &n->y, &n->z, &n->u, &n->v, &n->w);

            /* Sample distance to next collision */

            double d = sampleDistanceToCollision(n);

            if (d < 0.0)
            {
                /* Neutron is outside the geometry */
                
                n->status = NEUTRON_DEAD_LEAKAGE;

                return -1.0;
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

    return -1.0;
}
