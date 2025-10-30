#include "header.h"

long sampleInitialSource(void) {

    /* If bank capacity is zero, initialize it to a size corresponding 
       to the number of particles sized by a fixed factor */

    if (DATA.bank_cap == 0) 
    {
        DATA.bank_cap = (size_t)fmax(GLOB.n_particles * GLOB.nbuf_factor, (double)MIN_BANK_SIZE);
        if (VERBOSITY >= 1)
            fprintf(stdout, "Allocating neutron bank with capacity for %.0E neutrons.\n", (double)DATA.bank_cap);
        DATA.bank = (Neutron*)calloc(DATA.bank_cap, sizeof(Neutron));
        if (!DATA.bank) 
        {
            fprintf(stderr,"[ERROR] Memory allocation failed.\n");
            return -1;
        }
    }

    /* Reset number of neutrons in bank */

    DATA.n_bank = 0;

    /* Get geometry bounding surface */

    Surface *bb_surf = &DATA.surfs[DATA.outside_surf_idx];

    /* Sample from a user defined neutron source */

    if (GLOB.mode == RUNMODE_EXTERNAL_SOURCE) 
    {
        for (long i = 0; i < GLOB.n_particles; i++)
        {
            /* Add new neutron to bank */

            Neutron *n = &DATA.bank[DATA.n_bank++];

            /* Initialize misc. parameters */

            n->status = NEUTRON_ALIVE;
            n->id = (DATA.cur_gen - 1) * GLOB.n_particles + i;
            n->mat_idx = -1;
            n->fission_yield = 0;
            n->path_length = 0.0;
            n->fast_path_length = 0.0;
            n->time = 0.0;
            n->slowing_down_time = 0.0;
            n->genc = 0;

            /* Initialize particle private seeding */

            n->seed = GLOB.seed + n->id;
            xoshiro256ss_seed(&n->state, n->seed);

            /* Initial direction is isotropic */

            sampleNeutronDirection(n);

            /* Get properties from source */

            if (DATA.src_type == SRC_MONO_POINT) 
            {
                n->E = DATA.src->mono.E;
                n->x = DATA.src->mono.x;
                n->y = DATA.src->mono.y;
                n->z = DATA.src->mono.z;
            }
            else if (DATA.src_type == SRC_FISS_MAT)
            {
                FissileMaterialSource *src = &DATA.src->fmat;

                if (src->mat_idx < 0)
                {
                    /* Find material corresponding to mat_idx and save it for future calls */

                    for (long m = 0; m < (long)DATA.n_mats; m++)
                    {
                        Material *mat = &DATA.mats[m];

                        if (!strcmp(mat->name, src->mat_name))
                        {
                            src->mat_idx = m;
                            break;
                        }
                    }

                    if (src->mat_idx < 0)
                    {
                        fprintf(stderr, "Unable to link source to material \"%s\". Material not found.\n", src->mat_name);
                        return EXIT_FAILURE;
                    }
                }

                /* Get bounds to sample in */

                double xmin = isnan(src->xmin) ? DATA.x_min : src->xmin;
                double xmax = isnan(src->xmax) ? DATA.x_max : src->xmax;
                double ymin = isnan(src->ymin) ? DATA.y_min : src->ymin;
                double ymax = isnan(src->ymax) ? DATA.y_max : src->ymax;
                double zmin = isnan(src->zmin) ? DATA.z_min : src->zmin;
                double zmax = isnan(src->zmax) ? DATA.z_max : src->zmax;

                double dx = xmax - xmin;
                double dy = ymax - ymin;
                double dz = zmax - zmin;

                /* Sample points until valid found */

                int err;
                long counter = 0, mat_idx;
                double x, y, z, f;
                do 
                {
                    x = xmin + randd(&n->state) * dx;
                    y = ymin + randd(&n->state) * dy;
                    z = zmin + randd(&n->state) * dz;

                    f = surfaceTest(bb_surf->type, bb_surf->params, bb_surf->n_params, x, y, z);

                    mat_idx = getMaterialAtPosition(x, y, z, &err);

                    counter++;
                }
                while ((counter < 1000) && (f > 0 || mat_idx != src->mat_idx));

                if (counter >= 1e3)
                    return EXIT_FAILURE;

                n->x = x;
                n->y = y;
                n->z = z;
                n->E = sampleMaxwellianEnergy(n, TNUC_FISSION);
            }
        }   
    }

    /* Sample from fissile material (criticality simulation) */
    else
    {
        fprintf(stderr, "[ERROR] Initial source sampling for criticality mode not implemented.\n");
        return EXIT_FAILURE;
    }

    if (VERBOSITY >= 1)
        fprintf(stdout, "Added %zu source neutrons into bank.\n", DATA.n_bank);

    return 0;
}
