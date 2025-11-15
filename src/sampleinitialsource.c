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

    Surface *bb_surf = &DATA.surfs[DATA.boundary_surf_idx];

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
                NeutronPointSource *src = &DATA.src->npoint;

                n->E = src->E;
                n->x = src->x;
                n->y = src->y;
                n->z = src->z;
            }
            else if (DATA.src_type == SRC_FISS_MAT)
            {
                NeutronFissionSource *src = &DATA.src->fission;

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

                /* Sample fission energy from Maxwellian distribution */

                n->E = sampleMaxwellianEnergy(n, TNUC_FISSION);
            }
        }   
    }

    /* Sample from fissile material (criticality simulation) */
    
    else
    {
        long *mat_idxs = (long *)calloc(DATA.n_mats, sizeof(long));
        long n_mat_idxs = 0;

        /* Find fissile material(s) */

        for (long m = 0; m < (long)DATA.n_mats; m++)
        {
            Material *mat = &DATA.mats[m];
            bool fissile = false;

            for (size_t n = 0; n < mat->n_nucs; n++)
            {
                MaterialNuclide *nuc = &mat->nucs[n];

                if (nuc->nuc_data.has_nubar)
                {
                    fissile = true;
                    break;
                }
            }

            if (fissile)
                mat_idxs[n_mat_idxs++] = m;
            
        }

        if (n_mat_idxs <= 0)
        {
            fprintf(stderr, "Unable to sample fission source. No material with fission data found.\n");
            return EXIT_FAILURE;
        }

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


            /* Get bounds to sample in */

            double xmin = DATA.x_min;
            double xmax = DATA.x_max;
            double ymin = DATA.y_min;
            double ymax = DATA.y_max;
            double zmin = DATA.z_min;
            double zmax = DATA.z_max;

            double dx = xmax - xmin;
            double dy = ymax - ymin;
            double dz = zmax - zmin;

            /* Sample points until valid found */

            int err;
            long counter = 0;
            bool fissile;
            double x, y, z, f;
            do 
            {
                x = xmin + randd(&n->state) * dx;
                y = ymin + randd(&n->state) * dy;
                z = zmin + randd(&n->state) * dz;

                f = surfaceTest(bb_surf->type, bb_surf->params, bb_surf->n_params, x, y, z);

                long mat_idx = getMaterialAtPosition(x, y, z, &err);

                fissile = false;

                for (long m = 0; m < n_mat_idxs; m++)
                {
                    if (mat_idx == mat_idxs[m])
                    {
                        fissile = true;
                        break;
                    }
                }

                counter++;
            }
            while ((counter < 1000) && (f > 0 || !fissile));

            if (counter >= 1e3)
                return EXIT_FAILURE;

            n->x = x;
            n->y = y;
            n->z = z;
            n->E = sampleMaxwellianEnergy(n, TNUC_FISSION);
        }
    }

    if (VERBOSITY >= 1)
        fprintf(stdout, "Added %zu source neutrons into bank.\n", DATA.n_bank);

    return 0;
}
