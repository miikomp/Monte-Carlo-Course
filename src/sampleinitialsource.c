#include "header.h"

long sampleInitialSource(void) {

    /* If bank capacity is zero, initialize it to a size corresponding 
       to the number of particles sized by a fixed factor */

    if (DATA.bank_cap == 0) 
    {
        DATA.bank_cap = (size_t)fmax(GLOB.n_particles * 1.2 * GLOB.nbuf_factor, (double)MIN_BANK_SIZE);
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

    /* Sample from a user defined neutron source */

    if (GLOB.mode == RUNMODE_EXTERNAL_SOURCE) 
    {

        /* For now only Monoenergetic point source is included */

        for (long i = 0; i < GLOB.n_particles; i++)
        {
            /* Add new neutron to bank */

            Neutron *n = &DATA.bank[DATA.n_bank++];

            /* Get properties from source */

            if (DATA.src_type == SRC_MONO_POINT) {
                n->E = DATA.src->mono.E;
                n->x = DATA.src->mono.x;
                n->y = DATA.src->mono.y;
                n->z = DATA.src->mono.z;
            }

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
