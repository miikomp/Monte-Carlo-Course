#include "header.h"

long sampleInitialSource(void) {

    fprintf(stdout, "\nSampling initial neutron source...\n");

    /* If bank capacity is zero, initialize it to value slightly more than the current 
    number of particles to avoid having to realloc during a run */

    if (DATA.bank_cap == 0) {
        DATA.bank_cap = (size_t)fmax(GLOB.n_particles * 1.1, (double)MIN_BANK_SIZE);
        DATA.bank = (Neutron*)calloc(DATA.bank_cap, sizeof(Neutron));
        if (!DATA.bank) {
            fprintf(stderr,"[ERROR] Memory allocation failed.\n");
            return -1;
        }
    }

    /* Build a bank of neutrons by sampling a source */

    if (DATA.src != NULL) 
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
            n->id = (DATA.generation - 1) * GLOB.n_particles + i;
            n->mat_idx = -1;
            n->fission_yield = 0;
            n->path_length = 0.0;
            n->fast_path_length = 0.0;
            n->time = 0.0;
            n->slowing_down_time = 0.0;

            /* Initialize particle private seeding */

            n->seed = GLOB.seed + n->id;
            xoshiro256ss_seed(&n->state, n->seed);

            /* Initial direction is isotropic */

            sampleNeutronDirection(n);
        }   
    }

    if (VERBOSITY >= 1)
        fprintf(stdout, "Added %zu source neutrons into bank.\n", DATA.n_bank);

    /* Calculate and report total memory allocated for the bank */

    size_t total_bytes = 0;
    for (size_t idx = 0; idx < DATA.bank_cap; ++idx)
        total_bytes += sizeof(Neutron);

    GLOB.mem_nbank = total_bytes;
    fprintf(stdout, "Memory allocated for neutron bank: %.2f MB\n", (double)total_bytes / (1024.0 * 1024.0));

    fprintf(stdout, "DONE.\n");

    return 0;
}
