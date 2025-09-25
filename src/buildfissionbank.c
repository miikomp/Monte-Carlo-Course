#include "header.h"

int buildFissionBank()
{   
    /* Count fission neutrons */

    size_t n_fission_neutrons = 0;
    size_t n_fission_sites = 0;
    for (size_t i = 0; i < DATA.n_bank; i++)
    {
        Neutron *n = &DATA.bank[i];
        if (n->status == NEUTRON_DEAD_FISSION)
        {
            n_fission_neutrons += n->fission_yield;
            n_fission_sites++;
        }
    }

    /* No fissions revive all neutrons in the bank and continue */

    if (n_fission_neutrons == 0)
    {
        for (size_t i = 0; i < DATA.n_bank; i++)
        {
            Neutron *n = &DATA.bank[i];
            n->status = NEUTRON_ALIVE;
        }
        return EXIT_SUCCESS;
    }

    /* Build fission bank */

    Neutron *new_bank = (Neutron *)calloc(n_fission_neutrons, sizeof(Neutron));
    if (!new_bank)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    if (VERBOSITY >= 2)
        fprintf(stdout, "Building fission bank with %zu neutrons from %zu fission sites...\n", 
            n_fission_neutrons, n_fission_sites);

    /* Fill fission bank by creating new neutrons according to sampled fission yield at each fission site */

    size_t new_idx = 0;
    for (size_t parent_idx = 0; parent_idx < DATA.n_bank; parent_idx++)
    {
        Neutron *parent = &DATA.bank[parent_idx];
        if (parent->status != NEUTRON_DEAD_FISSION || parent->fission_yield <= 0)
            continue;

        for (int k = 0; k < parent->fission_yield; ++k)
        {
            if (new_idx >= n_fission_neutrons)
                break;

            Neutron *new_neutron = &new_bank[new_idx];

            new_neutron->status = NEUTRON_ALIVE;
            new_neutron->id = (DATA.generation - 1) * GLOB.n_particles + new_idx;
            new_neutron->mat_idx = -1;
            new_neutron->x = parent->x;
            new_neutron->y = parent->y;
            new_neutron->z = parent->z;
            new_neutron->fission_yield = 0;
            new_neutron->path_length = 0.0;

            new_neutron->seed = GLOB.seed + new_neutron->id;
            xoshiro256ss_seed(&new_neutron->state, new_neutron->seed);

            sampleNeutronDirection(new_neutron);

            new_neutron->E = sampleMaxwellianEnergy(new_neutron, TNUC_FISSION);

            ++new_idx;
        }

        if (new_idx >= n_fission_neutrons)
            break;
    }

    /* Reallocate bank if needed */

    if (n_fission_neutrons > DATA.bank_cap)
    {
        Neutron *tmp = (Neutron *)realloc(DATA.bank, n_fission_neutrons * sizeof(Neutron));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            free(new_bank);
            return 1;
        }
        DATA.bank = tmp;
        DATA.bank_cap = n_fission_neutrons;
    }

    /* Copy fission bank into neutron bank in DATA */

    memcpy(DATA.bank, new_bank, n_fission_neutrons * sizeof(Neutron));

    DATA.n_bank = n_fission_neutrons;

    /* For subcritical system expand the bank via duplication */

    if (DATA.n_bank < (size_t)GLOB.n_particles)
    {
        size_t base_count = GLOB.n_particles;
        size_t missing = 1 + base_count - DATA.n_bank;

        size_t expanded_size = base_count + missing;

        if (expanded_size > DATA.bank_cap)
        {
            Neutron *tmp = (Neutron *)realloc(DATA.bank, expanded_size * sizeof(Neutron));
            if (!tmp)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                return 1;
            }
            DATA.bank = tmp;
            DATA.bank_cap = expanded_size;
        }

        DATA.n_bank = expanded_size;

        /* Randomly sample neutrons to duplicate */

        for (size_t i = 0; i < missing; ++i)
        {
            size_t src_idx = randd(&GLOB.rng_state) * base_count;
            DATA.bank[base_count + i] = DATA.bank[src_idx];
        }

        if (VERBOSITY >= 2)
            fprintf(stdout, "Expanded neutron bank to %zu neutrons\n", DATA.n_bank);
    }

    if (VERBOSITY >= 2)
        fprintf(stdout, "DONE.\n");

    /* Free temporary bank */

    free(new_bank);

    return 0;
}
