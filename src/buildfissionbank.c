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

    /* No fissions, sample initial source again */

    if (n_fission_neutrons == 0)
    {
        if (sampleInitialSource() < 0)
        {
            fprintf(stderr, "[ERROR] Failed to re-sample initial source.\n");
            return EXIT_FAILURE;
        }

        fprintf(stdout, "\n[WARNING] No fission events. Re-sampling initial source\n");
        return EXIT_SUCCESS;
    }

    /* Build fission bank large enough to hold the next cur_gen */

    size_t new_bank_size = GLOB.n_particles * 1.2;
    Neutron *new_bank = (Neutron *)calloc(new_bank_size, sizeof(Neutron));
    if (!new_bank)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return EXIT_FAILURE;
    }

    size_t new_idx = 0;

    if (VERBOSITY >= 2)
        fprintf(stdout, "Building fission bank with %zu neutrons from %zu fission sites...\n", 
            n_fission_neutrons, n_fission_sites);

    /* Fill fission bank by walking the fission sites. Depending on if the system is 
       sub or supercritical, we may need to duplicate neutrons or skip sites altogether
       This is done by calculating the number of guaranteed repeats and the chance for 
       an additional repeat. For a subcritical system this means an N number of repeats and
       a possible additional repeat. For a supercritical system this means 0 guaranteed repeats
       and a chance to include that site */

    double factor = (double)GLOB.n_particles / (double)n_fission_neutrons;
    long repeats = (int)factor; // cast to int to floor
    double a = factor - repeats;

    /* Loop over sites */
    for (size_t i = 0; i < DATA.n_bank; i++)
    {
        /* Get neutron from bank */

        Neutron n = DATA.bank[i];

        /* Check if neutron caused a fission */

        if (n.status != NEUTRON_DEAD_FISSION)
            continue;

        /* Sample if an additional repeat is performed */

        double xi = randd(&n.state);
        long reps = repeats;
        if (xi < a)
            reps++;

        /* Loop over guaranteed repeats */

        for (long j = 0; j < reps; j++)
        {
            /* Loop over fission yield */

            for (int k = 0; k < n.fission_yield; k++)
            {
                /* Check if bank is full */

                if (new_idx >= new_bank_size)
                {
                    size_t expanded_size = new_bank_size * 1.2;
                    Neutron *tmp = (Neutron *)realloc(new_bank, expanded_size * sizeof(Neutron));
                    if (!tmp)
                    {
                        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                        free(new_bank);
                        return EXIT_FAILURE;
                    }
                    new_bank = tmp;
                    new_bank_size = expanded_size;
                }

                /* Create a secondary fission neutron */

                Neutron *new_neutron = &new_bank[new_idx];

                /* Initialize the neutron */

                initFissionNeutron(&n, new_neutron, new_idx);

                new_idx++;
            }

        }
    }
    if (VERBOSITY >= 2)
        fprintf(stdout, "Sampled bank size: %zu neutrons.\n", new_idx);

    /* Update global bank */
    if (new_idx > DATA.bank_cap)
    {
        Neutron *tmp = (Neutron *)realloc(DATA.bank, new_idx * sizeof(Neutron));
        if (!tmp)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            free(new_bank);
            return EXIT_FAILURE;
        }
        DATA.bank = tmp;
        DATA.bank_cap = new_idx;
    }

    DATA.n_bank = new_idx;
    memcpy(DATA.bank, new_bank, new_idx * sizeof(Neutron));

    if (VERBOSITY >= 2)
        fprintf(stdout, "DONE.\n");

    /* Free temporary bank */

    free(new_bank);

    return 0;
}
