#include "header.h"

void combNeutronBank();

int runTransport(void) 
{
    /*
    1) On all but the first generation either:
        a) Build a fission neutron bank from fission sites of last generation
        b) Re-use last generation's bank
    2) Initialize a GenerationScores struct to hold scores for this generation
    3) Loop over all neutrons in the bank by comb sampling #pragma omp parallel for
        a) While neutron is alive repeatedly:
            i) Sample distance to next collision
            ii) Move neutron to collision site
            iii) Sample reaction type
            iv) Score path length, collisions, fission yield
            v) If fission, sample number of neutrons produced and add to fission bank
            vi) If capture kill neutron
            vii) If scatter sample new energy and direction
                - Target energy is either ignored, or free gas model is used (E < 200eV)
            viii) if maximum number of collisions or energy below cutoff kill neutron
        b) Use #pragma reduction to sum scores into GenerationScores struct
    4) After all neutrons are done, add generation scores to RES.avg_scores in some valid way (single threaded?)
    5) Repeat until all generations are done
    */


    /* ########################################################################################## */

    /* Loop over generations (single-threaded) */

    for (long g = 1; g <= GLOB.n_generations + GLOB.n_inactive; g++) 
    {
        DATA.generation = (uint64_t)g;

        /* Initialize generation score struct to reduce during parallel run */

        GenerationScores gen_scores;
        memset(&gen_scores, 0, sizeof(GenerationScores));

        /* On all but the first run, build a fission neutron bank from fission sites of last generation */
        
        if (g > 1) 
        {
            if (buildFissionBank() != 0) 
            {
                fprintf(stderr, "[ERROR] Failed to build fission bank for generation %zu.\n", g);
                return -1;
            }
        }

        combNeutronBank();

        if (g > GLOB.n_inactive)
        {
            fprintf(stdout, "\n--- Generation %ld ---\n"
                            "Simulating %zu neutrons.\n", g - GLOB.n_inactive, DATA.n_bank);
        }
        else
        {
            fprintf(stdout, "Inactive cycle %2ld/%2ld -", g, GLOB.n_inactive);
        }

        /* Loop over all neutrons in bank (parallel) */

        #pragma omp parallel default(none)\
                shared(GLOB, DATA, RES, stdout, stderr) \
                reduction(+:gen_scores)
        {
            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < DATA.n_bank; i++)
            {
                /* Get neutron from bank */

                Neutron *n = &DATA.bank[i];
                if (n->status != NEUTRON_ALIVE) 
                    continue;

                gen_scores.n_histories++;

                /* Sample collisions until dead */

                size_t k = 0;
                while (n->status == NEUTRON_ALIVE) 
                {
                    /* Sample distance to next collision */

                    double d = sampleDistanceToCollision(n);
                    if (d < 0.0)
                    {
                        /* Neutron is outside the geometry */

                        n->status = NEUTRON_DEAD_LEAKAGE;
                        continue;
                    }

                    /* Move neutron to collision site */

                    n->x += d * n->u;
                    n->y += d * n->v;
                    n->z += d * n->w;
                    n->path_length += d;
                    gen_scores.total_path_length += d;

                    /* Sample collision nuclide in active material */

                    int nuc_idx = sampleCollisionNuclide(n);
                    if (nuc_idx < 0)
                    {
                        n->status = NEUTRON_DEAD_LEAKAGE;
                        continue;
                    }

                    /* Sample interaction type */

                    int mt = sampleInteractionType(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);
                    if (mt < 0)
                    {
                        n->status = NEUTRON_DEAD_LEAKAGE;
                        continue;
                    }

                    /* Score collision */

                    gen_scores.total_collisions++;

                    /* For the first set number of collision score energy */
                    if (k < MAX_COLLISION_BINS) 
                    {
                        gen_scores.collision_energy_sum[k] += n->E;
                        gen_scores.collision_energy_count[k] += 1;
                        k++;
                        if (k > gen_scores.max_collision_bin)
                            gen_scores.max_collision_bin = k;
                    }

                    /* Handle interaction */

                    if (MT_IS_ELASTIC_SCATTER(mt))
                    {
                        gen_scores.total_elastic_scatters++;
                        
                        handleElasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);
                    }
                    else if (MT_IS_FISSION(mt))
                    {
                        gen_scores.total_fissions++;

                        handleFission(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);

                        gen_scores.total_fission_yield += n->fission_yield;
                        
                    }
                    else if (MT_IS_INELASTIC_SCATTER(mt))
                    {
                        gen_scores.total_inelastic_scatters++;
                        
                    }
                    else if (MT_IS_CAPTURE(mt))
                    {
                        gen_scores.total_captures++;
                        n->status = NEUTRON_DEAD_CAPTURE;
                    }
                    else
                    {
                        gen_scores.total_unknowns++;
        
                    }

                    /* Check if cut-off has been reached */
                    if (n->E < GLOB.energy_cutoff)
                    {
                        n->status = NEUTRON_DEAD_TERMINATED;
                    }

                }
            }
        }
        
        /* Store generation scores if in active cycle */
        if (g > GLOB.n_inactive)
        {
            if (g - 1 - GLOB.n_inactive < RES.n_generations && RES.avg_scores)
                RES.avg_scores[g - 1 - GLOB.n_inactive] = gen_scores;
        }

        /* Calculate k-eff */
        double k_eff = (gen_scores.n_histories > 0) ? ((double)gen_scores.total_fission_yield / (double)gen_scores.n_histories) : 0.0;
        
        if (g > GLOB.n_inactive)
            fprintf(stdout, "Generation %ld k-eff: %.6f\n", g - GLOB.n_inactive, k_eff);
        else
            fprintf(stdout, " keff: %.6lf\n", k_eff);
        
    }
    return 0;
}

/**
 * @brief Walks the built fission bank and comb samples it down to target_count neutrons.
 * 
 */
void combNeutronBank()
{
    size_t target_count = GLOB.n_particles;
    size_t bank_size = DATA.n_bank;

    if (bank_size == 0 || target_count == 0)
    {
        DATA.n_bank = 0;
        return;
    }

    if (bank_size <= target_count)
        return;

    const double stride = (double)bank_size / (double)target_count;

    double position = randd(&GLOB.rng_state) * stride;

    for (size_t write_idx = 0; write_idx < target_count; write_idx++)
    {
        size_t selected_idx = (size_t)position;
        if (selected_idx >= bank_size)
            selected_idx = bank_size - 1;

        if (selected_idx != write_idx)
            DATA.bank[write_idx] = DATA.bank[selected_idx];

        position += stride;
    }

    DATA.n_bank = target_count;
}
