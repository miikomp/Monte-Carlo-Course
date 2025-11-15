#include "header.h"

int runExternalSourceSimulation(void)
{
    /* Allocate thread-local buffers*/

    Neutron **thread_bufs = (Neutron **)calloc((size_t)GLOB.n_threads, sizeof(Neutron *));
    size_t  *thread_buf_count = (size_t *)calloc((size_t)GLOB.n_threads, sizeof(size_t));
    size_t  *thread_buf_cap   = (size_t *)calloc((size_t)GLOB.n_threads, sizeof(size_t));

    if (!thread_bufs || !thread_buf_count || !thread_buf_cap)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed for neutron buffers.\n");
        free(thread_bufs);
        free(thread_buf_count);
        free(thread_buf_cap);
        return EXIT_FAILURE;
    }
    
    size_t base_cap = 2.0 * (double)GLOB.n_particles / (double)GLOB.n_threads;

    for (int t = 0; t < GLOB.n_threads; ++t)
    {
        thread_buf_cap[t] = base_cap;
        thread_bufs[t] = (Neutron *)calloc(thread_buf_cap[t], sizeof(Neutron));
        if (!thread_bufs[t])
        {
            fprintf(stderr, "[ERROR] Memory allocation failed for neutron buffer.\n");
            for (int j = 0; j < t; ++j)
                free(thread_bufs[j]);
            free(thread_bufs);
            free(thread_buf_count);
            free(thread_buf_cap);
            return EXIT_FAILURE;
        }
    }

    /* If doing track plots grab ptrs to arrays */

    double *track_points = NULL;
    size_t *track_counts = NULL;
    bool do_tracks = false;

    if (GLOB.trackplotmode && DATA.tracks != NULL && DATA.track_counts != NULL)
    {
        track_points = DATA.tracks;
        track_counts = DATA.track_counts;
        do_tracks = true;
    }


    for (long c = 1; c <= GLOB.n_cycles + GLOB.n_inactive; c++)
    {
        DATA.cur_cycle = (uint64_t)c;
        DATA.cur_gen = 1;

        TransportRunScores cycle_scores;
        memset(&cycle_scores, 0, sizeof(TransportRunScores));

        if (c > 1 && sampleInitialSource() != 0)
        {
            fprintf(stderr, "[ERROR] Failed to sample initial source for cycle %ld.\n", c);
            exit(EXIT_FAILURE);
            break;
        }

        while (DATA.n_bank > 0)
        {
            memset(thread_buf_count, 0, (size_t)GLOB.n_threads * sizeof(size_t));

            #pragma omp parallel default(none) \
                    shared(GLOB, DATA, RES, stdout, stderr, thread_bufs, thread_buf_cap, thread_buf_count, base_cap, c, track_counts, track_points, do_tracks) \
                    reduction(+:cycle_scores)
            {
                /* Get pointer to thread-local buffers */

                const int tid = omp_get_thread_num();

                Neutron *local_buf = thread_bufs[tid];
                size_t local_cap = thread_buf_cap[tid];
                size_t local_count = 0;

                /* Reset thread-local buffer count */

                thread_buf_count[tid] = 0;

                /* Do not allow threads to progress until all thread buffer counts are reset */

                #pragma omp barrier

                /* Loop through the neutrons in the bank */

                #pragma omp for schedule(dynamic)
                for (size_t i = 0; i < DATA.n_bank; i++)
                {
                    double *points = NULL;
                    size_t *count = NULL;

                    if (do_tracks && i < (size_t)GLOB.n_tracks)
                    {
                        points = track_points + (size_t)i * (MAX_COLLISION_BINS + 1u) * 3u;
                        count  = &track_counts[i];
                    }

                    /* Get neutron from bank */

                    Neutron *n = &DATA.bank[i];
                    if (n->status != NEUTRON_ALIVE)
                        continue;

                    cycle_scores.n_histories++;

                    /* Transport neutron until death */

                    while (n->status == NEUTRON_ALIVE)
                    {
                        /* Add point to track segment array if plotting tracks */

                        if (do_tracks && count && *count < (MAX_COLLISION_BINS + 1))
                        {
                            double *slot = points + (*count) * 3u;
                            slot[0] = n->x;
                            slot[1] = n->y;
                            slot[2] = n->z;
                            (*count)++;
                        }

                        /* Do boundary conditions */
                        
                        applyBoundaryConditions(&n->x, &n->y, &n->z, &n->u, &n->v, &n->w);
                    
                        /* Sample distance to next collision within current material */

                        double d = sampleDistanceToCollision(n);

                        /* Check distance to nearest boundary */

                        double d0 = distanceToNearestBoundary(n->x, n->y, n->z, n->u, n->v, n->w);

                        /* Check if collision is after boundary cross */

                        if (isfinite(d0) && d0 < d)
                        {
                            /* Move over  boundary*/

                            n->x += n->u * (d0 + STEP_INTPL);
                            n->y += n->v * (d0 + STEP_INTPL);
                            n->z += n->w * (d0 + STEP_INTPL);

                            /* Continue loop */

                            continue;
                        }

                        if (d < 0.0)
                        {
                            n->status = NEUTRON_DEAD_LEAKAGE;
                            continue;
                        }
                        

                        /* Move neutron */

                        n->x += d * n->u;
                        n->y += d * n->v;
                        n->z += d * n->w;

                        n->path_length += d;
                        cycle_scores.total_path_length += d;

                        if (n->E > E_THERMAL)
                        {
                            n->fast_path_length += d;
                            cycle_scores.total_fast_path_length += d;
                        }

                        double v = getVelocityCmPerS(n->E);
                        if (v > 0.0)
                        {
                            double dt = d / v;
                            n->time += dt;
                            cycle_scores.total_time += dt;
                            if (n->E > E_THERMAL)
                            {
                                n->slowing_down_time += dt;
                                cycle_scores.total_slowing_down_time += dt;
                            }
                        }

                        /* Get material at current position */

                        int err;
                        n->mat_idx = getMaterialAtPosition(n->x, n->y, n->z, &err);

                        if (n->mat_idx < 0)
                            continue;

                        /* Sample collision nuclide */

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

                        cycle_scores.total_collisions++;

                        /* Handle interaction */

                        if (MT_IS_ELASTIC_SCATTER(mt))
                        {
                            cycle_scores.total_elastic_scatters++;
                            handleElasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);
                        }
                        else if (MT_IS_FISSION(mt))
                        {
                            cycle_scores.total_fissions++;

                            if (n->E > E_EPITHERMAL)
                                cycle_scores.total_fast_fissions++;
                            else if (n->E <= E_THERMAL)
                                cycle_scores.total_thermal_fissions++;

                            handleFission(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);

                            if (n->fission_yield > 0)
                            {
                                cycle_scores.total_fission_yield += n->fission_yield;

                                size_t required = local_count + (size_t)n->fission_yield;
                                if (required > local_cap)
                                {
                                    size_t new_cap = local_cap ? local_cap : base_cap;
                                    while (new_cap < required)
                                    {
                                        if (new_cap == 0)
                                            new_cap = base_cap;
                                        else
                                            new_cap *= 2;
                                        if (new_cap < base_cap)
                                            new_cap = base_cap;
                                    }

                                    Neutron *tmp = (Neutron *)realloc(local_buf, new_cap * sizeof(Neutron));
                                    if (!tmp)
                                        exit(EXIT_FAILURE);
                                    local_buf = tmp;
                                    local_cap = new_cap;
                                    thread_bufs[tid] = tmp;
                                    
                                }
                                for (int f = 0; f < n->fission_yield; ++f)
                                {
                                    Neutron *new_neutron = &local_buf[local_count++];
                                    initFissionNeutron(n, new_neutron, (long)(DATA.cur_gen + c + f));
                                }
                                
                            }
                        }
                        else if (MT_IS_INELASTIC_SCATTER(mt))
                        {
                            cycle_scores.total_inelastic_scatters++;
                            handleInelasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data, mt);
                        }
                        else if (MT_IS_CAPTURE(mt))
                        {
                            cycle_scores.total_captures++;
                            n->status = NEUTRON_DEAD_CAPTURE;
                        }
                        else
                        {
                            cycle_scores.total_unknowns++;
                        }

                        /* If neutron is alive, check for cut-off */

                        if (n->status == NEUTRON_ALIVE)
                        {
                            if (checkNeutronCutoff(n) != 0)
                                cycle_scores.total_terminated++;
                        }

                    }
                }

                /* Update thread-local buffer counts to main array */

                thread_buf_count[tid] = local_count;
                thread_buf_cap[tid] = local_cap;
            }

            /* If in trackplotter mode, we can exit now. Secondary neutrons are not tracked. */

            if (GLOB.trackplotmode)
                return EXIT_SUCCESS;

            /* Combine buffers into neutron bank */

            size_t new_bank_size = 0;
            for (int t = 0; t < GLOB.n_threads; ++t)
                new_bank_size += thread_buf_count[t];

            if (new_bank_size == 0)
            {
                DATA.n_bank = 0;
                break;
            }

            /* Check if need to reallocate */

            if (new_bank_size > DATA.bank_cap)
            {
                size_t expanded_size = (size_t)(new_bank_size * 1.2);
                if (expanded_size < new_bank_size)
                    expanded_size = new_bank_size;
                Neutron *tmp = (Neutron *)realloc(DATA.bank, expanded_size * sizeof(Neutron));
                if (!tmp)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed for neutron bank.\n");
                    exit(EXIT_FAILURE);
                }
                DATA.bank = tmp;
                DATA.bank_cap = expanded_size;
            }

            /* Copy neutrons from thread buffers to neutron bank */

            size_t bank_index = 0;
            for (int t = 0; t < GLOB.n_threads; ++t)
            {
                size_t count = thread_buf_count[t];
                if (count == 0)
                    continue;
                memcpy(&DATA.bank[bank_index], thread_bufs[t], count * sizeof(Neutron));
                bank_index += count;
            }

            DATA.n_bank = bank_index;
            DATA.cur_gen += 1;

            if (VERBOSITY >= 2)
                fprintf(stdout, " Cycle %ld - Gen %ld: %zu neutrons left.\n", c, DATA.cur_gen, DATA.n_bank);
        }

        if (c > GLOB.n_inactive)
        {
            long idx = c - 1 - GLOB.n_inactive;
            if (idx >= 0 && idx < RES.n_iterations && RES.avg_scores)
                RES.avg_scores[idx] = cycle_scores;
        }

        /* Print cycle summary*/

        uint64_t active_histories = 0;
        if (cycle_scores.n_histories > cycle_scores.total_terminated)
            active_histories = cycle_scores.n_histories - cycle_scores.total_terminated;

        double keff = active_histories > 0 ? (double)cycle_scores.total_fission_yield / (double)active_histories : 0.0;
        fprintf(stdout, "---- Cycle %2ld/%2ld - keff: %.6f\n",
                c, GLOB.n_cycles + GLOB.n_inactive, keff);
    }

    /* Prepare for exit */

    for (int i = 0; i < GLOB.n_threads; ++i)
        free(thread_bufs[i]);
    free(thread_bufs);
    free(thread_buf_count);
    free(thread_buf_cap);

    return EXIT_SUCCESS;
}
