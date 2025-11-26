#include "header.h"
#include "detectorbuffer.h"

int runCriticalitySimulation(void) 
{
    /*
    1) On all but the first generation either:
        a) Build a fission neutron bank from fission sites of last generation
        b) Re-use last generation's bank
    2) Initialize a TransportRunScores struct to hold scores for this generation
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
        b) Use #pragma reduction to sum scores into TransportRunScores struct
    4) After all neutrons are done, add generation scores to RES.avg_scores in some valid way (single threaded?)
    5) Repeat until all generations are done
    */


    /* ########################################################################################## */

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

    /* Loop over generations (single-threaded) */

    for (long g = 1; g <= GLOB.n_generations + GLOB.n_inactive; g++) 
    {
        DATA.cur_gen = (uint64_t)g;

        /* Initialize generation score struct to reduce during parallel run */

        TransportRunScores gen_scores;
        memset(&gen_scores, 0, sizeof(TransportRunScores));

        /* On all but the first run, build a fission neutron bank from fission sites of last generation */
        
        if (g > 1) 
        {
            if (buildFissionBank() != 0) 
            {
                fprintf(stderr, "[ERROR] Failed to build fission bank for generation %zu.\n", g);
                return EXIT_FAILURE;
            }
        }
        if (!GLOB.trackplotmode)
        {
            if (g > GLOB.n_inactive)
            {
                fprintf(stdout, "\n--- Generation %ld ---\n"
                                "Simulating %zu neutrons.\n", g - GLOB.n_inactive, DATA.n_bank);
            }
            else
            {
                fprintf(stdout, "Inactive cycle %2ld/%2ld -", g, GLOB.n_inactive);
            }
        }
        /* Loop over all neutrons in bank (parallel) */

        #pragma omp parallel default(none)\
                shared(GLOB, DATA, RES, stdout, stderr, g, \
                       do_tracks, track_points, track_counts) \
                reduction(+:gen_scores)
        {
            const size_t n_detectors = DATA.n_detectors;
            DetectorHistoryBuffer *det_buffers = createDetectorHistoryBuffers(n_detectors);

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

                gen_scores.n_histories++;

                /* Sample collisions until dead */
                
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

                    n->last_mt = 0;
                    n->last_nuc = 0;

                    /* Do tracking to next collision */

                    uint64_t dt_count, dt_virtual_count, total_count;
                    double d = trackingRoutine(n, &dt_count, &dt_virtual_count, &total_count);

                    gen_scores.dt_count += dt_count;
                    gen_scores.dt_virtual_count += dt_virtual_count;
                    gen_scores.total_count += total_count;

                    if (d < 0.0)
                        break;

                    /* Score path length */

                    n->path_length += d;
                    gen_scores.total_path_length += d;

                    if (n->E >= E_THERMAL)
                    {
                        n->fast_path_length += d;
                        gen_scores.total_fast_path_length += d;
                    }

                    /* Update time */
        
                    double v = getVelocityCmPerS(n->E);
                    if (v > 0.0)
                    {
                        double dt = d / v;
                        n->time += dt;
                        gen_scores.total_time += dt;
                        if (n->E > E_THERMAL)
                        {
                            n->slowing_down_time += dt;
                            gen_scores.total_slowing_down_time += dt;
                        }
                    }

                    /* Get material at position moved to */

                    int err;
                    n->mat_idx = getMaterialAtPosition(n->x, n->y, n->z, &err);

                    /* Sample collision nuclide in active material */

                    int nuc_idx = sampleCollisionNuclide(n);
                    if (nuc_idx < 0)
                    {
                        n->status = NEUTRON_DEAD_LEAKAGE;
                        continue;
                    }

                    n->last_nuc = nuc_idx;

                    /* Sample interaction type */

                    int mt = sampleInteractionType(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);
                    if (mt < 0)
                    {
                        n->status = NEUTRON_DEAD_LEAKAGE;
                        continue;
                    }

                    n->last_mt = mt;

                    /* Score collision */

                    gen_scores.total_collisions++;

                    /* Handle interaction */

                    if (MT_IS_ELASTIC_SCATTER(mt))
                    {
                        gen_scores.total_elastic_scatters++;
                        
                        handleElasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);
                    }
                    else if (MT_IS_FISSION(mt))
                    {
                        gen_scores.total_fissions++;

                        if (n->E > E_EPITHERMAL)
                            gen_scores.total_fast_fissions++;
                        else if (n->E <= E_THERMAL)
                            gen_scores.total_thermal_fissions++;

                        handleFission(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data);

                        gen_scores.total_fission_yield += n->fission_yield;
                    }
                    else if (MT_IS_INELASTIC_SCATTER(mt))
                    {
                        gen_scores.total_inelastic_scatters++;

                        handleInelasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data, mt);
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

                    /* Score detectors using per-history buffers (active cycles only) */

                    if (det_buffers && g > GLOB.n_inactive)
                    {
                        for (size_t d = 0; d < n_detectors; ++d)
                        {
                            long bin = computeDetectorBin(n, d);
                            if (bin >= 0)
                                detectorHistoryBufferAccumulate(&det_buffers[d], (size_t)bin, 1.0);
                        }
                    }

                    /* Check if cut-off has been reached */
                    if (n->E < GLOB.energy_cutoff)
                    {
                        n->status = NEUTRON_DEAD_TERMINATED;
                    }

                }

                if (det_buffers)
                    flushDetectorHistoryBuffers(det_buffers, n_detectors);
            }

            if (det_buffers)
                destroyDetectorHistoryBuffers(det_buffers, n_detectors);
        }

        /* If in trackplotter mode, we can exit now. Secondary neutrons are not tracked. */

        if (GLOB.trackplotmode)
        {
            return EXIT_SUCCESS;
        }
            
        /* ###################################################################################### */
        /* Store generation scores if in active cycle */

        if (g > GLOB.n_inactive)
        {
            if (g - 1 - GLOB.n_inactive < RES.n_iterations && RES.avg_scores)
                RES.avg_scores[g - 1 - GLOB.n_inactive] = gen_scores;

            uint64_t histories_this_gen = gen_scores.n_histories;
            for (size_t d = 0; d < DATA.n_detectors; ++d)
                DATA.detectors[d].n_histories += histories_this_gen;
        }

        /* Print generation summary */
        
        double k_eff = (gen_scores.n_histories > 0) ? ((double)gen_scores.total_fission_yield / (double)gen_scores.n_histories) : 0.0;
        if (g > GLOB.n_inactive)
            fprintf(stdout, "Generation %ld k-eff: %.6f\n", g - GLOB.n_inactive, k_eff);
        else
            fprintf(stdout, " keff: %.6lf\n", k_eff);
    }

    /* Return succesfully */

    return EXIT_SUCCESS;
}
