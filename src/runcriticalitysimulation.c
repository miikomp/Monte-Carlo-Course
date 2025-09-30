#include "header.h"

inline size_t energyGridBinIndex(const EnergyGrid *grid, double energy);

int runCriticalitySimulation(void) 
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

    /* Get pointer to detector array */

    const size_t n_detectors = DATA.n_detectors;
    Detector **detectors = DATA.detectors;

    /* Initialize detector scoring related variables */

    size_t *detector_offsets = NULL;
    size_t *detector_bin_counts = NULL;
    size_t total_detector_bins = 0;
    int max_threads = omp_get_max_threads();
    double *detector_thread_sums = NULL;
    double *detector_thread_sums_sq = NULL;
    double *detector_history_counts = NULL;

    /* Create thread local storage for detector scores */

    if (n_detectors > 0)
    {
        detector_offsets = (size_t*)malloc(n_detectors * sizeof(size_t));
        detector_bin_counts = (size_t*)calloc(n_detectors, sizeof(size_t));
        if (!detector_offsets || !detector_bin_counts)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            free(detector_offsets);
            free(detector_bin_counts);
            free(detector_thread_sums);
            free(detector_thread_sums_sq);
            free(detector_history_counts);
            return EXIT_FAILURE;
        }

        for (size_t d = 0; d < n_detectors; ++d)
        {
            detector_offsets[d] = total_detector_bins;
            Detector *det = detectors[d];
            size_t bins = 0;
            if (det)
            {
                if (det->type == DETECTOR_TYPE_REACTION_RATE)
                {
                    size_t energy_bins = det->data.reaction_rate.n_energy_bins;
                    if (energy_bins == 0)
                        energy_bins = 1;
                    bins = det->data.reaction_rate.n_nuclides * energy_bins * RRDET_MODE_COUNT;
                }
                else if (det->type == DETECTOR_TYPE_ENERGY_SPECTRUM && det->data.energy_spectrum.grid.enabled)
                    bins = det->data.energy_spectrum.grid.n_bins;
            }
            detector_bin_counts[d] = bins;
            total_detector_bins += bins;
        }

        if (total_detector_bins > 0)
        {
            size_t per_thread_bytes = (size_t)max_threads * total_detector_bins;
            detector_thread_sums = (double*)calloc(per_thread_bytes, sizeof(double));
            detector_thread_sums_sq = (double*)calloc(per_thread_bytes, sizeof(double));
            detector_history_counts = (double*)calloc(per_thread_bytes, sizeof(double));
            if (!detector_thread_sums || !detector_thread_sums_sq || !detector_history_counts)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                free(detector_offsets);
                free(detector_bin_counts);
                free(detector_thread_sums);
                free(detector_thread_sums_sq);
                free(detector_history_counts);
                return EXIT_FAILURE;
            }
        }
    }

    /* Loop over generations (single-threaded) */

    for (long g = 1; g <= GLOB.n_generations + GLOB.n_inactive; g++) 
    {
        DATA.generation = (uint64_t)g;

        /* Initialize generation score struct to reduce during parallel run */

        GenerationScores gen_scores;
        memset(&gen_scores, 0, sizeof(GenerationScores));

        /* Clear thread-local detector scoring arrays */

        const int do_detectors = (total_detector_bins > 0) && (g > GLOB.n_inactive);
        if (do_detectors)
        {
            size_t total_bins_bytes = (size_t)max_threads * total_detector_bins * sizeof(double);
            memset(detector_thread_sums, 0, total_bins_bytes);
            memset(detector_thread_sums_sq, 0, total_bins_bytes);
        }

        /* On all but the first run, build a fission neutron bank from fission sites of last generation */
        
        if (g > 1) 
        {
            if (buildFissionBank() != 0) 
            {
                fprintf(stderr, "[ERROR] Failed to build fission bank for generation %zu.\n", g);
                free(detector_offsets);
                free(detector_bin_counts);
                free(detector_thread_sums);
                free(detector_thread_sums_sq);
                free(detector_history_counts);
                return EXIT_FAILURE;
            }
        }

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
                shared(GLOB, DATA, RES, stdout, stderr, \
                       detectors, n_detectors, detector_offsets, detector_bin_counts, \
                       detector_thread_sums, detector_thread_sums_sq, \
                       detector_history_counts, total_detector_bins, do_detectors) \
                reduction(+:gen_scores)
        {
            ReactionRateMode r_mode = RRDET_MODE_COUNT;
            const int tid = omp_get_thread_num();
            double *history_counts = NULL;
            double *local_sum = NULL;
            double *local_sq = NULL;

            if (do_detectors)
            {
                history_counts = detector_history_counts + (size_t)tid * total_detector_bins;
                local_sum = detector_thread_sums + (size_t)tid * total_detector_bins;
                local_sq  = detector_thread_sums_sq + (size_t)tid * total_detector_bins;
            }

            #pragma omp for schedule(dynamic)
            for (size_t i = 0; i < DATA.n_bank; i++)
            {
                /* Get neutron from bank */

                Neutron *n = &DATA.bank[i];
                if (n->status != NEUTRON_ALIVE) 
                    continue;

                gen_scores.n_histories++;

                if (do_detectors)
                    memset(history_counts, 0, total_detector_bins * sizeof(double));

                /* Sample collisions until dead */

                size_t k = 0;
                while (n->status == NEUTRON_ALIVE) 
                {
                    double segment_energy = n->E;
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

                    /* Update path length */

                    n->path_length += d;
                    gen_scores.total_path_length += d;

                    if (n->E > E_THERMAL)
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

                    n->mat_idx = getMaterialAtPosition(n->x, n->y, n->z);

                    if (do_detectors && d > 0.0)
                    {
                        if (n->mat_idx >= 0)
                        {
                            for (size_t det_idx = 0; det_idx < n_detectors; ++det_idx)
                            {
                                Detector *detector = detectors[det_idx];
                                if (!detector || detector->type != DETECTOR_TYPE_ENERGY_SPECTRUM)
                                    continue;

                                size_t bin_count = detector_bin_counts[det_idx];
                                if (bin_count == 0)
                                    continue;

                                EnergySpectrumDetector *es = &detector->data.energy_spectrum;
                                if (!es->grid.enabled)
                                    continue;
                                if (es->has_material_filter && es->material_index != n->mat_idx)
                                    continue;

                                size_t energy_bin = energyGridBinIndex(&es->grid, segment_energy);
                                if (energy_bin >= es->grid.n_bins)
                                    continue;

                                size_t offset = detector_offsets[det_idx];
                                history_counts[offset + energy_bin] += d;
                            }
                        }
                    }

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

                        r_mode = RRDET_MODE_ELASTIC;
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

                        r_mode = RRDET_MODE_FISSION;

                        /* Score fission time */
                        size_t time_bin = (size_t)(n->time / TIME_BIN_WIDTH);
                        if (time_bin < MAX_TIME_BINS && n->fission_yield > 0)
                        {
                            gen_scores.fission_time_yield[time_bin] += (double)n->fission_yield;
                            gen_scores.fission_time_events[time_bin] += 1;
                        }

                    }
                    else if (MT_IS_INELASTIC_SCATTER(mt))
                    {
                        gen_scores.total_inelastic_scatters++;

                        handleInelasticScatter(n, &DATA.mats[n->mat_idx].nucs[nuc_idx].nuc_data, mt);

                        r_mode = RRDET_MODE_INELASTIC;
                    }
                    else if (MT_IS_CAPTURE(mt))
                    {
                        gen_scores.total_captures++;
                        n->status = NEUTRON_DEAD_CAPTURE;

                        r_mode = RRDET_MODE_CAPTURE;
                    }
                    else
                    {
                        gen_scores.total_unknowns++;
                    }

                    /* Score detectors if applicable */

                    if (do_detectors && r_mode != RRDET_MODE_COUNT)
                    {
                        for (size_t d = 0; d < n_detectors; ++d)
                        {
                            Detector *det = detectors[d];
                            if (!det || det->type != DETECTOR_TYPE_REACTION_RATE)
                                continue;

                            ReactionRateDetector *rr = &det->data.reaction_rate;
                            if (rr->material_index == n->mat_idx && (size_t)nuc_idx < rr->n_nuclides)
                            {
                                size_t energy_bins = (rr->n_energy_bins > 0) ? rr->n_energy_bins : 1;
                                size_t energy_bin = 0;
                                if (rr->energy_grid.enabled && rr->n_energy_bins > 0)
                                {
                                    energy_bin = energyGridBinIndex(&rr->energy_grid, segment_energy);
                                    if (energy_bin >= energy_bins)
                                        continue;
                                }
                                size_t flatten_stride = energy_bins * RRDET_MODE_COUNT;
                                size_t bin = detector_offsets[d] + (size_t)nuc_idx * flatten_stride + energy_bin * RRDET_MODE_COUNT + (size_t)r_mode;
                                history_counts[bin] += 1.0;
                            }
                        }
                    }

                    /* Check if cut-off has been reached */
                    if (n->E < GLOB.energy_cutoff)
                    {
                        n->status = NEUTRON_DEAD_TERMINATED;
                    }

                }

                if (do_detectors)
                {
                    for (size_t b = 0; b < total_detector_bins; ++b)
                    {
                        double val = history_counts[b];
                        local_sum[b] += val;
                        local_sq[b]  += val * val;
                    }
                }
            }
        }
        /* ###################################################################################### */
        /* Store generation scores if in active cycle */

        if (g > GLOB.n_inactive)
        {
            if (g - 1 - GLOB.n_inactive < RES.n_generations && RES.avg_scores)
                RES.avg_scores[g - 1 - GLOB.n_inactive] = gen_scores;
        }
        
        /* Score detectors from thread private buffers if in active cycle */

        if (g > GLOB.n_inactive && n_detectors > 0)
        {
            uint64_t histories_this_gen = gen_scores.n_histories;
            for (size_t d = 0; d < n_detectors; ++d)
                detectors[d]->n_histories += histories_this_gen;
        }

        if (g > GLOB.n_inactive && total_detector_bins > 0)
        {
            for (size_t d = 0; d < n_detectors; ++d)
            {
                Detector *det = detectors[d];
                size_t bins = detector_bin_counts ? detector_bin_counts[d] : 0;
                if (!det || bins == 0)
                    continue;

                size_t offset = detector_offsets[d];
                if (det->type == DETECTOR_TYPE_REACTION_RATE)
                {
                    ReactionRateDetector *rr = &det->data.reaction_rate;
                    size_t energy_bins = (rr->n_energy_bins > 0) ? rr->n_energy_bins : 1;
                    size_t flatten_stride = energy_bins * RRDET_MODE_COUNT;
                    for (size_t nuc = 0; nuc < rr->n_nuclides; ++nuc)
                    {
                        for (size_t ebin = 0; ebin < energy_bins; ++ebin)
                        {
                            for (size_t mode = 0; mode < RRDET_MODE_COUNT; ++mode)
                            {
                                size_t bin_index = offset + nuc * flatten_stride + ebin * RRDET_MODE_COUNT + mode;
                                double sum_val = 0.0;
                                double sum_sq_val = 0.0;
                                for (int tid = 0; tid < max_threads; ++tid)
                                {
                                    size_t base = (size_t)tid * total_detector_bins + bin_index;
                                    sum_val    += detector_thread_sums[base];
                                    sum_sq_val += detector_thread_sums_sq[base];
                                }

                                size_t tally_index = ebin * RRDET_MODE_COUNT + mode;
                                ReactionRateTally *tally = &rr->nuclides[nuc].tallies[tally_index];
                                tally->sum    += sum_val;
                                tally->sum_sq += sum_sq_val;
                            }
                        }
                    }
                }
                else if (det->type == DETECTOR_TYPE_ENERGY_SPECTRUM)
                {
                    EnergySpectrumDetector *es = &det->data.energy_spectrum;
                    for (size_t bin = 0; bin < es->grid.n_bins; ++bin)
                    {
                        size_t bin_index = offset + bin;
                        double sum_val = 0.0;
                        double sum_sq_val = 0.0;
                        for (int tid = 0; tid < max_threads; ++tid)
                        {
                            size_t base = (size_t)tid * total_detector_bins + bin_index;
                            sum_val    += detector_thread_sums[base];
                            sum_sq_val += detector_thread_sums_sq[base];
                        }

                        es->track_length_sum[bin]    += sum_val;
                        es->track_length_sum_sq[bin] += sum_sq_val;
                    }
                }
            }
        }

        /* Print generation summary */
        
        double k_eff = (gen_scores.n_histories > 0) ? ((double)gen_scores.total_fission_yield / (double)gen_scores.n_histories) : 0.0;
        if (g > GLOB.n_inactive)
            fprintf(stdout, "Generation %ld k-eff: %.6f\n", g - GLOB.n_inactive, k_eff);
        else
            fprintf(stdout, " keff: %.6lf\n", k_eff);
    }


    /* Return succesfully */

    free(detector_offsets);
    free(detector_bin_counts);
    free(detector_thread_sums);
    free(detector_thread_sums_sq);
    free(detector_history_counts);
    return EXIT_SUCCESS;
}

inline size_t energyGridBinIndex(const EnergyGrid *grid, double energy)
{
    if (!grid)
        return 0;

    if (!grid->enabled || !grid->edges || grid->n_bins == 0)
        return grid->n_bins;

    if (energy < grid->E_min)
        return grid->n_bins;

    if (energy > grid->E_max)
        return grid->n_bins;

    if (energy == grid->E_max)
        return (grid->n_bins > 0) ? grid->n_bins - 1 : grid->n_bins;

    const double *edges = grid->edges;
    size_t left = 0;
    size_t right = grid->n_bins;
    while (left < right)
    {
        size_t mid = (left + right) >> 1;
        double low = edges[mid];
        double high = edges[mid + 1];
        if (energy < low)
        {
            right = mid;
        }
        else if (energy >= high)
        {
            left = mid + 1;
        }
        else
        {
            return mid;
        }
    }

    return grid->n_bins;
}
