#include "header.h"

void processDetectorResults(void)
{
    /* Mode labels for stdout */

    static const char *mode_labels[RRDET_MODE_COUNT] = {
    "Elastic scatter",
    "Fission",
    "Inelastic scatter",
    "Capture"
    };

    const size_t n_detectors = DATA.n_detectors;

    if (n_detectors == 0)
    {
        if (VERBOSITY >= 1)
            fprintf(stdout, "\n[INFO] No detectors defined.\n");
        return;
    }

    /* Process detector results while writing to stdout and mfile */

    fprintf(stdout, "\nProcessing detector results...\n");

    char mpath[MAX_PATH];
    const char *inputf = GLOB.inputfname;
    snprintf(mpath, sizeof(mpath), "%s_det.m", inputf);

    FILE *mfile = fopen(mpath, "w");
    if (!mfile)
    {
        fprintf(stderr, "[ERROR] Could not open '%s' for detector output.\n", mpath);
    }
    else
    {
        fprintf(mfile, "%%%% Detector results auto-generated\n");
        fprintf(mfile, "%% Input file: %s\n\n", inputf);
    }

    for (size_t d = 0; d < n_detectors; ++d)
    {
        Detector *det = DATA.detectors[d];
        if (!det)
            continue;

        switch (det->type)
        {
            case DETECTOR_TYPE_REACTION_RATE:
            {
                ReactionRateDetector *rr = &det->data.reaction_rate;
                const uint64_t n_histories = det->n_histories;

                fprintf(stdout, "\nDetector %zu: %s (material: %s)\n", d + 1, det->name, rr->material_name);
                fprintf(stdout, "  Tallied histories: %lu\n", n_histories);

                if (n_histories == 0)
                {
                    fprintf(stdout, "  No active histories recorded for this detector.\n");
                    if (mfile)
                        fprintf(mfile, "DET_%zu_HISTORIES = 0;\n\n", d + 1);
                    break;
                }

                const size_t n_nuclides = rr->n_nuclides;
                const size_t n_modes = RRDET_MODE_COUNT;
                const size_t n_energy_bins = (rr->n_energy_bins > 0) ? rr->n_energy_bins : 1;

                if (n_nuclides == 0)
                {
                    fprintf(stdout, "  [WARNING] No nuclides linked to detector.\n");
                    if (mfile)
                        fprintf(mfile, "DET_%zu_HISTORIES = %lu; %% No nuclides linked\n\n", d + 1, n_histories);
                    break;
                }

                size_t aggregated_cells = n_nuclides * n_modes;
                size_t total_cells = n_nuclides * n_energy_bins * n_modes;
                size_t total_bin_cells = n_nuclides * n_energy_bins;

                double *mode_sum = (double*)calloc(aggregated_cells, sizeof(double));
                double *mode_sum_sq = (double*)calloc(aggregated_cells, sizeof(double));
                double *mode_mean = (double*)calloc(aggregated_cells, sizeof(double));
                double *mode_ci = (double*)calloc(aggregated_cells, sizeof(double));
                double *mode_frac = (double*)calloc(aggregated_cells, sizeof(double));
                double *mode_frac_ci = (double*)calloc(aggregated_cells, sizeof(double));
                double *nuclide_sum = (double*)calloc(n_nuclides, sizeof(double));
                double *nuclide_sum_sq = (double*)calloc(n_nuclides, sizeof(double));
                double *nuclide_frac = (double*)calloc(n_nuclides, sizeof(double));
                double *nuclide_frac_ci = (double*)calloc(n_nuclides, sizeof(double));
                double *bin_sum = (double*)calloc(total_bin_cells, sizeof(double));
                double *bin_sum_sq = (double*)calloc(total_bin_cells, sizeof(double));
                double *bin_mean = (double*)calloc(total_bin_cells, sizeof(double));
                double *bin_ci = (double*)calloc(total_bin_cells, sizeof(double));
                double *energy_mode_sum = (double*)calloc(total_cells, sizeof(double));
                double *energy_mode_sum_sq = (double*)calloc(total_cells, sizeof(double));
                double *det_bin_sum = (double*)calloc(n_energy_bins, sizeof(double));
                if (!mode_sum || !mode_sum_sq || !mode_mean || !mode_ci ||
                    !mode_frac || !mode_frac_ci || !nuclide_sum || !nuclide_sum_sq ||
                    !nuclide_frac || !nuclide_frac_ci || !bin_sum || !bin_sum_sq ||
                    !bin_mean || !bin_ci || !energy_mode_sum || !energy_mode_sum_sq ||
                    !det_bin_sum)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                    free(mode_sum); free(mode_sum_sq); free(mode_mean); free(mode_ci);
                    free(mode_frac); free(mode_frac_ci);
                    free(nuclide_sum); free(nuclide_sum_sq);
                    free(nuclide_frac); free(nuclide_frac_ci);
                    free(bin_sum); free(bin_sum_sq); free(bin_mean); free(bin_ci);
                    free(energy_mode_sum); free(energy_mode_sum_sq);
                    free(det_bin_sum);
                    if (mfile)
                        fclose(mfile);
                    exit(EXIT_FAILURE);
                }

                double det_sum = 0.0;
                double det_sum_sq = 0.0;
                bool warned_bin_mismatch = false;

                for (size_t n = 0; n < n_nuclides; ++n)
                {
                    ReactionRateNuclide *dn = &rr->nuclides[n];
                    size_t bins_this = (dn->n_energy_bins > 0) ? dn->n_energy_bins : n_energy_bins;
                    if (bins_this != n_energy_bins && !warned_bin_mismatch)
                    {
                        fprintf(stderr, "[WARN] Detector %s: nuclide %s energy bin count mismatch (expected %zu, got %zu). Using minimum.\n",
                                det->name, dn->nuclide_name, n_energy_bins, bins_this);
                        warned_bin_mismatch = true;
                    }
                    size_t loop_bins = (bins_this < n_energy_bins) ? bins_this : n_energy_bins;

                    for (size_t e = 0; e < loop_bins; ++e)
                    {
                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t tally_idx = e * n_modes + m;
                            ReactionRateTally *tally = &dn->tallies[tally_idx];
                            double sum = tally->sum;
                            double sum_sq = tally->sum_sq;

                            size_t agg_idx = n * n_modes + m;
                            mode_sum[agg_idx] += sum;
                            mode_sum_sq[agg_idx] += sum_sq;

                            size_t bin_idx = n * n_energy_bins + e;
                            bin_sum[bin_idx] += sum;
                            bin_sum_sq[bin_idx] += sum_sq;

                            size_t energy_idx = ((n * n_energy_bins) + e) * n_modes + m;
                            energy_mode_sum[energy_idx] = sum;
                            energy_mode_sum_sq[energy_idx] = sum_sq;

                            nuclide_sum[n] += sum;
                            nuclide_sum_sq[n] += sum_sq;
                            det_sum += sum;
                            det_sum_sq += sum_sq;
                            det_bin_sum[e] += sum;
                        }
                    }
                }

                double det_mean = (n_histories > 0) ? det_sum / (double)n_histories : 0.0;
                double det_std = 0.0;
                if (n_histories > 1)
                {
                    double numerator = det_sum_sq - (det_sum * det_sum) / (double)n_histories;
                    if (numerator < 0.0)
                        numerator = 0.0;
                    det_std = sqrt(numerator / (double)(n_histories - 1));
                }
                double det_se = (n_histories > 1) ? det_std / sqrt((double)n_histories) : 0.0;

                for (size_t idx = 0; idx < total_bin_cells; ++idx)
                {
                    double sum = bin_sum[idx];
                    double sum_sq = bin_sum_sq[idx];
                    double mean = (n_histories > 0) ? sum / (double)n_histories : 0.0;
                    double ci = 0.0;
                    if (n_histories > 1)
                    {
                        double numerator = sum_sq - (sum * sum) / (double)n_histories;
                        if (numerator < 0.0)
                            numerator = 0.0;
                        double sigma = sqrt(numerator / (double)(n_histories - 1));
                        ci = CI95_FACTOR * sigma / sqrt((double)n_histories);
                    }
                    bin_mean[idx] = mean;
                    bin_ci[idx] = ci;
                }

                for (size_t n = 0; n < n_nuclides; ++n)
                {
                    double mean_total = (n_histories > 0) ? nuclide_sum[n] / (double)n_histories : 0.0;
                    double std_total = 0.0;
                    if (n_histories > 1)
                    {
                        double numerator = nuclide_sum_sq[n] - (nuclide_sum[n] * nuclide_sum[n]) / (double)n_histories;
                        if (numerator < 0.0)
                            numerator = 0.0;
                        std_total = sqrt(numerator / (double)(n_histories - 1));
                    }
                    double se_total = (n_histories > 1) ? std_total / sqrt((double)n_histories) : 0.0;
                    double frac_det = (det_sum > 0.0) ? nuclide_sum[n] / det_sum : 0.0;
                    double frac_det_ci = 0.0;
                    if (n_histories > 1 && det_mean > 0.0)
                    {
                        double se_frac = sqrt((se_total * se_total) / (det_mean * det_mean) +
                                              (frac_det * frac_det) * (det_se * det_se) / (det_mean * det_mean));
                        frac_det_ci = CI95_FACTOR * se_frac;
                    }

                    nuclide_frac[n] = frac_det;
                    nuclide_frac_ci[n] = frac_det_ci;

                    ReactionRateNuclide *dn = &rr->nuclides[n];
                    fprintf(stdout, "    idx %d: Nuclide %5s  frac(total) = %12.6e +/- %9.3e\n",
                            dn->nuclide_index, dn->nuclide_name, frac_det, frac_det_ci);

                    for (size_t m = 0; m < n_modes; ++m)
                    {
                        size_t agg_idx = n * n_modes + m;
                        double sum = mode_sum[agg_idx];
                        double sum_sq = mode_sum_sq[agg_idx];
                        double mean = (n_histories > 0) ? sum / (double)n_histories : 0.0;
                        double std = 0.0;
                        if (n_histories > 1)
                        {
                            double numerator = sum_sq - (sum * sum) / (double)n_histories;
                            if (numerator < 0.0)
                                numerator = 0.0;
                            std = sqrt(numerator / (double)(n_histories - 1));
                        }
                        double se = (n_histories > 1) ? std / sqrt((double)n_histories) : 0.0;
                        double ci = (n_histories > 1) ? CI95_FACTOR * se : 0.0;

                        mode_mean[agg_idx] = mean;
                        mode_ci[agg_idx] = ci;

                        double frac = (nuclide_sum[n] > 0.0) ? sum / nuclide_sum[n] : 0.0;
                        double frac_ci = 0.0;
                        if (n_histories > 1 && mean_total > 0.0)
                        {
                            double se_frac = sqrt((se * se) / (mean_total * mean_total) +
                                                  (frac * frac) * (se_total * se_total) / (mean_total * mean_total));
                            frac_ci = CI95_FACTOR * se_frac;
                        }
                        mode_frac[agg_idx] = frac;
                        mode_frac_ci[agg_idx] = frac_ci;

                        fprintf(stdout,
                                "      %-18s     frac(nucl) = %12.6e +/- %9.3e   mean = %12.6e +/- %9.3e\n",
                                mode_labels[m], frac, frac_ci, mean, ci);
                    }
                }

                const double *edges = (rr->energy_grid.enabled && rr->energy_grid.edges) ? rr->energy_grid.edges : NULL;
                fprintf(stdout, "\n      Energy-bin breakdown (bin means per history)\n");
                for (size_t e = 0; e < n_energy_bins; ++e)
                {
                    double e_low = (edges && rr->energy_grid.n_bins == n_energy_bins) ? edges[e] : 0.0;
                    double e_high = (edges && rr->energy_grid.n_bins == n_energy_bins) ? edges[e + 1] : 0.0;
                    if (edges && rr->energy_grid.n_bins == n_energy_bins)
                        fprintf(stdout, "      Bin %2zu: [%11.4e, %11.4e)\n", e, e_low, e_high);
                    else
                        fprintf(stdout, "      Bin %2zu: (no explicit energy bounds)\n", e);

                    double det_bin_mean = (n_histories > 0) ? det_bin_sum[e] / (double)n_histories : 0.0;
                    fprintf(stdout, "        Total detector mean in bin: %12.6e\n", det_bin_mean);

                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        ReactionRateNuclide *dn = &rr->nuclides[n];
                        size_t bin_idx = n * n_energy_bins + e;
                        double mean_total = bin_mean[bin_idx];
                        double ci_total = bin_ci[bin_idx];
                        fprintf(stdout, "        Nuclide %5s: mean = %12.6e +/- %9.3e\n",
                                dn->nuclide_name, mean_total, ci_total);

                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t energy_idx = ((n * n_energy_bins) + e) * n_modes + m;
                            double sum = energy_mode_sum[energy_idx];
                            double sum_sq = energy_mode_sum_sq[energy_idx];
                            double mean = (n_histories > 0) ? sum / (double)n_histories : 0.0;
                            double ci = 0.0;
                            if (n_histories > 1)
                            {
                                double numerator = sum_sq - (sum * sum) / (double)n_histories;
                                if (numerator < 0.0)
                                    numerator = 0.0;
                                double sigma = sqrt(numerator / (double)(n_histories - 1));
                                ci = CI95_FACTOR * sigma / sqrt((double)n_histories);
                            }
                            fprintf(stdout,
                                    "          %-18s mean = %12.6e +/- %9.3e\n",
                                    mode_labels[m], mean, ci);
                        }
                    }
                }

                if (mfile)
                {
                    fprintf(mfile, "%% Detector %zu: %s (material: %s)\n",
                            d + 1, det->name, rr->material_name);
                    fprintf(mfile, "DET_%zu_HISTORIES = %lu\n", d + 1, n_histories);

                    fprintf(mfile, "DET_%zu_NUCLIDE_IDX = [", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        ReactionRateNuclide *dn = &rr->nuclides[n];
                        fprintf(mfile, "%s%d", n == 0 ? "" : " ", dn->nuclide_index);
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_MEAN = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t idx = n * n_modes + m;
                            fprintf(mfile, "%s%.12e", m == 0 ? "  " : " ", mode_mean[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_CI95 = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t idx = n * n_modes + m;
                            fprintf(mfile, "%s%.12e", m == 0 ? "  " : " ", mode_ci[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_NUCLIDE_FRAC = [", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                        fprintf(mfile, "%s%.12e", n == 0 ? "" : " ", nuclide_frac[n]);
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_NUCLIDE_FRAC_CI = [", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                        fprintf(mfile, "%s%.12e", n == 0 ? "" : " ", nuclide_frac_ci[n]);
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_MODE_FRAC = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t idx = n * n_modes + m;
                            fprintf(mfile, "%s%.12e", m == 0 ? "  " : " ", mode_frac[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_MODE_FRAC_CI = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t m = 0; m < n_modes; ++m)
                        {
                            size_t idx = n * n_modes + m;
                            fprintf(mfile, "%s%.12e", m == 0 ? "  " : " ", mode_frac_ci[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    if (edges && rr->energy_grid.n_bins == n_energy_bins)
                    {
                        fprintf(mfile, "DET_%zu_ENERGY_BOUNDS = [\n", d + 1);
                        for (size_t e = 0; e < n_energy_bins; ++e)
                            fprintf(mfile, "  %.12e %.12e\n", edges[e], edges[e + 1]);
                        fprintf(mfile, "];\n");
                    }
                    else
                    {
                        fprintf(mfile, "DET_%zu_ENERGY_BOUNDS = [];\n", d + 1);
                    }

                    fprintf(mfile, "DET_%zu_BIN_MEAN = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t e = 0; e < n_energy_bins; ++e)
                        {
                            size_t idx = n * n_energy_bins + e;
                            fprintf(mfile, "%s%.12e", e == 0 ? "  " : " ", bin_mean[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_BIN_CI95 = [\n", d + 1);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t e = 0; e < n_energy_bins; ++e)
                        {
                            size_t idx = n * n_energy_bins + e;
                            fprintf(mfile, "%s%.12e", e == 0 ? "  " : " ", bin_ci[idx]);
                        }
                        fprintf(mfile, n + 1 < n_nuclides ? "\n" : "\n");
                    }
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_MODE_ENERGY_MEAN = zeros(%zu, %zu);\n", d + 1, n_nuclides, n_energy_bins * n_modes);
                    fprintf(mfile, "DET_%zu_MODE_ENERGY_CI95 = zeros(%zu, %zu);\n", d + 1, n_nuclides, n_energy_bins * n_modes);
                    for (size_t n = 0; n < n_nuclides; ++n)
                    {
                        for (size_t e = 0; e < n_energy_bins; ++e)
                        {
                            for (size_t m = 0; m < n_modes; ++m)
                            {
                                size_t energy_idx = ((n * n_energy_bins) + e) * n_modes + m;
                                double sum = energy_mode_sum[energy_idx];
                                double sum_sq = energy_mode_sum_sq[energy_idx];
                                double mean = (n_histories > 0) ? sum / (double)n_histories : 0.0;
                                double ci = 0.0;
                                if (n_histories > 1)
                                {
                                    double numerator = sum_sq - (sum * sum) / (double)n_histories;
                                    if (numerator < 0.0)
                                        numerator = 0.0;
                                    double sigma = sqrt(numerator / (double)(n_histories - 1));
                                    ci = CI95_FACTOR * sigma / sqrt((double)n_histories);
                                }
                                size_t col = e * n_modes + m;
                                fprintf(mfile, "DET_%zu_MODE_ENERGY_MEAN(%zu,%zu) = %.12e;\n",
                                        d + 1, n + 1, col + 1, mean);
                                fprintf(mfile, "DET_%zu_MODE_ENERGY_CI95(%zu,%zu) = %.12e;\n",
                                        d + 1, n + 1, col + 1, ci);
                            }
                        }
                    }
                    fprintf(mfile, "\n");
                }

                free(mode_sum);
                free(mode_sum_sq);
                free(mode_mean);
                free(mode_ci);
                free(mode_frac);
                free(mode_frac_ci);
                free(nuclide_sum);
                free(nuclide_sum_sq);
                free(nuclide_frac);
                free(nuclide_frac_ci);
                free(bin_sum);
                free(bin_sum_sq);
                free(bin_mean);
                free(bin_ci);
                free(energy_mode_sum);
                free(energy_mode_sum_sq);
                free(det_bin_sum);
                break;
            }
            case DETECTOR_TYPE_ENERGY_SPECTRUM:
            {
                EnergySpectrumDetector *es = &det->data.energy_spectrum;
                const uint64_t n_histories = det->n_histories;
                const size_t n_bins = es->grid.n_bins;
                const char *mat_desc = es->has_material_filter ? es->material_name : "all materials";

                fprintf(stdout, "\nDetector %zu: %s (energy spectrum, material: %s, bins: %zu)\n",
                        d + 1, det->name, mat_desc, n_bins);
                fprintf(stdout, "  Tallied histories: %lu\n", n_histories);

                if (n_histories == 0)
                {
                    fprintf(stdout, "  No active histories recorded for this detector.\n");
                    if (mfile)
                        fprintf(mfile, "DET_%zu_HISTORIES = 0;\n\n", d + 1);
                    break;
                }

                const double *edges = es->grid.edges;
                if (n_bins == 0 || !edges)
                {
                    fprintf(stdout, "  [WARNING] Energy grid missing for detector. Skipping.\n");
                    if (mfile)
                        fprintf(mfile, "DET_%zu_HISTORIES = %lu; %% Missing energy grid\n\n", d + 1, n_histories);
                    break;
                }

                double *bin_mean = (double*)calloc(n_bins, sizeof(double));
                double *bin_ci = (double*)calloc(n_bins, sizeof(double));
                if (!bin_mean || !bin_ci)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                    free(bin_mean);
                    free(bin_ci);
                    if (mfile)
                        fclose(mfile);
                    exit(EXIT_FAILURE);
                }

                fprintf(stdout, "  Bin edges (MeV) and track-length tallies (cm per history):\n");
                for (size_t b = 0; b < n_bins; ++b)
                {
                    double sum = es->track_length_sum ? es->track_length_sum[b] : 0.0;
                    double sum_sq = es->track_length_sum_sq ? es->track_length_sum_sq[b] : 0.0;
                    double mean = sum / (double)n_histories;
                    double ci = 0.0;
                    if (n_histories > 1)
                    {
                        double numerator = sum_sq - (sum * sum) / (double)n_histories;
                        if (numerator < 0.0)
                            numerator = 0.0;
                        double sigma = sqrt(numerator / (double)(n_histories - 1));
                        ci = CI95_FACTOR * sigma / sqrt((double)n_histories);
                    }
                    bin_mean[b] = mean;
                    bin_ci[b] = ci;

                    fprintf(stdout,
                            "    [%11.4e, %11.4e)  mean = %12.6e +/- %9.3e\n",
                            edges[b], edges[b + 1], mean, ci);
                }

                if (mfile)
                {
                    fprintf(mfile, "%% Detector %zu: %s (energy spectrum, material: %s)\n",
                            d + 1, det->name, mat_desc);
                    fprintf(mfile, "DET_%zu_HISTORIES = %lu;\n", d + 1, n_histories);

                    fprintf(mfile, "DET_%zu_BIN_EDGES = [", d + 1);
                    for (size_t b = 0; b <= n_bins; ++b)
                        fprintf(mfile, "%s%.12e", b == 0 ? "" : " ", edges[b]);
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_TRACK_MEAN = [", d + 1);
                    for (size_t b = 0; b < n_bins; ++b)
                        fprintf(mfile, "%s%.12e", b == 0 ? "" : " ", bin_mean[b]);
                    fprintf(mfile, "];\n");

                    fprintf(mfile, "DET_%zu_TRACK_CI95 = [", d + 1);
                    for (size_t b = 0; b < n_bins; ++b)
                        fprintf(mfile, "%s%.12e", b == 0 ? "" : " ", bin_ci[b]);
                    fprintf(mfile, "];\n\n");
                }

                free(bin_mean);
                free(bin_ci);
                break;
            }
            default:
            {
                if (VERBOSITY >= 1)
                    fprintf(stdout, "\nDetector %zu: %s (unsupported type). Skipping.\n", d + 1, det->name);
                break;
            }
        }
    }

    if (mfile)
    {
        fclose(mfile);
        fprintf(stdout, "\nWrote detector results to '%s'.\n", mpath);
    }
}
