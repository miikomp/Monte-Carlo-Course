#include "header.h"

static inline void compute_statistics(const ReactionRateTally *tally, uint64_t n_histories, double *mean, double *std, double *ci95);

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
    const char *fname = GLOB.fname ? GLOB.fname : "results";
    snprintf(mpath, sizeof(mpath), "%s_det.m", fname);

    FILE *mfile = fopen(mpath, "w");
    if (!mfile)
    {
        fprintf(stderr, "[ERROR] Could not open '%s' for detector output.\n", mpath);
    }
    else
    {
        fprintf(mfile, "%%%% Detector results auto-generated\n");
        fprintf(mfile, "%% Input file: %s\n\n", fname);
    }

    for (size_t d = 0; d < n_detectors; ++d)
    {
        ReactionRateDetector *det = DATA.detectors[d];
        const uint64_t n_histories = det->n_histories;

        fprintf(stdout, "\nDetector %zu: %s (material: %s)\n", d + 1, det->name, det->material_name);
        fprintf(stdout, "  Tallied histories: %lu\n", n_histories);

        if (n_histories == 0)
        {
            fprintf(stdout, "  No active histories recorded for this detector.\n");
            if (mfile)
                fprintf(mfile, "DET_%zu_HISTORIES = 0;\n\n", d + 1);
            continue;
        }

        const size_t n_nuclides = det->n_nuclides;
        const size_t n_modes = RRDET_MODE_COUNT;

        double *mode_sum = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *mode_sum_sq = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *mode_mean = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *mode_ci = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *mode_frac = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *mode_frac_ci = (double*)calloc(n_nuclides * n_modes, sizeof(double));
        double *nuclide_sum = (double*)calloc(n_nuclides, sizeof(double));
        double *nuclide_sum_sq = (double*)calloc(n_nuclides, sizeof(double));
        double *nuclide_frac = (double*)calloc(n_nuclides, sizeof(double));
        double *nuclide_frac_ci = (double*)calloc(n_nuclides, sizeof(double));
        if (!mode_sum || !mode_sum_sq || !mode_mean || !mode_ci ||
            !mode_frac || !mode_frac_ci || !nuclide_sum || !nuclide_sum_sq ||
            !nuclide_frac || !nuclide_frac_ci)
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            free(mode_sum); free(mode_sum_sq); free(mode_mean); free(mode_ci);
            free(mode_frac); free(mode_frac_ci);
            free(nuclide_sum); free(nuclide_sum_sq);
            free(nuclide_frac); free(nuclide_frac_ci);
            if (mfile)
                fclose(mfile);
            exit(EXIT_FAILURE);
        }

        double det_sum = 0.0;
        double det_sum_sq = 0.0;

        for (size_t n = 0; n < n_nuclides; ++n)
        {
            ReactionRateNuclide *dn = &det->nuclides[n];
            for (size_t m = 0; m < n_modes; ++m)
            {
                size_t idx = n * n_modes + m;
                mode_sum[idx] = dn->tallies[m].sum;
                mode_sum_sq[idx] = dn->tallies[m].sum_sq;
                nuclide_sum[n] += mode_sum[idx];
                nuclide_sum_sq[n] += mode_sum_sq[idx];
            }
            det_sum += nuclide_sum[n];
            det_sum_sq += nuclide_sum_sq[n];
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

        for (size_t n = 0; n < n_nuclides; ++n)
        {
            ReactionRateNuclide *dn = &det->nuclides[n];
            fprintf(stdout, "    idx %d: Nuclide %5s", dn->nuclide_index, dn->nuclide_name);

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

            fprintf(stdout, "  frac = %12.6e +/- %9.3e\n", frac_det, frac_det_ci);

            for (size_t m = 0; m < n_modes; ++m)
            {
                size_t idx = n * n_modes + m;
                double mean = (n_histories > 0) ? mode_sum[idx] / (double)n_histories : 0.0;
                double std = 0.0;
                if (n_histories > 1)
                {
                    double numerator = mode_sum_sq[idx] - (mode_sum[idx] * mode_sum[idx]) / (double)n_histories;
                    if (numerator < 0.0)
                        numerator = 0.0;
                    std = sqrt(numerator / (double)(n_histories - 1));
                }
                double se = (n_histories > 1) ? std / sqrt((double)n_histories) : 0.0;
                double ci = (n_histories > 1) ? CI95_FACTOR * se : 0.0;

                mode_mean[idx] = mean;
                mode_ci[idx] = ci;

                double frac = (nuclide_sum[n] > 0.0) ? mode_sum[idx] / nuclide_sum[n] : 0.0;
                double frac_ci = 0.0;
                if (n_histories > 1 && mean_total > 0.0)
                {
                    double se_frac = sqrt((se * se) / (mean_total * mean_total) +
                                          (frac * frac) * (se_total * se_total) / (mean_total * mean_total));
                    frac_ci = CI95_FACTOR * se_frac;
                }
                mode_frac[idx] = frac;
                mode_frac_ci[idx] = frac_ci;

                fprintf(stdout,
                        "      %-18s     frac = %12.6e +/- %9.3e   mean = %12.6e +/- %9.3e\n",
                        mode_labels[m], frac, frac_ci, mean, ci);
            }
        }

        if (mfile)
        {
            fprintf(mfile, "%% Detector %zu: %s (material: %s)\n",
                    d + 1, det->name, det->material_name);
            fprintf(mfile, "DET_%zu_HISTORIES = %lu\n", d + 1, n_histories);

            fprintf(mfile, "DET_%zu_NUCLIDE_IDX = [", d + 1);
            for (size_t n = 0; n < n_nuclides; ++n)
            {
                ReactionRateNuclide *dn = &det->nuclides[n];
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
            fprintf(mfile, "];\n\n");
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
    }

    if (mfile)
    {
        fclose(mfile);
        fprintf(stdout, "\nWrote detector results to '%s'.\n", mpath);
    }
}

void compute_statistics(const ReactionRateTally *tally,
                                      uint64_t n_histories,
                                      double *mean,
                                      double *std,
                                      double *ci95)
{
    double mu = 0.0;
    double sigma = 0.0;
    double ci = 0.0;

    if (n_histories > 0)
        mu = tally->sum / (double)n_histories;

    if (n_histories > 1)
    {
        double numerator = tally->sum_sq - (tally->sum * tally->sum) / (double)n_histories;
        if (numerator < 0.0)
            numerator = 0.0;
        sigma = sqrt(numerator / (double)(n_histories - 1));
        ci = CI95_FACTOR * sigma / sqrt((double)n_histories);
    }

    if (mean) *mean = mu;
    if (std)  *std  = sigma;
    if (ci95) *ci95 = ci;
}
