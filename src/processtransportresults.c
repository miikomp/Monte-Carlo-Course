#include <inttypes.h>

#include "header.h"

typedef struct {
    const char *perhist_label;
    const char *perhist_name;
    const char *rate_label;
    const char *rate_name;
    double (*extract)(const TransportRunScores *gs);
    bool    apply_norm;  // true if metric should be scaled by normalization factor
} MetricSpec;

static double safe_norm(uint64_t count, double value) 
{
    return (count > 0) ? value / (double)count : 0.0;
}

static inline uint64_t active_histories(const TransportRunScores *gs)
{
    if (!gs)
        return 0;
    if (gs->n_histories <= gs->total_terminated)
        return 0;
    return gs->n_histories - gs->total_terminated;
}

static double safe_norm_active(const TransportRunScores *gs, double value)
{
    return safe_norm(active_histories(gs), value);
}

static double metric_path_length(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, gs->total_path_length);
}

static double metric_fast_path_length(const TransportRunScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_fast_path_length);
}

static double metric_time(const TransportRunScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_time);
}

static double metric_time_fast(const TransportRunScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_slowing_down_time);
}

static double metric_fission_yield(const TransportRunScores *gs) 
{
    return safe_norm_active(gs, (double)gs->total_fission_yield);
}

static double metric_collisions(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_collisions);
}

static double metric_captures(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_captures);
}

static double metric_elastic(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_elastic_scatters);
}

static double metric_inelastic(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_inelastic_scatters);
}

static double metric_fissions(const TransportRunScores *gs) 
{
    return safe_norm_active(gs, (double)gs->total_fissions);
}

static double metric_fast_fissions(const TransportRunScores *gs) 
{
    return safe_norm_active(gs, (double)gs->total_fast_fissions);
}

static double metric_thermal_fissions(const TransportRunScores *gs) 
{
    return safe_norm_active(gs, (double)gs->total_thermal_fissions);
}

static double metric_leakages(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_leakages);
}

static double metric_unknowns(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_unknowns);
}

static double metric_keff(const TransportRunScores *gs)
{
    double den = (double)active_histories(gs);
    return (den > 0.0) ? (double)gs->total_fission_yield / den : 0.0;
}
static double metric_terminated(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_terminated);
}

static void compute_stats(const double *values, size_t n, double *mean, double *std, double *ci) 
{
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i)
        sum += values[i];
    double mu = (n > 0) ? sum / (double)n : 0.0;

    double var = 0.0;
    for (size_t i = 0; i < n; ++i) 
    {
        double d = values[i] - mu;
        var += d * d;
    }
    double sigma = (n > 1) ? sqrt(var / (double)(n - 1)) : 0.0;
    double ci_hw = (n > 1) ? CI95_FACTOR * sigma / sqrt((double)n) : 0.0;

    if (mean) *mean = mu;
    if (std)  *std  = sigma;
    if (ci)   *ci   = ci_hw;
}

void processTransportResults(void) 
{
    if (RES.n_iterations <= 0 || !RES.avg_scores) 
    {
        fprintf(stdout, "\n[INFO] No transport results to process.\n");
        return;
    }

    const size_t n_iter = (size_t)RES.n_iterations;

    const MetricSpec metrics[] = {
        {"Path length / history",        "PATH_LENGTH",        NULL,                 NULL,                 metric_path_length,       false},
        {"Fast path length / history",   "FAST_PATH_LENGTH",   NULL,                 NULL,                 metric_fast_path_length,  false},
        {"Flight time / history",        "FLIGHT_TIME",        NULL,                 NULL,                 metric_time,              false},
        {"Fast flight time / history",   "FAST_FLIGHT_TIME",   NULL,                 NULL,                 metric_time_fast,         false},
        {"Fission yield / history",      "FISSION_YIELD",      "Fission rate [1/s]",      "FISSION_RATE",      metric_fission_yield,     true},
        {"Collisions / history",         "COLLISIONS",         "Collision rate [1/s]",    "COLLISION_RATE",    metric_collisions,        true},
        {"Captures / history",           "CAPTURES",           "Capture rate [1/s]",      "CAPTURE_RATE",      metric_captures,          true},
        {"Elastic scatters / history",   "ELASTIC_SCATTERS",   "Elastic rate [1/s]",      "ELASTIC_RATE",      metric_elastic,           true},
        {"Inelastic scatters / history", "INELASTIC_SCATTERS", "Inelastic rate [1/s]",    "INELASTIC_RATE",    metric_inelastic,         true},
        {"Fissions / history",           "FISSIONS",           "Fission rate [1/s]",      "TOTAL_FISSION_RATE",metric_fissions,          true},
        {"Thermal fissions / history",   "THERMAL_FISSIONS",   "Thermal fission rate [1/s]", "THERMAL_FISSION_RATE", metric_thermal_fissions,  true},
        {"Fast fissions / history",      "FAST_FISSIONS",      "Fast fission rate [1/s]", "FAST_FISSION_RATE", metric_fast_fissions,     true},
        {"Leakages / history",           "LEAKAGES",           "Leakage rate [1/s]",      "LEAKAGE_RATE",      metric_leakages,          true},
        {"Terminated / history",         "TERMINATED",         "Termination rate [1/s]",  "TERMINATION_RATE",  metric_terminated,        true},
        {"Unknown outcomes / history",   "UNKNOWN_OUTCOMES",   "Unknown rate [1/s]",      "UNKNOWN_RATE",      metric_unknowns,          true},
        {"k-eff",                          "KEFF",              NULL,                       NULL,                metric_keff,              false}
    };

    const size_t nmetrics = sizeof(metrics) / sizeof(metrics[0]);

    double *histories = (double*)calloc(n_iter, sizeof(double));
    if (!histories) 
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return;
    }

    double **metric_values = (double**)calloc(nmetrics, sizeof(double*));
    double *summary_mean = (double*)calloc(nmetrics, sizeof(double));
    double *summary_ci = (double*)calloc(nmetrics, sizeof(double));
    if (!metric_values || !summary_mean || !summary_ci) 
    {
        free(histories);
        if (metric_values) free(metric_values);
        if (summary_mean) free(summary_mean);
        if (summary_ci) free(summary_ci);
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return;
    }
    for (size_t m = 0; m < nmetrics; ++m) 
    {
        metric_values[m] = (double*)calloc(n_iter, sizeof(double));
        if (!metric_values[m]) 
        {
            for (size_t j = 0; j < m; ++j) free(metric_values[j]);
            for (size_t j = 0; j < m; ++j)
                free(metric_values[j]);
            free(metric_values);
            free(histories);
            free(summary_mean);
            free(summary_ci);
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            return;
        }
    }

    for (size_t g = 0; g < n_iter; ++g)
    {
        const TransportRunScores *gs = &RES.avg_scores[g];
        histories[g] = (double)gs->n_histories;

        for (size_t m = 0; m < nmetrics; ++m)
            metric_values[m][g] = metrics[m].extract(gs);
    }

    if (VERBOSITY >= 1)
    {
        fprintf(stdout, "\nTransport per-generation statistics:\n");
        fprintf(stdout, "  %s  Histories\n", (GLOB.mode == RUNMODE_EXTERNAL_SOURCE) ? "Cycle" : "Generation");
        for (size_t g = 0; g < n_iter; ++g)
            fprintf(stdout, "  %-10zu  %-10.0f\n", g + 1, histories[g]);
    }

    fprintf(stdout, "\nTransport per-history statistics (95%% Confidence-interval):\n");
    for (size_t m = 0; m < nmetrics; ++m) 
    {
        double mean, std, ci;
        compute_stats(metric_values[m], n_iter, &mean, &std, &ci);
        summary_mean[m] = mean;
        summary_ci[m] = ci;
        double rel = (mean != 0.0) ? fabs(ci / mean) * 100.0 : 0.0;
        fprintf(stdout, "  %-28s  mean = %12.6e +/- %9.3e (%.6lf%%)\n",
                metrics[m].perhist_label, mean, ci, rel);
    }

    const double norm_factor = getNormalizationFactor();
    const char *norm_desc = getNormalizationModeDescription();

    fprintf(stdout, "\nTotal reaction rates (normalized to %.3e %s):\n",
            norm_factor, norm_desc);
    for (size_t m = 0; m < nmetrics; ++m)
    {
        if (!metrics[m].apply_norm)
            continue;

        double scaled_mean = summary_mean[m] * norm_factor;
        double scaled_ci = summary_ci[m] * norm_factor;
        fprintf(stdout, "  %-28s  %12.6e +/- %9.3e\n",
                metrics[m].rate_label ? metrics[m].rate_label : metrics[m].perhist_label,
                scaled_mean, scaled_ci);
    }

    /* ########################################################################################## */
    /* Write results to MATLAB-compatible .m file */

    char mpath[MAX_PATH];
    snprintf(mpath, sizeof(mpath), "%s_res.m", GLOB.inputfname);

    FILE *mfile = fopen(mpath, "w");
    if (!mfile) 
    {
        fprintf(stderr, "[ERROR] Could not open '%s' for writing.\n", mpath);
    } 
    else 
    {
        fprintf(mfile, "%%%% Transport results auto-generated\n");
        fprintf(mfile, "%% Input file: %s\n", GLOB.inputf);
        fprintf(mfile, "%% Normalisation: %.8E (%s)\n\n", norm_factor, norm_desc);
        fprintf(mfile, "N_GENERATIONS = %zu;\n\n", n_iter);

        fprintf(mfile, "%% Per-history averages\n");
        for (size_t m = 0; m < nmetrics; ++m) 
        {
            fprintf(mfile, "%-30s = [ %.8E %.8E ];\n",
                    metrics[m].perhist_name,
                    summary_mean[m], summary_ci[m]);
        }
        fprintf(mfile, "\n%% Reaction rates (normalized)\n");
        for (size_t m = 0; m < nmetrics; ++m)
        {
            if (!metrics[m].apply_norm)
                continue;

            const char *rate_name = metrics[m].rate_name ? metrics[m].rate_name
                                                         : metrics[m].perhist_name;
            fprintf(mfile, "%-30s = [ %.8E %.8E ];\n",
                    rate_name,
                    summary_mean[m] * norm_factor,
                    summary_ci[m] * norm_factor);
        }
        fprintf(mfile, "\n");

        fclose(mfile);
        fprintf(stdout, "\nWrote transport results to '%s'.\n", mpath);
    }

    for (size_t m = 0; m < nmetrics; ++m)
        free(metric_values[m]);
    free(metric_values);
    free(histories);
    free(summary_mean);
    free(summary_ci);
}
