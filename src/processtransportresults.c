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

static double metric_flux(const TransportRunScores *gs)
{
    /* Track-length estimator of total flux (path length per history) */
    return metric_path_length(gs);
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

static double metric_ana_keff(const TransportRunScores *gs)
{
    double den = (double)active_histories(gs);
    return (den > 0.0) ? (double)gs->total_fission_yield / den : 0.0;
}
static double metric_imp_keff(const TransportRunScores *gs)
{
    double den = gs->total_collisions - (gs->total_elastic_scatters + gs->total_inelastic_scatters) + gs->total_leakages;
    return (den > 0.0) ? gs->total_fission_yield / den : 0.0;
}
static double metric_terminated(const TransportRunScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_terminated);
}

static double metric_dt_usefrac(const TransportRunScores *gs)
{
    if (!gs || gs->total_count == 0)
        return 0.0;
    return (double)gs->dt_count / (double)gs->total_count;
}

static double metric_surf_usefrac(const TransportRunScores *gs)
{
    double dt_frac = metric_dt_usefrac(gs);
    double surf_frac = 1.0 - dt_frac;
    return (surf_frac < 0.0) ? 0.0 : surf_frac;
}

static double metric_dt_eff(const TransportRunScores *gs)
{
    if (!gs || gs->dt_count == 0)
        return 0.0;
    double virt = (double)gs->dt_virtual_count;
    double attempts = (double)gs->dt_count;
    if (virt > attempts)
        return 0.0;
    return 1.0 - virt / attempts;
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
    size_t n_hist = 0;

    const MetricSpec metrics[] = {
        {"Path length / history",        "PATH_LENGTH",        NULL,                 NULL,                       metric_path_length,       false},
        {"Fast path length / history",   "FAST_PATH_LENGTH",   NULL,                 NULL,                       metric_fast_path_length,  false},
        {"Flight time / history",        "FLIGHT_TIME",        NULL,                 NULL,                       metric_time,              false},
        {"Fast flight time / history",   "FAST_FLIGHT_TIME",   NULL,                 NULL,                       metric_time_fast,         false},
        {"Collisions / history",         "COLLISIONS",         "Total reaction rate", "COLLISION_RATE",    metric_collisions,        true},
        {"Fission prod. / history",      "FISSION_PROD",       "Production rate",     "PRODUCTION_RATE",      metric_fission_yield,     true},
        {"Captures / history",           "CAPTURES",           "Capture rate",        "CAPTURE_RATE",      metric_captures,          true},
        {"Elastic scatters / history",   "ELASTIC_SCATTERS",   "Elastic scat. rate",        "ELASTIC_RATE",      metric_elastic,           true},
        {"Inelastic scatters / history", "INELASTIC_SCATTERS", "Inelastic scat. rate",      "INELASTIC_RATE",    metric_inelastic,         true},
        {"Fissions / history",           "FISSIONS",           "Fission rate",        "TOTAL_FISSION_RATE",metric_fissions,          true},
        {"Thermal fissions / history",   "THERMAL_FISSIONS",   "Thermal fission rate","THERMAL_FISSION_RATE", metric_thermal_fissions,  true},
        {"Fast fissions / history",      "FAST_FISSIONS",      "Fast fission rate",   "FAST_FISSION_RATE", metric_fast_fissions,     true},
        {"Flux / history",               "FLUX",               "Total flux",                "TOTAL_FLUX",        metric_flux,              true},
        {"Leakages / history",           "LEAKAGES",           "Leakage rate",        "LEAKAGE_RATE",      metric_leakages,          true},
        {"Terminated / history",         "TERMINATED",         "Termination rate",    "TERMINATION_RATE",  metric_terminated,        true},
        {"Unknown outcomes / history",   "UNKNOWN_OUTCOMES",   "Unknown rate",        "UNKNOWN_RATE",      metric_unknowns,          true},
        {"k-eff analog",                 "ANA_KEFF",           "k-eff analog",        "ANA_KEFF",               metric_ana_keff,          false},
        {"k-eff implicit",               "IMP_KEFF",           "k-eff implicit",      "IMP_KEFF",           metric_imp_keff,            false},
        {"Delta tracking use frac",      "DT_USEFRAC",         NULL,                  NULL,                 metric_dt_usefrac,         false},
        {"Surface tracking use frac",    "SURF_USEFRAC",       NULL,                  NULL,                 metric_surf_usefrac,       false},
        {"Delta tracking efficiency",    "DT_EFF",             NULL,                  NULL,                 metric_dt_eff,             false}
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
        n_hist += gs->n_histories;
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

    fprintf(stdout, "\nNon-normalised transport statistics (95%% Confidence-interval):\n");
    for (size_t m = 0; m < nmetrics; ++m) 
    {
        double mean, std, ci;
        compute_stats(metric_values[m], n_iter, &mean, &std, &ci);
        summary_mean[m] = mean;
        summary_ci[m] = ci;
        double rel = (mean != 0.0) ? fabs(ci / mean) * 100.0 : 0.0;
        fprintf(stdout, "  %-28s  %12.6E +/- %9.3E (%.6lf%%)\n",
                metrics[m].perhist_label, mean, ci, rel);
    }

    const double norm_factor = (GLOB.norm_factor > 0.0) ? GLOB.norm_factor : 1.0;
    double effective_norm = norm_factor;

    /* For power normalization, convert target fission rate to source strength (histories/s)
       using the simulated fissions per history. */
    if (GLOB.norm_mode == NORM_POWER)
    {
        size_t fission_idx = SIZE_MAX;
        for (size_t m = 0; m < nmetrics; ++m)
        {
            if (metrics[m].perhist_name && strcmp(metrics[m].perhist_name, "FISSIONS") == 0)
            {
                fission_idx = m;
                break;
            }
        }

        double fissions_per_history = (fission_idx < nmetrics) ? summary_mean[fission_idx] : 0.0;
        if (fissions_per_history > 0.0)
            effective_norm = norm_factor / fissions_per_history;
        else
            effective_norm = 0.0;
    }
    char norm_desc[MAX_STR_LEN];
    getNormalizationModeDescription(norm_desc);

    fprintf(stdout, "\nReaction rates (normalized to %s):\n", norm_desc);
    for (size_t m = 0; m < nmetrics; ++m)
    {
        if (!metrics[m].apply_norm)
            continue;

        double scaled_mean = summary_mean[m] * effective_norm;
        double scaled_ci = summary_ci[m] * effective_norm;
        double rel = (scaled_mean != 0.0) ? fabs(scaled_ci / scaled_mean) * 100.0 : 0.0;
        fprintf(stdout, "  %-28s  %12.6E +/- %9.3E (%.6lf%%)\n",
                metrics[m].rate_label ? metrics[m].rate_label : metrics[m].perhist_label,
                scaled_mean, scaled_ci, rel);
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
        fprintf(mfile, "%% Normalisation: %.8E (%s)\n\n", effective_norm, norm_desc);
        fprintf(mfile, "N_GENERATIONS = %zu;\n", n_iter);
        fprintf(mfile, "N_HISTORIES = %zu;\n\n", n_hist);

        fprintf(mfile, "%% Non-normalised statistics\n");
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
                    summary_mean[m] * effective_norm,
                    summary_ci[m] * effective_norm);
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
