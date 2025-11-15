#include <inttypes.h>

#include "header.h"

typedef struct {
    const char *label;
    const char *matlab_name;
    double (*extract)(const TransportRunScores *gs);
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
        {"Path length",        "PATH_LENGTH",          metric_path_length},
        {"Fast path length",   "FAST_PATH_LENGTH",     metric_fast_path_length},
        {"Flight time",        "FLIGHT_TIME",          metric_time},
        {"Fast flight time",   "FAST_FLIGHT_TIME",     metric_time_fast},
        {"Fission yield",      "FISSION_YIELD",        metric_fission_yield},
        {"Collisions",         "COLLISIONS",           metric_collisions},
        {"Captures",           "CAPTURES",             metric_captures},
        {"Elastic scatters",   "ELASTIC_SCATTERS",     metric_elastic},
        {"Inelastic scatters", "INELASTIC_SCATTERS",   metric_inelastic},
        {"Fissions",           "FISSIONS",             metric_fissions},
        {"Thermal Fissions",   "THERMAL_FISSIONS",     metric_thermal_fissions},
        {"Fast Fissions",      "FAST_FISSIONS",        metric_fast_fissions},
        {"Leakages",           "LEAKAGES",             metric_leakages},
        {"Terminated",         "TERMINATED",           metric_terminated},
        {"Unknown outcomes",   "UNKNOWN_OUTCOMES",     metric_unknowns}
    };

    const size_t nmetrics = sizeof(metrics) / sizeof(metrics[0]);

    double *histories = (double*)calloc(n_iter, sizeof(double));
    if (!histories) 
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return;
    }

    double **metric_values = (double**)calloc(nmetrics, sizeof(double*));
    if (!metric_values) 
    {
        free(histories);
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        return;
    }
    for (size_t m = 0; m < nmetrics; ++m) 
    {
        metric_values[m] = (double*)calloc(n_iter, sizeof(double));
        if (!metric_values[m]) 
        {
            for (size_t j = 0; j < m; ++j) free(metric_values[j]);
            free(metric_values);
            free(histories);
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
        fprintf(stdout, "  %-20s  mean = %12.6e +/- %9.3e (%.6lf%%)\n",
                metrics[m].label, mean, ci, (ci / mean) * 100.0);
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
        fprintf(mfile, "%% Input file: %s\n\n", GLOB.inputf);
        fprintf(mfile, "N_GENERATIONS = %zu;\n\n", n_iter);

        fprintf(mfile, "HISTORIES = [\n");
        for (size_t g = 0; g < n_iter; ++g)
            fprintf(mfile, "  %.0f", histories[g]);
        fprintf(mfile, "];\n\n");

        for (size_t m = 0; m < nmetrics; ++m) 
        {
            fprintf(mfile, "%s = [\n", metrics[m].matlab_name);
            for (size_t g = 0; g < n_iter; ++g)
                fprintf(mfile, "  %.12e", metric_values[m][g]);
            fprintf(mfile, "];\n\n");
        }

        fclose(mfile);
        fprintf(stdout, "\nWrote transport results to '%s'.\n", mpath);
    }

    for (size_t m = 0; m < nmetrics; ++m)
        free(metric_values[m]);
    free(metric_values);
    free(histories);
}
