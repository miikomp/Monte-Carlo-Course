#include <inttypes.h>

#include "header.h"

typedef struct {
    const char *label;
    const char *matlab_name;
    double (*extract)(const GenerationScores *gs);
} MetricSpec;

static double safe_norm(uint64_t count, double value) 
{
    return (count > 0) ? value / (double)count : 0.0;
}

static double metric_path_length(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, gs->total_path_length);
}

static double metric_fast_path_length(const GenerationScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_fast_path_length);
}

static double metric_time(const GenerationScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_time);
}

static double metric_time_fast(const GenerationScores *gs)
{
    return safe_norm(gs->n_histories, gs->total_time_fast);
}

static double metric_fission_yield(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, gs->total_fission_yield);
}

static double metric_collisions(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_collisions);
}

static double metric_captures(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_captures);
}

static double metric_elastic(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_elastic_scatters);
}

static double metric_inelastic(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_inelastic_scatters);
}

static double metric_fissions(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_fissions);
}

static double metric_fast_fissions(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_fast_fissions);
}

static double metric_leakages(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_leakages);
}

static double metric_unknowns(const GenerationScores *gs) 
{
    return safe_norm(gs->n_histories, (double)gs->total_unknowns);
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
    if (RES.n_generations <= 0 || !RES.avg_scores) 
    {
        fprintf(stdout, "\n[INFO] No transport results to process.\n");
        return;
    }

    const size_t ngen = (size_t)RES.n_generations;

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
        {"Fast Fissions",      "FAST_FISSIONS",        metric_fast_fissions},
        {"Leakages",           "LEAKAGES",             metric_leakages},
        {"Unknown outcomes",   "UNKNOWN_OUTCOMES",     metric_unknowns}
    };

    const size_t nmetrics = sizeof(metrics) / sizeof(metrics[0]);

    double *histories = (double*)calloc(ngen, sizeof(double));
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
        metric_values[m] = (double*)calloc(ngen, sizeof(double));
        if (!metric_values[m]) 
        {
            for (size_t j = 0; j < m; ++j) free(metric_values[j]);
            free(metric_values);
            free(histories);
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            return;
        }
    }

    size_t max_bin = 0;

    size_t max_time_bin = 0;
    double *time_yield = (double*)calloc(MAX_TIME_BINS, sizeof(double));
    uint64_t *time_events = (uint64_t*)calloc(MAX_TIME_BINS, sizeof(uint64_t));
    if (!time_yield || !time_events)
    {
        fprintf(stderr, "[ERROR] Memory allocation failed.\n");
        free(time_yield);
        free(time_events);
        for (size_t m = 0; m < nmetrics; ++m)
            free(metric_values[m]);
        free(metric_values);
        free(histories);
        return;
    }

    const GenerationScores *ref_gen = NULL;

    for (size_t g = 0; g < ngen; ++g)
    {
        const GenerationScores *gs = &RES.avg_scores[g];
        histories[g] = (double)gs->n_histories;

        for (size_t m = 0; m < nmetrics; ++m)
            metric_values[m][g] = metrics[m].extract(gs);

        for (size_t b = 0; b < MAX_COLLISION_BINS; ++b)
        {
            if (gs->collision_energy_count[b] == 0)
                continue;
            if (b + 1 > max_bin)
                max_bin = b + 1;
        }

        for (size_t b = 0; b < MAX_TIME_BINS; ++b)
        {
            double yield = gs->fission_time_yield[b];
            uint64_t events = gs->fission_time_events[b];
            if (yield <= 0.0 && events == 0)
                continue;
            time_yield[b] += yield;
            time_events[b] += events;
            if (b + 1 > max_time_bin)
                max_time_bin = b + 1;
        }

        if (!ref_gen && gs->collision_energy_count[0] > 0)
            ref_gen = gs;
    }
    if (ref_gen)
        max_bin = ref_gen->max_collision_bin;

    if (VERBOSITY >= 1)
    {
        fprintf(stdout, "\nTransport per-generation statistics:\n");
        fprintf(stdout, "  Generation  Histories\n");
        for (size_t g = 0; g < ngen; ++g)
            fprintf(stdout, "  %-10zu  %-10.0f\n", g + 1, histories[g]);
    }

    fprintf(stdout, "\nTransport per-history statistics (95%% Confidence-interval):\n");
    for (size_t m = 0; m < nmetrics; ++m) 
    {
        double mean, std, ci;
        compute_stats(metric_values[m], ngen, &mean, &std, &ci);
        fprintf(stdout, "  %-20s  mean = %12.6e +/- %9.3e (%.6lf%%)\n",
                metrics[m].label, mean, ci, (ci / mean) * 100.0);
    }
    

    if (VERBOSITY >= 2)
    {
        if (ref_gen && max_bin > 0) 
        {
            fprintf(stdout, "\nAverage neutron energy by collision index:\n");
            for (size_t b = 0; b < max_bin; ++b) 
            {
                uint64_t cnt = ref_gen->collision_energy_count[b];
                double avgE = (cnt > 0)
                    ? ref_gen->collision_energy_sum[b] / (double)cnt
                    : 0.0;
                fprintf(stdout, "  Collision %-3zu  <E> = %12.6e MeV  samples = %lu\n",
                        b + 1, avgE, cnt);
            }
        } 
        else {
            fprintf(stdout, "\nNo collision energy tallies were recorded.\n");
        }
    }

    /* ########################################################################################## */
    /* Write results to MATLAB-compatible .m file */

    char mpath[MAX_PATH];
    snprintf(mpath, sizeof(mpath), "%s_res.m", GLOB.fname ? GLOB.fname : "results");

    FILE *mfile = fopen(mpath, "w");
    if (!mfile) 
    {
        fprintf(stderr, "[ERROR] Could not open '%s' for writing.\n", mpath);
    } 
    else 
    {
        fprintf(mfile, "%%%% Transport results auto-generated\n");
        fprintf(mfile, "%% Input file: %s\n\n", GLOB.fname ? GLOB.fname : "<unknown>");
        fprintf(mfile, "N_GENERATIONS = %zu;\n\n", ngen);

        fprintf(mfile, "HISTORIES = [\n");
        for (size_t g = 0; g < ngen; ++g)
            fprintf(mfile, "  %.0f", histories[g]);
        fprintf(mfile, "];\n\n");

        for (size_t m = 0; m < nmetrics; ++m) 
        {
            fprintf(mfile, "%s = [\n", metrics[m].matlab_name);
            for (size_t g = 0; g < ngen; ++g)
                fprintf(mfile, "  %.12e", metric_values[m][g]);
            fprintf(mfile, "];\n\n");
        }

        if (ref_gen && max_bin > 0) {
            fprintf(mfile, "MEAN_ENERGY_PER_COLLISION = [\n");
            for (size_t b = 0; b < max_bin; ++b) 
            {
                uint64_t cnt = ref_gen->collision_energy_count[b];
                double avgE = (cnt > 0)
                    ? ref_gen->collision_energy_sum[b] / (double)cnt
                    : 0.0;
                fprintf(mfile, "  %.12e", avgE);
            }
            fprintf(mfile, "];\n\n");

            fprintf(mfile, "ENERGY_SAMPLES_PER_COLLISION = [\n");
            for (size_t b = 0; b < max_bin; ++b)
                fprintf(mfile, "  %lu", ref_gen->collision_energy_count[b]);
            fprintf(mfile, "];\n\n");
        }

        if (max_time_bin > 0)
        {
            fprintf(mfile, "FISSION_TIME_BIN_WIDTH = %.12e;\n", (double)TIME_BIN_WIDTH);
            fprintf(mfile, "FISSION_TIME_BIN_COUNT = %zu;\n\n", max_time_bin);

            fprintf(mfile, "FISSION_TIME_BIN_CENTERS = [\n");
            for (size_t b = 0; b < max_time_bin; ++b)
            {
                double center = (b + 0.5) * (double)TIME_BIN_WIDTH;
                fprintf(mfile, "  %.12e\n", center);
            }
            fprintf(mfile, "];\n\n");

            fprintf(mfile, "FISSION_NEUTRONS_PER_BIN = [\n");
            for (size_t b = 0; b < max_time_bin; ++b)
                fprintf(mfile, "  %.12e\n", time_yield[b]);
            fprintf(mfile, "];\n\n");

            fprintf(mfile, "FISSION_EVENTS_PER_BIN = [\n");
            for (size_t b = 0; b < max_time_bin; ++b)
                fprintf(mfile, "  %" PRIu64 "\n", time_events[b]);
            fprintf(mfile, "];\n\n");
        }

        fclose(mfile);
        fprintf(stdout, "\nWrote transport results to '%s'.\n", mpath);
    }

    for (size_t m = 0; m < nmetrics; ++m)
        free(metric_values[m]);
    free(metric_values);
    free(histories);
    free(time_yield);
    free(time_events);
}
