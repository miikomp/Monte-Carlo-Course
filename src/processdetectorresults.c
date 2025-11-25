#include "header.h"

static void makeDetectorVariableName(const char *src, char *dst, size_t len)
{
    if (!dst || len == 0)
        return;

    size_t j = 0;
    while (src && *src && j < len - 1)
    {
        char c = (char)toupper((unsigned char)*src++);
        if (!isalnum((unsigned char)c))
            c = '_';
        dst[j++] = c;
    }
    dst[j] = '\0';
}

static void writeBoundsTable(FILE *fp, const char *var, const char *suffix, const double *edges, size_t n_bins)
{
    if (!fp || !var || !suffix || !edges || n_bins == 0)
        return;

    fprintf(fp, "%s%s = [\n", var, suffix);
    for (size_t i = 0; i < n_bins; ++i)
    {
        double emin = edges[i];
        double emax = edges[i + 1];
        double emid = 0.5 * (emin + emax);
        fprintf(fp, "  % .12E % .12E % .12E\n", emin, emax, emid);
    }
    fprintf(fp, "];\n\n");
}

void processDetectorResults(void)
{
    const size_t n_detectors = DATA.n_detectors;

    if (n_detectors == 0)
    {
        if (VERBOSITY >= 1)
            fprintf(stdout, "\n[INFO] No detectors defined.\n");
        return;
    }

    const double norm_factor = (GLOB.norm_factor > 0.0) ? GLOB.norm_factor : 1.0;

    char norm_desc[MAX_STR_LEN];
    getNormalizationModeDescription(norm_desc, norm_factor);

    fprintf(stdout, "\nProcessing detector results (normalized to %s):\n", norm_desc);

    char mpath[MAX_PATH];
    snprintf(mpath, sizeof(mpath), "%s_det.m", GLOB.inputfname);

    FILE *mfile = fopen(mpath, "w");
    if (!mfile)
    {
        fprintf(stderr, "[ERROR] Could not open '%s' for detector output.\n", mpath);
    }
    else
    {
        fprintf(mfile, "%%%% Detector results auto-generated\n");
        fprintf(mfile, "%% Input file: %s\n", GLOB.inputfname);
        fprintf(mfile, "%% Normalisation: %.12E (%s)\n\n", norm_factor, norm_desc);
    }

    for (size_t d = 0; d < n_detectors; ++d)
    {
        Detector *det = &DATA.detectors[d];
        uint64_t histories = det->n_histories;

        if (!mfile)
            continue;

        char det_var[2 * MAX_STR_LEN];
        makeDetectorVariableName(det->name, det_var, sizeof(det_var));

        size_t n_time = det->time.active ? det->time.n_bins : 1u;
        size_t n_energy = det->energy.active ? det->energy.n_bins : 1u;
        size_t n_x = det->mesh_x.active ? det->mesh_x.n_bins : 1u;
        size_t n_y = det->mesh_y.active ? det->mesh_y.n_bins : 1u;
        size_t n_z = det->mesh_z.active ? det->mesh_z.n_bins : 1u;
        size_t n_resp = det->responses.n_bins;

        fprintf(mfile, "%s = [\n", det_var);

        const DetectorAxisLayout *layout = det->layout.axes;

        for (size_t t = 0; t < n_time; ++t)
        {
            size_t t_offset = det->time.active ? t * layout[DETECTOR_AXIS_TIME].stride : 0u;
            unsigned long t_col = det->time.active ? (unsigned long)(t + 1u) : 1u;

            for (size_t e = 0; e < n_energy; ++e)
            {
                size_t e_offset = det->energy.active ? e * layout[DETECTOR_AXIS_ENERGY].stride : 0u;
                unsigned long e_col = det->energy.active ? (unsigned long)(e + 1u) : 1u;

                for (size_t x = 0; x < n_x; ++x)
                {
                    size_t x_offset = det->mesh_x.active ? x * layout[DETECTOR_AXIS_MESH_X].stride : 0u;
                    unsigned long x_col = det->mesh_x.active ? (unsigned long)(x + 1u) : 1u;

                    for (size_t y = 0; y < n_y; ++y)
                    {
                        size_t y_offset = det->mesh_y.active ? y * layout[DETECTOR_AXIS_MESH_Y].stride : 0u;
                        unsigned long y_col = det->mesh_y.active ? (unsigned long)(y + 1u) : 1u;

                        for (size_t z = 0; z < n_z; ++z)
                        {
                            size_t z_offset = det->mesh_z.active ? z * layout[DETECTOR_AXIS_MESH_Z].stride : 0u;
                            unsigned long z_col = det->mesh_z.active ? (unsigned long)(z + 1u) : 1u;

                            for (size_t r = 0; r < n_resp; ++r)
                            {
                                size_t r_offset = r * layout[DETECTOR_AXIS_RESPONSE].stride;
                                size_t idx = t_offset + e_offset + x_offset + y_offset + z_offset + r_offset;

                                double sum = det->scores ? det->scores[idx] : 0.0;
                                double mean = (histories > 0) ? sum / (double)histories : 0.0;
                                double err = 0.0;
                                if (det->scores_sq && histories > 1)
                                {
                                    double sum_sq = det->scores_sq[idx];
                                    double num = sum_sq - (sum * sum) / (double)histories;
                                    if (num < 0.0)
                                        num = 0.0;
                                    double sigma = sqrt(num / (double)(histories - 1));
                                    err = sigma / sqrt((double)histories);
                                }

                                double rel_err = 0.0;
                                if (fabs(mean) > 0.0)
                                    rel_err = err / fabs(mean);

                                double scaled_mean = mean * norm_factor;

                                fprintf(mfile, "  %3lu %3lu %3lu %3lu %3lu %3lu % .8E % .5lf\n",
                                        t_col, e_col, (unsigned long)(r + 1u),
                                        x_col, y_col, z_col, scaled_mean, rel_err);
                            }
                        }
                    }
                }
            }
        }

        fprintf(mfile, "];\n\n");

        if (det->time.active)
            writeBoundsTable(mfile, det_var, "T", det->time.edges, det->time.n_bins);
        if (det->energy.active)
            writeBoundsTable(mfile, det_var, "E", det->energy.edges, det->energy.n_bins);
        if (det->mesh_x.active)
            writeBoundsTable(mfile, det_var, "X", det->mesh_x.edges, det->mesh_x.n_bins);
        if (det->mesh_y.active)
            writeBoundsTable(mfile, det_var, "Y", det->mesh_y.edges, det->mesh_y.n_bins);
        if (det->mesh_z.active)
            writeBoundsTable(mfile, det_var, "Z", det->mesh_z.edges, det->mesh_z.n_bins);
    }

    if (mfile)
    {
        fprintf(stdout, "\nWrote detector results to '%s'.\n", mpath);
        fclose(mfile);
    }
}
