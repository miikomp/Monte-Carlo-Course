#include "header.h"

int resolveDetectors(void)
{
    fprintf(stdout, "\nResolving detectors...\n");

    for (size_t i = 0; i < DATA.n_detectors; ++i)
    {
        Detector *det = &DATA.detectors[i];
        det->n_histories = 0u;
        if (det->has_material_filter)
        {
            long idx = -1;
            for (size_t m = 0; m < DATA.n_mats; ++m)
            {
                if (!strcmp(DATA.mats[m].name, det->material_filter_name))
                {
                    idx = (long)m;
                    break;
                }
            }

            if (idx < 0)
            {
                fprintf(stderr, "[ERROR] Material filter \"%s\" for detector \"%s\" not found.\n",
                        det->material_filter_name, det->name);
                return EXIT_FAILURE;
            }

            det->material_filter_index = idx;
        }
        else
        {
            det->material_filter_index = -1;
        }

        if (!det->responses.active || det->responses.n_bins == 0 || det->responses.bins == NULL)
        {
            fprintf(stderr, "[ERROR] Detector \"%s\" has no response bins defined.\n", det->name);
            return EXIT_FAILURE;
        }

        /* Resolve response targets */
        for (size_t r = 0; r < det->responses.n_bins; ++r)
        {
            DetectorResponseBin *bin = &det->responses.bins[r];
            if (bin->response_id < 0)
            {
                if (bin->rmat_name[0] == '\0')
                {
                    fprintf(stderr, "[ERROR] Detector \"%s\" reaction %ld requires a material name.\n",
                            det->name, bin->response_id);
                    return EXIT_FAILURE;
                }

                long mat_idx = -1;
                for (size_t m = 0; m < DATA.n_mats; ++m)
                {
                    if (!strcmp(DATA.mats[m].name, bin->rmat_name))
                    {
                        mat_idx = (long)m;
                        break;
                    }
                }

                if (mat_idx < 0)
                {
                    fprintf(stderr, "[ERROR] Material \"%s\" referenced by detector \"%s\" not found.\n",
                            bin->rmat_name, det->name);
                    return EXIT_FAILURE;
                }

                bin->rmat_id = mat_idx;
            }
            else if (bin->response_id > 0)
            {
                if (bin->rmat_name[0] == '\0')
                {
                    fprintf(stderr, "[ERROR] Detector \"%s\" microscopic reaction %ld requires a nuclide ZA identifier.\n",
                            det->name, bin->response_id);
                    return EXIT_FAILURE;
                }

                char *endptr = NULL;
                long za = strtol(bin->rmat_name, &endptr, 10);
                if (endptr == bin->rmat_name || *endptr != '\0' || za <= 0)
                {
                    fprintf(stderr, "[ERROR] Invalid nuclide identifier \"%s\" in detector \"%s\".\n",
                            bin->rmat_name, det->name);
                    return EXIT_FAILURE;
                }

                bool has_nuclide = false;
                for (size_t m = 0; m < DATA.n_mats && !has_nuclide; ++m)
                {
                    Material *mat = &DATA.mats[m];
                    for (size_t n = 0; n < mat->n_nucs; ++n)
                    {
                        if (mat->nucs[n].nuc_data.ZA == za)
                        {
                            has_nuclide = true;
                            break;
                        }
                    }
                }

                if (!has_nuclide)
                {
                    fprintf(stderr, "[ERROR] Nuclide ZA %ld referenced by detector \"%s\" not found in any material.\n",
                            za, det->name);
                    return EXIT_FAILURE;
                }

                bin->rmat_id = za;
            }
            else
            {
                bin->rmat_id = -1;
            }
        }

        /* Resolve time axis */
        if (det->time.active)
        {
            if (det->time.n_bins == 0 || det->time.max <= det->time.min)
            {
                fprintf(stderr, "[ERROR] Invalid time bin specification for detector \"%s\".\n", det->name);
                return EXIT_FAILURE;
            }

            size_t n = det->time.n_bins;
            double step = (det->time.max - det->time.min) / (double)n;
            det->time.edges = (double*)malloc((n + 1u) * sizeof(double));
            if (!det->time.edges)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed for detector \"%s\" time bins.\n", det->name);
                return EXIT_FAILURE;
            }

            det->time.edges[0] = det->time.min;
            for (size_t b = 1; b <= n; ++b)
                det->time.edges[b] = det->time.edges[b - 1] + step;
            det->time.edges[n] = det->time.max;
        }
        else
        {
            det->time.edges = NULL;
        }

        /* Resolve energy axis */
        if (det->energy.active)
        {
            if (det->energy.n_bins == 0 || det->energy.max <= det->energy.min)
            {
                fprintf(stderr, "[ERROR] Invalid energy bin specification for detector \"%s\".\n", det->name);
                return EXIT_FAILURE;
            }
            if (det->energy.type == ENERGY_GRID_LOG && det->energy.min <= 0.0)
            {
                fprintf(stderr, "[ERROR] Logarithmic energy bins for detector \"%s\" require Emin > 0.\n", det->name);
                return EXIT_FAILURE;
            }

            size_t n = det->energy.n_bins;
            det->energy.edges = (double*)malloc((n + 1u) * sizeof(double));
            if (!det->energy.edges)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed for detector \"%s\" energy bins.\n", det->name);
                return EXIT_FAILURE;
            }

            if (det->energy.type == ENERGY_GRID_LOG)
            {
                double ratio = pow(det->energy.max / det->energy.min, 1.0 / (double)n);
                det->energy.edges[0] = det->energy.min;
                for (size_t b = 1; b <= n; ++b)
                    det->energy.edges[b] = det->energy.edges[b - 1] * ratio;
                det->energy.edges[n] = det->energy.max;
            }
            else
            {
                double step = (det->energy.max - det->energy.min) / (double)n;
                for (size_t b = 0; b <= n; ++b)
                    det->energy.edges[b] = det->energy.min + step * (double)b;
                det->energy.edges[n] = det->energy.max;
            }
        }
        else
        {
            det->energy.edges = NULL;
        }

        /* Resolve spatial meshes */
        DetectorCartesianAxis *mesh_axes[3] = { &det->mesh_x, &det->mesh_y, &det->mesh_z };
        const char *mesh_labels[3] = { "dx", "dy", "dz" };

        for (size_t m = 0; m < 3; ++m)
        {
            DetectorCartesianAxis *axis = mesh_axes[m];
            if (!axis->active)
            {
                axis->edges = NULL;
                continue;
            }

            if (axis->n_bins == 0 || axis->max <= axis->min)
            {
                fprintf(stderr, "[ERROR] Invalid '%s' specification for detector \"%s\".\n",
                        mesh_labels[m], det->name);
                return EXIT_FAILURE;
            }

            if (axis->type != DETECTOR_MESH_AXIS_UNIFORM)
            {
                fprintf(stderr, "[ERROR] Detector \"%s\" uses unsupported mesh axis type for '%s'.\n",
                        det->name, mesh_labels[m]);
                return EXIT_FAILURE;
            }

            double step = (axis->max - axis->min) / (double)axis->n_bins;
            axis->spacing = step;
            axis->origin = axis->min;
            axis->edges = (double*)malloc((axis->n_bins + 1u) * sizeof(double));
            if (!axis->edges)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed for detector \"%s\" '%s' bins.\n",
                        det->name, mesh_labels[m]);
                return EXIT_FAILURE;
            }

            axis->edges[0] = axis->min;
            for (size_t b = 1; b <= axis->n_bins; ++b)
                axis->edges[b] = axis->edges[b - 1] + step;
            axis->edges[axis->n_bins] = axis->max;
        }

        /* Build bin layout and allocate tallies */
        bool axis_active[DETECTOR_AXIS_COUNT] = {
            det->time.active,
            det->energy.active,
            det->mesh_x.active,
            det->mesh_y.active,
            det->mesh_z.active,
            true
        };

        size_t axis_bins[DETECTOR_AXIS_COUNT] = {
            det->time.active   ? det->time.n_bins   : 1u,
            det->energy.active ? det->energy.n_bins : 1u,
            det->mesh_x.active ? det->mesh_x.n_bins : 1u,
            det->mesh_y.active ? det->mesh_y.n_bins : 1u,
            det->mesh_z.active ? det->mesh_z.n_bins : 1u,
            det->responses.n_bins
        };

        size_t total_bins = 1u;
        for (size_t k = 0; k < DETECTOR_AXIS_COUNT; ++k)
        {
            size_t count = axis_bins[k];
            if (count == 0)
            {
                fprintf(stderr, "[ERROR] Detector \"%s\" axis %zu has zero bins.\n", det->name, k);
                return EXIT_FAILURE;
            }

            DetectorAxisLayout *layout = &det->layout.axes[k];
            layout->active = axis_active[k];
            layout->n_bins = count;
            layout->stride = total_bins;

            if (count > 0 && total_bins > SIZE_MAX / count)
            {
                fprintf(stderr, "[ERROR] Detector \"%s\" bin count exceeds supported size.\n", det->name);
                return EXIT_FAILURE;
            }
            total_bins *= count;
        }
        det->layout.total_bins = total_bins;

        det->scores = (double*)calloc(total_bins, sizeof(double));
        det->scores_sq = (double*)calloc(total_bins, sizeof(double));
        if (!det->scores || !det->scores_sq)
        {
            free(det->scores);
            free(det->scores_sq);
            det->scores = NULL;
            det->scores_sq = NULL;
            fprintf(stderr, "[ERROR] Memory allocation failed for detector \"%s\" tallies.\n", det->name);
            return EXIT_FAILURE;
        }

        if (VERBOSITY >= 2)
            fprintf(stdout, "Resolved detector \"%s\" with %zu total bins.\n", det->name, total_bins);
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}
