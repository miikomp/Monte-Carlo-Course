#include "header.h"

bool createEnergyGrid(EnergyGrid *grid, double emin, double emax, size_t n_bins, EnergyBinSpacing spacing);

int processDetectors(void)
{
    /* Process detectors */

    fprintf(stdout, "\nProcessing detectors...\n");

    for (size_t i = 0; i < DATA.n_detectors; i++) 
    {
        if (DATA.detectors[i] == NULL) 
        {
            fprintf(stderr, "[ERROR] NULL detector definition found.\n");
            return EXIT_FAILURE;
        }

        /* Get pointer to detector */

        Detector *det = DATA.detectors[i];

        /* Process based on detector type */

        if (det->type == DETECTOR_TYPE_REACTION_RATE)
        {
            ReactionRateDetector *rr = &det->data.reaction_rate;

            /* Find material index based on provided material name */

            int mat_idx = -1;
            for (size_t j = 0; j < DATA.n_mats; ++j) 
            {
                if (!strcmp(rr->material_name, DATA.mats[j].name)) 
                {
                    mat_idx = (int)j;
                    break;
                }
            }

            if (mat_idx < 0) 
            {
                fprintf(stderr, "[ERROR] Material \"%s\" for detector %zu not found.\n", rr->material_name, i);
                return EXIT_FAILURE;
            }

            rr->material_index = mat_idx;

            size_t n_energy_bins = 1;
            if (rr->energy_grid.enabled)
            {
                if (!createEnergyGrid(&rr->energy_grid,
                                      rr->energy_grid.E_min,
                                      rr->energy_grid.E_max,
                                      rr->energy_grid.n_bins,
                                      rr->energy_grid.spacing))
                {
                    fprintf(stderr, "[ERROR] Failed to initialize energy grid for detector \"%s\".\n", det->name);
                    return EXIT_FAILURE;
                }

                if (rr->energy_grid.n_bins == 0)
                {
                    fprintf(stderr, "[ERROR] Detector \"%s\" energy grid has zero bins.\n", det->name);
                    return EXIT_FAILURE;
                }

                n_energy_bins = rr->energy_grid.n_bins;
            }
            rr->n_energy_bins = n_energy_bins;

            /* Create all nuclide bins */
            size_t n_nucs = DATA.mats[mat_idx].n_nucs;
            rr->n_nuclides = n_nucs;
            rr->nuclides = (ReactionRateNuclide*)calloc(n_nucs, sizeof(ReactionRateNuclide));
            if (!rr->nuclides) 
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                return EXIT_FAILURE;
            }

            /* Add data in nuclide bins */
            for (size_t j = 0; j < n_nucs; j++)
            {
                ReactionRateNuclide *det_nuc = &rr->nuclides[j];
                det_nuc->nuclide_index = j;
                snprintf(det_nuc->nuclide_name, sizeof(det_nuc->nuclide_name), "%s", DATA.mats[mat_idx].nucs[j].nuc_data.name);
                det_nuc->n_energy_bins = n_energy_bins;
                det_nuc->tallies = (ReactionRateTally*)calloc(n_energy_bins * RRDET_MODE_COUNT, sizeof(ReactionRateTally));
                if (!det_nuc->tallies)
                {
                    fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                    for (size_t k = 0; k < j; ++k)
                        free(rr->nuclides[k].tallies);
                    free(rr->nuclides);
                    rr->nuclides = NULL;
                    return EXIT_FAILURE;
                }
            }

            if (VERBOSITY >= 2)
            {
                fprintf(stdout, "Processed detector \"%s\" linked to material \"%s\" with %zu nuclide(s).\n", 
                        det->name, 
                        rr->material_name, 
                        rr->n_nuclides
                    );
            }
        }
        else if (det->type == DETECTOR_TYPE_ENERGY_SPECTRUM)
        {
            EnergySpectrumDetector *es = &det->data.energy_spectrum;

            /* Find material index based on provided material name if applicable */

            if (es->has_material_filter) 
            {
                int mat_idx = -1;
                for (size_t j = 0; j < DATA.n_mats; ++j) 
                {
                    if (!strcmp(es->material_name, DATA.mats[j].name)) 
                    {
                        mat_idx = (int)j;
                        break;
                    }
                }

                if (mat_idx < 0) 
                {
                    fprintf(stderr, "[ERROR] Material \"%s\" for detector %zu not found.\n", es->material_name, i);
                    return EXIT_FAILURE;
                }

                es->material_index = mat_idx;
            } 
            else 
            {
                es->material_index = -1;
            }

            /* Create energy bins */

            if (!es->grid.enabled)
            {
                fprintf(stderr, "[ERROR] Energy spectrum detector \"%s\" has no energy grid defined.\n", det->name);
                return EXIT_FAILURE;
            }

            if (!createEnergyGrid(&es->grid, es->grid.E_min, es->grid.E_max, es->grid.n_bins, es->grid.spacing))
            {
                fprintf(stderr, "[ERROR] Failed to initialize energy grid for detector \"%s\".\n", det->name);
                return EXIT_FAILURE;
            }

            size_t n_bins = es->grid.n_bins;
            es->track_length_sum = (double*)calloc(n_bins, sizeof(double));
            es->track_length_sum_sq = (double*)calloc(n_bins, sizeof(double));
            if (!es->track_length_sum || !es->track_length_sum_sq)
            {
                fprintf(stderr, "[ERROR] Memory allocation failed.\n");
                free(es->track_length_sum);
                free(es->track_length_sum_sq);
                free(es->grid.edges);
                es->grid.edges = NULL;
                return EXIT_FAILURE;
            }

            if (VERBOSITY >= 2)
            {
                const char *mat_desc = es->has_material_filter ? es->material_name : "all materials";
                fprintf(stdout, "Processed detector \"%s\" linked to %s with %zu bins between %e - %e\n", 
                        det->name, 
                        mat_desc, 
                        es->grid.n_bins,
                        es->grid.E_min,
                        es->grid.E_max
                    );
            }
        }
        else
        {
            fprintf(stderr, "[ERROR] Unsupported detector type %d for detector %zu.\n", det->type, i);
            return EXIT_FAILURE;
        }
        
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}

bool createEnergyGrid(EnergyGrid *grid, double emin, double emax, size_t n_bins, EnergyBinSpacing spacing)
{
    if (!grid || n_bins == 0)
        return false;

    if (emax <= emin)
        return false;

    if (spacing == ENERGY_BIN_SPACING_LOG && (emin <= 0.0 || emax <= 0.0))
        return false;

    double *edges = (double*)malloc((n_bins + 1) * sizeof(double));
    if (!edges)
        return false;

    if (spacing == ENERGY_BIN_SPACING_LINEAR)
    {
        double step = (emax - emin) / (double)n_bins;
        for (size_t idx = 0; idx <= n_bins; ++idx)
            edges[idx] = emin + step * (double)idx;
        edges[n_bins] = emax;
    }
    else
    {
        double ratio = pow(emax / emin, 1.0 / (double)n_bins);
        edges[0] = emin;
        for (size_t idx = 1; idx <= n_bins; ++idx)
            edges[idx] = edges[idx - 1] * ratio;
        edges[n_bins] = emax;
    }

    if (grid->edges)
    {
        free(grid->edges);
    }

    grid->edges = edges;
    grid->enabled = true;
    grid->n_bins = n_bins;
    grid->E_min = emin;
    grid->E_max = emax;
    grid->spacing = spacing;

    return true;
}
