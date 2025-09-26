#include "header.h"

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

        ReactionRateDetector *det = DATA.detectors[i];

        /* Find material index based on provided material name */

        int mat_idx = -1;
        for (size_t j = 0; j < DATA.n_mats; ++j) 
        {
            if (!strcmp(det->material_name, DATA.mats[j].name)) 
            {
                mat_idx = (int)j;
                break;
            }
        }

        if (mat_idx < 0) 
        {
            fprintf(stderr, "[ERROR] Material \"%s\" for detector %zu not found.\n", det->material_name, i);
            return EXIT_FAILURE;
        }

        det->material_index = mat_idx;

        /* Create all nuclide bins */
        size_t n_nucs = DATA.mats[mat_idx].n_nucs;
        det->n_nuclides = n_nucs;
        det->nuclides = (ReactionRateNuclide*)calloc(n_nucs, sizeof(ReactionRateNuclide));
        if (!det->nuclides) 
        {
            fprintf(stderr, "[ERROR] Memory allocation failed.\n");
            return EXIT_FAILURE;
        }

        /* Add data in nuclide bins */
        for (size_t j = 0; j < n_nucs; j++)
        {
            ReactionRateNuclide *det_nuc = &det->nuclides[j];
            det_nuc->nuclide_index = j;
            snprintf(det_nuc->nuclide_name, sizeof(det_nuc->nuclide_name), "%s", DATA.mats[mat_idx].nucs[j].nuc_data.name);
        }

        if (VERBOSITY >= 2)
        {
            fprintf(stdout, "Processed detector \"%s\" linked to material \"%s\" with %zu nuclide(s).\n", 
                    det->name, 
                    det->material_name, 
                    det->n_nuclides
                );
        }
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}