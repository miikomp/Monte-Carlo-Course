#include "header.h"

int resolveCells() {

    fprintf(stdout, "\nProcessing geometry cells...\n");

    for (size_t c = 0; c < DATA.n_cells; c++)
    {
        Cell *cell = &DATA.cells[c];

        /* Cell can be filled with universe or material */
        if (cell->unifilled)
        {
            /* Check if given universe exists */

            for (size_t u = 0; u < DATA.n_unis; u++)
            {
                if (!strcmp(cell->filluni_name, DATA.unis[u].name))
                {
                    cell->filluni_idx = (int)u;
                    break;
                }
            }
            if (cell->filluni_idx == -1)
            {
                fprintf(stdout, "[ERROR] Universe \"%s\" for cell \"%s\" not found.\n", cell->filluni_name, cell->name);
                return EXIT_FAILURE;
            }
        }
        else
        {
            /* Check if given material exists */

            for (size_t m = 0; m < DATA.n_mats; m++)
            {
                if (!strcmp(cell->mat_name, DATA.mats[m].name))
                {
                    cell->mat_idx = (int)m;
                    break;
                }
            }
            if (cell->mat_idx == -1)
            {
                if (!strcmp(cell->mat_name, "outside"))
                    cell->mat_idx = -1; // outside cell
                else
                {
                    fprintf(stdout, "[ERROR] Material \"%s\" for cell \"%s\" not found.\n", cell->mat_name, cell->name);
                    return EXIT_FAILURE;
                }
            }
        }
        /* Check that all surfaces are valid */
        for (size_t s = 0; s < cell->n_surfs; s++)
        {
            /* Move pointer to correct surface name */

            char *surf_name = cell->surf_names + s * MAX_STR_LEN;

            /* Find corresponding surface */

            for (size_t s0 = 0; s0 < DATA.n_surf; s0++)
            {
                if (!strcmp(surf_name, DATA.surfs[s0].name))
                {
                    cell->surf_idxs[s] = (int)s0;
                    break;
                }
            }
            if (cell->surf_idxs[s] == -1)
            {
                fprintf(stdout, "[ERROR] Surface \"%s\" for cell \"%s\" not found.\n", surf_name, cell->name);
                return EXIT_FAILURE;
            }
        }
    }

    fprintf(stdout, "DONE.\n");

    return EXIT_SUCCESS;
}
