#include "header.h"

int resolveTransformations() {

    fprintf(stdout, "\nProcessing transformations...\n");

    /* Loop over all transformations */

    for (size_t t = 0; t < DATA.n_transforms; t++)
    {
        Transform* T = &DATA.transforms[t];

        switch (T->type)
        {
            case TRA_SURFACE:
            {
                /* Find targeted surface */

                bool found = false;
                for (size_t s = 0; s < DATA.n_surf; s++)
                {
                    Surface* S = &DATA.surfs[s];

                    if (!strcmp(T->target_name, S->name))
                    {
                        if (S->t_idx >= 0)
                            fprintf(stdout, "[WARNING] Multiple transformations applied to surface \"%s\".\n", S->name);
                    
                        S->t_idx = t;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    fprintf(stderr, "[ERROR] Transformation references unknown surface \"%s\".\n", T->target_name);
                    return EXIT_FAILURE;
                }

                break;
            }
            case TRA_UNIVERSE:
            {
                /* Find targeted universe */

                bool found = false;
                for (size_t u = 0; u < DATA.n_unis; u++)
                {
                    Universe *U = &DATA.unis[u];

                    if (!strcmp(T->target_name, U->name))
                    {
                        if (U->t_idx >= 0)
                            fprintf(stdout, "[WARNING] Multiple transformation applied to universe \"%s\"\n", U->name);
                        
                        if (u == 0)
                        {
                            fprintf(stdout, "[ERROR] Universe transformation cannot be applied to root universe 0\n");
                            return EXIT_FAILURE;
                        }

                        U->t_idx = t;
                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    fprintf(stderr, "[ERROR] Transformation references unknown universe \"%s\".\n", T->target_name);
                    return EXIT_FAILURE;
                }

                break;
            }
            default:
            {
                fprintf(stderr, "[ERROR] Unknown transformation type (%u)\n", T->type);
                return EXIT_FAILURE;
            }
        }
    }

    fprintf(stdout, "DONE.\n");
    return EXIT_SUCCESS;
}